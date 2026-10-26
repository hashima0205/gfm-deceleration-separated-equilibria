# 4gfm_simulation_RL.jl
using PowerDynamics
using PowerDynamics: Library
using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t
using PythonCall

include("compornent.jl")

const SIM_STATE = Dict{String, Any}()

# =======================================================
# 物理パラメータの定数定義
# =======================================================
const f0 = 60.0
const ω0_val = 2π * f0
const Lf_val = 0.0196
const C_val  = 0.1086
const Lg_val = 0.0209
const Kq_val = 0.05
const tau_p  = 0.02
const tau_q  = 0.02

# =======================================================
# コンポーネント定義 (名前被りによるバグを修正)
# =======================================================
# 指定時刻で事故を切り替える。MTKの記号式でも使えるifelseで記述する。
# 有効区間は [t_fault, t_clear): 投入時刻にON、解消時刻にOFF。
fault_conductance(t, G, t_fault, t_clear) =
    ifelse(t < t_fault, 0.0, ifelse(t < t_clear, G, 0.0))

@component function FaultLoad(; name, P_load=0.0, Q_load=0.0, G_fault_val=0.0, t_clear=1.0, defaults...)
    @named terminal = Terminal()
    @parameters begin
        # MTKのバグを防ぐため、引数名とパラメータ名を明確に分ける
        P_load_p = P_load
        Q_load_p = Q_load
        t_fault_p = 1.0     
        t_clear_p = t_clear   
        G_fault_p = G_fault_val
        # -1: 時刻のifelse式。数値積分時はイベントで0/1を切り替える。
        fault_enabled_p = -1.0
    end
    vars = @variables begin
        V_mag(t)
        G_fault_var(t)
    end
    eqs = [
        V_mag ~ sqrt(terminal.u_r^2 + terminal.u_i^2 + 1e-6)
        G_fault_var ~ ifelse(fault_enabled_p < 0,
                            fault_conductance(t, G_fault_p, t_fault_p, t_clear_p),
                            ifelse(fault_enabled_p > 0, G_fault_p, 0.0))
        terminal.i_r ~ (P_load_p * terminal.u_r + Q_load_p * terminal.u_i) / V_mag^2 + G_fault_var * terminal.u_r
        terminal.i_i ~ (P_load_p * terminal.u_i - Q_load_p * terminal.u_r) / V_mag^2 + G_fault_var * terminal.u_i
    ]
    sys = System(eqs, t; name, systems=[terminal])
    set_mtk_defaults!(sys, defaults)
    return sys
end

# 不連続点の直前の積分には旧状態を使い、到達後に事故状態を切り替える。
# NetworkDynamicsのイベント処理がパラメータ履歴を保存し、DAEの整合性を更新する。
# これにより、強い事故でもtanhによる平滑化や時間のずらしを必要としない。
function solve_fault_switched(prob, alg; kwargs...)
    s = NWState(prob.f, copy(prob.u0), copy(prob.p), prob.tspan[1])
    enabled = :fault_load₊fault_enabled_p
    onset = :fault_load₊t_fault_p
    clearing = :fault_load₊t_clear_p
    t_on = s.p[VIndex(4, onset)]
    t_off = s.p[VIndex(4, clearing)]
    s.p[VIndex(4, enabled)] = if t_on <= prob.tspan[1] < t_off
        1.0
    else
        0.0
    end
    events = sort(unique(filter(tt -> prob.tspan[1] < tt <= prob.tspan[2], [t_on, t_off])))
    if isempty(events)
        switched_prob = ODEProblem(s.nw, s, prob.tspan)
    else
        affect = NetworkDynamics.ComponentAffect([], [enabled, onset, clearing]) do u, p, ctx
            p[enabled] = if p[onset] <= ctx.t < p[clearing]
                1.0
            else
                0.0
            end
        end
        cb = NetworkDynamics.PresetTimeComponentCallback(events, affect)
        switched_prob = ODEProblem(s.nw, s, prob.tspan; add_comp_cb=VIndex(4)=>cb)
    end
    return solve(switched_prob, alg; kwargs...)
end

@component function SimpleLoad(; name, P_load=0.0, Q_load=0.0, defaults...)
    @named terminal = Terminal()
    @parameters begin
        P_load_p = P_load
        Q_load_p = Q_load
    end
    vars = @variables begin
        V_mag(t)
    end
    eqs = [
        V_mag ~ sqrt(terminal.u_r^2 + terminal.u_i^2 + 1e-6)
        terminal.i_r ~ -(P_load_p * terminal.u_r + Q_load_p * terminal.u_i) / V_mag^2
        terminal.i_i ~ -(P_load_p * terminal.u_i - Q_load_p * terminal.u_r) / V_mag^2
    ]
    sys = System(eqs, t; name, systems=[terminal])
    set_mtk_defaults!(sys, defaults)
    return sys
end

# =======================================================
# バス構築用ヘルパー関数
# =======================================================
function make_gfm(idx, P_set, V_set, Kp_val)
    name_sym = Symbol("inv", idx)
    @named inv_model = DroopInverter_addlimit(
        name = name_sym, filter_type = :LCL,
        vsrc₊Lf = Lf_val, vsrc₊C = C_val, vsrc₊Lg = Lg_val, vsrc₊ω0 = ω0_val,
        droop₊Kp = Kp_val, droop₊Kq = Kq_val, droop₊τ_p = tau_p, droop₊τ_q = tau_q, droop₊ω0 = ω0_val, droop₊Vset = V_set
    )
    bus = compile_bus(MTKBus(inv_model); vidx=idx, name=Symbol("bus", idx), assume_io_coupling=true)
    set_pfmodel!(bus, pfPV(P=P_set, V=1.0)) 
    return bus
end

function make_load(idx, P_load, Q_load)
    name_sym = Symbol("load", idx)
    @named PQLoad = SimpleLoad(; P_load=P_load, Q_load=Q_load)
    bus = compile_bus(MTKBus(PQLoad); vidx=idx, name=Symbol("bus", idx), assume_io_coupling=true)
    set_pfmodel!(bus, pfPQ(P=-P_load, Q=-Q_load)) 
    return bus
end

function make_fault_node(idx, G_fault, t_cle)
    name_sym = Symbol("fault", idx)
    @named fault_load = FaultLoad(; P_load=0.0, Q_load=0.0, G_fault_val=G_fault, t_clear=t_cle)
    bus = compile_bus(MTKBus(fault_load); vidx=idx, name=Symbol("bus", idx))
    set_pfmodel!(bus, pfPQ(P=0.0, Q=0.0))
    return bus
end

function make_slack(idx, V, δ)
    name_sym = Symbol("slack", idx)
    @named symbolic_slack = Library.VδConstraint(; V=V, δ=δ)
    bus = compile_bus(MTKBus(symbolic_slack); vidx=idx, name=Symbol("bus", idx))
    set_pfmodel!(bus, pfSlack(V=V))  
    return bus
end

# =======================================================
# 初期化関数
# =======================================================
function init_sim(p_set, p_set_near_bus, p_set_far_bus; Kp_gain=0.02, t_cle=1.6, G_fault=0.0, P_load=0.4, max_t=3.0)
    p_set_list = Dict{Int, Any}(
        1 => p_set_near_bus,
        2 => p_set,
        3 => p_set,
        7 => p_set_far_bus
    )
    Kp_val = Kp_gain * ω0_val

    buses = Any[]
    push!(buses, make_gfm(1, p_set_list[1], 1.0, Kp_val))
    push!(buses, make_gfm(2, p_set_list[2], 1.0, Kp_val))
    push!(buses, make_gfm(3, p_set_list[3], 1.0, Kp_val))
    push!(buses, make_fault_node(4, G_fault, t_cle))
    push!(buses, make_load(5, P_load, 0.0))
    push!(buses, make_slack(6, 1.0, 0.0))
    push!(buses, make_gfm(7, p_set_list[7], 1.0, Kp_val))
    push!(buses, make_load(8, 0.0, 0.0))

    y_line = [
        6.39, 1.69, 1.91, 2.14, 2.17, 2.21, 8.96, 1.80, 0.68, 1.50,
        1.90, 1.48, 2.90, 2.14, 3.43, 4.47, 1.40, 1.96, 1.89, 1.08
    ]
    lines = Any[]
    theta_rad = 75.0 * π / 180.0
    line_data = [
        (1,4,y_line[14]), (4,6,y_line[15]), (5,6,y_line[9]), (3,5,y_line[6]), (2,3,y_line[3]),
        (2,5,y_line[4]), (4,5,y_line[8]), (2,8,y_line[5]), (2,7,y_line[1]), (7,8,y_line[2]), (5,8,y_line[7])
    ]

    for (idx, (src_node, dst_node, y_mag)) in enumerate(line_data)
        z_mag = 1.0 / y_mag
        r_val = z_mag * cos(theta_rad)
        x_val = z_mag * sin(theta_rad)
        name_sym = Symbol("line", idx)
        @named line_sym = Library.PiLine(R=r_val, X=x_val)
        l = compile_line(MTKLine(line_sym); src=src_node, dst=dst_node, name=name_sym)
        push!(lines, l)
    end

    nw = Network(buses, lines; warn_order=false)
       
    pfnw = PowerDynamics.powerflow_model(nw)
    pfs0 = NWState(pfnw)
    pfs = PowerDynamics.solve_powerflow(nw; pfnw=pfnw, pfs0=pfs0, verbose=false)
    interf = NetworkDynamics.interface_values(pfs)
    pfconstraints = PowerDynamics.specialize_pfinitconstraints(nw, pfs)
    pfformulas = PowerDynamics.specialize_pfinitformulas(nw, pfs)

    s0 = NetworkDynamics.initialize_componentwise(
        nw;
        default_overrides = interf,
        additional_initconstraint = pfconstraints,
        additional_initformula = pfformulas,
        verbose = false,
        subverbose = false,
        tol = 1e-5,   
        nwtol = 1e-5, 
        t = 0.0
    )

    # ===========================================================
    # 修正: 初期化の失敗を検知してPython側に例外（Error）を投げる
    # ===========================================================
  

    # 極短時間(1e-5秒)だけテスト実行し、初期状態が破綻していないか確認
    prob_test = ODEProblem(nw, s0, (0.0, 1e-5))
    sol_test = solve(prob_test, Rodas5P())
    if string(sol_test.retcode) != "Success"
        error("Initialization failed: ODE unstable at t=0.")
    end

    if any(isnan, prob_test.u0)
        error("Initialization failed: NaN detected in initial conditions.")
    end

    prob = ODEProblem(nw, s0, (0.0, max_t))
    SIM_STATE["prob"] = prob
    SIM_STATE["u0"] = copy(prob.u0)
    SIM_STATE["t_current"] = 0.0
    SIM_STATE["t_cle"] = t_cle
    SIM_STATE["max_t"] = max_t

    return true
end

# =======================================================
# ステップ関数
# =======================================================
function run_step(dt, dummy_action, obs_requests_py)
    prob = SIM_STATE["prob"]
    u0 = SIM_STATE["u0"]
    t_start = SIM_STATE["t_current"]
    t_cle = SIM_STATE["t_cle"]
    max_t = SIM_STATE["max_t"]
    t_end = t_start + dt

    prob_step = ODEProblem(prob.f, u0, (t_start, t_end), prob.p)

    tstops = Float64[]
    if t_start < 1.0 <= t_end
        push!(tstops, 1.0)
    end
    if t_start < t_cle <= t_end
        push!(tstops, t_cle)
    end

    # ステップ終端と一致する場合も、投入・解消時刻を不連続点として通知する。
    sol = solve_fault_switched(prob_step, Rodas5P(), tstops=tstops, d_discontinuities=tstops)

    SIM_STATE["u0"] = copy(sol.u[end])
    SIM_STATE["t_current"] = t_end

    obs_requests = pyconvert(Dict, obs_requests_py)
    obs = Dict{String, Any}()
    
    for (b_idx_any, keys_any) in obs_requests
        b_idx = Int(b_idx_any)
        keys = [String(k) for k in keys_any]
        
        t_last = sol.t[end]
        for k in keys
            if k == "V_mag"
                obs["$(b_idx)_V_mag"] = sol(t_last; continuity=:right, idxs=VIndex(b_idx, :busbar₊u_mag))
            elseif k == "I"
                obs["$(b_idx)_I"] = sol(t_last; continuity=:right, idxs=VIndex(b_idx, :busbar₊i_mag))
            elseif b_idx in [1, 2, 3,7]
                inv_sym = Symbol("inv$(b_idx)")
                if k == "V_out"
                    obs["$(b_idx)_V_out"] = sol(t_last; continuity=:right, idxs=VIndex(b_idx, Symbol(inv_sym, "₊droop₊V_out₊u")))
                elseif k == "delta"
                    obs["$(b_idx)_delta"] = sol(t_last; continuity=:right, idxs=VIndex(b_idx, Symbol(inv_sym, "₊droop₊δ")))
                elseif k == "P"
                    obs["$(b_idx)_P"] = sol(t_last; continuity=:right, idxs=VIndex(b_idx, Symbol(inv_sym, "₊droop₊Pfilt")))
                elseif k == "Q"
                    obs["$(b_idx)_Q"] = sol(t_last; continuity=:right, idxs=VIndex(b_idx, Symbol(inv_sym, "₊droop₊Qfilt")))
                end
            else
                obs["$(b_idx)_$(k)"] = NaN
            end
        end
    end

    obs["time"] = t_end
    obs["retcode"] = string(sol.retcode)
    obs["max_t"] = max_t

    return obs
end
