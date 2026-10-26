using PowerDynamics
using PowerDynamics: Library
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using PowerDynamics.Library.ComposableInverter: CC1, VC, LCLFilter, DroopOuter, _ri_to_dq, _dq_to_ri
using ModelingToolkit: ModelingToolkit, @named, t_nounits as t, D_nounits as Dt, @component
using ModelingToolkit: @mtkmodel, @variables, @parameters, @unpack, Num, System, Equation, connect, setmetadata

@component function DroopInverter_addlimit(; name=:droop_inv, filter_type, defaults...)
    @named droop = DroopOuter(; pq_input=false)
    @named vsrc = VoltageSource_addlimit(; Vset_input=true, filter_type)
    @named terminal = Terminal()

    eqs = [
        connect(vsrc.Vset_in, droop.V_out)
        connect(vsrc.δset_in, droop.δ_out)
        connect(vsrc.terminal, terminal)
        droop.P_meas_in.u ~ vsrc.filter.V_C_r * vsrc.filter.i_g_r + vsrc.filter.V_C_i * vsrc.filter.i_g_i
        droop.Q_meas_in.u ~ vsrc.filter.V_C_i * vsrc.filter.i_g_r - vsrc.filter.V_C_r * vsrc.filter.i_g_i
    ]

    sys = System(eqs, t; name, systems=[droop, vsrc, terminal])
    set_mtk_defaults!(sys, defaults)
    return sys
end

@component function VoltageSource_addlimit(; name, Vset_input=false, filter_type=:LC, defaults...)
    @parameters begin
        Rf  = 0.0139; Lf = 0.0196; C = 0.1086; ω0  = 2π*60; I_max = 1.2
        CC1_KP = 0.9817; CC1_KI = 0.0018; CC1_F  = 0; CC1_Fcoupl = 1
        VC_KP = 1.4476; VC_KI = 0.0273; VC_F = 0; VC_Fcoupl = 1
        R_virt = 0; X_virt = 0
    end
    if filter_type == :LCL
        @parameters begin
            Rg  = 0.01; Lg = 0.007
        end
    end
    if !Vset_input
        @parameters begin
            Vset, [guess=1]; δset, [guess=0]
        end
    end

    if filter_type == :LC
        @named filter = LCFilter(; ω0, Rf, Lf, C)
    else  
        @named filter = LCLFilter(; ω0, Rf, Lf, C, Rg, Lg)
    end

    @named cc1 = CC1(; Lf=Lf, F=CC1_F, Fcoupl=CC1_Fcoupl, KP=CC1_KP, KI=CC1_KI)
    systems = @named begin
        terminal = Terminal()
        vc = VC_addlimit(; C=C, F=VC_F, Fcoupl=VC_Fcoupl, KP=VC_KP, KI=VC_KI, I_max=I_max)
    end
    push!(systems, filter)
    push!(systems, cc1)
    if Vset_input
        @named Vset_in = RealInput()
        @named δset_in = RealInput()
        append!(systems, [Vset_in, δset_in])
    end
    _Vset = Vset_input ? Vset_in.u : Vset
    _δ = Vset_input ? δset_in.u : δset

    eqs = [
        connect(terminal, filter.terminal)
        [cc1.i_f_d, cc1.i_f_q] .~ _ri_to_dq(filter.i_f_r, filter.i_f_i, _δ)
        [cc1.V_C_d, cc1.V_C_q] .~ _ri_to_dq(filter.V_C_r, filter.V_C_i, _δ)
        [vc.V_C_d, vc.V_C_q] .~ _ri_to_dq(filter.V_C_r, filter.V_C_i, _δ)
        [vc.i_g_d, vc.i_g_q] .~ _ri_to_dq(filter.i_g_r, filter.i_g_i, _δ)
        [filter.V_I_r, filter.V_I_i] .~ _dq_to_ri(cc1.V_I_d, cc1.V_I_q, _δ)
        cc1.i_f_ref_d ~ vc.i_f_ref_d
        cc1.i_f_ref_q ~ vc.i_f_ref_q
        vc.V_C_ref_d ~ _Vset - (vc.i_g_d * R_virt - vc.i_g_q * X_virt)
        vc.V_C_ref_q ~       - (vc.i_g_q * R_virt + vc.i_g_d * X_virt)
    ]
    sys = System(eqs, t; name, systems)
    set_mtk_defaults!(sys, defaults)
    return sys
end

@component function VC_addlimit(; name, C, F, Fcoupl=1, KP, KI, I_max, defaults...)
    vars = @variables begin
        γ_d(t), [guess=0]
        γ_q(t), [guess=0]
        i_f_ref_d(t); i_f_ref_q(t)
        V_C_d(t); V_C_q(t)
        V_C_ref_d(t); V_C_ref_q(t)
        i_g_d(t); i_g_q(t)
    end
    
    unlim_d = -Fcoupl*C*V_C_q + KP*(V_C_ref_d - V_C_d) + KI*γ_d + F*i_g_d
    unlim_q =  Fcoupl*C*V_C_d + KP*(V_C_ref_q - V_C_q) + KI*γ_q + F*i_g_q
    
    I_mag = sqrt(unlim_d^2 + unlim_q^2 + 1e-6)
    
    eqs = [
        Dt(γ_d) ~ V_C_ref_d - V_C_d
        Dt(γ_q) ~ V_C_ref_q - V_C_q
        
        i_f_ref_d ~ ifelse(I_mag > I_max, unlim_d * (I_max / I_mag), unlim_d)
        i_f_ref_q ~ ifelse(I_mag > I_max, unlim_q * (I_max / I_mag), unlim_q)
    ]
    sys = System(eqs, t; name)
    set_mtk_defaults!(sys, defaults)
    return sys
end

# compornent.jl の末尾を以下に差し替えてください

# 1. 内部で利用するための定電力負荷コンポーネントを定義
@component function ConstantPQLoad(; name, P_load=0.0, Q_load=0.0, defaults...)
    @named terminal = Terminal()
    @parameters begin
        P_load = P_load
        Q_load = Q_load
    end
    vars = @variables begin
        V_mag(t)
    end
    eqs = [
        V_mag ~ sqrt(terminal.u_r^2 + terminal.u_i^2 + 1e-6)
        # 電力を消費するため、系統へ注入する電流はマイナス
        terminal.i_r ~ -(P_load * terminal.u_r + Q_load * terminal.u_i) / V_mag^2
        terminal.i_i ~ -(P_load * terminal.u_i - Q_load * terminal.u_r) / V_mag^2
    ]
    sys = System(eqs, t; name, systems=[terminal])
    set_mtk_defaults!(sys, defaults)
    return sys
end

# 2. GFMインバータと負荷の複合コンポーネント
@component function DroopInverter_WithLoad(; name=:droop_inv_load, filter_type, P_load=0.0, Q_load=0.0, defaults...)
    @named droop = DroopOuter(; pq_input=false)
    @named vsrc = VoltageSource_addlimit(; Vset_input=true, filter_type)
    @named int_load = ConstantPQLoad(; P_load=P_load, Q_load=Q_load) # 内部に負荷を組み込む
    @named terminal = Terminal()

    eqs = [
        connect(vsrc.Vset_in, droop.V_out)
        connect(vsrc.δset_in, droop.δ_out)
        
        # MTKのconnectマクロを使って、インバータ、負荷、外部ターミナルを1点に並列接続
        # (これで重複エラーを回避しつつ、自動的に電流の足し合わせが行われます)
        connect(vsrc.terminal, int_load.terminal, terminal)
        
        # ドループ制御用の測定値 (インバータ自身の出力電流のみを使用)
        droop.P_meas_in.u ~ vsrc.filter.V_C_r * vsrc.filter.i_g_r + vsrc.filter.V_C_i * vsrc.filter.i_g_i
        droop.Q_meas_in.u ~ vsrc.filter.V_C_i * vsrc.filter.i_g_r - vsrc.filter.V_C_r * vsrc.filter.i_g_i
    ]

    sys = System(eqs, t; name, systems=[droop, vsrc, int_load, terminal])
    set_mtk_defaults!(sys, defaults)
    return sys
end