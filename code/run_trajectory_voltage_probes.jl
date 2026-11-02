# One model build for trajectory-centered sections and a controlled load-level probe.
include(joinpath(@__DIR__,"model.jl"))
using SHA,Serialization,Printf
const OUTP=joinpath(PACK,"basin_recalculation/trajectory_voltage_probes")
const STP=vec(readdlm(joinpath(DATA,"seeds/tau010_T015_h150_stable.csv"),',',Float64))
const PREP=vec(readdlm(joinpath(DATA,"seeds/nf02_1_pre0.csv"),',',Float64))
const POSTP=vec(readdlm(joinpath(CHECK,"tau010_T015_h150/post150_parameters.csv"),',',Float64))
const BEFOREP=vec(readdlm(joinpath(CHECK,"tau010_T015_h150/pre_parameters.csv"),',',Float64))
const CLEARP=vec(readdlm(joinpath(CHECK,"tau010_T015_h150/post030_initial.csv"),',',Float64))
const IX=NWState(NW,Float64.(1:length(STP)),copy(POSTP),0.)
const DI010=[Int(IX[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]) for b in ALLG]
const PI010=[Int(IX[VIndex(b,Symbol("inv$(b)₊droop₊Pfilt"))]) for b in ALLG]
const KP010=[NWState(NW,copy(STP),copy(POSTP),0.).p[VIndex(b,Symbol("inv$(b)₊droop₊Kp"))] for b in ALLG]
const PS010=[.5,.4,.4,.4]
const HELDP=setdiff(DYN,vcat(DELTA,PF))
const ANCHOR_TIMES=[0.,1.,3.,10.]
const ANGLE_OFFSETS=[-.15,-.075,0.,.075,.15]
const FREQ_OFFSETS=[-.05,-.025,0.,.025,.05]
const LOAD_FACTORS=[1.,.97,.90,.80]
const VOLTAGE_TAUS=[.1,1.2]
include(joinpath(@__DIR__,"adaptive_probe_metrics.jl"))

function tablep(file,header,rows)
    path=joinpath(OUTP,file);mkpath(dirname(path))
    open(path*".tmp","w") do io
        println(io,header)
        for row in rows;println(io,join(row,','));end
    end
    mv(path*".tmp",path;force=true)
end
phasep(x)=sum(x[DELTA])/3
freqp(x)=sum(KP010[2:4].*(PS010[2:4]-x[PI010[2:4]]))/(6pi)
function busvp(x,b)
    b==6 && return 1.
    hypot(x[Int(IX[VIndex(b,:busbar₊u_r)])],x[Int(IX[VIndex(b,:busbar₊u_i)])])
end
function endpointp(x,p)
    y=copy(x);before=copy(y[DYN]);er=algfix!(y,p)
    @assert er<1e-8 && y[DYN]==before
    y,er
end
function integratep(tag,x,p,a,b;step=.05)
    path=joinpath(OUTP,tag*".bin");mkpath(dirname(path))
    if isfile(path);return deserialize(path);end
    y,er=endpointp(x,p)
    times=sort(unique(vcat(collect(a:step:b),collect(a:.0001:min(a+.02,b)),[a,b])))
    tick=time()
    sol=solve(ODEProblem(NW,y,(a,b),copy(p)),Rodas5P();reltol=1e-9,abstol=1e-11,dtmax=.005,
              saveat=times,save_everystep=false,dense=false,maxiters=10^8)
    @assert string(sol.retcode)=="Success" "Solver failed in $tag"
    result=(t=collect(sol.t),u=[copy(v) for v in sol.u],solve_s=time()-tick)
    serialize(path*".tmp",result);mv(path*".tmp",path;force=true)
    result
end
function appendp!(t,u,result)
    for j in eachindex(result.t)
        result.t[j]>t[end] || continue
        push!(t,result.t[j]);push!(u,result.u[j])
    end
end
function trajectoryp(tag,t,u)
    tablep(tag*".csv","t_post_s,theta_infinite_bus_rad,df_mean_Hz,min_bus_voltage_pu,min_GFM_voltage_pu",
           ([tt,phasep(x),freqp(x),minimum(busvp(x,b) for b in 1:8),minimum(busvp(x,b) for b in ALLG)] for (tt,x) in zip(t,u)))
end

function neighborhoodp()
    # Reintegrate once to retain full states at actual trajectory anchors.
    ref=integratep("reference_post160",CLEARP,POSTP,0.,160.;step=.02)
    at(t)=begin
        j=argmin(abs.(ref.t.-t));@assert abs(ref.t[j]-t)<1e-8;copy(ref.u[j])
    end
    expected=vec(readdlm(joinpath(CHECK,"tau010_T015_h150/post150_end.csv"),',',Float64))
    reference_error=maximum(abs.(at(150.)-expected))
    @assert reference_error<2e-7
    tablep("reference_validation.csv","point,max_full_state_difference",[["post150",reference_error]])
    trajectoryp("reference_post160",ref.t,ref.u)
    rows=[];hist=[];audits=[];initials=[];centers=[];machine_rows=[];anchors=[]
    for (ai,anchor_time) in enumerate(ANCHOR_TIMES)
        base=at(anchor_time);theta0=phasep(base);frequency0=freqp(base)
        push!(anchors,vcat([ai,anchor_time,theta0,frequency0],base))
        for (point_index,(da,df)) in enumerate(Iterators.product(ANGLE_OFFSETS,FREQ_OFFSETS))
            tag=@sprintf("near_a%d_p%02d",ai,point_index)
            x=copy(base);x[DELTA].+=da;x[PF].-=df*2pi ./ KP010[2:4]
            x,er=endpointp(x,POSTP)
            @assert x[HELDP]==base[HELDP] && abs(phasep(x)-theta0-da)<1e-9 && abs(freqp(x)-frequency0-df)<1e-9
            push!(initials,vcat([tag,ai,anchor_time,da,df],x))
            push!(audits,[tag,er,phasep(x)-theta0-da,freqp(x)-frequency0-df,0.])
            t=[0.];u=[copy(x)];last=0.;label="unresolved";turns=0.;machines=[];solve_s=0.
            for stop in [30.,60.,150.]
                result=integratep(tag*"_h"*string(Int(stop)),x,POSTP,last,stop)
                appendp!(t,u,result);x=copy(u[end]);last=t[end];solve_s+=result.solve_s
                label,turns,machines=classify_probe(t,u,STP)
                push!(hist,[tag,ai,anchor_time,stop,label,result.solve_s])
                label in ["unresolved","other"] || break
            end
            f=zeros(length(x));NW(f,x,POSTP,0.);res=maximum(abs.(f[ALG]));@assert res<1e-7
            branches=round.(Int,(x[DI010]-STP[DI010])/(2pi))
            if da==df==0
                error=maximum(abs.(x-at(anchor_time+last)))
                @assert error<2e-7
                push!(centers,[tag,anchor_time,last,error])
            end
            push!(rows,vcat([tag,ai,anchor_time,da,df,theta0+da,frequency0+df,label,last,turns,res,solve_s],branches))
            for m in machines;push!(machine_rows,vcat([tag,last],collect(m)));end
            savevec(joinpath(OUTP,tag*"_end.csv"),x)
            trajectoryp(tag,t,u)
            tablep("neighborhood_samples.csv","case,anchor_id,anchor_post_s,angle_offset_rad,df_offset_Hz,theta_infinite_bus_rad,df_Hz,classification,horizon_s,remote_internal_turns,max_algebraic_residual,solve_s,branch_b1,branch_b2,branch_b3,branch_b7",rows)
            tablep("neighborhood_history.csv","case,anchor_id,anchor_post_s,horizon_s,classification,solve_s",hist)
            tablep("neighborhood_initialization.csv","case,algebraic_residual,theta_error_rad,df_error_Hz,held_change",audits)
            tablep("neighborhood_initial_states.csv",join(vcat(["case","anchor_id","anchor_post_s","angle_offset_rad","df_offset_Hz"],["x$j" for j in 1:66]),','),initials)
            tablep("neighborhood_machines.csv","case,horizon_s,bus,classification,internal_turns,tail_df_max_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",machine_rows)
            tablep("center_validation.csv","case,anchor_post_s,additional_horizon_s,max_full_state_difference",centers)
            println("NEAR ",length(rows),"/100 anchor=",anchor_time," ",label," horizon=",last);flush(stdout)
        end
    end
    tablep("anchors.csv",join(vcat(["anchor_id","anchor_post_s","theta_infinite_bus_rad","df_Hz"],["x$j" for j in 1:66]),','),anchors)
end

function voltagep()
    rows=[];roots=[];parameter_rows=[]
    previous_pre=copy(PREP);previous_post=copy(STP);previous_factor=1.
    for (li,factor) in enumerate(LOAD_FACTORS)
        # Small continuation steps keep the root search on the connected branch.
        for r in range(previous_factor,factor;length=max(2,ceil(Int,abs(factor-previous_factor)/.01)+1))
            pp=copy(BEFOREP);pq=copy(POSTP)
            NWState(NW,copy(PREP),pp,0.).p[VIndex(5,:PQLoad₊P_load_p)]*=r
            NWState(NW,copy(STP),pq,0.).p[VIndex(5,:PQLoad₊P_load_p)]*=r
            previous_pre,ep=stationary(NW,previous_pre,pp)
            previous_post,eq=stationary(NW,previous_post,pq)
            @assert ep<1e-8 && eq<1e-8 "Operating-point continuation failed at factor=$r"
        end
        previous_factor=factor
        for tau in VOLTAGE_TAUS
            tag=@sprintf("load%03d_tau%03d",round(Int,100factor),round(Int,100tau))
            pre=tauall!(copy(BEFOREP),tau);post=tauall!(copy(POSTP),tau)
            NWState(NW,copy(PREP),pre,0.).p[VIndex(5,:PQLoad₊P_load_p)]*=factor
            NWState(NW,copy(STP),post,0.).p[VIndex(5,:PQLoad₊P_load_p)]*=factor
            on=copy(pre);NWState(NW,copy(PREP),on,0.).p[VIndex(4,:fault_load₊fault_enabled_p)]=1.
            for (stage,target,p) in [("pre",previous_pre,pre),("post",previous_post,post)]
                er,ev,n=eqcheck(target,p);@assert er<1e-8 && ev<0
                push!(roots,[tag,stage,er,ev,n,busvp(target,5),minimum(busvp(target,b) for b in ALLG),phasep(target)])
            end
            for (stage,p) in [("pre",pre),("fault",on),("post",post)],(j,b) in enumerate(ALLG)
                state=NWState(NW,copy(previous_pre),p,0.)
                tp=state.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))];tq=state.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_q"))]
                ps=state.p[VIndex(b,Symbol("inv$(b)₊droop₊Pset"))]
                @assert tp==tau && tq==.02 && ps==PS010[j]
                push!(parameter_rows,[tag,stage,b,tp,tq,ps])
            end
            a=integratep(tag*"_pre",previous_pre,pre,0.,1.;step=.01)
            b=integratep(tag*"_fault",a.u[end],on,1.,1.15;step=.0005)
            x=copy(b.u[end]);t=[0.];u=[endpointp(x,post)[1]];last=0.;solve_s=a.solve_s+b.solve_s
            for stop in [30.,60.,150.]
                result=integratep(tag*"_post"*string(Int(stop)),u[end],post,last,stop)
                appendp!(t,u,result);last=stop;solve_s+=result.solve_s
            end
            label,turns,machines=classify_probe(t,u,previous_post)
            minbus=minimum(busvp(xx,bb) for xx in u for bb in 1:8)
            mingfm=minimum(busvp(xx,bb) for xx in u for bb in ALLG)
            push!(rows,[tag,factor,tau,busvp(previous_pre,5),minimum(busvp(previous_pre,bb) for bb in ALLG),
                        busvp(previous_post,5),minimum(busvp(previous_post,bb) for bb in ALLG),label,turns,minbus,mingfm,solve_s])
            savevec(joinpath(OUTP,tag*"_pre_root.csv"),previous_pre)
            savevec(joinpath(OUTP,tag*"_post_root.csv"),previous_post)
            savevec(joinpath(OUTP,tag*"_clear.csv"),u[1])
            savevec(joinpath(OUTP,tag*"_post_end.csv"),u[end])
            trajectoryp(tag,t,u)
            tablep("voltage_cases.csv","case,load_factor,taup,pre_V5_pu,pre_min_GFM_V_pu,post_eq_V5_pu,post_eq_min_GFM_V_pu,classification,remote_internal_turns,post_min_bus_voltage_pu,post_min_GFM_voltage_pu,solve_s",rows)
            tablep("voltage_equilibria.csv","case,stage,residual,max_real_eigen,unstable_modes,V5_pu,min_GFM_V_pu,mean_internal_phase_rad",roots)
            tablep("voltage_parameters.csv","case,stage,bus,taup,tauq,Pset",parameter_rows)
            println("VOLTAGE ",tag," preV5=",busvp(previous_pre,5)," postV5=",busvp(previous_post,5)," ",label," turns=",turns);flush(stdout)
        end
    end
end

function mainp()
    mkpath(OUTP)
    files=[@__FILE__,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(files,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    signature=bytes2hex(sha256(string(VERSION,STP,PREP,POSTP,BEFOREP,CLEARP,ANCHOR_TIMES,ANGLE_OFFSETS,FREQ_OFFSETS,LOAD_FACTORS,VOLTAGE_TAUS)*join([bytes2hex(sha256(read(f))) for f in files])))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature "Inputs/code changed: use a new output directory"
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    @assert length(DYN)==52 && length(HELDP)==46 && length(ALG)==14
    @assert eqcheck(STP,POSTP)[1]<1e-8
    tick=time();neighborhoodp();voltagep()
    println("DONE neighborhood + voltage probes; batch_wall_s=",time()-tick);flush(stdout)
end
mainp()
