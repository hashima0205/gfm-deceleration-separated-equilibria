# B=1, fixed k near the previously tracked fold; identical chord grid for two taus.
const B100_BASE=joinpath(@__DIR__,"run_capacitor_k_strength.jl")
b100source=first(split(read(B100_BASE,String),"\nfunction kmain()"))
@assert occursin("const KB=.5561806200995012",b100source)
b100source=replace(b100source,"const KB=.5561806200995012"=>"const KB=1.0",
    "basin_recalculation/capacitor_k_strength\\\""=>"basin_recalculation/B100_tau_comparison\\\"")
include_string(Main,b100source,B100_BASE)
@assert KB==1. && endswith(OUTP,"B100_tau_comparison")
const BK=.6812730259441179
const BTAUS=[.5,1.]
const BGLOBAL_THETA=[-.4,-.15,0.,.4,.8,1.2,1.6]
const BGLOBAL_FREQ=[-.15,-.10,-.05,0.,.05]
const BLOCAL_ALPHA=[-.5,0.,.5,.9,1.,1.1,1.5]
const BLOCAL_FREQ=[-.01,-.005,-.0025,0.,.0025,.005,.01]
const BHORIZONS=[30.,60.,150.,300.]
function bgrid(xs,xu)
    gap=phasep(xu)-phasep(xs);@assert abs(gap)>1e-4
    points=Tuple{Float64,Float64,String}[]
    for df in BGLOBAL_FREQ,theta in BGLOBAL_THETA;push!(points,(theta/gap,df,"global"));end
    for df in BLOCAL_FREQ,alpha in BLOCAL_ALPHA
        if !any(abs(a-alpha)<1e-12 && f==df for (a,f,g) in points);push!(points,(alpha,df,"local"));end
    end
    @assert length(points)==83
    points
end
function bbasin(xs,xu,tau,points)
    p=kparams(POSTP,BK,KB,tau);rows=[];hist=[];initials=[];audits=[];machines=[]
    for (i,(alpha,df,grid)) in enumerate(points)
        tag=@sprintf("B100_tau%03d_b%03d",round(Int,100tau),i)
        x=copy(xs);x[DYN]=xs[DYN]+alpha*(xu[DYN]-xs[DYN]);x[PF].-=df*2pi ./ KP010[2:4]
        before=copy(x[DYN]);x,er=endpointp(x,p);@assert x[DYN]==before
        push!(initials,vcat([tag,grid,alpha,df],x));push!(audits,[tag,er,maximum(abs.(x[DYN]-before))])
        t=[0.];u=[copy(x)];last=0.;label="unresolved";turns=0.;ms=[]
        for stop in BHORIZONS
            result=integratep(tag*"_h"*string(Int(stop)),u[end],p,last,stop)
            appendp!(t,u,result);last=stop;label,turns,ms=classify_probe(t,u,xs)
            push!(hist,[tag,stop,label,result.solve_s])
            label in ["unresolved","other"] || break
        end
        branches=round.(Int,(u[end][DI010]-xs[DI010])/(2pi))
        if alpha==1. && df==0.;label="saddle_initial";end
        residual=0.;ff=similar(x)
        for (tt,xx) in zip(t,u);NW(ff,xx,p,tt);residual=max(residual,maximum(abs.(ff[ALG])));end
        @assert residual<1e-7
        push!(rows,vcat([tag,grid,BK,KB,tau,alpha,df,phasep(x),freqp(x),label,last,turns,residual],branches))
        for m in ms;push!(machines,vcat([tag],collect(m)));end
        savevec(joinpath(OUTP,tag*"_end.csv"),u[end]);trajectoryp(tag,t,u)
        suffix=@sprintf("tau%03d",round(Int,100tau))
        tablep("basin_"*suffix*".csv","case,grid,k,B_pu,taup,alpha,df_offset_Hz,theta_infinite_bus_rad,df_Hz,classification,horizon_s,remote_turns,max_algebraic_residual,branch_b1,branch_b2,branch_b3,branch_b7",rows)
        tablep("basin_"*suffix*"_history.csv","case,horizon_s,classification,solve_s",hist)
        tablep("basin_"*suffix*"_initial.csv",join(vcat(["case","grid","alpha","df_offset_Hz"],["x$j" for j in 1:66]),','),initials)
        tablep("basin_"*suffix*"_audit.csv","case,algebraic_residual,differential_change",audits)
        tablep("basin_"*suffix*"_machines.csv","case,bus,classification,turns,tail_max_abs_df_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",machines)
        println("BASIN ",suffix," ",i,"/",length(points)," ",grid," ",label," t=",last);flush(stdout)
    end
end
function bmain()
    mkpath(OUTP)
    files=[@__FILE__,B100_BASE,K_SOURCE,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(files,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    seed=joinpath(PACK,"basin_recalculation/fast_tau_search/runs")
    xp=vec(readdlm(joinpath(seed,"fast_12_pre_root.csv"),',',Float64))
    xs=vec(readdlm(joinpath(seed,"fast_12_post_root.csv"),',',Float64))
    xu0=vec(readdlm(joinpath(PACK,"basin_recalculation/capacitor_k_strength/k080_saddle.csv"),',',Float64))
    signature=bytes2hex(sha256(string(VERSION,xp,xs,xu0,BK,KB,BTAUS,BGLOBAL_THETA,BGLOBAL_FREQ,BLOCAL_ALPHA,BLOCAL_FREQ,BHORIZONS)*join(bytes2hex(sha256(read(f))) for f in files)))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    tick=time()
    # Route below both folds before approaching the narrow stable/saddle pair.
    xu=kfollow(xu0,POSTP,.8,.5561806200995012,.6,.5561806200995012,"saddle_strengthen")
    xu=kfollow(xu,POSTP,.6,.5561806200995012,.6,KB,"saddle_compensate")
    xu=kfollow(xu,POSTP,.6,KB,BK,KB,"saddle_approach")
    for (name,x) in [("pre",xp),("stable",xs),("saddle",xu)];savevec(joinpath(OUTP,name*".csv"),x);end
    points=bgrid(xs,xu)
    tablep("grid.csv","point,grid,alpha,df_offset_Hz",[[i,g,a,f] for (i,(a,f,g)) in enumerate(points)])
    rows=[];events=[];metrics=[]
    for tau in BTAUS
        for (name,x,base) in [("pre",xp,BEFOREP),("stable",xs,POSTP),("saddle",xu,POSTP)]
            p=kparams(base,BK,KB,tau);er,ev,n=eqcheck(x,p);@assert er<1e-8
            @assert name=="saddle" ? n>=1 : ev<0
            push!(rows,[tau,name,er,ev,n,busvp(x,5),phasep(x)])
        end
        tablep("equilibria.csv","taup,root,residual,max_real_eigen,unstable_modes,V5_pu,mean_phase_rad",rows)
        tag=@sprintf("B100_tau%03d",round(Int,100tau))
        label,turns,ms=kphysical(tag,xp,xs,BK,tau)
        push!(events,[tag,BK,KB,tau,busvp(xp,5),busvp(xs,5),2.2682301152333055*busvp(xp,5)^2,KB*busvp(xp,5)^2,label,turns])
        for m in ms;push!(metrics,vcat([tag],collect(m)));end
        tablep("events.csv","case,k,B_pu,taup,pre_V5_pu,post_eq_V5_pu,pre_load_P_pu,pre_cap_Q_pu,classification,remote_turns",events)
        tablep("event_machines.csv","case,bus,classification,turns,tail_max_abs_df_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",metrics)
    end
    for tau in BTAUS;bbasin(xs,xu,tau,points);end
    println("DONE B100 tau comparison; batch_wall_s=",time()-tick);flush(stdout)
end
bmain()
