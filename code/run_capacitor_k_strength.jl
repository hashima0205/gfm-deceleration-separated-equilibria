# One model build: fixed capacitor/load, stronger grid, equilibria/events/basins.
const K_SOURCE=joinpath(@__DIR__,"run_trajectory_voltage_probes.jl")
ksrc=first(split(read(K_SOURCE,String),"\nfunction mainp()"))
ksrc=replace(ksrc,"basin_recalculation/trajectory_voltage_probes\""=>"basin_recalculation/capacitor_k_strength\"")
include_string(Main,ksrc,K_SOURCE)
const KB=.5561806200995012
const KVALUES=[.8,.75,.7]
const KTAUS=[.1,1.2]
const KALPHA=[-5.,-3.,-1.,-.5,0.,.5,.9,1.,1.1,1.5,3.]
const KFREQ=[-.15,-.10,-.05,0.,.05]
function kparams(base,k,B,tau)
    p=tauall!(copy(base),tau);s=NWState(NW,copy(PREP),p,0.)
    s.p[VIndex(5,:PQLoad₊Q_load_p)]=-B*s.p[VIndex(5,:PQLoad₊vref2)]
    for (j,y) in [(2,3.43),(3,.68)]
        s.p[EIndex(Symbol("line",j),:line_sym₊R)]=k/y*cosd(75.)
        s.p[EIndex(Symbol("line",j),:line_sym₊X)]=k/y*sind(75.)
    end
    p
end
function kfollow(x,base,fromk,fromB,tok,toB,tag)
    distance=max(abs(tok-fromk),abs(toB-fromB));distance==0 && return copy(x)
    w=0.;dw=min(1.,.005/distance);rows=[];count=0
    while w<1-1e-12
        trial=min(1.,w+dw)
        pp(q)=kparams(base,fromk+q*(tok-fromk),fromB+q*(toB-fromB),.1)
        f1=similar(x);f2=similar(x);h=1e-5
        NW(f1,x,pp(w+h),0.);NW(f2,x,pp(w-h),0.)
        tangent=-(fd(NW,x,pp(w))\((f1-f2)/(2h)))
        y=x+(trial-w)*tangent;p=pp(trial);f=similar(x);er=Inf
        for correction in 1:25
            NW(f,y,p,0.);er=maximum(abs.(f));er<1e-9 && break
            all(isfinite,y) || break
            dy=fd(NW,y,p)\f;norm(dy,Inf)<100 || break
            y-=dy
        end
        NW(f,y,p,0.);er=maximum(abs.(f));accept=er<1e-8
        push!(rows,[w,trial,fromk+trial*(tok-fromk),fromB+trial*(toB-fromB),er,accept])
        tablep(tag*"_continuation.csv","from_w,trial_w,k,B,residual,accepted",rows)
        if accept
            x=y;w=trial;dw=min(.01/distance,dw*1.3);count+=1
            if count%25==0;println("CONTINUE ",tag," k=",fromk+w*(tok-fromk)," B=",fromB+w*(toB-fromB));flush(stdout);end
        else
            dw/=2;@assert dw*distance>1e-8 "Equilibrium continuation unresolved: $tag at k=$(fromk+trial*(tok-fromk)) B=$(fromB+trial*(toB-fromB)) residual=$er"
        end
    end
    x
end
function kphysical(tag,xp,xs,k,tau)
    pre=kparams(BEFOREP,k,KB,tau);post=kparams(POSTP,k,KB,tau)
    on=copy(pre);NWState(NW,copy(xp),on,0.).p[VIndex(4,:fault_load₊fault_enabled_p)]=1.
    for (name,p) in [("pre",pre),("fault",on),("post",post)];savevec(joinpath(OUTP,tag*"_"*name*"_parameters.csv"),p);end
    a=integratep(tag*"_pre",xp,pre,0.,1.;step=.01)
    f=integratep(tag*"_fault",a.u[end],on,1.,1.15;step=.0005)
    t=[0.];u=[endpointp(f.u[end],post)[1]];last=0.
    for stop in [30.,60.,150.]
        result=integratep(tag*"_post"*string(Int(stop)),u[end],post,last,stop)
        appendp!(t,u,result);last=stop
    end
    label,turns,ms=classify_probe(t,u,xs)
    savevec(joinpath(OUTP,tag*"_clear.csv"),u[1]);savevec(joinpath(OUTP,tag*"_end.csv"),u[end])
    trajectoryp(tag,t,u)
    tablep(tag*"_fault.csv","t_absolute_s,theta_infinite_bus_rad,df_mean_Hz,V5_pu",([tt,phasep(x),freqp(x),busvp(x,5)] for (tt,x) in zip(f.t,f.u)))
    println("PHYSICAL ",tag," ",label," turns=",turns);flush(stdout)
    label,turns,ms
end
function kbasin(xs,xu,k,tau)
    p=kparams(POSTP,k,KB,tau);rows=[];hist=[];initials=[];audits=[];machines=[]
    for (i,(alpha,df)) in enumerate(Iterators.product(KALPHA,KFREQ))
        tag=@sprintf("k080_tau%03d_b%03d",round(Int,100tau),i)
        x=copy(xs);x[DYN]=xs[DYN]+alpha*(xu[DYN]-xs[DYN]);x[PF].-=df*2pi ./ KP010[2:4]
        before=copy(x[DYN]);x,er=endpointp(x,p);@assert x[DYN]==before
        push!(initials,vcat([tag,alpha,df],x));push!(audits,[tag,er,maximum(abs.(x[DYN]-before))])
        t=[0.];u=[copy(x)];last=0.;label="unresolved";turns=0.;ms=[]
        for stop in [30.,60.,150.]
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
        push!(rows,vcat([tag,k,tau,alpha,df,phasep(x),freqp(x),label,last,turns,residual],branches))
        for m in ms;push!(machines,vcat([tag],collect(m)));end
        savevec(joinpath(OUTP,tag*"_end.csv"),u[end]);trajectoryp(tag,t,u)
        suffix=@sprintf("tau%03d",round(Int,100tau))
        tablep("basin_"*suffix*".csv","case,k,taup,alpha,df_offset_Hz,theta_infinite_bus_rad,df_Hz,classification,horizon_s,remote_turns,max_algebraic_residual,branch_b1,branch_b2,branch_b3,branch_b7",rows)
        tablep("basin_"*suffix*"_history.csv","case,horizon_s,classification,solve_s",hist)
        tablep("basin_"*suffix*"_initial.csv",join(vcat(["case","alpha","df_offset_Hz"],["x$j" for j in 1:66]),','),initials)
        tablep("basin_"*suffix*"_audit.csv","case,algebraic_residual,differential_change",audits)
        tablep("basin_"*suffix*"_machines.csv","case,bus,classification,turns,tail_max_abs_df_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",machines)
        println("BASIN ",suffix," ",i,"/",length(KALPHA)*length(KFREQ)," ",label," t=",last);flush(stdout)
    end
end
function kmain()
    mkpath(OUTP)
    files=[@__FILE__,K_SOURCE,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(files,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    signature=bytes2hex(sha256(string(VERSION,PREP,STP,POSTP,BEFOREP,KB,KVALUES,KTAUS,KALPHA,KFREQ)*join(bytes2hex(sha256(read(f))) for f in files)))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    tick=time();roots=[];events=[];metrics=[]
    xp0=vec(readdlm(joinpath(PACK,"basin_recalculation/capacitor_probe_predictor/v095_tau010_pre_root.csv"),',',Float64))
    xp=kfollow(xp0,BEFOREP,1.0299,KB,.8,KB,"pre_to_k080")
    xs0=kfollow(STP,POSTP,1.0299,0.,.8,0.,"stable_grid")
    xu0=kfollow(vec(readdlm(joinpath(DATA,"seeds/nf02_1_post1.csv"),',',Float64)),POSTP,1.0299,0.,.8,0.,"saddle_grid")
    xs=kfollow(xs0,POSTP,.8,0.,.8,KB,"stable_cap")
    xu=kfollow(xu0,POSTP,.8,0.,.8,KB,"saddle_cap")
    lastk=.8
    for k in KVALUES
        xp=kfollow(xp,BEFOREP,lastk,KB,k,KB,@sprintf("pre_k%03d",100k))
        xs=kfollow(xs,POSTP,lastk,KB,k,KB,@sprintf("stable_k%03d",100k))
        xu=kfollow(xu,POSTP,lastk,KB,k,KB,@sprintf("saddle_k%03d",100k));lastk=k
        for (name,x) in [("pre",xp),("stable",xs),("saddle",xu)];savevec(joinpath(OUTP,@sprintf("k%03d_%s.csv",100k,name)),x);end
        for tau in KTAUS
            for (name,x,base) in [("pre",xp,BEFOREP),("stable",xs,POSTP),("saddle",xu,POSTP)]
                p=kparams(base,k,KB,tau);er,ev,n=eqcheck(x,p);@assert er<1e-8
                @assert name=="saddle" ? n>=1 : ev<0
                push!(roots,[k,tau,name,er,ev,n,busvp(x,5),phasep(x)])
            end
            tablep("equilibria.csv","k,taup,root,residual,max_real_eigen,unstable_modes,V5_pu,mean_phase_rad",roots)
            tag=@sprintf("k%03d_tau%03d",100k,round(Int,100tau))
            label,turns,ms=kphysical(tag,xp,xs,k,tau)
            G=2.2682301152333055
            push!(events,[tag,k,tau,KB,busvp(xp,5),busvp(xs,5),G*busvp(xp,5)^2,G*busvp(xs,5)^2,label,turns])
            for m in ms;push!(metrics,vcat([tag],collect(m)));end
            tablep("events.csv","case,k,taup,B_pu,pre_V5_pu,post_eq_V5_pu,pre_load_P_pu,post_eq_load_P_pu,classification,remote_turns",events)
            tablep("event_machines.csv","case,bus,classification,turns,tail_max_abs_df_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",metrics)
        end
    end
    # Both roots lie on the same 52-dimensional chord section. Include the
    # phase range near the actual clearing projection, not only root vicinity.
    xs=vec(readdlm(joinpath(OUTP,"k080_stable.csv"),',',Float64));xu=vec(readdlm(joinpath(OUTP,"k080_saddle.csv"),',',Float64))
    for tau in KTAUS;kbasin(xs,xu,.8,tau);end
    println("DONE stronger-grid capacitor study; batch_wall_s=",time()-tick);flush(stdout)
end
kmain()
