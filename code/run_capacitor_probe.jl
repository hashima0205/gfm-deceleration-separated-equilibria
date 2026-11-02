# Constant shunt susceptance at Bus5, present before/during/after the fault.
const CAP_SOURCE=joinpath(@__DIR__,"run_trajectory_voltage_probes.jl")
definitions=first(split(read(CAP_SOURCE,String),"\nfunction mainp()"))
definitions=replace(definitions,"basin_recalculation/trajectory_voltage_probes\""=>"basin_recalculation/capacitor_probe_predictor\"")
include_string(Main,definitions,CAP_SOURCE)
const TARGETS_CAP=[busvp(PREP,5),.93,.95,.98]
function capparams(base,B,tau)
    p=tauall!(copy(base),tau);s=NWState(NW,copy(PREP),p,0.)
    @assert s.p[VIndex(5,:PQLoad₊z_fraction)]==1.
    @assert s.p[VIndex(5,:PQLoad₊Q_load_p)]==0.
    s.p[VIndex(5,:PQLoad₊Q_load_p)]=-B*s.p[VIndex(5,:PQLoad₊vref2)]
    p
end
function capfollow(x,from,to,base,tag)
    r=from;step=min(.001,max(to-from,0.));attempts=[]
    while r<to-1e-12
        trial=min(to,r+step);p=capparams(base,trial,.1)
        # Tangent predictor avoids a poorly scaled residual line search
        # forcing microunit parameter increments on an otherwise smooth branch.
        p0=capparams(base,r,.1);h=1e-5
        fp=similar(x);fm=similar(x)
        NW(fp,x,capparams(base,r+h,.1),0.)
        NW(fm,x,capparams(base,r-h,.1),0.)
        tangent=-(fd(NW,x,p0)\((fp-fm)/(2h)))
        y=x+(trial-r)*tangent;f=similar(x);er=Inf
        for correction in 1:20
            NW(f,y,p,0.);er=maximum(abs.(f))
            er<1e-8 && break
            all(isfinite,y) || break
            dy=fd(NW,y,p)\f
            norm(dy,Inf)<100 || break
            y-=dy
        end
        NW(f,y,p,0.);er=maximum(abs.(f))
        ev=NaN
        if er<1e-8;_,ev,_=eqcheck(y,p);end
        accepted=er<1e-8 && ev<0
        push!(attempts,[r,trial,step,er,ev,accepted])
        tablep(tag*"_continuation.csv","last_B_pu,trial_B_pu,step,residual,max_real_eigen,accepted",attempts)
        if accepted;x=y;r=trial;step=min(.01,step*1.4)
        else
            step/=2
            if step<1e-7
                println("ROOT_UNRESOLVED ",tag," B=",trial," residual=",er);flush(stdout)
                return nothing
            end
        end
    end
    x
end
function captarget(x,B,target)
    for iteration in 1:1000
        V=busvp(x,5);abs(V-target)<1e-6 && return x,B
        # Use a conservative voltage-sensitivity estimate, then update from
        # accepted roots; cap increments and retry root failures locally.
        step=min(.005,max(1e-6,(target-V)/2))
        for retry in 1:18
            y,er=stationary(NW,x,capparams(BEFOREP,B+step,.1))
            if er<1e-8 && busvp(y,5)>V && busvp(y,5)<=target+1e-6
                x=y;B+=step;break
            end
            step/=2
            @assert retry<18 "Pre-fault voltage targeting unresolved"
        end
        @assert B<3 "Target not reached within B=3 pu"
    end
    error("Pre-fault voltage target iteration limit")
end
function capmain()
    mkpath(OUTP)
    files=[@__FILE__,CAP_SOURCE,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(files,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    signature=bytes2hex(sha256(string(VERSION,TARGETS_CAP,PREP,STP,BEFOREP,POSTP)*join([bytes2hex(sha256(read(f))) for f in files])))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    rows=[];roots=[];params=[];machine_rows=[]
    xp=copy(PREP);xq=copy(STP);B=0.;lastB=0.;tick=time()
    for (i,target) in enumerate(TARGETS_CAP)
        xp,B=captarget(xp,B,target)
        xq=capfollow(xq,lastB,B,POSTP,"target$(i)_post")
        post_found=xq!==nothing
        println("TARGET ",target," B=",B," preV=",busvp(xp,5)," post_found=",post_found);flush(stdout)
        for tau in [.1,1.2]
            tag=@sprintf("v%03d_tau%03d",round(Int,100target),round(Int,100tau))
            pre=capparams(BEFOREP,B,tau);post=capparams(POSTP,B,tau)
            on=copy(pre);NWState(NW,copy(xp),on,0.).p[VIndex(4,:fault_load₊fault_enabled_p)]=1.
            for (stage,x,p) in [("pre",xp,pre),("post",xq,post)]
                x===nothing && continue
                er,ev,n=eqcheck(x,p);@assert er<1e-8 && ev<0
                push!(roots,[tag,stage,er,ev,n,busvp(x,5),phasep(x)])
                savevec(joinpath(OUTP,tag*"_"*stage*"_root.csv"),x)
            end
            for (stage,p) in [("pre",pre),("fault",on),("post",post)]
                savevec(joinpath(OUTP,tag*"_"*stage*"_parameters.csv"),p)
                s=NWState(NW,copy(xp),p,0.)
                G=s.p[VIndex(5,:PQLoad₊P_load_p)]/s.p[VIndex(5,:PQLoad₊vref2)]
                bc=-s.p[VIndex(5,:PQLoad₊Q_load_p)]/s.p[VIndex(5,:PQLoad₊vref2)]
                @assert abs(bc-B)<1e-12
                for (j,b) in enumerate(ALLG)
                    tp=s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))];tq=s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_q"))];ps=s.p[VIndex(b,Symbol("inv$(b)₊droop₊Pset"))]
                    @assert tp==tau && tq==.02 && ps==PS010[j]
                    push!(params,[tag,stage,b,tp,tq,ps,G,bc])
                end
            end
            a=integratep(tag*"_pre",xp,pre,0.,1.;step=.01)
            f=integratep(tag*"_fault",a.u[end],on,1.,1.15;step=.0005)
            t=[0.];u=[endpointp(f.u[end],post)[1]];last=0.;solve_s=a.solve_s+f.solve_s
            for stop in [30.,60.,150.]
                result=integratep(tag*"_post"*string(Int(stop)),u[end],post,last,stop)
                appendp!(t,u,result);last=stop;solve_s+=result.solve_s
            end
            targetroot=post_found ? xq : fill(NaN,length(xp))
            label,turns,ms=classify_probe(t,u,targetroot)
            for m in ms;push!(machine_rows,vcat([tag],collect(m)));end
            G=NWState(NW,copy(xp),pre,0.).p[VIndex(5,:PQLoad₊P_load_p)]/NWState(NW,copy(xp),pre,0.).p[VIndex(5,:PQLoad₊vref2)]
            postV=post_found ? busvp(xq,5) : NaN
            push!(rows,[tag,target,B,tau,busvp(xp,5),postV,G*busvp(xp,5)^2,G*postV^2,B*busvp(xp,5)^2,B*postV^2,post_found,label,turns,minimum(busvp(x,5) for x in f.u),minimum(busvp(x,5) for x in u),solve_s])
            savevec(joinpath(OUTP,tag*"_clear.csv"),u[1]);savevec(joinpath(OUTP,tag*"_end.csv"),u[end])
            trajectoryp(tag,t,u)
            tablep(tag*"_fault.csv","t_absolute_s,theta_infinite_bus_rad,df_mean_Hz,V5_pu",([tt,phasep(x),freqp(x),busvp(x,5)] for (tt,x) in zip(f.t,f.u)))
            tablep("cases.csv","case,target_pre_V5_pu,B_pu,taup,pre_V5_pu,post_eq_V5_pu,pre_load_P_pu,post_load_P_pu,pre_cap_Q_pu,post_cap_Q_pu,post_root_found,classification,remote_turns,fault_min_V5_pu,post_min_V5_pu,solve_s",rows)
            tablep("equilibria.csv","case,stage,residual,max_real_eigen,unstable_modes,V5_pu,theta_rad",roots)
            tablep("parameters.csv","case,stage,bus,taup,tauq,Pset,load_G_pu,cap_B_pu",params)
            tablep("machines.csv","case,bus,classification,turns,tail_max_abs_df_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",machine_rows)
            println("CASE ",tag," B=",B," preV=",busvp(xp,5)," postV=",postV," ",label," turns=",turns);flush(stdout)
        end
        if !post_found; xq=copy(STP);lastB=0.;else;lastB=B;end
    end
    println("DONE capacitor probe; batch_wall_s=",time()-tick);flush(stdout)
end
capmain()
