# Preserve the completed 100-point study; refine only the load continuation.
const PROBE_SOURCE=joinpath(@__DIR__,"run_trajectory_voltage_probes.jl")
definitions=first(split(read(PROBE_SOURCE,String),"\nfunction mainp()"))
definitions=replace(definitions,"basin_recalculation/trajectory_voltage_probes\""=>"basin_recalculation/trajectory_voltage_probes/voltage_refined\"")
coarse=r"(?s)        for r in range\(previous_factor,factor;.*?\n        end"
@assert occursin(coarse,definitions)
definitions=replace(definitions,coarse=>"        previous_pre,previous_post=refine_roots(previous_pre,previous_post,previous_factor,factor)")
include_string(Main,definitions,PROBE_SOURCE)

function refine_roots(pre,post,from,to)
    from==to && return copy(pre),copy(post)
    r=from;step=min(.001,from-to);attempts=[];tries=0
    while r>to+1e-12
        tries+=1;@assert tries<5000
        trial=max(to,r-step)
        pp=copy(BEFOREP);pq=copy(POSTP)
        NWState(NW,copy(PREP),pp,0.).p[VIndex(5,:PQLoad₊P_load_p)]*=trial
        NWState(NW,copy(STP),pq,0.).p[VIndex(5,:PQLoad₊P_load_p)]*=trial
        xp,ep=stationary(NW,pre,pp);xq,eq=stationary(NW,post,pq)
        evp=NaN;evq=NaN
        if ep<1e-8 && eq<1e-8
            _,evp,_=eqcheck(xp,pp);_,evq,_=eqcheck(xq,pq)
        end
        accepted=ep<1e-8 && eq<1e-8 && evp<0 && evq<0
        push!(attempts,[from,to,r,trial,step,ep,eq,evp,evq,accepted])
        tablep(@sprintf("continuation_%03d.csv",round(Int,100to)),
               "from_factor,target_factor,last_accepted_factor,trial_factor,step,pre_residual,post_residual,pre_max_real_eigen,post_max_real_eigen,accepted",attempts)
        if accepted
            pre=xp;post=xq;r=trial;step=min(.002,step*1.4)
        else
            println("CONTINUATION retry factor=",trial," pre_res=",ep," post_res=",eq," pre_ev=",evp," post_ev=",evq);flush(stdout)
            step/=2
            @assert step>=1e-7 "Could not continue the locally stable branch; this is not a nonexistence proof."
        end
    end
    pre,post
end

function refined_main()
    mkpath(OUTP)
    files=[@__FILE__,PROBE_SOURCE,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(files,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    signature=bytes2hex(sha256(string(VERSION,STP,PREP,POSTP,BEFOREP,LOAD_FACTORS,VOLTAGE_TAUS)*join([bytes2hex(sha256(read(f))) for f in files])))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    tick=time();voltagep()
    println("DONE refined voltage probe; batch_wall_s=",time()-tick);flush(stdout)
end
refined_main()
