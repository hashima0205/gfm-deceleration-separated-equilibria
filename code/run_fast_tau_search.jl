# Persistent single Julia process: fixed-B near-fold search plus compensation sweep.
const FAST_BASE=joinpath(@__DIR__,"run_capacitor_k_strength.jl")
fastsource=first(split(read(FAST_BASE,String),"\nfunction kmain()"))
fastsource=replace(fastsource,"basin_recalculation/capacitor_k_strength\\\""=>"basin_recalculation/fast_tau_search/runs\\\"")
include_string(Main,fastsource,FAST_BASE)
@assert endswith(OUTP,"fast_tau_search/runs")
const FASTDIR=dirname(OUTP)
function fastmove(x,base,fromk,fromB,k,B,tag)
    if B==fromB;return kfollow(x,base,fromk,fromB,k,B,tag);end
    safe=.9*min(fromk,k)
    x=kfollow(x,base,fromk,fromB,safe,fromB,tag*"_strengthen")
    x=kfollow(x,base,safe,fromB,safe,B,tag*"_compensate")
    kfollow(x,base,safe,B,k,B,tag*"_approach")
end
function fastmain()
    mkpath(OUTP)
    files=[@__FILE__,FAST_BASE,K_SOURCE,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(files,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    signature=bytes2hex(sha256(string(VERSION,POSTP,BEFOREP)*join(bytes2hex(sha256(read(f))) for f in files)))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    # Case files can be appended while this process remains alive. A unique
    # input signature per tag prevents reuse with different physical settings.
    processed=Set{String}();rows=[];roots=[];machines=[]
    xp=vec(readdlm(joinpath(PACK,"basin_recalculation/capacitor_k_strength/k080_pre.csv"),',',Float64))
    xs=vec(readdlm(joinpath(PACK,"basin_recalculation/capacitor_k_strength/k080_stable.csv"),',',Float64))
    lastk=.8;lastB=KB;tick=time()
    while true
        cases=readdlm(joinpath(FASTDIR,"candidates.csv"),',',Any;header=true)
        data,header=cases;names=vec(String.(header))
        col(name)=findfirst(==(name),names)
        for i in axes(data,1)
            tag=String(data[i,col("case")]);tag in processed && continue
            B=Float64(data[i,col("B_pu")]);k=Float64(data[i,col("k")]);critical=Float64(data[i,col("critical_k")])
            fg=col("fault_G_pu")===nothing ? 20. : Float64(data[i,col("fault_G_pu")])
            T=col("fault_duration_s")===nothing ? .15 : Float64(data[i,col("fault_duration_s")])
            expectedV=Float64(data[i,col("pre_V5_pu")])
            inputhash=bytes2hex(sha256(string(signature,k,B,fg,T,.1)))
            marker=joinpath(OUTP,tag*"_inputs.sha256")
            if isfile(marker);@assert strip(read(marker,String))==inputhash;else;write(marker,inputhash*"\n");end
            println("START ",tag," k=",k," B=",B," Gfault=",fg);flush(stdout)
            xp=fastmove(xp,BEFOREP,lastk,lastB,k,B,tag*"_pre")
            xs=fastmove(xs,POSTP,lastk,lastB,k,B,tag*"_stable");lastk=k;lastB=B
            pre=kparams(BEFOREP,k,B,.1);post=kparams(POSTP,k,B,.1)
            @assert abs(busvp(xp,5)-expectedV)<1e-6
            stable=true
            for (name,x,p) in [("pre",xp,pre),("post",xs,post)]
                er,ev,n=eqcheck(x,p);@assert er<1e-8
                stable &= ev<0
                push!(roots,[tag,name,er,ev,n,busvp(x,5),phasep(x)])
                savevec(joinpath(OUTP,tag*"_"*name*"_root.csv"),x)
                savevec(joinpath(OUTP,tag*"_"*name*"_parameters.csv"),p)
            end
            @assert stable "Candidate has no locally stable pre/post pair: $tag"
            on=copy(pre);ss=NWState(NW,copy(xp),on,0.)
            ss.p[VIndex(4,:fault_load₊fault_enabled_p)]=1.
            ss.p[VIndex(4,:fault_load₊G_fault_p)]=-fg
            savevec(joinpath(OUTP,tag*"_fault_parameters.csv"),on)
            a=integratep(tag*"_pre",xp,pre,0.,1.;step=.01)
            f=integratep(tag*"_fault",a.u[end],on,1.,1.0 + T;step=.0005)
            t=[0.];u=[endpointp(f.u[end],post)[1]];last=0.;label="unresolved";turns=0.;ms=[]
            for stop in [30.,60.,150.,300.]
                result=integratep(tag*"_post"*string(Int(stop)),u[end],post,last,stop)
                appendp!(t,u,result);last=stop;label,turns,ms=classify_probe(t,u,xs)
                label in ["unresolved","other"] || break
            end
            beforeV=busvp(xp,5);postV=busvp(xs,5)
            push!(rows,[tag,k,B,critical,(critical-k)/critical,fg,T,beforeV,postV,2.2682301152333055*beforeV^2,B*beforeV^2,label,last,turns])
            for m in ms;push!(machines,vcat([tag],collect(m)));end
            savevec(joinpath(OUTP,tag*"_clear.csv"),u[1]);savevec(joinpath(OUTP,tag*"_end.csv"),u[end]);trajectoryp(tag,t,u)
            tablep(tag*"_fault.csv","t_absolute_s,theta_infinite_bus_rad,df_mean_Hz,V5_pu",([tt,phasep(x),freqp(x),busvp(x,5)] for (tt,x) in zip(f.t,f.u)))
            tablep("results.csv","case,k,B_pu,critical_k,k_margin_fraction,fault_G_pu,fault_duration_s,pre_V5_pu,post_eq_V5_pu,pre_load_P_pu,pre_cap_Q_pu,classification,horizon_s,remote_turns",rows)
            tablep("equilibria.csv","case,stage,residual,max_real_eigen,unstable_modes,V5_pu,mean_phase_rad",roots)
            tablep("machines.csv","case,bus,classification,turns,tail_max_abs_df_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",machines)
            push!(processed,tag)
            println("RESULT ",tag," ",label," horizon=",last," turns=",turns," preV=",beforeV);flush(stdout)
        end
        write(joinpath(FASTDIR,"ready.txt"),string(length(processed))*"\n")
        if isfile(joinpath(FASTDIR,"finish.flag"));break;end
        sleep(.5)
    end
    println("DONE fast-tau search; cases=",length(processed)," batch_wall_s=",time()-tick);flush(stdout)
end
fastmain()
