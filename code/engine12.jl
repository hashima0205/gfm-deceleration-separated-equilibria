const FAULTBUS=MODE=="line" ? 9 : 4
const OPENED=MODE=="line" ? [7,12] : [7]
function tauall!(p,tau)
    s=NWState(NW,copy(PRE),p,0.)
    for b in ALLG
        s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]=tau
        @assert s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_q"))]==.02
    end
    p
end
function algfix!(x,p)
    f=similar(x)
    for k in 1:12
        NW(f,x,p,0.);maximum(abs.(f[ALG]))<1e-10&&return maximum(abs.(f[ALG]))
        J=fd(NW,x,p);x[ALG]-=J[ALG,ALG]\f[ALG]
    end
    NW(f,x,p,0.);maximum(abs.(f[ALG]))
end
function eqcheck(x,p)
    f=similar(x);NW(f,x,p,0.);J=fd(NW,x,p)
    e=eigvals(Diagonal(1 ./ MASS[DYN])*(J[DYN,DYN]-J[DYN,ALG]*(J[ALG,ALG]\J[ALG,DYN])))
    maximum(abs.(f)),maximum(real.(e)),sum(real.(e).>1e-7)
end
function appendrow(path,header,row)
    fresh=!isfile(path);open(path,"a") do io
        fresh&&println(io,header);println(io,join(row,","));flush(io)
    end
end
function savevec(path,x)
    writedlm(path*".tmp",x,',');mv(path*".tmp",path;force=true)
end
function stages(tau,event)
    pre=tauall!(copy(P_PRE),tau);s=NWState(NW,copy(PRE),pre,0.)
    s.p[VIndex(FAULTBUS,:fault_load₊fault_enabled_p)]=0.;s.p[VIndex(FAULTBUS,:fault_load₊G_fault_p)]=-20.
    on=copy(pre);s=NWState(NW,copy(PRE),on,0.)
    s.p[VIndex(FAULTBUS,:fault_load₊fault_enabled_p)]=event in ["main","fault_only"] ? 1. : 0.
    post=copy(pre);s=NWState(NW,copy(PRE),post,0.)
    if event in ["main","opening_only"]
        for j in OPENED;s.p[EIndex(Symbol("line",j),:line_sym₊active)]=0.;end
        if MODE=="line";s.p[VIndex(FAULTBUS,:fault_load₊fault_enabled_p)]=1.;end
    end
    pre,on,post
end
const KEYNAMES=["delta_rad","omega_rad_s","Pf_pu","Qf_pu","Pcap_pu","Qcap_pu","Vc_r_pu","Vc_i_pu","Ig_r_pu","Ig_i_pu","If_r_pu","If_i_pu","U_r_pu","U_i_pu","Pterminal_pu","Qterminal_pu","ref_in_d_pu","ref_in_q_pu","ref_out_d_pu","ref_out_q_pu","limit_pu","saturated","VI_r_pu","VI_i_pu"]
function obs(x,p,t)
    s=NWState(NW,x,p,t);values=Float64[]
    for b in ALLG
        v(k)=s[VIndex(b,Symbol("inv$(b)₊",k))];a(k)=s.p[VIndex(b,Symbol("inv$(b)₊vsrc₊",k))]
        vc=complex(v("vsrc₊filter₊V_C_r"),v("vsrc₊filter₊V_C_i"));ig=complex(v("vsrc₊filter₊i_g_r"),v("vsrc₊filter₊i_g_i"))
        u=complex(s[VIndex(b,:busbar₊u_r)],s[VIndex(b,:busbar₊u_i)]);pc=vc*conj(ig);pt=u*conj(ig)
        ud=-a("VC_Fcoupl")*a("C")*v("vsrc₊vc₊V_C_q")+a("VC_KP")*(v("vsrc₊vc₊V_C_ref_d")-v("vsrc₊vc₊V_C_d"))+a("VC_KI")*v("vsrc₊vc₊γ_d")+a("VC_F")*v("vsrc₊vc₊i_g_d")
        uq=a("VC_Fcoupl")*a("C")*v("vsrc₊vc₊V_C_d")+a("VC_KP")*(v("vsrc₊vc₊V_C_ref_q")-v("vsrc₊vc₊V_C_q"))+a("VC_KI")*v("vsrc₊vc₊γ_q")+a("VC_F")*v("vsrc₊vc₊i_g_q")
        od=v("vsrc₊vc₊i_f_ref_d");oq=v("vsrc₊vc₊i_f_ref_q");fac=min(1.,a("I_max")/sqrt(ud^2+uq^2+1e-6))
        @assert max(abs(od-fac*ud),abs(oq-fac*uq))<1e-8
        append!(values,[v("droop₊δ"),v("droop₊ω"),v("droop₊Pfilt"),v("droop₊Qfilt"),real(pc),imag(pc),real(vc),imag(vc),real(ig),imag(ig),v("vsrc₊filter₊i_f_r"),v("vsrc₊filter₊i_f_i"),real(u),imag(u),real(pt),imag(pt),ud,uq,od,oq,a("I_max"),fac<1.,v("vsrc₊filter₊V_I_r"),v("vsrc₊filter₊V_I_i")])
    end
    f=similar(x);NW(f,x,p,t);push!(values,maximum(abs.(f[ALG])));values
end
function auditparams(tag,stage,p,tau)
    s=NWState(NW,copy(PRE),p,0.)
    for b in ALLG
        vals=[s.p[VIndex(b,Symbol("inv$(b)₊",k))] for k in ["droop₊τ_p","droop₊τ_q","droop₊Pset","droop₊Qset","droop₊Kp","droop₊Kq","vsrc₊VC_F","vsrc₊CC1_F","vsrc₊VC_KP","vsrc₊VC_KI","vsrc₊CC1_KP","vsrc₊CC1_KI","vsrc₊I_max"]]
        @assert vals[1]==tau&&vals[2]==.02&&vals[3]>0&&vals[7]==1&&vals[8]==1
        appendrow(joinpath(DATA,"parameter_audit.csv"),"case,stage,bus,tau_p_s,tau_q_s,Pset_pu,Qset_pu,Kp,Kq,VC_F,CC1_F,VC_KP,VC_KI,CC1_KP,CC1_KI,Imax_pu",vcat([tag,stage,b],vals))
    end
end
function segment(tag,label,x,p,a,b,clear;level="A",solver="Rodas5P",step=.02)
    dir=joinpath(CHECK,tag);mkpath(dir);file=joinpath(dir,label*".csv");state=joinpath(dir,label*"_end.csv")
    if isfile(file)&&isfile(state);return vec(readdlm(state,',',Float64)),"Success";end
    x=copy(x);before=copy(x[DYN]);err=algfix!(x,p);@assert err<1e-8&&x[DYN]==before
    savevec(joinpath(dir,label*"_initial.csv"),x);savevec(joinpath(dir,label*"_parameters.csv"),p)
    rt,at,dt=level=="A" ? (1e-9,1e-11,.005) : (1e-10,1e-12,.0025)
    dense_times=collect(a:.0001:min(b,a+.02));samples=sort(unique(vcat(collect(a:step:b),dense_times,[a,b])))
    tick=time();sol=solve(ODEProblem(NW,x,(a,b),copy(p)),solver=="Rodas4P" ? Rodas4P() : Rodas5P();reltol=rt,abstol=at,dtmax=dt,maxiters=10^8,saveat=samples,save_everystep=false,dense=false)
    open(file*".tmp","w") do io
        println(io,join(vcat(["t_absolute_s","t_post_s","stage"],["b$(b)_$(k)" for b in ALLG for k in KEYNAMES],"alg_residual"),","))
        for (t,u) in zip(sol.t,sol.u);println(io,join(vcat([t,t-clear,label],obs(u,p,t)),","));end
    end
    mv(file*".tmp",file;force=true)
    appendrow(joinpath(dir,"progress.csv"),"stage,end_time,retcode,runtime_s,solver,reltol,abstol,dtmax,regular_save_s",[label,sol.t[end],sol.retcode,time()-tick,solver,rt,at,dt,step])
    string(sol.retcode)=="Success"&&savevec(state,sol.u[end])
    copy(sol.u[end]),string(sol.retcode)
end
function recovercheck(x,p,target)
    s=NWState(NW,x,p,0.);ss=NWState(NW,target,p,0.)
    df=maximum(abs(s[VIndex(b,Symbol("inv$(b)₊droop₊ω"))]/(2pi)-60) for b in ALLG)
    er=maximum(abs(atan(sin(s[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]-ss[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]),cos(s[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]-ss[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]))) for b in ALLG)
    df,er
end
function physical(tag,tau,T,event;level="A",solver="Rodas5P",force_horizon=30.)
    println("START ",tag);flush(stdout)
    try
        pre,on,post=stages(tau,event)
        for (stage,p) in [("pre",pre),("fault",on),("post",post)];auditparams(tag,stage,p,tau);end
        target=event in ["main","opening_only"] ? STABLE : PRE
        for (stage,x,p) in [("pre",PRE,pre),("post",target,post)]
            er,ev,n=eqcheck(x,p);appendrow(joinpath(DATA,"equilibria.csv"),"case,stage,residual,max_real_eigen,unstable_modes",[tag,stage,er,ev,n]);@assert er<1e-8&&ev<0
        end
        clear=1+T
        x,rc=segment(tag,"pre",PRE,pre,0.,1.,clear;level=level,solver=solver,step=.01);@assert rc=="Success"
        x,rc=segment(tag,"fault",x,on,1.,clear,clear;level=level,solver=solver,step=.0005);@assert rc=="Success"
        x,rc=segment(tag,"post030",x,post,clear,clear+30,clear;level=level,solver=solver);@assert rc=="Success"
        df,er=recovercheck(x,post,target);horizon=30.
        if df>=1e-5||er>=.001||force_horizon>30
            for stop in [100.,200.,300.]
                x,rc=segment(tag,"post"*lpad(Int(stop),3,'0'),x,post,clear+horizon,clear+stop,clear;level=level,solver=solver,step=.05);horizon=stop;rc=="Success"||break
            end
        end
        appendrow(joinpath(DATA,"cases.csv"),"case,model,tau_s,T_s,event,horizon_s,level,solver,retcode",[tag,MODE=="line" ? "line" : "bus",tau,T,event,horizon,level,solver,rc])
        println("DONE ",tag," horizon=",horizon," final=",recovercheck(x,post,target));flush(stdout)
    catch e
        appendrow(joinpath(DATA,"cases.csv"),"case,model,tau_s,T_s,event,horizon_s,level,solver,retcode",[tag,MODE,tau,T,event,NaN,level,solver,"exception"])
        println("ERROR ",tag," ",sprint(showerror,e,catch_backtrace()));flush(stdout)
    end
end
function basin_batch()
    p=stages(.02,"main")[3]
    for (name,x) in [("stable",STABLE),("saddle",UNSTABLE)]
        e,v,n=eqcheck(x,p);appendrow(joinpath(DATA,"basin_equilibria.csv"),"name,residual,max_real_eigen,unstable_modes",[name,e,v,n])
    end
    for row in eachrow(readdlm(joinpath(DATA,"basin_manifest.csv"),',',Any;skipstart=1))
        tag=string(row[1]);theta=Float64(row[2]);df=Float64(row[3]);kind=string(row[4]);println("START ",tag);flush(stdout)
        try
            x=vec(readdlm(joinpath(DATA,"seeds",tag*".csv"),',',Float64));old=copy(x[DYN]);er=algfix!(x,p);s=NWState(NW,x,p,0.)
            th=sum(s[VIndex(b,Symbol("inv$(b)₊droop₊δ"))] for b in REMOTE)/3
            fq=sum(s[VIndex(b,Symbol("inv$(b)₊droop₊ω"))]/(2pi)-60 for b in REMOTE)/3
            @assert er<1e-8&&maximum(abs.(x[DYN]-old))==0&&abs(th-theta)<1e-9&&abs(fq-df)<1e-9
            appendrow(joinpath(DATA,"basin_initialization.csv"),"case,theta_rad,df_Hz,kind,alg_residual,theta_error_rad,df_error_Hz,dynamic_change",[tag,theta,df,kind,er,th-theta,fq-df,maximum(abs.(x[DYN]-old))])
            y,rc=segment(tag,"post030",x,p,0.,30.,0.;step=.05);h=30.;fq,pe=recovercheck(y,p,STABLE)
            if fq>=1e-5||pe>=.001||kind!="coarse"
                y,rc=segment(tag,"post300",y,p,30.,300.,0.;step=.1);h=300.
            end
            fq,pe=recovercheck(y,p,STABLE)
            if tag in ["basin_005","basin_006","basin_007","basin_008"] || (kind!="coarse" && (fq>=1e-5||pe>=.001))
                for stop in [600.,1000.]
                    y,rc=segment(tag,"post"*string(Int(stop)),y,p,h,stop,0.;step=.1);h=stop
                end
            end
            appendrow(joinpath(DATA,"basin_runs.csv"),"case,theta_rad,df_Hz,kind,horizon_s,retcode",[tag,theta,df,kind,h,rc])
            println("DONE ",tag," ",h);flush(stdout)
        catch e
            println("ERROR ",tag," ",sprint(showerror,e));flush(stdout)
            appendrow(joinpath(DATA,"basin_runs.csv"),"case,theta_rad,df_Hz,kind,horizon_s,retcode",[tag,theta,df,kind,NaN,"initialization_or_solver_failure"])
        end
    end
end
