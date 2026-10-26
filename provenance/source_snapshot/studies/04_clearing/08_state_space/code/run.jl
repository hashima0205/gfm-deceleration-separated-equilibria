# Reuse model definitions, without executing or overwriting the previous study.
const BASIN=normpath(joinpath(@__DIR__,"..","output"))
basefile=normpath(joinpath(@__DIR__,"../../07_grid_strength/code/run.jl"))
definitions=first(split(read(basefile,String),"\ninclude(\"seed_support.jl\")"))
include_string(Main,definitions,basefile)
const BB,PRE,P_PRE,POST,P_POST,BASEERR=baseop
const REMOTE=[2,3,7]
const MASS=diag(Matrix(BB.f.mass_matrix))
const ALG=findall(iszero,MASS)
const IDX=NWState(BB.f,Float64.(1:length(POST)),copy(P_POST),0.)
const DELTA=[Int(IDX[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]) for b in REMOTE]
const PF=[Int(IDX[VIndex(b,Symbol("inv$(b)₊droop₊Pfilt"))]) for b in REMOTE]
const FAST=setdiff(1:length(POST),vcat(DELTA,PF))
const OBSIDS=Any[]
const OBSLABELS=String[]
for b in [1,2,3,7]
    for key in ["droop₊δ","droop₊ω","droop₊Pfilt","droop₊Qfilt","droop₊Pmeas","vsrc₊filter₊i_f_r","vsrc₊filter₊i_f_i","vsrc₊filter₊V_C_r","vsrc₊filter₊V_C_i","vsrc₊vc₊γ_d","vsrc₊vc₊γ_q"]
        push!(OBSIDS,VIndex(b,Symbol("inv$(b)₊",key)));push!(OBSLABELS,"b$(b)_$key")
    end
end
const SEEDROWS=[split(l,",") for l in readlines(joinpath(BASIN,"section_states.csv"))[2:end]]
function params(tau)
    p=copy(P_POST);s=NWState(BB.f,copy(POST),p,0.)
    for b in REMOTE;s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]=tau;end
    @assert s.p[EIndex(:line7,:line_sym₊active)]==0
    @assert s.p[VIndex(4,:fault_load₊G_fault_p)]==0
    p
end
function section(tag)
    x=copy(POST);s=NWState(BB.f,x,copy(P_POST),0.)
    for row in SEEDROWS
        row[1]==tag||continue
        s[VIndex(parse(Int,row[2]),Symbol(row[3]))]=parse(Float64,row[4])
    end
    x
end
function algebraic!(x,p)
    f=similar(x)
    for iter in 1:8
        BB.f(f,x,p,0.)
        maximum(abs.(f[ALG]))<1e-10&&return maximum(abs.(f[ALG]))
        j=fd(BB.f,x,p);x[ALG]-=j[ALG,ALG]\f[ALG]
    end
    BB.f(f,x,p,0.);maximum(abs.(f[ALG]))
end
function trial(tag,x,p,horizon;verify=false)
    err=algebraic!(x,p);@assert err<1e-8
    f=similar(x);BB.f(f,x,p,0.);fast=maximum(abs.(f[FAST]))
    theta0=sum(x[DELTA])/3
    println("START ",tag," theta=",theta0);flush(stdout)
    started=time()
    sol=solve(ODEProblem(BB.f,x,(0.,horizon),copy(p)),verify ? Rodas4P() : Rodas5P();
        reltol=1e-9,abstol=1e-11,dtmax=verify ? .0005 : .01,saveat=.1,save_everystep=false,dense=false)
    vals=zeros(length(sol.t),length(OBSIDS))
    for (j,(t,u)) in enumerate(zip(sol.t,sol.u))
        s=NWState(BB.f,copy(u),copy(p),t)
        vals[j,:]=[s[id] for id in OBSIDS]
    end
    open(joinpath(BASIN,"runs",tag*".csv"),"w") do io
        println(io,join(vcat("t",OBSLABELS),","))
        for j in eachindex(sol.t);println(io,join(vcat(sol.t[j],vals[j,:]),","));end
    end
    open(joinpath(BASIN,"runs",tag*"_initial.csv"),"w") do io
        println(io,"index,value");for j in eachindex(x);println(io,"$j,$(x[j])");end
    end
    thetaf=sum(sol.u[end][DELTA])/3
    ka=findall(>=(max(0.,sol.t[end]-10)),sol.t)
    frequencies=hcat([vals[:,findfirst(==("b$(b)_droop₊ω"),OBSLABELS)]/(2pi).-60 for b in REMOTE]...)
    phase=hcat([vals[:,findfirst(==("b$(b)_droop₊δ"),OBSLABELS)] for b in REMOTE]...)
    target=POST[DELTA]
    phaseerr=maximum(abs.(atan.(sin.(phase[end,:]-target),cos.(phase[end,:]-target))))
    maxdf=maximum(abs.(frequencies[ka,:]));meandf=sum(frequencies[ka,:])/length(frequencies[ka,:])
    turns=(thetaf-theta0)/(2pi)
    label=string(sol.retcode)!="Success" ? "numerical_failure" :
        (maxdf<1e-4 && phaseerr<.002 ? "recovered" :
        (turns < -2 && maximum(frequencies[ka,:])<-.001 ? "negative_running" : "unresolved"))
    coherence=maximum(maximum(phase,dims=2)-minimum(phase,dims=2))
    maxI=maximum(maximum(hypot.(vals[:,findfirst(==("b$(b)_vsrc₊filter₊i_f_r"),OBSLABELS)],vals[:,findfirst(==("b$(b)_vsrc₊filter₊i_f_i"),OBSLABELS)])) for b in REMOTE)
    maxalg=0.
    for j in 1:10:length(sol.t);BB.f(f,sol.u[j],p,sol.t[j]);maxalg=max(maxalg,maximum(abs.(f[ALG])));end
    println("DONE ",tag," ",label," runtime=",time()-started);flush(stdout)
    return [label,string(sol.retcode),sol.t[end],theta0,frequencies[1,1],thetaf,turns,meandf,maxdf,phaseerr,coherence,maxI,err,maxalg,fast,time()-started]
end
function batchfile(manifest,status;verify=false)
    rows=readlines(manifest);cols=split(rows[1],",")
    open(status,"w") do io
        println(io,rows[1]*",classification,retcode,last_t,actual_theta0,actual_df2_0,theta_end,turns,tail_df,tail_max_df,phase_error,phase_spread_max,remote_max_current,initial_alg_residual,max_alg_residual,initial_fast_residual,runtime_s");flush(io)
        for row in rows[2:end]
            a=split(row,",");tag=String(a[1]);tau=parse(Float64,a[3]);horizon=parse(Float64,a[8])
            try
                x=section(tag);p=params(tau)
                result=trial(tag,x,p,horizon;verify=verify)
                println(io,row*","*join(result,","));flush(io)
            catch e
                println("ERROR ",tag," ",sprint(showerror,e));flush(stdout)
                println(io,row*",error,"*join(fill("NaN",15),","));flush(io)
            end
        end
    end
end
if get(ENV,"BASIN_DEFINITIONS_ONLY","0")!="1"
    batchfile(joinpath(BASIN,"manifest.csv"),joinpath(BASIN,"full_status.csv"))
end
