# Stable equilibrium existence is checked for EVERY controller setting.
using DelimitedFiles
const HERE=dirname(@__DIR__)
const OUT=joinpath(HERE,"output")
old=joinpath(HERE,"../05_bus1_remote/code/equilibrium/bootstrap.jl")
include_string(Main,first(split(read(old,String),"\nfor pl in")),old)
function fd(nw,x,p;h=1e-6)
    n=length(x);J=zeros(n,n);a=zeros(n);b=zeros(n)
    for j in 1:n
        step=h*max(1.,abs(x[j]));xp=copy(x);xm=copy(x);xp[j]+=step;xm[j]-=step
        nw(a,xp,p,0.);nw(b,xm,p,0.);J[:,j]=(a-b)/(2step)
    end
    J
end
function stationary(nw,x,p)
    x=copy(x);f=zeros(length(x))
    for iter in 1:40
        nw(f,x,p,0.);err=maximum(abs.(f));err<1e-8 && return x,err
        step=fd(nw,x,p)\f;ok=false
        for k in 0:14
            q=x-step*2.0^(-k);ff=similar(f);nw(ff,q,p,0.)
            if maximum(abs.(ff))<err;x=q;ok=true;break;end
        end
        ok||break
    end
    nw(f,x,p,0.);return x,maximum(abs.(f))
end
function stability(bb,x,p)
    J=fd(bb.f,x,p);m=diag(Matrix(bb.f.mass_matrix));di=findall(!iszero,m);al=findall(iszero,m)
    maximum(real.(eigvals(J[di,di]-J[di,al]*(J[al,al]\J[al,di]))))
end
function operating(pl)
    bb,u,p,s,vref=prepare_eq(pl)
    pre,err=stationary(bb.f,u,p);err<1e-8||error("prefault equilibrium failed")
    post=copy(pre);pp=copy(p);q=NWState(bb.f,post,pp,0.)
    q.p[EIndex(:line7,:line_sym₊active)]=0.
    # Seed from the independently verified low-load stable equilibrium.
    # A line-admittance homotopy need not follow the desired endpoint branch.
    known=joinpath(HERE,"../05_bus1_remote/output/equilibrium/full_equilibrium_states.csv")
    for row in [split(l,",") for l in readlines(known)[2:end]]
        row[1]=="L0.8_branch0"||continue
        q[VIndex(parse(Int,row[3]),Symbol(row[4]))]=parse(Float64,row[5])
    end
    xx,e=stationary(bb.f,post,pp);post.=xx;e<1e-8||error("postfault equilibrium failed residual $e")
    return bb,pre,p,post,pp,err
end
function simulate(bb,pre,p,pl,tau,dur,horizon,tag;verify=false,perturb=false,posteq=nothing)
    uu=copy(pre);pp=copy(p);s=NWState(bb.f,uu,pp,0.)
    for b in [2,3,7];s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]=tau;end
    s.p[VIndex(4,:fault_load₊G_fault_p)]=perturb ? 0. : -20.
    s.p[VIndex(4,:fault_load₊t_clear_p)]=1+dur
    if perturb
        uu.=posteq;s.p[EIndex(:line7,:line_sym₊active)]=0.
        for b in [2,3,7];s[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]+=.01;end
    end
    ids=Any[];labels=String[]
    for b in [1,2,3,7]
        for k in ["droop₊δ","droop₊ω","droop₊Pfilt","vsrc₊filter₊i_f_r","vsrc₊filter₊i_f_i","vsrc₊cc1₊i_f_ref_d","vsrc₊cc1₊i_f_ref_q"]
            push!(ids,VIndex(b,Symbol("inv$(b)₊",k)));push!(labels,"b$(b)_$k")
        end
    end
    for b in 1:8, k in ["u_r","u_i","i_r","i_i"]
        push!(ids,VIndex(b,Symbol("busbar₊",k)));push!(labels,"net_b$(b)_$k")
    end
    clear=1+dur;last=0.;status="Success"
    open(joinpath(OUT,"runs",tag*".csv"),"w") do io
        println(io,join(vcat("t",labels),","))
        intervals=perturb ? [(0.,horizon)] : [(0.,clear),(clear,horizon)]
        for (j,(a,b)) in enumerate(intervals)
            if j==2;s.p[EIndex(:line7,:line_sym₊active)]=0.;end
            stops=filter(t->a<t<=b,[1.,clear])
            sample=sort(unique(vcat(collect(a:.02:b),[a,b],filter(t->a<=t<=b,collect(.995:.0002:1.01)),filter(t->a<=t<=b,collect(clear-.005:.0002:clear+.02)))))
            sol=solve_fault_switched(ODEProblem(bb.f,uu,(a,b),copy(pp)),verify ? Rodas4P() : Rodas5P();reltol=1e-10,abstol=1e-12,dtmax=verify ? .0005 : .005,saveat=sample,save_everystep=false,dense=false,tstops=stops)
            for t in sol.t
                vals=sol(t;idxs=ids,continuity=:right);println(io,join(vcat(t,vals),","));last=t
            end
            status=string(sol.retcode);status=="Success"||break
            uu=copy(sol.u[end])
        end
    end
    return status,last
end
track=length(ARGS)>0 ? ARGS[1] : "screen"
manifest=length(ARGS)>1 ? ARGS[2] : joinpath(HERE,"code/screen.csv")
specs=[parse.(Float64,split(x,",")) for x in readlines(manifest)[2:end] if !isempty(x)]
cache=Dict{Float64,Any}()
open(joinpath(OUT,track*"_status.csv"),"w") do io
    println(io,"case,load,tau,duration,horizon,max_real_eigenvalue,eq_residual,prefault_residual,status,last_t,perturbation")
    for (i,(pl,tau,dur,horizon,pert)) in enumerate(specs)
        tag="$(track)_$(lpad(i,3,'0'))_L$(pl)_tau$(tau)_T$(dur)";println("START ",tag);flush(stdout)
        try
            op=get!(()->operating(pl),cache,pl);bb,pre,p,post,pp,preerr=op
            eqp=copy(pp);q=NWState(bb.f,copy(post),eqp,0.)
            for b in [2,3,7];q.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]=tau;end
            f=zeros(length(post));bb.f(f,post,eqp,0.);eqerr=maximum(abs.(f));ev=stability(bb,post,eqp)
            eqerr<1e-8 && ev<0||error("Not a stable equilibrium: residual=$eqerr ev=$ev")
            status,last=simulate(bb,pre,p,pl,tau,dur,horizon,tag;verify=startswith(track,"verify"),perturb=pert==1,posteq=post)
            println(io,join([tag,pl,tau,dur,horizon,ev,eqerr,preerr,status,last,pert],","));flush(io)
            println("DONE ",tag," ",status," eigenvalue=",ev);flush(stdout)
        catch e
            println(io,join([tag,pl,tau,dur,horizon,NaN,NaN,NaN,"error",0,pert],","));flush(io)
            println("ERROR ",tag," ",sprint(showerror,e));flush(stdout)
        end
    end
end
