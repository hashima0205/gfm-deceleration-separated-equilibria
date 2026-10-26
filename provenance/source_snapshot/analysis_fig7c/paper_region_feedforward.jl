# Preserve the paper's constant-power load and operating point while changing FF.
# Default integration is retained as a numerical diagnostic, not a stability proof.
using PowerDynamics, ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
root = dirname(@__DIR__)
src = replace(read(joinpath(root, "4gfm_simulation_RL.jl"), String), "using PythonCall" => "")
src = replace(src, "include(\"compornent.jl\")" => "include(" * repr(joinpath(root, "compornent.jl")) * ")")
include_string(Main, src, joinpath(root, "4gfm_simulation_RL.jl"))
mode = isempty(ARGS) ? "center" : ARGS[1]
out = joinpath(@__DIR__, "paper_region_feedforward"); mkpath(out)
Y = zeros(ComplexF64,8,8)
for (a,b,m) in [(1,4,2.14),(4,6,3.43),(5,6,.68),(3,5,2.21),(2,3,1.91),(2,5,2.14),(4,5,1.8),(2,8,2.17),(2,7,6.39),(7,8,1.69),(5,8,8.96)]
    y=m*exp(-im*75pi/180); Y[a,a]+=y;Y[b,b]+=y;Y[a,b]-=y;Y[b,a]-=y
end
Z = inv(Y[[1,2,3,4,5,7,8],[1,2,3,4,5,7,8]])
names = ["V5", "fault_ur", "fault_ui", "fault_ir", "fault_ii"]
ids = Any[VIndex(5,:busbar₊u_mag)]
append!(ids,[VIndex(4,Symbol("busbar₊",k)) for k in ["u_r","u_i","i_r","i_i"]])
for b in [1,2,3,7]
    for k in ["u_r", "u_i"]
        push!(names,"b$(b)_terminal_$(k)"); push!(ids,VIndex(b,Symbol("busbar₊",k)))
    end
    for k in ["droop₊δ","droop₊ω","droop₊Pfilt","vsrc₊cc1₊i_f_d","vsrc₊cc1₊i_f_q","vsrc₊cc1₊i_f_ref_d","vsrc₊cc1₊i_f_ref_q"]
        push!(names,"b$(b)_$(k)"); push!(ids,VIndex(b,Symbol("inv$(b)₊",k)))
    end
end

push!(names,"fault_G");push!(ids,VIndex(4,:fault_load₊G_fault_var))

function transfer(base,vf,cf)
    u=copy(base.u0); p=copy(base.p); s=NWState(base.f,u,p,0.)
    before=zeros(length(u)); base.f(before,u,p,0.)
    for b in [1,2,3,7]
        # Preserve both controller outputs at the original electrical equilibrium.
        for axis in ["d","q"]
            vg=VIndex(b,Symbol("inv$(b)₊vsrc₊vc₊γ_",axis))
            cg=VIndex(b,Symbol("inv$(b)₊vsrc₊cc1₊γ_",axis))
            s[vg] -= vf*s[VIndex(b,Symbol("inv$(b)₊vsrc₊vc₊i_g_",axis))]/.0273
            s[cg] -= cf*s[VIndex(b,Symbol("inv$(b)₊vsrc₊cc1₊V_C_",axis))]/.0018
        end
        s.p[VIndex(b,Symbol("inv$(b)₊vsrc₊VC_F"))]=vf
        s.p[VIndex(b,Symbol("inv$(b)₊vsrc₊CC1_F"))]=cf
    end
    after=zeros(length(u)); base.f(after,u,p,0.)
    err=maximum(abs.(after-before))
    err<1e-7 || error("Equilibrium transfer error $err")
    return u,p,err
end

function jacdet(tt,v,g,pl)
    gt=fault_conductance(tt, g, 1., 1.2)
    z=Z[5,5]+gt*Z[5,4]*Z[4,5]/(1-gt*Z[4,4])
    v2=v^2; e=1e-6
    return 1+2real(z)*pl*e/(v2+e)^2+abs2(z)*pl^2*(e-v2)/(v2+e)^3
end

summary=open(joinpath(out,mode*"_status.csv"),"w")
println(summary,"case,P1,Pload,VC_F,CC1_F,G_fault,method,status,last_t,first_branch_limit_t,transfer_error")
function run_case(base,p1,pl,vf,cf,g,method)
    u,p,err=transfer(base,vf,cf)
    s=NWState(base.f,u,p,0.)
    s.p[VIndex(4,:fault_load₊G_fault_p)]=g
    label="P$(p1)_L$(pl)_VC$(vf)_CC$(cf)_G$(Int(g))_$(method)"
    println("START ",label);flush(stdout)
    ts=0.;status="finished";first_limit=NaN
    strict=method!="default"; verify=method=="verify"
    horizon=strict ? 20. : 5.
    open(joinpath(out,label*".csv"),"w") do io
        println(io,join(vcat("t","success",names,"jacdet"),","))
        q=[s[id] for id in ids]
        println(io,join(vcat(0.,1,q,jacdet(0.,q[1],g,pl)),","))
        while ts<horizon
            fine=.99<=ts<1.21
            te=ts==0. ? 1e-5 : min(horizon,ts+(strict ? (fine ? .0001 : .005) : .01))
            ps=ODEProblem(base.f,u,(ts,te),p)
            stops=filter(x->ts<x<=te,[1.,1.2])
            sol=if strict
                solve_fault_switched(ps,verify ? Rodas4P() : Rodas5P();reltol=verify ? 1e-10 : 1e-9,abstol=verify ? 1e-12 : 1e-11,dtmax=fine ? (verify ? .00005 : .0001) : .001,tstops=stops,d_discontinuities=stops)
            else
                solve_fault_switched(ps,Rodas5P();tstops=stops,d_discontinuities=stops)
            end
            tt=sol.t[end];q=sol(tt;continuity=:right,idxs=ids);ok=string(sol.retcode)=="Success"
            j=jacdet(tt,q[1],g,pl)
            println(io,join(vcat(tt,Int(ok),q,j),","))
            if ok && j<=1e-6 && isnan(first_limit); first_limit=tt; end
            ts=tt
            if !ok;status=string(sol.retcode);break;end
            if strict && j<=1e-6;status="branch_guard";break;end
            u=copy(sol.u[end])
        end
    end
    println(summary,join([label,p1,pl,vf,cf,g,method,status,ts,first_limit,err],","));flush(summary)
    println("DONE ",label," ",status," t=",ts," first_limit=",first_limit);flush(stdout)
end

points=mode=="nearby" ? [(pl,p1) for pl in [1.16,1.2,1.24] for p1 in [.6,.7,.8] if (pl,p1)!=(1.2,.7)] : mode=="verify_corner" ? [(1.24,.8)] : [(1.2,.7)]
for (pl,p1) in points
    far=round(.3+round((.9-p1)/3,digits=3),digits=3)
    init_sim(far,p1,far;P_load=pl,G_fault=20.,t_cle=1.2,max_t=20.)
    base=SIM_STATE["prob"]
    for vf in [0,1], cf in [0,1]
        mode=="verify_corner" && (vf,cf)!=(1,1) && continue
        if startswith(mode,"verify")
            for g in [20.,-20.];run_case(base,p1,pl,vf,cf,g,"verify");end
        else
            run_case(base,p1,pl,vf,cf,20.,"default")
            run_case(base,p1,pl,vf,cf,20.,"strict")
            run_case(base,p1,pl,vf,cf,-20.,"strict")
        end
    end
end
close(summary)
