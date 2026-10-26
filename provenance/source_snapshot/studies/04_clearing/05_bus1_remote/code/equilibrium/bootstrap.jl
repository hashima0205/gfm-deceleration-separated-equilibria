using Serialization
const EQROOT=normpath(joinpath(@__DIR__,"../../../../.."))
const EQOUT=joinpath(EQROOT,"studies/04_clearing/05_bus1_remote/output/equilibrium")
source=joinpath(EQROOT,"analysis_fig7c/postfault_deceleration_search.jl")
bootstrap=first(split(read(source,String),"\nsummary=open("))
include_string(Main,bootstrap,source)
function prepare_eq(pl)
    init_sim(.1,.7,.1;P_load=pl,G_fault=0.,t_cle=1.2,max_t=120.)
    base=SIM_STATE["prob"];u,p,err=transfer(base,1,1);s=NWState(base.f,u,p,0.)
    v02=s[VIndex(5,:busbar₊u_r)]^2+s[VIndex(5,:busbar₊u_i)]^2+1e-6
    s.p[VIndex(5,:PQLoad₊z_fraction)]=1.
    s.p[VIndex(5,:PQLoad₊vref2)]=v02
    s.p[VIndex(4,:fault_load₊G_fault_p)]=0.
    s.p[VIndex(4,:fault_load₊fault_enabled_p)]=0.
    return base,u,p,s,v02
end
for pl in [.8,.9,1.]
    base,u,p,s,v02=prepare_eq(pl)
    open(joinpath(EQOUT,"initial_L$(pl).csv"),"w") do io
        println(io,"bus,Pset,Qset,Vset,Kq,Kp,Rg,Lg,C,Vcr,Vci,ur,ui,load_conductance")
        for b in [1,2,3,7]
            values=[s.p[VIndex(b,Symbol("inv$(b)₊",k))] for k in ["droop₊Pset","droop₊Qset","droop₊Vset","droop₊Kq","droop₊Kp","vsrc₊Rg","vsrc₊Lg","vsrc₊C"]]
            append!(values,[s[VIndex(b,Symbol("inv$(b)₊vsrc₊filter₊",k))] for k in ["V_C_r","V_C_i"]])
            append!(values,[s[VIndex(b,Symbol("busbar₊",k))] for k in ["u_r","u_i"]])
            println(io,join(vcat(b,values,pl/v02),","))
        end
    end
end
base,u,p,s,v02=prepare_eq(1.)
s.p[EIndex(:line7,:line_sym₊active)]=0.
println("EQUILIBRIUM_READY states=",length(u));flush(stdout)
