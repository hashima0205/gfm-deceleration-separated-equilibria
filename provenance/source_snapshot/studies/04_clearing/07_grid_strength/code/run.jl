# Vary the physical external-grid impedance, leaving all setpoints unchanged.
src=first(split(read(joinpath(@__DIR__,"../../06_equilibrium_basin/code/scan.jl"),String),"\ntrack=length(ARGS)"))
# HERE/OUT resolve into this scenario; the shared bootstrap remains next door.
src=replace(src,"\"vsrc₊filter₊i_f_r\",\"vsrc₊filter₊i_f_i\""=>"\"vsrc₊filter₊i_f_r\",\"vsrc₊filter₊i_f_i\",\"vsrc₊filter₊V_C_r\",\"vsrc₊filter₊V_C_i\",\"vsrc₊filter₊i_g_r\",\"vsrc₊filter₊i_g_i\"")
include_string(Main,src,@__FILE__)
const branches=[(1,4,2.14),(4,6,3.43),(5,6,.68),(3,5,2.21),(2,3,1.91),(2,5,2.14),(4,5,1.8),(2,8,2.17),(2,7,6.39),(7,8,1.69),(5,8,8.96)]
baseop=operating(.85)
function gridparams(bb,u,p,k,angle,mode)
    q=copy(p);s=NWState(bb.f,copy(u),q,0.)
    for (j,(a,b,y)) in enumerate(branches)
        selected=mode=="all" || j in [2,3]
        z=(selected ? k : 1.)/y;th=deg2rad(selected ? angle : 75.)
        s.p[EIndex(Symbol("line",j),:line_sym₊R)]=z*cos(th)
        s.p[EIndex(Symbol("line",j),:line_sym₊X)]=z*sin(th)
    end
    q
end
function continuation_rootpath(bb,start,p,k,angle,mode)
    x=copy(start);err=0.;q=copy(p)
    n=max(1,ceil(Int,abs(k-1)/.05),ceil(Int,abs(angle-75)/2))
    for w in range(0,1,length=n+1)[2:end]
        q=gridparams(bb,x,p,1+w*(k-1),75+w*(angle-75),mode)
        xx,err=stationary(bb.f,x,q)
        if err>=1e-8
            # Do not mistake a failed continuation step for nonexistence.
            qt=gridparams(bb,start,p,k,angle,mode);xt,et=stationary(bb.f,start,qt)
            return xt,qt,et
        end
        x=xx
    end
    x,q,err
end
function maxcurrent(bb,x,p)
    s=NWState(bb.f,copy(x),p,0.)
    maximum(hypot(s[VIndex(b,Symbol("inv$(b)₊vsrc₊filter₊i_f_r"))],s[VIndex(b,Symbol("inv$(b)₊vsrc₊filter₊i_f_i"))]) for b in [1,2,3,7])
end
function export_settings()
    bb,pre,p,post,pp,er=baseop;s=NWState(bb.f,copy(pre),p,0.)
    open(joinpath(OUT,"fixed_settings.csv"),"w") do io
        println(io,"bus,Pset,Qset,Vset,Kq,Kp,Rg,Lg,C,Vcr,Vci,ur,ui,load_conductance")
        conductance=s.p[VIndex(5,:PQLoad₊P_load_p)]/s.p[VIndex(5,:PQLoad₊vref2)]
        for b in [1,2,3,7]
            values=[s.p[VIndex(b,Symbol("inv$(b)₊",key))] for key in ["droop₊Pset","droop₊Qset","droop₊Vset","droop₊Kq","droop₊Kp","vsrc₊Rg","vsrc₊Lg","vsrc₊C"]]
            append!(values,[s[VIndex(b,Symbol("inv$(b)₊vsrc₊filter₊",key))] for key in ["V_C_r","V_C_i"]])
            append!(values,[s[VIndex(b,Symbol("busbar₊",key))] for key in ["u_r","u_i"]])
            println(io,join(vcat(b,values,conductance),","))
        end
    end
end
function batch(track,manifest)
    cache=Dict{Any,Any}();bb,pre0,p0,post0,pp0,er=baseop
    open(joinpath(OUT,track*"_status.csv"),"w") do io
        println(io,"case,mode,k,angle,tau,duration,horizon,pre_residual,pre_eigen,pre_current,post_residual,post_eigen,post_current,status,last_t");flush(io)
        for (i,line) in enumerate(readlines(manifest)[2:end])
            isempty(line)&&continue
            a=split(line,",");mode=String(a[1]);k,angle,tau,dur,horizon=parse.(Float64,a[2:end])
            tag="$(track)_$(lpad(i,3,'0'))";println("START ",tag," ",line);flush(stdout)
            try
                pre,p,pe,post,pp,qe=get!(cache,(mode,k,angle)) do
                    x,pr,e=rootpath(bb,pre0,p0,k,angle,mode)
                    y,po,f=rootpath(bb,post0,pp0,k,angle,mode)
                    (x,pr,e,y,po,f)
                end
                p=copy(p);pp=copy(pp)
                for (x,q) in [(pre,p),(post,pp)]
                    s=NWState(bb.f,copy(x),q,0.)
                    for b in [2,3,7];s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]=tau;end
                end
                evpre=pe<1e-8 ? stability(bb,pre,p) : NaN
                evpost=qe<1e-8 ? stability(bb,post,pp) : NaN
                ipre=pe<1e-8 ? maxcurrent(bb,pre,p) : NaN
                ipost=qe<1e-8 ? maxcurrent(bb,post,pp) : NaN
                if pe>=1e-8
                    status,last="prefault_root_unresolved",0.
                elseif evpre>=0 || ipre>1.2+1e-6
                    status,last="prefault_invalid",0.
                else
                    status,last=simulate(bb,pre,p,.85,tau,dur,horizon,tag;verify=startswith(track,"verify"))
                end
                println(io,join([tag,mode,k,angle,tau,dur,horizon,pe,evpre,ipre,qe,evpost,ipost,status,last],","));flush(io)
                if qe<1e-8
                    open(joinpath(OUT,"runs",tag*"_equilibrium.csv"),"w") do eqio
                        println(eqio,"state_index,pre,post")
                        for j in eachindex(pre);println(eqio,join([j,pre[j],post[j]],","));end
                    end
                end
                println("DONE ",tag," ",status," pre=",evpre," post=",evpost," post_residual=",qe);flush(stdout)
            catch e
                println(io,join([tag,mode,k,angle,tau,dur,horizon,NaN,NaN,NaN,NaN,NaN,NaN,"error",0],","));flush(io)
                println("ERROR ",tag," ",sprint(showerror,e));flush(stdout)
            end
        end
    end
end
include("seed_support.jl")
export_settings()
track=length(ARGS)>0 ? ARGS[1] : "coarse"
manifest=length(ARGS)>1 ? ARGS[2] : joinpath(@__DIR__,"coarse.csv")
batch(track,manifest)
