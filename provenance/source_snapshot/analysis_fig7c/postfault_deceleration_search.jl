# Independent diagnostic: passive fault, VC=CC1=1, independently varied dispatch/load.
# Original model files are read but never modified. Optional ZIP part matches pre-fault P.
bootstrap=first(split(read(joinpath(@__DIR__,"paper_region_feedforward.jl"),String),"summary=open("))
zip_patch="""
src=replace(src,"Q_load_p = Q_load\\n    end"=>"Q_load_p = Q_load\\n        z_fraction = 0.0\\n        vref2 = 1.0\\n    end")
src=replace(src,"-(P_load_p * terminal.u_r + Q_load_p * terminal.u_i) / V_mag^2"=>"-(P_load_p * terminal.u_r + Q_load_p * terminal.u_i) * ((1-z_fraction)/V_mag^2 + z_fraction/vref2)")
src=replace(src,"-(P_load_p * terminal.u_i - Q_load_p * terminal.u_r) / V_mag^2"=>"-(P_load_p * terminal.u_i - Q_load_p * terminal.u_r) * ((1-z_fraction)/V_mag^2 + z_fraction/vref2)")
"""
bootstrap=replace(bootstrap,"include_string(Main, src,"=>zip_patch*"\ninclude_string(Main, src,")
include_string(Main,bootstrap,joinpath(@__DIR__,"paper_region_feedforward.jl"))
out=joinpath(@__DIR__,"postfault_deceleration_search");mkpath(out)
track=isempty(ARGS) ? "generation" : ARGS[1]
verify=startswith(track,"verify")
for b in [1,2,3,7], k in ["vsrc₊vc₊γ_d","vsrc₊vc₊γ_q","vsrc₊vc₊V_C_d","vsrc₊vc₊V_C_q","vsrc₊cc1₊V_I_d","vsrc₊cc1₊V_I_q"]
    push!(names,"b$(b)_$(k)");push!(ids,VIndex(b,Symbol("inv$(b)₊",k)))
end
if length(ARGS)>=3
    out=abspath(ARGS[3]);mkpath(out)
end
summary=open(joinpath(out,track*"_status.csv"),"w")
println(summary,"case,P1,Pother,Pload,z_fraction,duration,G_fault,method,status,last_t,min_jacdet,transfer_error,prefault_frequency_error,prefault_current_max,post_line_factor,current_limit,fault_switch");flush(summary)
function study_case(base,p1,po,pl,alpha,dur;g=-20.,horizon=15.,linefactor=1.,imax=1.2)
    clear=1+dur
    label="$(track)_P$(p1)_O$(po)_L$(pl)_Z$(alpha)_T$(dur)_G$(g)"
    linefactor!=1. && (label*="_line$(linefactor)")
    imax!=1.2 && (label*="_I$(imax)")
    println("START ",label);flush(stdout)
    u,p,err=transfer(base,1,1);s=NWState(base.f,u,p,0.)
    for b in [1,2,3,7];s.p[VIndex(b,Symbol("inv$(b)₊vsrc₊I_max"))]=imax;end
    s.p[VIndex(4,:fault_load₊G_fault_p)]=g
    s.p[VIndex(4,:fault_load₊t_clear_p)]=clear
    v02=s[VIndex(5,:busbar₊u_r)]^2+s[VIndex(5,:busbar₊u_i)]^2+1e-6
    s.p[VIndex(5,:PQLoad₊z_fraction)]=alpha;s.p[VIndex(5,:PQLoad₊vref2)]=v02
    yy=copy(Y);yy[5,5]+=alpha*pl/v02;zz=inv(yy[[1,2,3,4,5,7,8],[1,2,3,4,5,7,8]])
    yypost=copy(yy);dy=(linefactor-1)*3.43*exp(-im*75pi/180)
    yypost[4,4]+=dy;yypost[6,6]+=dy;yypost[4,6]-=dy;yypost[6,4]-=dy
    zzpost=inv(yypost[[1,2,3,4,5,7,8],[1,2,3,4,5,7,8]])
    pc=(1-alpha)*pl
    function guard(tt,v)
        gt=fault_conductance(tt, g, 1., clear)
        zmat=tt>clear ? zzpost : zz
        z=zmat[5,5]+gt*zmat[5,4]*zmat[4,5]/(1-gt*zmat[4,4]);x=v^2;e=1e-6
        1+2real(z)*pc*e/(x+e)^2+abs2(z)*pc^2*(e-x)/(x+e)^3
    end
    ts=0.;status="finished";minj=Inf;pre_f=0.;pre_i=0.
    open(joinpath(out,label*".csv"),"w") do io
        println(io,join(vcat("t","success",names,"jacdet"),","))
        q=[s[id] for id in ids];println(io,join(vcat(0.,1,q,guard(0.,q[1])),","))
        while ts<horizon
            if linefactor!=1. && ts>=clear
                # Network edge positions are reordered; select the original line by name.
                s.p[EIndex(:line2,:line_sym₊active)]=linefactor
            end
            fine=(.99<=ts<1.04)||(clear-.01<=ts<clear+.04)
            interval=fine ? (verify ? .0001 : .0002) : (ts<clear+.1 && ts>.99 ? .001 : (verify ? .002 : .01))
            te=ts==0. ? 1e-5 : min(horizon,ts+interval)
            if ts<clear<te;te=clear;end
            ps=ODEProblem(base.f,u,(ts,te),p);stops=filter(x->ts<x<=te,[1.,clear])
            sol=solve_fault_switched(ps,verify ? Rodas4P() : Rodas5P();reltol=verify ? 1e-10 : 1e-9,abstol=verify ? 1e-12 : 1e-11,dtmax=verify ? .0002 : (fine ? .0001 : .001),tstops=stops,d_discontinuities=stops)
            tt=sol.t[end];q=sol(tt;continuity=:right,idxs=ids);ok=string(sol.retcode)=="Success";j=guard(tt,q[1])
            println(io,join(vcat(tt,Int(ok),q,j),","));ts=tt
            if !ok;status=string(sol.retcode);break;end
            minj=min(minj,j)
            if j<=1e-6;status="branch_guard";break;end
            if .5<tt<.99
                for b in [1,2,3,7]
                    om=q[findfirst(==("b$(b)_droop₊ω"),names)]
                    a=q[findfirst(==("b$(b)_vsrc₊cc1₊i_f_d"),names)];c=q[findfirst(==("b$(b)_vsrc₊cc1₊i_f_q"),names)]
                    pre_f=max(pre_f,abs(om/(2pi)-60));pre_i=max(pre_i,hypot(a,c))
                end
                if pre_f>.001 || pre_i>imax+.0001;status="invalid_prefault";break;end
            end
            u=copy(sol.u[end])
        end
    end
    println(summary,join([label,p1,po,pl,alpha,dur,g,verify ? "verify" : "strict",status,ts,minj,err,pre_f,pre_i,linefactor,imax,"if"],","));flush(summary)
    println("DONE ",label," ",status," t=",ts);flush(stdout)
end

# CSV manifests allow adaptive follow-up and independent verification without source edits.
if length(ARGS)>1
    lines=readlines(ARGS[2]); specs=[parse.(Float64,split(l,",")) for l in lines[2:end] if !isempty(strip(l))]
else
    pvals=track=="charging" ? [-.3,-.6,-.9] : [.1,.4,.7]
    specs=[[p1,po,pl,a,d,-20.,15.] for p1 in pvals for po in [.2,.5] for pl in [.6,1.2,1.8] for a in [0.] for d in [.2,1.]]
end
lastpoint=nothing;base=nothing
for spec in specs
    p1,po,pl,alpha,dur,g,horizon=spec[1:7]
    linefactor=length(spec)>7 ? spec[8] : 1.
    imax=length(spec)>8 ? spec[9] : 1.2
    point=(p1,po,pl)
    try
        if point!=lastpoint
            init_sim(po,p1,po;P_load=pl,G_fault=0.,t_cle=1+dur,max_t=horizon)
            global base=SIM_STATE["prob"];global lastpoint=point
        end
        study_case(base,p1,po,pl,alpha,dur;g,horizon,linefactor,imax)
    catch e
        println("ERROR ",spec," ",sprint(showerror,e));flush(stdout)
        println(summary,join(["init_or_case_error",p1,po,pl,alpha,dur,g,verify ? "verify" : "strict","error",0.,NaN,NaN,NaN,NaN,linefactor,imax,"if"],","));flush(summary)
        global lastpoint=nothing
    end
end
close(summary)
