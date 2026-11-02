# One Julia process; classify at 30/60/150 s and extend only unresolved points.
include(joinpath(@__DIR__, "model.jl"))
using SHA, Serialization, Printf

const OUT010 = joinpath(PACK, "basin_recalculation", "tau010_adaptive")
mkpath(OUT010)
const ST010 = vec(readdlm(joinpath(DATA,"seeds/tau010_T015_h150_stable.csv"), ',', Float64))
const SA010 = vec(readdlm(joinpath(DATA,"seeds/tau010_T015_h150_saddle.csv"), ',', Float64))
const PA010 = vec(readdlm(joinpath(CHECK,"tau010_T015_h150/post150_parameters.csv"), ',', Float64))
const XGRID010 = [-.4,-.15,0.,.4,.8,1.2,1.6]
const YGRID010 = [-.15,-.10,-.05,0.,.05]
const HORIZONS010 = [30.,60.,150.]
const INDEX010 = NWState(NW, Float64.(1:length(ST010)), copy(PA010), 0.)
const DI010 = [Int(INDEX010[VIndex(b,Symbol("inv$(b)₊droop₊δ"))]) for b in ALLG]
const PI010 = [Int(INDEX010[VIndex(b,Symbol("inv$(b)₊droop₊Pfilt"))]) for b in ALLG]
const STATE010 = NWState(NW, copy(ST010), copy(PA010), 0.)
const KP010 = [STATE010.p[VIndex(b,Symbol("inv$(b)₊droop₊Kp"))] for b in ALLG]
const PS010 = [STATE010.p[VIndex(b,Symbol("inv$(b)₊droop₊Pset"))] for b in ALLG]
const HELD010 = setdiff(DYN, vcat(DELTA,PF))
const ORIGIN010 = sum(ST010[DELTA])/3
const GAP010 = abs.(atan.(sin.(SA010[DI010]-ST010[DI010]), cos.(SA010[DI010]-ST010[DI010])))
const TOL010 = [gap>1e-8 ? min(.002,.1*gap) : .002 for gap in GAP010]

function table010(file, header, rows)
    tmp=joinpath(OUT010,file*".tmp")
    open(tmp,"w") do io
        println(io,header)
        for row in rows; println(io,join(row,',')); end
    end
    mv(tmp,joinpath(OUT010,file);force=true)
end

function interp010(t,a,q)
    j=clamp(searchsortedlast(t,q),1,length(t)-1)
    a[j]+(a[j+1]-a[j])*(q-t[j])/(t[j+1]-t[j])
end

function classify010(t,u)
    finish=t[end];tail=findall(t .>= finish-min(50.,finish/3))
    machines=[]
    for j in eachindex(ALLG)
        angle=[x[DI010[j]] for x in u]
        freq=[KP010[j]*(PS010[j]-x[PI010[j]])/(2pi) for x in u]
        shift=angle .- angle[1];turns=shift[end]/(2pi)
        error=abs(atan(sin(angle[end]-ST010[DI010[j]]),cos(angle[end]-ST010[DI010[j]])))
        width=maximum(angle[tail])-minimum(angle[tail])
        windows=diff([interp010(t,angle,q) for q in range(.4*finish,finish;length=4)])/(2pi)
        dfmax=maximum(abs.(freq[tail]))
        label=if dfmax<1e-4 && error<TOL010[j] && width<.005
            maximum(abs.(shift))>=2pi ? "recovery_after_slip" : "recovery"
        elseif turns< -2 && maximum(windows)< -.2 && sum(windows)/(.6*finish)< -1e-4
            "negative_running"
        elseif turns>2 && minimum(windows)>.2 && sum(windows)/(.6*finish)>1e-4
            "positive_running"
        else
            "unresolved"
        end
        push!(machines,(ALLG[j],label,turns,dfmax,error,width,minimum(windows),maximum(windows)))
    end
    labels=[m[2] for m in machines]
    label=if all(l in ["recovery","recovery_after_slip"] for l in labels)
        "recovery_after_slip" in labels ? "recovery_after_slip" : "recovery"
    elseif "negative_running" in labels && "positive_running" in labels
        "other"
    elseif "negative_running" in labels
        "negative_running"
    elseif "positive_running" in labels
        "positive_running"
    else
        "unresolved"
    end
    turns=sum(u[end][DELTA]-u[1][DELTA])/(6pi)
    label,turns,machines
end

function main010()
    started=time()
    sources=[@__FILE__,joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(sources,sort([joinpath(root,f) for (root,dirs,files) in walkdir(SNAP) for f in files if endswith(f,".jl") || endswith(f,".csv")]))
    fingerprint=bytes2hex(sha256(string(VERSION,PA010,ST010,SA010,XGRID010,YGRID010,HORIZONS010)*join([bytes2hex(sha256(read(f))) for f in sources])))
    sig=joinpath(OUT010,"configuration.sha256")
    if isfile(sig)
        @assert strip(read(sig,String))==fingerprint "Inputs/code changed: preserve the old output and use a new output directory."
    else
        @assert isempty(readdir(OUT010))
        write(sig,fingerprint*"\n")
    end
    @assert length(DYN)==52 && length(ALG)==14 && length(HELD010)==46
    @assert PS010==[.5,.4,.4,.4]
    parameter_rows=[]
    for (j,b) in enumerate(ALLG)
        tp=STATE010.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]
        tq=STATE010.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_q"))]
        vf=STATE010.p[VIndex(b,Symbol("inv$(b)₊vsrc₊VC_F"))]
        cf=STATE010.p[VIndex(b,Symbol("inv$(b)₊vsrc₊CC1_F"))]
        @assert tp==.1 && tq==.02 && vf==cf==1
        push!(parameter_rows,[b,tp,tq,PS010[j],KP010[j],vf,cf])
    end
    table010("parameters.csv","bus,taup,tauq,Pset,Kp,VC_F,CC1_F",parameter_rows)
    rootrows=[]
    for (name,x) in [("stable",ST010),("saddle",SA010)]
        er,ev,n=eqcheck(x,PA010)
        @assert er<1e-8 && (name=="stable" ? ev<0 : n==1)
        push!(rootrows,[name,er,ev,n])
    end
    table010("equilibria.csv","root,residual,max_real_eigen,unstable_modes",rootrows)
    savevec(joinpath(OUT010,"post_parameters.csv"),PA010)
    rows=[];history=[];audit=[];initials=[];machine_rows=[]
    for (i,(theta,frequency)) in enumerate(Iterators.product(XGRID010,YGRID010))
        # Julia product visits theta first; retain this order in the triangle table.
        tag=@sprintf("tau010_grid_%03d",i)
        pointstart=time();x=copy(ST010)
        x[DELTA].+=theta
        x[PF].-=frequency*2pi ./ KP010[2:4]
        before=copy(x[DYN]);er=algfix!(x,PA010)
        s=NWState(NW,x,copy(PA010),0.)
        actual_theta=sum(x[DELTA])/3-ORIGIN010
        actual_df=sum(s[VIndex(b,Symbol("inv$(b)₊droop₊ω"))]/(2pi)-60 for b in REMOTE)/3
        @assert er<1e-8 && x[DYN]==before && x[HELD010]==ST010[HELD010]
        @assert abs(actual_theta-theta)<1e-9 && abs(actual_df-frequency)<1e-9
        push!(initials,vcat([tag],x))
        push!(audit,[tag,er,actual_theta-theta,actual_df-frequency,0.])
        t=Float64[0.];u=[copy(x)];last=0.;label="unresolved";turns=0.;machines=[]
        for stop in HORIZONS010
            checkpoint=joinpath(OUT010,tag*"_h"*string(Int(stop))*".bin")
            runtime=0.;retcode="Success"
            if isfile(checkpoint)
                tt,uu,runtime=deserialize(checkpoint)
            else
                before=copy(x[DYN]);er=algfix!(x,PA010);@assert er<1e-8 && x[DYN]==before
                samples=sort(unique(vcat(collect(last:.05:stop),collect(last:.0001:min(stop,last+.02)),[last,stop])))
                tick=time()
                sol=solve(ODEProblem(NW,copy(x),(last,stop),copy(PA010)),Rodas5P();
                    reltol=1e-9,abstol=1e-11,dtmax=.005,maxiters=10^8,
                    saveat=samples,save_everystep=false,dense=false)
                runtime=time()-tick;retcode=string(sol.retcode)
                tt=collect(sol.t);uu=[copy(v) for v in sol.u]
                if retcode=="Success"
                    serialize(checkpoint*".tmp",(tt,uu,runtime));mv(checkpoint*".tmp",checkpoint;force=true)
                end
            end
            for j in eachindex(tt)
                tt[j]>t[end] || continue
                push!(t,tt[j]);push!(u,uu[j])
            end
            x=copy(u[end]);last=t[end]
            label,turns,machines=classify010(t,u)
            retcode=="Success" || (label="numerical_failure")
            push!(history,[tag,theta,frequency,stop,last,label,turns,runtime,retcode])
            println("POINT ",i,"/35 ",tag," t=",last," ",label," solve_s=",round(runtime;digits=3));flush(stdout)
            label=="unresolved" || break
        end
        f=zeros(length(x));maxres=0.
        for (tt,xx) in zip(t,u)
            NW(f,xx,PA010,tt);maxres=max(maxres,maximum(abs.(f[ALG])))
        end
        @assert maxres<1e-7
        push!(rows,[tag,.1,theta,frequency,label,last,turns,maxres,time()-pointstart])
        for m in machines;push!(machine_rows,vcat([tag,last],collect(m)));end
        savevec(joinpath(OUT010,tag*"_end.csv"),x)
        table010(tag*"_trajectory.csv","t_s,theta_relative_rad,df_mean_Hz",
            ([tt,sum(xx[DELTA])/3-ORIGIN010,sum(KP010[2:4].*(PS010[2:4]-xx[PI010[2:4]]))/(6pi)] for (tt,xx) in zip(t,u)))
        table010("plot_samples_tau0.10.csv","case,taup,theta_relative_rad,df_Hz,classification,completed_horizon_s,remote_internal_turns,max_algebraic_residual,wall_s",rows)
        table010("history.csv","case,theta_relative_rad,df_Hz,requested_horizon_s,completed_horizon_s,classification,remote_internal_turns,solve_s,retcode",history)
        table010("initial_states.csv",join(vcat(["case"],["x$j" for j in 1:length(x)]),','),initials)
        table010("initialization.csv","case,algebraic_residual,theta_error_rad,df_error_Hz,held_change",audit)
        table010("machines.csv","case,horizon_s,bus,classification,internal_turns,tail_df_max_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",machine_rows)
    end
    faces=[];nx=length(XGRID010)
    for j in 0:length(YGRID010)-2, i in 0:nx-2
        a=j*nx+i;b=a+1;c=a+nx+1;d=a+nx
        push!(faces,[a,b,c]);push!(faces,[a,c,d])
    end
    table010("plot_triangles_tau0.10.csv","vertex0,vertex1,vertex2",faces)
    println("DONE 35 points; batch_wall_s=",time()-started);flush(stdout)
end
main010()
