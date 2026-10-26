# Finite rerun of selected saved initial states; original results remain immutable.
include(joinpath(@__DIR__,"model.jl"))
using SHA,Printf
const BASIN_INPUT=joinpath(PACK,"results/conditional_basins/data")
const BASIN_OUTPUT=joinpath(PACK,"basin_recalculation")
mkpath(BASIN_OUTPUT);mkpath(joinpath(DATA,"results"))
const BASIN_STABLE=vec(readdlm(joinpath(BASIN_INPUT,"references/section_stable.csv"),',',Float64))
const BASIN_SADDLE=vec(readdlm(joinpath(BASIN_INPUT,"references/section_saddle.csv"),',',Float64))
const BASIN_HELD=setdiff(DYN,vcat(DELTA,PF))
const BASIN_ORIGIN=sum(BASIN_STABLE[DELTA])/3
const BASIN_JOBS=[split(line,',') for line in readlines(joinpath(BASIN_INPUT,"jobs.csv"))[2:end]]
const BASIN_INITIAL=Dict(a[1]=>parse.(Float64,a[2:end]) for a in [split(line,',') for line in readlines(joinpath(BASIN_INPUT,"initial_states.csv"))[2:end]])
selected=isempty(ARGS) ? ["refine_00247","focus_00148"] : ARGS==["--all"] ? [a[1] for a in BASIN_JOBS] : ARGS
@assert all(tag in [a[1] for a in BASIN_JOBS] for tag in selected) "Unknown saved basin case. See results/conditional_basins/data/jobs.csv."
python=get(ENV,"GFM_PYTHON","python")
for source in selected
    a=only(filter(a->a[1]==source,BASIN_JOBS))
    tau=parse(Float64,a[2]);theta=parse(Float64,a[3]);frequency=parse(Float64,a[4]);section=a[5]
    horizon=parse(Float64,a[7]);level=a[8];step=parse(Float64,a[9]);ts=@sprintf("%.2f",tau)
    tag="saved_basin_"*source
    println("START ",tag);flush(stdout)
    p=vec(readdlm(joinpath(BASIN_INPUT,"references/parameters_"*ts*".csv"),',',Float64))
    x=copy(BASIN_INITIAL[source]);initial=copy(x[DYN]);err=algfix!(x,p)
    base=section=="equilibrium" ? BASIN_STABLE : vec(readdlm(joinpath(BASIN_INPUT,"references/clearing_"*ts*".csv"),',',Float64))
    st=NWState(NW,x,p,0.)
    angle=sum(x[DELTA])/3-BASIN_ORIGIN
    freq=sum(st[VIndex(b,Symbol("inv$(b)₊droop₊ω"))]/(2pi)-60 for b in REMOTE)/3
    @assert err<1e-8&&x[DYN]==initial&&maximum(abs.(x[BASIN_HELD]-base[BASIN_HELD]))<1e-12
    @assert abs(angle-theta)<1e-9&&abs(freq-frequency)<1e-9
    for (i,b) in enumerate(ALLG)
        @assert st.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]==tau
        @assert st.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_q"))]==.02
        @assert st.p[VIndex(b,Symbol("inv$(b)₊droop₊Pset"))]==[.5,.4,.4,.4][i]
    end
    for (name,root) in [("stable",BASIN_STABLE),("saddle",BASIN_SADDLE)]
        er,ev,n=eqcheck(root,p);@assert er<1e-8
        @assert name=="stable" ? ev<0 : n==1
        appendrow(joinpath(BASIN_OUTPUT,"equilibria.csv"),"case,root,residual,max_real_eigen,unstable_modes",[source,name,er,ev,n])
    end
    appendrow(joinpath(BASIN_OUTPUT,"initialization.csv"),"case,algebraic_residual,theta_error_rad,df_error_Hz,held_change",[source,err,angle-theta,freq-frequency,maximum(abs.(x[BASIN_HELD]-base[BASIN_HELD]))])
    codepaths=vcat([joinpath(@__DIR__,f) for f in ["run_basin.jl","model.jl","engine12.jl"]],sort([joinpath(root,f) for (root,dirs,files) in walkdir(SNAP) for f in files if endswith(f,".jl")||endswith(f,".csv")]))
    fingerprint=bytes2hex(sha256(join(a,",")*join(string.(BASIN_INITIAL[source]))*join(string.(p))*join([bytes2hex(sha256(read(f))) for f in codepaths])))
    directory=joinpath(CHECK,tag);mkpath(directory);signature=joinpath(directory,"configuration.sha256")
    if isfile(signature)
        @assert strip(read(signature,String))==fingerprint "Changed saved basin input or code: use a fresh checkout/workspace."
    else
        @assert isempty(readdir(directory))
        write(signature,fingerprint*"\n")
    end
    last=0.
    for stop in sort(unique(vcat([30.,100.,300.],collect(600.:300.:horizon),[horizon])))
        stop>horizon&&continue
        x,rc=segment(tag,"post"*lpad(Int(stop),3,'0'),x,p,last,stop,0.;level=level,step=step)
        @assert rc=="Success";last=stop
    end
    run(`$python $(joinpath(@__DIR__,"classify_basin_run.py")) $source $tag`)
    println("DONE ",tag);flush(stdout)
end
