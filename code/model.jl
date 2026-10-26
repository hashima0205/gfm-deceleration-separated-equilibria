using DelimitedFiles,LinearAlgebra
const PACK=normpath(joinpath(@__DIR__,".."));const SNAP=joinpath(PACK,"provenance/source_snapshot")
const DATA=joinpath(PACK,"data");const CHECK=joinpath(PACK,"checkpoints");const MODE="scr"
ENV["BASIN_DEFINITIONS_ONLY"]="1"
println("BOOT Julia ",VERSION);flush(stdout)
include(joinpath(SNAP,"studies/04_clearing/08_state_space/code/run.jl"))
const NW=BB.f;const ALLG=[1,2,3,7];const DYN=findall(!iszero,MASS)
const STABLE=vec(readdlm(joinpath(DATA,"references/stable.csv"),',',Float64))
const UNSTABLE=vec(readdlm(joinpath(DATA,"references/saddle.csv"),',',Float64))
const PARAM=vec(readdlm(joinpath(DATA,"references/post_parameters.csv"),',',Float64))
include(joinpath(@__DIR__,"engine12.jl"))
const THETA0=sum(STABLE[DELTA])/3
for x in [STABLE,UNSTABLE]
    e,v,n=eqcheck(x,PARAM);@assert e<1e-8;println("ROOT ",e," ",v," ",n)
end
s=NWState(NW,copy(STABLE),PARAM,0.)
for b in ALLG
    @assert s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_p"))]==1.
    @assert s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_q"))]==.02
end
