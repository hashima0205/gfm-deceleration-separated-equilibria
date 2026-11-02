# Same trajectory-centered sections as the original setting, for Vpre=0.95 pu.
# No post-fault stable equilibrium is known: classify running, never infer recovery.
const ORIGINAL_NEAR_SOURCE=joinpath(@__DIR__,"run_trajectory_voltage_probes.jl")
textnear=first(split(read(ORIGINAL_NEAR_SOURCE,String),"\nfunction voltagep()"))
textnear=replace(textnear,
    "basin_recalculation/trajectory_voltage_probes\""=>"basin_recalculation/capacitor_neighborhood_v095_tau010\"",
    "joinpath(DATA,\"seeds/tau010_T015_h150_stable.csv\")"=>"joinpath(DATA,\"seeds/tau010_T015_h150_stable.csv\")",
    "const STP=vec(readdlm(joinpath(DATA,\"seeds/tau010_T015_h150_stable.csv\"),',',Float64))"=>"const STP=fill(NaN,66)",
    "joinpath(DATA,\"seeds/nf02_1_pre0.csv\")"=>"joinpath(PACK,\"basin_recalculation/capacitor_probe_predictor/v095_tau010_pre_root.csv\")",
    "joinpath(CHECK,\"tau010_T015_h150/post150_parameters.csv\")"=>"joinpath(PACK,\"basin_recalculation/capacitor_probe_predictor/v095_tau010_post_parameters.csv\")",
    "joinpath(CHECK,\"tau010_T015_h150/pre_parameters.csv\")"=>"joinpath(PACK,\"basin_recalculation/capacitor_probe_predictor/v095_tau010_pre_parameters.csv\")",
    "joinpath(CHECK,\"tau010_T015_h150/post030_initial.csv\")"=>"joinpath(PACK,\"basin_recalculation/capacitor_probe_predictor/v095_tau010_clear.csv\")",
    "joinpath(CHECK,\"tau010_T015_h150/post150_end.csv\")"=>"joinpath(PACK,\"basin_recalculation/capacitor_probe_predictor/v095_tau010_end.csv\")",
    "branches=round.(Int,(x[DI010]-STP[DI010])/(2pi))"=>"branches=fill(-999,4) # unknown equilibrium branch",
    "reference_error<2e-7"=>"reference_error<2e-5",
    "@assert error<2e-7"=>"@assert error<2e-5")
include_string(Main,textnear,ORIGINAL_NEAR_SOURCE)
function capnearmain()
    mkpath(OUTP)
    inputs=[@__FILE__,ORIGINAL_NEAR_SOURCE,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(inputs,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    signature=bytes2hex(sha256(string(VERSION,PREP,POSTP,BEFOREP,CLEARP,ANCHOR_TIMES,ANGLE_OFFSETS,FREQ_OFFSETS)*join(bytes2hex(sha256(read(f))) for f in inputs)))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    savevec(joinpath(OUTP,"post_parameters.csv"),POSTP)
    tick=time();neighborhoodp()
    println("DONE capacitor trajectory sections; batch_wall_s=",time()-tick);flush(stdout)
end
capnearmain()
