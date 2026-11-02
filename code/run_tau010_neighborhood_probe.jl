# Run only the neighborhood stage, preserving its original integration code/cache signature.
source_file=joinpath(@__DIR__,"run_trajectory_voltage_probes.jl")
source_text=read(source_file,String)
combined="tick=time();neighborhoodp();voltagep()"
@assert occursin(combined,source_text)
source_text=replace(source_text,combined=>"tick=time();neighborhoodp()";count=1)
include_string(Main,source_text,source_file)
