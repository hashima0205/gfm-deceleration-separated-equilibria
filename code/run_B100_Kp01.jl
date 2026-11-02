# Gain comparison: common-base mp=1% vs 2%, fixed B/k, one model build.
const KP01_BASE=joinpath(@__DIR__,"run_capacitor_k_strength.jl")
const KP01_GRID=joinpath(@__DIR__,"run_B100_tau_comparison.jl")
const GAIN_MULTIPLIER=Ref(.5)
prefix=first(split(read(KP01_BASE,String),"\nfunction kmain()"))
marker="    for (j,y) in [(2,3.43),(3,.68)]"
@assert occursin(marker,prefix)
gainpatch=raw"""    for b in ALLG
        s.p[VIndex(b,Symbol("inv$(b)₊droop₊Kp"))]*=GAIN_MULTIPLIER[]
    end
"""
prefix=replace(prefix,marker=>gainpatch*marker,"const KB=.5561806200995012"=>"const KB=1.0",
    "basin_recalculation/capacitor_k_strength\\\""=>"basin_recalculation/B100_Kp01_tau_comparison\\\"",
    "for stop in [30.,60.,150.]"=>"for stop in [30.,60.,150.,300.]")
include_string(Main,prefix,KP01_BASE)
@assert KB==1. && endswith(OUTP,"B100_Kp01_tau_comparison")
const BASE_KP_VALUES=copy(KP010)
function usegain!(factor)
    GAIN_MULTIPLIER[]=factor;KP010 .= BASE_KP_VALUES .* factor
end
# Reuse the exact original wide/local grids and basin integration method.
gridtext=read(KP01_GRID,String)
a=first(findfirst("const BK=",gridtext));b=first(findfirst("\nfunction bmain()",gridtext))
include_string(Main,gridtext[a:b-1],KP01_GRID)
function kp01main()
    mkpath(OUTP)
    seed=joinpath(PACK,"basin_recalculation/B100_tau_comparison")
    xp=vec(readdlm(joinpath(seed,"pre.csv"),',',Float64));xs=vec(readdlm(joinpath(seed,"stable.csv"),',',Float64));xu=vec(readdlm(joinpath(seed,"saddle.csv"),',',Float64))
    files=[@__FILE__,KP01_BASE,KP01_GRID,K_SOURCE,joinpath(@__DIR__,"adaptive_probe_metrics.jl"),joinpath(@__DIR__,"model.jl"),joinpath(@__DIR__,"engine12.jl")]
    append!(files,sort([joinpath(r,f) for (r,ds,fs) in walkdir(SNAP) for f in fs if endswith(f,".jl") || endswith(f,".csv")]))
    signature=bytes2hex(sha256(string(VERSION,xp,xs,xu,BK,KB,BTAUS,BGLOBAL_THETA,BGLOBAL_FREQ,BLOCAL_ALPHA,BLOCAL_FREQ,BHORIZONS)*join(bytes2hex(sha256(read(f))) for f in files)))
    sig=joinpath(OUTP,"configuration.sha256")
    if isfile(sig);@assert strip(read(sig,String))==signature
    else;@assert isempty(readdir(OUTP));write(sig,signature*"\n");end
    tick=time();eq=[];events=[];machines=[];checks=[]
    for (name,x) in [("pre",xp),("stable",xs),("saddle",xu)];savevec(joinpath(OUTP,name*".csv"),x);end
    for factor in [1.,.5]
        usegain!(factor)
        for tau in BTAUS
            for (name,x,base) in [("pre",xp,BEFOREP),("stable",xs,POSTP),("saddle",xu,POSTP)]
                p=kparams(base,BK,KB,tau);s=NWState(NW,copy(x),p,0.)
                for (j,bus) in enumerate(ALLG)
                    @assert abs(s.p[VIndex(bus,Symbol("inv$(bus)₊droop₊Kp"))]-KP010[j])<1e-12
                end
                er,ev,n=eqcheck(x,p);@assert er<1e-8;@assert name=="saddle" ? n>=1 : ev<0
                push!(eq,[factor,.02factor,tau,name,er,ev,n,busvp(x,5),phasep(x)])
            end
            tablep("equilibria.csv","gain_factor,mp_common,taup,root,residual,max_real_eigen,unstable_modes,V5_pu,mean_phase_rad",eq)
            tag=@sprintf("B100_Kp%02d_tau%03d",round(Int,2factor),round(Int,100tau))
            label,turns,ms=kphysical(tag,xp,xs,BK,tau)
            if factor==1.
                # Verify the replay against the already computed mp=2% case.
                p=kparams(POSTP,BK,KB,tau);result=deserialize(joinpath(OUTP,tag*"_post150.bin"))
                oldtag=@sprintf("B100_tau%03d",round(Int,100tau))
                expected=vec(readdlm(joinpath(seed,oldtag*"_end.csv"),',',Float64))
                er=maximum(abs.(result.u[end]-expected));@assert er<2e-7
                push!(checks,[tau,er])
            end
            push!(events,[tag,factor,.02factor,BASE_KP_VALUES[1]*factor,BK,KB,tau,300.,busvp(xp,5),busvp(xs,5),label,turns])
            for m in ms;push!(machines,vcat([tag],collect(m)));end
            tablep("events.csv","case,gain_factor,mp_common,Kp_rad_s_pu,k,B_pu,taup,horizon_s,pre_V5_pu,post_eq_V5_pu,classification,remote_turns",events)
            tablep("event_machines.csv","case,bus,classification,turns,tail_max_abs_df_Hz,phase_error_rad,tail_phase_width_rad,min_window_turns,max_window_turns",machines)
            tablep("baseline_validation.csv","taup,post150_max_full_state_difference",checks)
        end
    end
    usegain!(.5)
    points=bgrid(xs,xu)
    tablep("grid.csv","point,grid,alpha,df_offset_Hz",[[i,g,a,f] for (i,(a,f,g)) in enumerate(points)])
    for tau in BTAUS;bbasin(xs,xu,tau,points);end
    println("DONE mp=1% comparison; batch_wall_s=",time()-tick);flush(stdout)
end
kp01main()
