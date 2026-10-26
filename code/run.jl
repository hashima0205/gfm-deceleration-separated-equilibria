include(joinpath(@__DIR__,"model.jl"))
using SHA
for name in ["seeds","results"];mkpath(joinpath(DATA,name));end
mkpath(CHECK)
const RUN_FAILURES=Ref(0)
function params(load,ps,taup,tauq)
 tuple((begin
  p=gridparams(BB,PRE,p,1.0299,75.,"tie");s=NWState(NW,copy(PRE),p,0.);s.p[VIndex(5,:PQLoad₊P_load_p)]*=load
  for (i,b) in enumerate(ALLG)
   s.p[VIndex(b,Symbol("inv$(b)₊droop₊Pset"))]=ps[i];s.p[VIndex(b,Symbol("inv$(b)₊droop₊τ_q"))]=tauq
  end;p
 end for p in stages(taup,"main"))...)
end
manifest=readlines(joinpath(DATA,isempty(ARGS) ? "representative.csv" : ARGS[1]))
@assert manifest[1]=="case,op,load_scale,P1,P2,P3,P7,taup,tauq,T,horizon,level"
for row in manifest[2:end]
 a=split(row,',');tag,op=a[1:2];g=parse(Float64,a[3]);ps=parse.(Float64,a[4:7]);taup=parse(Float64,a[8]);tauq=parse(Float64,a[9]);T=parse(Float64,a[10]);horizon=parse(Float64,a[11]);level=a[12]
 @assert occursin(r"^[A-Za-z0-9_.-]+$",tag)&&!occursin("..",tag)
 @assert op=="nf02_1"&&g==2.461761102780045&&ps==[.5,.4,.4,.4] "New load/dispatch needs matching equilibrium seeds."
 @assert isfinite(taup)&&taup>0&&tauq==.02&&isfinite(T)&&T>0&&isfinite(horizon)&&horizon>0&&isinteger(horizon)&&level in ["A","B"]
 # Refuse reuse after a manifest, runner, numerical engine, or model change.
 sources=vcat([joinpath(@__DIR__,f) for f in ["run.jl","model.jl","engine12.jl"]],
              sort([joinpath(root,f) for (root,dirs,files) in walkdir(SNAP) for f in files if endswith(f,".jl")||endswith(f,".csv")]))
 signature=bytes2hex(sha256(row*join([bytes2hex(sha256(read(f))) for f in sources])))
 directory=joinpath(CHECK,tag);mkpath(directory);signature_file=joinpath(directory,"configuration.sha256")
 if isfile(signature_file)
  @assert strip(read(signature_file,String))==signature "Changed condition: use a new case name."
 else
  @assert isempty(readdir(directory)) "Existing checkpoints lack a configuration signature."
  write(signature_file,signature*"\n")
 end
 println("START ",tag);flush(stdout)
 try
  pre,on,post=params(g,ps,taup,tauq);x=nothing;target=nothing;saddle=nothing
  for (stage,p) in [("pre",pre),("post",post)]
   for j in 0:20
    file=joinpath(DATA,"seeds",op*"_"*stage*string(j)*".csv");isfile(file)||continue
    y=vec(readdlm(file,',',Float64));er,ev,n=eqcheck(y,p)
    appendrow(joinpath(DATA,"equilibria.csv"),"case,op,stage,branch,residual,max_real_eigen,unstable_modes",[tag,op,stage,j,er,ev,n]);@assert er<1e-8
    if ev<0
     stage=="pre" ? (x=y) : (target=y)
    elseif stage=="post"&&n==1;saddle=y;end
   end
  end
  if x===nothing||target===nothing
   RUN_FAILURES[]+=1
   appendrow(joinpath(DATA,"cases.csv"),"case,op,load_scale,P1,P2,P3,P7,taup,tauq,T,horizon,level,retcode",vcat([tag,op,g],ps,[taup,tauq,T,0.,level,x===nothing ? "prefault_stable_root_unresolved" : "postfault_stable_root_unresolved"]));continue
  end
  savevec(joinpath(DATA,"seeds",tag*"_stable.csv"),target);saddle!==nothing&&savevec(joinpath(DATA,"seeds",tag*"_saddle.csv"),saddle)
  for (stage,p) in [("pre",pre),("fault",on),("post",post)]
   st=NWState(NW,copy(x),p,0.)
   for (i,b) in enumerate(ALLG)
    vals=[st.p[VIndex(b,Symbol("inv$(b)₊",k))] for k in ["droop₊τ_p","droop₊τ_q","droop₊Pset","droop₊Qset","vsrc₊VC_F","vsrc₊CC1_F"]];@assert vals[1]==taup&&vals[2]==tauq&&vals[3]==ps[i]&&vals[3]>0&&vals[5]==vals[6]==1
    appendrow(joinpath(DATA,"parameter_audit.csv"),"case,stage,bus,taup,tauq,Pset,Qset,VC_F,CC1_F",vcat([tag,stage,b],vals))
   end
  end
  clear=1+T;x,rc=segment(tag,"pre",x,pre,0.,1.,clear;level=level,step=.01);@assert rc=="Success"
  x,rc=segment(tag,"fault",x,on,1.,clear,clear;level=level,step=.0005);@assert rc=="Success";last=0.
  for stop in sort(unique(vcat([30.,100.,300.],collect(600.:300.:horizon),[horizon])))
   stop>horizon&&continue
   x,rc=segment(tag,"post"*lpad(Int(stop),3,'0'),x,post,clear+last,clear+stop,clear;level=level,step=.02);last=stop;rc=="Success"||break
  end
  rc=="Success"||(RUN_FAILURES[]+=1)
  appendrow(joinpath(DATA,"cases.csv"),"case,op,load_scale,P1,P2,P3,P7,taup,tauq,T,horizon,level,retcode",vcat([tag,op,g],ps,[taup,tauq,T,last,level,rc]));println("DONE ",tag," ",last," ",rc);flush(stdout)
 catch e
  RUN_FAILURES[]+=1
  println("ERROR ",tag," ",sprint(showerror,e,catch_backtrace()));flush(stdout)
  appendrow(joinpath(DATA,"cases.csv"),"case,op,load_scale,P1,P2,P3,P7,taup,tauq,T,horizon,level,retcode",vcat([tag,op,g],ps,[taup,tauq,T,0.,level,"exception"]))
 end
end
RUN_FAILURES[]>0&&error("$(RUN_FAILURES[]) case(s) failed; inspect cases.csv and logs.")
