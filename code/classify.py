"""Finite-horizon, per-machine classification; no interpolation of uncomputed labels."""
from pathlib import Path
import sys,json,numpy as np,pandas as pd
import observations12 as o
P=Path(__file__).resolve().parents[1];D=P/'data';B=[1,2,3,7];R=B[1:]
def measure(tag,d,stable,saddle):
 TARGET=o.targets(stable)
 gap=abs(np.angle(np.exp(1j*(o.targets(saddle)-TARGET))));TOL=np.where(gap>1e-8,np.minimum(.002,.1*gap),.002)
 p=d[d.stage.str.startswith('post')].drop_duplicates('t_post_s',keep='last');t=p.t_post_s.to_numpy();end=t[-1];tail=p[t>=end-min(50,end/3)];ms=[]
 for j,b in enumerate(B):
  a=p[f'b{b}_delta_rad'].to_numpy();shift=a-a[0];df=tail[f'b{b}_df_Hz'].to_numpy();windows=np.diff(np.interp(np.linspace(.4*end,end,4),t,a))/(2*np.pi);err=abs(np.angle(np.exp(1j*(a[-1]-TARGET[j]))));turns=shift[-1]/(2*np.pi);width=np.ptp(tail[f'b{b}_delta_rad'])
  if max(abs(df))<1e-4 and err<TOL[j] and width<.005:label='recovery_after_slip' if max(abs(shift))>=2*np.pi else 'recovery'
  elif turns<-2 and max(windows)<-.2 and windows.sum()/(.6*end)<-1e-4:label='negative_running'
  elif turns>2 and min(windows)>.2 and windows.sum()/(.6*end)>1e-4:label='positive_running'
  else:label='unresolved'
  phase=p[f'b{b}_terminal_phase_unwrapped_rad'].to_numpy();valid=bool(p[f'b{b}_terminal_phase_reliable'].all())
  ms.append(dict(bus=b,classification=label,internal_turns=float(turns),terminal_turns=float((phase[-1]-phase[0])/(2*np.pi)) if valid else None,terminal_phase_reliable=valid,post_max_phase_increment_rad=float(np.max(abs(np.diff(phase)))),tail_df_mean_Hz=float(np.trapezoid(df,tail.t_post_s)/(tail.t_post_s.iloc[-1]-tail.t_post_s.iloc[0])),tail_df_max_Hz=float(max(abs(df))),phase_error_rad=float(err),phase_tolerance_rad=float(TOL[j]),tail_phase_width_rad=float(width),min_window_turns=float(min(windows)),max_window_turns=float(max(windows)),min_voltage_pu=float(p[f'b{b}_V_pu'].min()),max_voltage_pu=float(p[f'b{b}_V_pu'].max()),max_limiter_input_pu=float(np.hypot(p[f'b{b}_ref_in_d_pu'],p[f'b{b}_ref_in_q_pu']).max()),saturation_sample_duration_s=float(np.trapezoid(p[f'b{b}_saturated'],t))))
 labs=[x['classification'] for x in ms]
 if all(x in ['recovery','recovery_after_slip'] for x in labs):label='recovery_after_slip' if 'recovery_after_slip' in labs else 'recovery'
 elif 'negative_running' in labs and 'positive_running' in labs:label='other'
 elif 'negative_running' in labs:label='negative_running'
 elif 'positive_running' in labs:label='positive_running'
 else:label='unresolved'
 # This duration comparison classifies saved trajectories; no missing full-state endpoint is inferred.
 out=dict(case=tag,horizon_s=float(end),classification=label,remote_internal_turns=float((p.theta_mean_rad.iloc[-1]-p.theta_mean_rad.iloc[0])/(2*np.pi)),max_remote_angle_spread_rad=float(np.ptp(p[[f'b{b}_delta_rad' for b in R]].to_numpy(),axis=1).max()),max_remote_df_spread_Hz=float(np.ptp(p[[f'b{b}_df_Hz' for b in R]].to_numpy(),axis=1).max()),max_algebraic_residual=float(p.alg_residual.max()),machines=ms)
 path=D/'results'/f'{tag}.json';tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2,allow_nan=False));tmp.replace(path)
 (D/'results'/f'{tag}_h{int(end)}.json').write_text(json.dumps(out,indent=2,allow_nan=False))
 return out
