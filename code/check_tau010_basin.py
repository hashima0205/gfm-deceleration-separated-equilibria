"""Check the adaptive basin's initial-state section, histories, and saved endpoints."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--section',choices=['fixed','chord'],default='fixed')
args=parser.parse_args()
chord=args.section=='chord'
P=ROOT/'basin_recalculation'/('tau010_chord' if chord else 'tau010_adaptive')
D=ROOT/'results/conditional_basins/data'
samples=pd.read_csv(P/'plot_samples_tau0.10.csv')
history=pd.read_csv(P/'history.csv')
initial=pd.read_csv(P/'initial_states.csv',float_precision='round_trip').set_index('case')
machines=pd.read_csv(P/'machines.csv')
params=pd.read_csv(P/'parameters.csv').set_index('bus')
roots=pd.read_csv(P/'equilibria.csv').set_index('root')
meta=json.loads((D/'definition.json').read_text())
states=pd.read_csv(D/'states.csv')
stable=np.loadtxt(ROOT/'data/seeds/tau010_T015_h150_stable.csv',delimiter=',')
saddle=np.loadtxt(ROOT/'data/seeds/tau010_T015_h150_saddle.csv',delimiter=',')
assert len(samples)==35 and samples['case'].is_unique
assert set(initial.index)==set(samples['case'])
assert np.isfinite(initial.to_numpy()).all()
assert (params.taup==.1).all() and (params.tauq==.02).all()
assert params.Pset.to_dict()=={1:.5,2:.4,3:.4,7:.4}
assert (params.VC_F==1).all() and (params.CC1_F==1).all()
assert roots.residual.max()<1e-8 and roots.loc['stable','max_real_eigen']<0
assert roots.loc['saddle','unstable_modes']==1
assert history.retcode.eq('Success').all()
assert samples.max_algebraic_residual.max()<1e-7
held=np.array(meta['held_differential_indices'])-1
index=lambda bus,key:int(states.loc[(states.bus==bus)&(states.symbol==f'inv{bus}₊droop₊{key}'),'index'].iloc[0])-1
delta=np.array([index(b,'δ') for b in [2,3,7]])
pf=np.array([index(b,'Pfilt') for b in [2,3,7]])
saddle_theta=np.mean(saddle[delta]-stable[delta])
xgrid=sorted(np.array([-.5,0.,.5,.9,1.,1.1,1.5])*saddle_theta) if chord else meta['theta_grid']
ygrid=[-.05,-.025,0,.025,.05] if chord else meta['df_grid']
expected=sorted((x,y) for x in xgrid for y in ygrid)
actual=sorted(map(tuple,samples[['theta_relative_rad','df_Hz']].to_numpy()))
assert np.allclose(actual,expected,rtol=0,atol=1e-12)
if chord:
    anchors=pd.read_csv(P/'slice_anchors.csv').set_index('root')
    assert anchors.max_full_state_difference.max()<1e-10
    for target,theta in [(stable,0.),(saddle,saddle_theta)]:
        hit=samples[np.isclose(samples.theta_relative_rad,theta,rtol=0,atol=1e-12)&np.isclose(samples.df_Hz,0)]
        assert len(hit)==1 and np.max(abs(initial.loc[hit['case'].iloc[0]].to_numpy()-target))<1e-10
for row in samples.itertuples(index=False):
    x=initial.loc[row.case].to_numpy()
    base=stable+(row.theta_relative_rad/saddle_theta)*(saddle-stable) if chord else stable
    assert np.allclose(x[held],base[held],rtol=0,atol=1e-12)
    assert np.allclose(x[delta],base[delta] if chord else stable[delta]+row.theta_relative_rad,rtol=0,atol=1e-12)
    assert np.allclose(x[pf]-base[pf],-row.df_Hz*2*np.pi/params.loc[[2,3,7],'Kp'],rtol=0,atol=1e-12)
    h=history[history['case']==row.case]
    stops=h.requested_horizon_s.tolist()
    assert stops==[30,60,150][:len(stops)] and h.completed_horizon_s.tolist()==stops
    assert h.classification.iloc[:-1].eq('unresolved').all()
    assert row.completed_horizon_s==stops[-1] and row.classification==h.classification.iloc[-1]
    end=np.loadtxt(P/f'{row.case}_end.csv',delimiter=',')
    t=pd.read_csv(P/f'{row.case}_trajectory.csv')
    assert t.t_s.is_monotonic_increasing and not t.t_s.duplicated().any()
    assert t.t_s.iloc[0]==0 and t.t_s.iloc[-1]==row.completed_horizon_s
    assert np.isclose(t.theta_relative_rad.iloc[0],row.theta_relative_rad,rtol=0,atol=1e-9)
    assert np.isclose(t.df_mean_Hz.iloc[0],row.df_Hz,rtol=0,atol=1e-9)
    assert np.isclose(np.mean(end[delta]-x[delta])/(2*np.pi),row.remote_internal_turns,rtol=0,atol=1e-10)
    assert np.isclose((t.theta_relative_rad.iloc[-1]-t.theta_relative_rad.iloc[0])/(2*np.pi),row.remote_internal_turns,rtol=0,atol=1e-10)
    m=machines[machines['case']==row.case]
    assert len(m)==4 and set(m.bus)=={1,2,3,7}
    for r in m.itertuples(index=False):
        assert np.isclose((end[index(r.bus,'δ')]-x[index(r.bus,'δ')])/(2*np.pi),r.internal_turns,rtol=0,atol=1e-10)
        if r.classification.startswith('recovery'):
            assert r.tail_df_max_Hz<1e-4 and r.tail_phase_width_rad<.005
            assert abs(np.angle(np.exp(1j*(end[index(r.bus,'δ')]-stable[index(r.bus,'δ')]))))<.002
    if row.theta_relative_rad==row.df_Hz==0:
        assert row.classification=='recovery' and np.max(abs(end-stable))<1e-7
faces=pd.read_csv(P/'plot_triangles_tau0.10.csv').to_numpy(dtype=int)
assert faces.shape==(48,3) and faces.min()==0 and faces.max()==34
xy=samples[['theta_relative_rad','df_Hz']].to_numpy()[faces]
a=xy[:,1]-xy[:,0];b=xy[:,2]-xy[:,0]
assert np.isclose(.5*np.sum(abs(a[:,0]*b[:,1]-a[:,1]*b[:,0])),(max(xgrid)-min(xgrid))*(max(ygrid)-min(ygrid)))
report=dict(points=len(samples),classification_counts=samples.classification.value_counts().to_dict(),
            section=args.section,
            stopping_horizon_counts={str(int(k)):int(v) for k,v in samples.completed_horizon_s.value_counts().sort_index().items()},
            total_integrated_time_s=float(samples.completed_horizon_s.sum()),
            all150_integrated_time_s=35*150,solve_wall_s=float(history.solve_s.sum()),
            max_algebraic_residual=float(samples.max_algebraic_residual.max()),
            scope='35 sampled initial states; adaptive finite observation, not a global or infinite-time basin')
(P/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
print('PASS: grid, initial-state section, staged continuation, endpoints, and mesh')
