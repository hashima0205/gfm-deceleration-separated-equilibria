"""Check saved basin data, initial-state construction, and cross-file consistency."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/conditional_basins';D=P/'data'
meta=json.loads((D/'definition.json').read_text())
jobs=pd.read_csv(D/'jobs.csv').set_index('case')
classification=pd.read_csv(D/'classification.csv').set_index('case')
initial=pd.read_csv(D/'initial_states.csv').set_index('case')
outcomes={z['case']:z for z in map(json.loads,(D/'outcomes.jsonl').read_text().splitlines())}
assert len(jobs)==len(classification)==len(initial)==len(outcomes)==281
assert jobs.index.is_unique and classification.index.is_unique and initial.index.is_unique
assert set(jobs.index)==set(classification.index)==set(initial.index)==set(outcomes)
machines=pd.read_csv(D/'machine_metrics.csv')
assert len(machines)==4*len(jobs) and set(machines['case'])==set(jobs.index)
assert all(set(group.bus)=={1,2,3,7} for _,group in machines.groupby('case'))
states=pd.read_csv(D/'states.csv')
delta=np.array([int(states.loc[states.symbol==f'inv{b}₊droop₊δ','index'].iloc[0])-1 for b in [2,3,7]])
pf=np.array([int(states.loc[states.symbol==f'inv{b}₊droop₊Pfilt','index'].iloc[0])-1 for b in [2,3,7]])
held=np.array(meta['held_differential_indices'])-1
assert len(held)==46 and len(delta)==len(pf)==3
stable=np.loadtxt(D/'references/section_stable.csv',delimiter=',')
saddle=np.loadtxt(D/'references/section_saddle.csv',delimiter=',')
assert abs(stable[delta].mean()-meta['origin_rad'])<1e-14
assert abs(saddle[delta].mean()-stable[delta].mean()-meta['saddle_theta_relative_rad'])<1e-14
coord_errors=[];background=[]
for case,r in jobs.iterrows():
    z=outcomes[case];c=classification.loc[case]
    assert z['classification']==c.classification
    assert z['horizon_s']>=r.horizon_s-1e-6
    assert abs(z['remote_internal_turns']-c.remote_internal_turns)<1e-9
    x=initial.loc[case].to_numpy(dtype=float);assert x.shape==(66,)
    base=stable if r.section=='equilibrium' else np.loadtxt(D/f'references/clearing_{r.taup:.2f}.csv',delimiter=',')
    assert np.max(abs(x[held]-base[held]))<1e-12
    error=max(abs(x[delta].mean()-meta['origin_rad']-r.theta_relative_rad),abs(1.2*np.mean(.4-x[pf])-r.df_Hz))
    coord_errors.append(error);assert error<1e-12

expected=pd.read_csv(ROOT/'tests/reference/expected_summary.csv')
assert meta['k']==1.0299 and meta['T']==.15 and meta['tauq']==.02
assert meta['Pset']==[.5,.4,.4,.4]
assert np.allclose(expected.load_scale,meta['load_scale'],rtol=0,atol=1e-14)
for tau in meta['taus']:
    samples=pd.read_csv(D/f'plot_samples_tau{tau:.2f}.csv')
    faces=pd.read_csv(D/f'plot_triangles_tau{tau:.2f}.csv').to_numpy(dtype=int)
    assert len(samples)==88 and np.min(faces)>=0 and np.max(faces)<len(samples)
    assert (samples.section=='equilibrium').all() and (samples.level=='A').all()
    assert samples['case'].is_unique
    for r in samples.itertuples():assert r.classification==classification.loc[r.case,'classification']
    background.append(dict(taup_s=tau,points=len(samples),**samples.classification.value_counts().to_dict()))
    eq=pd.read_csv(D/f'equilibrium_{tau:.2f}.csv').set_index('root')
    assert (eq.residual<1e-8).all() and eq.loc['stable','max_real_eigen']<0 and eq.loc['saddle','unstable_modes']==1
    d=pd.read_csv(D/f'trajectories/fault_tau{tau:.2f}.csv.gz')
    assert d.t_post_s.is_unique and np.all(np.diff(d.t_post_s)>0)
    assert abs(d.theta_relative_rad.iloc[0]-meta['sections'][str(tau)]['clearing_theta_rad'])<1e-12
    angle=d[[f'b{b}_delta_rad' for b in [2,3,7]]].mean(axis=1)
    assert np.max(abs(angle-meta['origin_rad']-d.theta_relative_rad))<1e-12
    ref=expected[np.isclose(expected.taup_s,tau)].iloc[0]
    assert d.t_post_s.iloc[-1]==ref.horizon_s
    turns=(angle.iloc[-1]-angle.iloc[0])/(2*np.pi)
    assert abs(turns-ref.remote_turns)<1e-9

hidden=pd.read_csv(D/'hidden_state_comparison.csv')
pairs=hidden[~hidden.classification_matches]
assert len(pairs)==3
for r in pairs.itertuples():
    for case,turns in [(r.clearing_case,r.clearing_turns),(r.equilibrium_case,r.equilibrium_turns)]:
        d=pd.read_csv(D/f'hidden_trajectories/{case}.csv.gz')
        assert d.t_post_s.iloc[-1]==1000
        assert abs(d.theta_relative_rad.iloc[0]-r.theta_rad)<1e-12
        assert abs(d.df_mean_Hz.iloc[0]-r.df_Hz)<1e-12
        assert abs((d.theta_relative_rad.iloc[-1]-d.theta_relative_rad.iloc[0])/(2*np.pi)-turns)<1e-8

audit=pd.read_csv(D/'initialization_audit.csv')
assert audit.algebraic_residual.max()<1e-8
assert audit.held_state_change.max()==0 and audit.dynamic_change.max()==0
checks=dict(unique_cases=len(jobs),background_points=sum(r['points'] for r in background),
            max_reconstructed_coordinate_error=float(max(coord_errors)),
            historical_initialization_algebraic_residual_max=float(audit.algebraic_residual.max()),
            background_counts=background,hidden_pairs_with_different_outcomes=len(pairs),
            scope='Integrity of saved data and recorded audits; no new integration of the full map.')
manifest=P/'SHA256SUMS'
if manifest.exists():
    for line in manifest.read_text().splitlines():
        digest,name=line.split(maxsplit=1)
        assert hashlib.sha256((P/name).read_bytes()).hexdigest()==digest,name
    checks['sha256_manifest']='PASS'
print(json.dumps(checks,indent=2))
print('PASS: saved data, 264 map points, full initial states, trajectories, and recorded root checks')
