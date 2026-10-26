"""Compare freshly integrated endpoint states with the original study31 results."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import observations12 as observations
from classify import measure

P=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--smoke',action='store_true')
parser.add_argument('--all',action='store_true',help='Check all three demo conditions')
args=parser.parse_args()
assert not (args.smoke and args.all)
cases=pd.read_csv(P/'data/cases.csv').drop_duplicates('case',keep='last')
names=['smoke_tau120'] if args.smoke else ['nf02_1_p1.200_q0.02_T0.15']
if args.all:names=['nf02_1_p1.160_q0.02_T0.15','nf02_1_p1.180_q0.02_T0.15','nf02_1_p1.200_q0.02_T0.15']
rows=[]
for name in names:
    row=cases.set_index('case').loc[name]
    assert row.retcode=='Success'
    expected_horizon=30 if args.smoke else (1000 if 'p1.200' in name else 300)
    assert row.horizon==expected_horizon
    source='nf02_1_p1.200_q0.02_T0.15' if args.smoke else name
    for stage in ['pre','fault',f'post{expected_horizon:03d}']:
        a=np.loadtxt(P/'checkpoints'/name/f'{stage}_end.csv',delimiter=',')
        b=np.loadtxt(P/'tests/reference'/f'{source}__{stage}_end.csv',delimiter=',')
        assert a.shape==b.shape==(66,)
        assert np.allclose(a,b,rtol=1e-9,atol=2e-7),(name,stage,float(np.max(abs(a-b))))
        # Max absolute difference is a numerical diagnostic with mixed units.
        rows.append(dict(case=name,stage=stage,max_endpoint_difference_mixed_units=float(np.max(abs(a-b)))))
    audit=pd.read_csv(P/'data/parameter_audit.csv')
    audit=audit[audit['case']==name].drop_duplicates(['stage','bus'],keep='last')
    assert len(audit)==12 and set(audit.bus)=={1,2,3,7}
    assert (audit.taup==float(row.taup)).all() and (audit.tauq==.02).all()
    assert (audit.VC_F==1).all() and (audit.CC1_F==1).all()
    assert (audit.Pset>0).all()
    equilibria=pd.read_csv(P/'data/equilibria.csv')
    equilibria=equilibria[(equilibria['case']==name)&(equilibria.stage=='post')].drop_duplicates('branch',keep='last')
    assert (equilibria.residual<1e-8).all()
    assert (equilibria.max_real_eigen<0).any() and (equilibria.unstable_modes==1).any()
    if not args.smoke:
        stable=np.loadtxt(P/'data/seeds'/f'{name}_stable.csv',delimiter=',')
        saddle=np.loadtxt(P/'data/seeds'/f'{name}_saddle.csv',delimiter=',')
        metrics=measure(name,observations.load(name),stable,saddle)
        expected=pd.read_csv(P/'tests/reference/expected_summary.csv').set_index('case').loc[name]
        assert metrics['classification']==expected.classification
        assert abs(metrics['remote_internal_turns']-expected.remote_turns)<1e-6
        gap=abs(np.angle(np.exp(1j*(observations.targets(saddle)[1:].mean()-observations.targets(stable)[1:].mean()))))
        assert abs(gap-expected.remote_stable_saddle_gap_rad)<1e-10
        print(name,metrics['classification'],metrics['remote_internal_turns'],'turns; gap',gap,'rad')
result=pd.DataFrame(rows)
result.to_csv(P/'data/validation.csv',index=False)
print(result.to_string(index=False))
print('PASS: original 66-state endpoints and all-GFM parameter audit')
