from pathlib import Path
import numpy as np,pandas as pd
import json,sys
import observations12 as o
from classify import measure
P=Path(__file__).resolve().parents[1];D=P/'data';rows=[];machines=[]
for r in pd.read_csv(D/'cases.csv').drop_duplicates('case',keep='last').itertuples(index=False):
 if r.retcode!='Success':rows.append(dict(**r._asdict(),classification=r.retcode));continue
 cache=D/'results'/f'{r.case}_h{int(r.horizon)}.json'
 if cache.exists() and '--fresh' not in sys.argv:
  z=json.loads(cache.read_text());machines.extend([dict(case=r.case,**m) for m in z.pop('machines')]);rows.append(dict(**r._asdict(),**{k:v for k,v in z.items() if k!='case'},saddle_found=(D/f'seeds/{r.case}_saddle.csv').exists()));continue
 d=o.load(r.case);d=d[d.t_post_s<=r.horizon+1e-8];stable=np.loadtxt(D/f'seeds/{r.case}_stable.csv',delimiter=',');f=D/f'seeds/{r.case}_saddle.csv';saddle=np.loadtxt(f,delimiter=',') if f.exists() else stable
 z=measure(r.case,d,stable,saddle);machines.extend([dict(case=r.case,**m) for m in z.pop('machines')]);rows.append(dict(**r._asdict(),**{k:v for k,v in z.items() if k!='case'},saddle_found=f.exists()))
if rows:
 s=pd.DataFrame(rows);s.to_csv(D/'summary.csv',index=False);pd.DataFrame(machines).to_csv(D/'machines.csv',index=False);print(s.classification.value_counts().to_string());print(s[~s.classification.str.startswith('recovery')][['case','classification','remote_internal_turns']].to_string(index=False))
