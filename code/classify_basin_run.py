"""Compare one recalculated basin trajectory with the saved finite-horizon result."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import observations12 as o
from classify import measure

P=Path(__file__).resolve().parents[1];D=P/'results/conditional_basins/data'
source,tag=sys.argv[1:]
stable=np.loadtxt(D/'references/section_stable.csv',delimiter=',')
saddle=np.loadtxt(D/'references/section_saddle.csv',delimiter=',')
actual=measure(tag,o.load(tag),stable,saddle)
expected=pd.read_csv(D/'classification.csv').set_index('case').loc[source]
result=dict(source_case=source,rerun_case=tag,expected_classification=expected.classification,
            actual_classification=actual['classification'],horizon_s=actual['horizon_s'],
            internal_turns=actual['remote_internal_turns'],expected_internal_turns=expected.remote_internal_turns,
            absolute_turn_difference=abs(actual['remote_internal_turns']-expected.remote_internal_turns))
out=P/'basin_recalculation';out.mkdir(exist_ok=True)
(out/f'{source}.json').write_text(json.dumps(result,indent=2)+'\n')
assert actual['classification']==expected.classification,result
assert result['absolute_turn_difference']<1e-6,result
print(json.dumps(result,indent=2))
