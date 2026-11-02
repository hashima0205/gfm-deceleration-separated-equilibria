"""Export verified fault-on observations for the saved three-condition basin plots."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / 'results/conditional_basins/data'


def main():
    meta = json.loads((BUNDLE / 'definition.json').read_text())
    manifest = pd.read_csv(ROOT / 'data/fault_overlay.csv')
    completed = pd.read_csv(ROOT / 'data/cases.csv').drop_duplicates('case', keep='last').set_index('case')
    audit = pd.read_csv(ROOT / 'data/parameter_audit.csv')
    exported, checks = [], []
    for row in manifest.itertuples(index=False):
        result = completed.loc[row.case]
        assert result.retcode == 'Success', row.case
        for key in manifest.columns.drop('case'):
            assert result[key] == getattr(row, key), (row.case, key)
        assert row.tauq == meta['tauq'] and row.T == meta['T']
        assert row.load_scale == meta['load_scale']
        assert [row.P1, row.P2, row.P3, row.P7] == meta['Pset']
        params = audit[audit['case'] == row.case].drop_duplicates(['stage', 'bus'], keep='last')
        assert len(params) == 12 and set(params.bus) == {1, 2, 3, 7}
        assert (params.taup == row.taup).all() and (params.tauq == row.tauq).all()
        assert (params.VC_F == 1).all() and (params.CC1_F == 1).all()
        assert np.allclose(params.Pset, params.bus.map({1: row.P1, 2: row.P2, 3: row.P3, 7: row.P7}), rtol=0, atol=0)
        directory = ROOT / 'checkpoints' / row.case
        reference = f'nf02_1_p{row.taup:.3f}_q0.02_T0.15'
        for stage in ['pre', 'fault', 'post030']:
            actual = np.loadtxt(directory / f'{stage}_end.csv', delimiter=',')
            expected = np.loadtxt(ROOT / 'tests/reference' / f'{reference}__{stage}_end.csv', delimiter=',')
            assert actual.shape == expected.shape == (66,)
            assert np.allclose(actual, expected, rtol=1e-9, atol=2e-7), (row.case, stage)
            checks.append(dict(case=row.case, stage=stage,
                               max_endpoint_difference_mixed_units=float(np.max(abs(actual - expected)))))
        columns = ['t_absolute_s', 't_post_s', 'stage'] + [
            f'b{bus}_{key}' for bus in [2, 3, 7] for key in ['delta_rad', 'omega_rad_s']]
        d = pd.read_csv(directory / 'fault.csv', usecols=columns).copy()
        assert d.stage.eq('fault').all() and d.t_absolute_s.is_monotonic_increasing
        assert np.isclose(d.t_post_s.iloc[0], -row.T) and np.isclose(d.t_post_s.iloc[-1], 0.)
        d['theta_relative_rad'] = d[[f'b{b}_delta_rad' for b in [2, 3, 7]]].mean(axis=1) - meta['origin_rad']
        d['df_mean_Hz'] = d[[f'b{b}_omega_rad_s' for b in [2, 3, 7]]].mean(axis=1) / (2 * np.pi) - 60
        saved = pd.read_csv(BUNDLE / f'trajectories/fault_tau{row.taup:.2f}.csv.gz').iloc[0]
        for key in ['theta_relative_rad', 'df_mean_Hz']:
            assert np.isclose(d[key].iloc[-1], saved[key], rtol=0, atol=1e-7), (row.case, key)
        d['taup'] = row.taup
        d['source_case'] = row.case
        exported.append(d[['taup', 'source_case', 't_absolute_s', 't_post_s', 'stage',
                           'theta_relative_rad', 'df_mean_Hz']])
    output = ROOT / 'data/fault_trajectories.csv.gz'
    pd.concat(exported, ignore_index=True).to_csv(output, index=False, compression={'method': 'gzip', 'mtime': 0})
    pd.DataFrame(checks).to_csv(ROOT / 'logs/fault_overlay_validation.csv', index=False)
    print(pd.DataFrame(checks).to_string(index=False))
    print('PASS: fault trajectories, full 66-state endpoints, parameters, and clearing continuity')
    print('Exported:', output)


if __name__ == '__main__':
    main()
