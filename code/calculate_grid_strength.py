"""Nominal-voltage, network-only fault levels; no assumed converter rating."""
import ast
import csv
import json
import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
source = (ROOT/'provenance/source_snapshot/studies/04_clearing/07_grid_strength/code/run.jl').read_text()
branches = ast.literal_eval(re.search(r'^const branches=(\[.*\])$', source, re.M).group(1))
k = json.loads((ROOT/'results/conditional_basins/data/definition.json').read_text())['k']
nodes = [1, 2, 3, 4, 5, 7, 8]
rows = []
for stage in ['pre', 'post']:
    Y = np.zeros((8, 8), complex)
    for j, (a, b, magnitude) in enumerate(branches, 1):
        if stage == 'post' and j == 7:
            continue
        y = magnitude/(k if j in [2, 3] else 1)*np.exp(-1j*np.deg2rad(75))
        a -= 1; b -= 1
        Y[a, a] += y; Y[b, b] += y
        Y[a, b] -= y; Y[b, a] -= y
    # Bus6 is the only external voltage source. Other converters and loads
    # are omitted; no converter/filter short-circuit contribution is credited.
    reduced = Y[np.ix_(np.array(nodes)-1, np.array(nodes)-1)]
    Z = np.linalg.inv(reduced)
    assert np.max(abs(reduced @ Z - np.eye(7))) < 1e-12
    for bus in [1, 2, 3, 5, 7]:
        z = Z[nodes.index(bus), nodes.index(bus)]
        fault_level = 1/abs(z)  # V_nom=1 pu; SCR=fault_level/rating_on_system_base
        rows.append(dict(stage=stage, bus=bus, Rth_pu=z.real, Xth_pu=z.imag,
                         Ssc_nominal_pu=fault_level))
        if stage == 'post' and bus == 5:
            assert abs(z-k/.68*np.exp(1j*np.deg2rad(75))) < 1e-12
out = ROOT/'data/grid_strength'
out.mkdir(exist_ok=True)
with (out/'network_fault_levels.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader(); writer.writerows(rows)
for row in rows:
    print(row)
