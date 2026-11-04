"""Offline paired geometry study; no model or reference-conformer imports."""
import csv
import json
import os
import sys
from pathlib import Path
import numpy as np

TASK = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get('WATER_PORT_STUDY_OUTPUT', str(TASK))).resolve()
ROOT = TASK.parents[2]
sys.path.insert(0, str(ROOT))
OLD = ROOT / 'fbdd_efficiency_generalization_20261103/independent'
DIAG = ROOT / 'independent_diagnostics_20261104'
from denovo_quality.common import clean, sha256, write_json, load_config

CONFIG = OUT / 'fixed_parameters.json'
QC = ROOT / 'geometry_constraints_20261103/quality_config.json'
WATER_TESTS = {'minimum_distance_to_waters', 'volume_overlap_with_waters'}

def load(p):
    return json.loads(Path(p).read_text())

def jsonl(p):
    return [json.loads(l) for l in Path(p).read_text().splitlines() if l.strip()]

def key(r):
    m = r.get('metadata', r)
    return m['receptor_id'], str(m['seed']), int(m.get('candidate_rank', m.get('rank', 0)))

def name(r):
    p, s, k = key(r)
    return f'{p}_{s}_{k:02d}'

def aux_hq(r):
    tests = r['modules']['posebusters']['tests']
    req = load(QC)['quality']['required_checks']
    return all(v is True for k, v in tests.items() if k not in WATER_TESTS) and all(r['checks'].get(k) is True for k in req if k != 'posebusters_all_pass')

def savecsv(path, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(k for r in records for k in r))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for r in records:
            writer.writerow({k: json.dumps(clean(v), ensure_ascii=False) if isinstance(v, (dict, list, tuple)) else clean(v) for k, v in r.items()})

def save_sdf(path, mol):
    from rdkit import Chem
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(path)) as w:
        w.write(mol)

def waters_for(row):
    return [r for r in load(OUT / 'fixed_waters.json') if r['receptor_id'] == row['receptor_id']]
