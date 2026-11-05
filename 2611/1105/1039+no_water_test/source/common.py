import sys,json,csv,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parents[1]
OLD=ROOT/'2611/1104/1425+水约束与连接区域柔性'
sys.path.insert(0,str(ROOT))
import numpy as np
from denovo_quality.common import clean,write_json,sha256,load_config
QC=ROOT/'geometry_constraints_20261103/quality_config.json'
WATER={'HOH','WAT','DOD'}
WATER_TESTS={'minimum_distance_to_waters','volume_overlap_with_waters'}
def load(p):return json.loads(Path(p).read_text())
def jsonl(p):return [json.loads(s) for s in Path(p).read_text().splitlines() if s.strip()]
def key(r):
 m=r.get('metadata',r);return m['receptor_id'],int(m['seed']),int(m.get('candidate_rank',m.get('rank',0)))
def name(r):return '%s_%d_%02d'%key(r)
def savecsv(p,rows):
 if not rows:return
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader()
  for r in rows:w.writerow({k:json.dumps(clean(v),ensure_ascii=False) if isinstance(v,(dict,list,tuple)) else clean(v) for k,v in r.items()})
def sdf(p,m):
 from rdkit import Chem
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 with Chem.SDWriter(str(p)) as w:w.write(m)
def install_pb():
 import copy
 from functools import lru_cache
 from denovo_quality import posebusters_adapter as pb
 if getattr(pb,'_no_water_installed',False):return
 original=pb.effective_config
 @lru_cache(maxsize=8)
 def effective(mode,threads):
  c=copy.deepcopy(original(mode,threads));kept=[]
  for m in c['modules']:
   names=set(pb.test_names({'modules':[m]}))
   if names&WATER_TESTS:
    assert names<=WATER_TESTS
   else:kept.append(m)
  c['modules']=kept
  assert not (set(pb.test_names(c))&WATER_TESTS)
  return c
 pb.effective_config=effective;pb._no_water_installed=True
