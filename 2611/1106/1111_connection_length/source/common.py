"""Isolated GEOM connection/geometry experiment; upstream files are read-only."""
import os,sys,json,hashlib
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]
ROOT=OUT.parent
RELEASES=ROOT.parent
PROJECT=RELEASES/'oracle_rebuild_20261031'
BUNDLE=RELEASES/'cutoff_geom3684_pdbbind_20261024'
DATA=Path('/data/zhouwanyao/cutoff_geom3684_pdbbind_20261024/geom')
ASSETS=RELEASES/'geom_m12/assets'
BASELINE=PROJECT/'results/connection_robust_20261101/training_v2/best.pt'
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'cache/matplotlib'))
os.environ.setdefault('WANDB_MODE','disabled')
os.environ.setdefault('PYTHONDONTWRITEBYTECODE','1')
sys.dont_write_bytecode=True
sys.path.insert(0,str(PROJECT))
sys.path.append(str(BUNDLE/'src'))
import bootstrap
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def dump(path,value):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def load(path):return json.loads(Path(path).read_text())
def stable(*args):return int(hashlib.sha256('|'.join(map(str,args)).encode()).hexdigest()[:16],16)
