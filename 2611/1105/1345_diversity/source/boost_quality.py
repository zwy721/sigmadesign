"""Additional CPU QC queue. Same module/config/hash, independent output cache.
The original eight expensive jobs keep running. Avoid duplicate active jobs;
queue their known expensive seed/strategy group after other pending work.
"""
from common import *
import time
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from denovo_quality.common import object_hash,versions
from evaluate import work
def old_state():
 done=set()
 p=OUT/'main/quality/stream_metrics.jsonl'
 if p.exists():
  for line in p.read_text().splitlines():
   try:done.add(json.loads(line)['molecule_id'])
   except json.JSONDecodeError:pass
 active=set()
 for pid in range(2055502,2055510):
  for path in (Path('/proc')/str(pid)/'fd').glob('*'):
   try:
    dest=path.readlink()
    if dest.name=='evaluation.log':active.add(dest.parent.name.rsplit('_',1)[0])
   except OSError:pass
 return done,active
def main():
 out=OUT/'main/quality_boost';out.mkdir(exist_ok=True);cfg=load_config(QC);install_pb();deps=versions();src=Path(__file__).parent
 compute=object_hash({str(p):sha256(p) for p in list((ROOT/'denovo_quality').glob('*.py'))+[src/'common.py',src/'evaluate.py']})
 plan=load(OUT/'frozen_plan.json');expected=len(plan['cohort'])*len(plan['seeds'])*len(plan['strategies']);seen=set();pending={};finished=[];start=time.perf_counter()
 with ProcessPoolExecutor(32) as pool:
  while True:
   old,active=old_state();attempts=[load(p) for p in (OUT/'main/structures').glob('*/result.json')];candidates=[m for r in attempts for m in r['manifest'] if m['molecule_id'] not in seen|old|active]
   candidates.sort(key=lambda m:(m['seed']==41008 and m['strategy'] in ['joint_prior_iterative','alternating_categorical'],m['seed'],m['molecule_id']))
   for m in candidates:
    if len(pending)>=64:break
    seen.add(m['molecule_id']);pending[pool.submit(work,(m,cfg,str(out),compute,deps,False))]=m
   done,_=wait(list(pending),timeout=2,return_when=FIRST_COMPLETED) if pending else (set(),set())
   for future in done:
    r=future.result();finished.append(r);del pending[future]
    with (out/'stream_metrics.jsonl').open('a') as file:file.write(json.dumps(clean(r),ensure_ascii=False,allow_nan=False)+'\n')
    if len(finished)%24==0:print('BOOST',len(attempts),'attempts;',len(finished),'QC records;',len(pending),'pending;',len(active),'original active',flush=True)
   full={m['molecule_id'] for r in attempts for m in r['manifest']}
   if len(attempts)==expected and full<=old|{r['molecule_id'] for r in finished} and not pending:break
   if not pending:time.sleep(3)
 write_json(out/'timing.json',dict(seconds=time.perf_counter()-start,records=len(finished)));print('BOOST COMPLETE',flush=True)
if __name__=='__main__':main()

