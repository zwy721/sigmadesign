"""CPU quality evaluation as completed attempts arrive. Cache matches evaluate.py."""
from common import *
import time
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from denovo_quality.common import object_hash,versions
from evaluate import work
def main():
 out=OUT/'main/quality';out.mkdir(parents=True,exist_ok=True);cfg=load_config(QC);install_pb();deps=versions();source=Path(__file__).parent
 compute=object_hash({str(p):sha256(p) for p in list((ROOT/'denovo_quality').glob('*.py'))+[source/'common.py',source/'evaluate.py']})
 plan=load(OUT/'frozen_plan.json');expected=len(plan['cohort'])*len(plan['seeds'])*len(plan['strategies']);seen=set();pending={};finished=[];start=time.perf_counter();known=set()
 with ProcessPoolExecutor(8) as pool:
  while True:
   for path in (OUT/'main/structures').glob('*/result.json'):
    if str(path) in known:continue
    r=load(path);known.add(str(path))
    for m in r['manifest']:
     if m['molecule_id'] in seen:raise AssertionError('Duplicate quality id')
     seen.add(m['molecule_id']);f=pool.submit(work,(m,cfg,str(out),compute,deps,False));pending[f]=m
   done,_=wait(list(pending),timeout=3,return_when=FIRST_COMPLETED) if pending else (set(),set())
   for f in done:
    r=f.result();finished.append(r);del pending[f]
    with (out/'stream_metrics.jsonl').open('a') as file:file.write(json.dumps(clean(r),ensure_ascii=False,allow_nan=False)+'\n')
    if len(finished)%24==0:print('STREAM',len(known),'/',expected,'attempts;',len(finished),'QC records;',len(pending),'pending',flush=True)
   if len(known)==expected and not pending:break
   if not pending:time.sleep(3)
 write_json(out/'stream_timing.json',dict(seconds=time.perf_counter()-start,records=len(finished),attempts=len(known)));print('STREAM COMPLETE',flush=True)
if __name__=='__main__':main()

