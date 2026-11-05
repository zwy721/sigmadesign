from common import *
import argparse,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from denovo_quality.common import object_hash,versions

def work(args):
 install_pb()
 from denovo_quality.pipeline import task
 return task(args)
def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(OUT/'manifest.json'));p.add_argument('--output',default=str(OUT/'quality'));p.add_argument('--workers',type=int,default=8);a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
 rows=load(a.manifest);cfg=load_config(QC);install_pb();deps=versions();compute=object_hash({str(p):sha256(p) for p in list((ROOT/'denovo_quality').glob('*.py'))+list(Path(__file__).parent.glob('*.py'))});start=time.perf_counter()
 write_json(out/'protocol.json',dict(manifest_sha256=sha256(a.manifest),config=cfg,PB_config=load(OUT/'posebusters_no_water.json'),versions=deps,computation_sha256=compute,input_count=len(rows)))
 results=[]
 with ProcessPoolExecutor(a.workers) as pool:
  futures=[pool.submit(work,(r,cfg,str(out),compute,deps,False)) for r in rows]
  for f in as_completed(futures):
   r=f.result();results.append(r)
   if len(results)%24==0:print('QC',len(results),'/',len(rows),'HQ_no_water',sum(x['checks'].get('high_quality') is True for x in results),flush=True)
 by={r['molecule_id']:r for r in results};results=[by[r['molecule_id']] for r in rows]
 with (out/'molecule_metrics.jsonl').open('w') as f:
  for r in results:f.write(json.dumps(clean(r),ensure_ascii=False,allow_nan=False)+'\n')
 write_json(out/'timing.json',dict(wall_seconds=time.perf_counter()-start,records=len(results),workers=a.workers));print('COMPLETE',flush=True)
if __name__=='__main__':main()
