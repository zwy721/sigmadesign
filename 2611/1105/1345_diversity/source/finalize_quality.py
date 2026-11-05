"""Merge all frozen QC records; preserve fixed queue priority, never quality-rank."""
from common import *
import time,argparse
from denovo_quality.common import object_hash,versions
def read_records(path):
 records=[]
 if path.exists():
  for line in path.read_text().splitlines():
   try:records.append(json.loads(line))
   except json.JSONDecodeError:pass
 return records
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--snapshot',action='store_true');args=parser.parse_args()
 manifest=load(OUT/'main/manifest.json');old=read_records(OUT/'main/quality/stream_metrics.jsonl');boost=read_records(OUT/'main/quality_boost/stream_metrics.jsonl');by={};differences=[]
 for queue,records in [('original',old),('boost',boost)]:
  for r in records:
   mid=r['molecule_id']
   if mid in by:
    if by[mid]['checks']!=r['checks']:differences.append(dict(molecule_id=mid,original=by[mid]['checks'],duplicate=r['checks']))
   else:by[mid]=r
 missing=[m for m in manifest if m['molecule_id'] not in by];write_json(OUT/'quality_completion.json',dict(expected=len(manifest),original_records=len(old),boost_records=len(boost),unique_completed=len(by),missing=len(missing),duplicate_check_differences=differences));savecsv(OUT/'PENDING_QUALITY_RECORDS.csv',missing)
 if not missing:
  with (OUT/'PENDING_QUALITY_RECORDS.csv').open('w',newline='') as file:csv.writer(file).writerow(list(manifest[0]) if manifest else ['molecule_id'])
 if missing and not args.snapshot:print('QC pending',len(missing),'/',len(manifest));return
 assert not differences,'Duplicate frozen check values differ; inspect quality_completion.json'
 results=[by[m['molecule_id']] for m in manifest if m['molecule_id'] in by];out=OUT/'main/quality';cfg=load_config(QC);src=Path(__file__).parent;compute=object_hash({str(p):sha256(p) for p in list((ROOT/'denovo_quality').glob('*.py'))+[src/'common.py',src/'evaluate.py']})
 with (out/'molecule_metrics.jsonl').open('w') as file:
  for r in results:file.write(json.dumps(clean(r),ensure_ascii=False,allow_nan=False)+'\n')
 write_json(out/'protocol.json',dict(manifest_sha256=sha256(OUT/'main/manifest.json'),config=cfg,PB_config=load(OUT/'posebusters_no_water.json'),versions=versions(),computation_sha256=compute,input_count=len(manifest),completed_count=len(results),pending_count=len(missing),execution='Original 8 workers plus additional 32 CPU workers, same QC functions/config/hash; original queue priority independent of pass/fail',source_paths=['main/quality','main/quality_boost']))
 print('QC merged',len(results),'pending',len(missing),flush=True)
if __name__=='__main__':main()

