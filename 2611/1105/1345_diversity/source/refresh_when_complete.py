"""Complete auxiliary QC delivery automatically under the unchanged protocol."""
from common import *
import time,subprocess
def ids(path):
 result=set()
 if path.exists():
  for line in path.read_text().splitlines():
   try:result.add(json.loads(line)['molecule_id'])
   except json.JSONDecodeError:pass
 return result
def main():
 expected={m['molecule_id'] for m in load(OUT/'main/manifest.json')};last=None
 while True:
  complete=ids(OUT/'main/quality/stream_metrics.jsonl')|ids(OUT/'main/quality_boost/stream_metrics.jsonl');pending=len(expected-complete)
  if pending!=last:print('Auxiliary QC pending',pending,flush=True);last=pending
  if pending==0:break
  time.sleep(30)
 for module in ['finalize_quality','analyze','reports','plot_results','verify']:
  subprocess.run([sys.executable,'-B',str(OUT/'source'/f'{module}.py')],check=True)
 state=load(OUT/'RUN_STATE.json');state.update(optional_initial_QC_pending=0,all_frozen_QC_completed=len(expected),verification=load(OUT/'verification.json')['passed'],timeout_policy_question='No runtime cap introduced; original auxiliary checks completed.');write_json(OUT/'RUN_STATE.json',state)
 subprocess.run([sys.executable,'-B',str(OUT/'source/fingerprint_delivery.py')],check=True)
 print('All auxiliary QC and reports refreshed',flush=True)
if __name__=='__main__':main()

