"""Detach the already authorized small formal experiment and record its PID."""
from common import *
import ast,subprocess,argparse,time
from datetime import datetime,timezone

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true');args=parser.parse_args()
 if args.resume:
  from run_experiment import validate_resume
  recovery=validate_resume()
 else:assert not (OUT/'formal_start.json').exists(),'A formal run already exists; use --resume without resetting the 24-hour clock'
 if (OUT/'process.json').exists():
  previous=load(OUT/'process.json')
  try:os.kill(previous['pid'],0)
  except ProcessLookupError:pass
  else:raise RuntimeError(f"Recorded experiment PID {previous['pid']} is still alive; refusing duplicate GPU job")
 for name in ['preflight','prelaunch_verification','dataset_verification','short_workflow_verification']:
  assert load(OUT/f'{name}.json')['passed'],name
 for p in (OUT/'source').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
 inputs=load(OUT/'data/complete.json');assert inputs['complete']
 assert all(sha(OUT/'data'/f'{s}.pt')==v for s,v in inputs['split_pt_sha256'].items())
 assert all(sha(p)==v['sha256'] for p,v in load(OUT/'audit_initial.json')['source_fingerprints'].items())
 env=dict(os.environ,CUDA_VISIBLE_DEVICES='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='2',PYTHONDONTWRITEBYTECODE='1',WANDB_MODE='disabled')
 argv=[sys.executable,'-B',str(OUT/'source/run_experiment.py')]
 if args.resume:
  argv.append('--resume');dump(OUT/'process_history'/f'{time.time_ns()}.json',load(OUT/'process.json'))
 with (OUT/'training.log').open('ab',buffering=0) as log:
  proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 dump(OUT/'process.json',dict(pid=proc.pid,UTC=datetime.now(timezone.utc).isoformat(),argv=argv,cwd=str(ROOT),log=str(OUT/'training.log'),status=str(OUT/'status.json'),physical_GPU=0,seeds=[20261106,20261107,20261108],wall_budget_hours=24,protocol_sha256=sha(OUT/'protocol.json'),resumed=args.resume,original_start_UTC=load(OUT/'formal_start.json')['UTC'] if args.resume else None,deadline_UTC=recovery['deadline_UTC'] if args.resume else None))
 print(json.dumps(load(OUT/'process.json'),ensure_ascii=False),flush=True)
if __name__=='__main__':main()
