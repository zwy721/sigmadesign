"""One-GPU, 24-hour, preregistered small-scale experiment supervisor."""
from common import *
import os,time,traceback,copy,signal,argparse,math
from datetime import datetime,timezone,timedelta
import numpy as np
import torch
import train as trainer
from train import *
from benchmarks import mmff,conformers
from inference import audit_native
from robust_connection_data import make_official_diffuser

GRID=[0.,.02,.1,.3,.6,.9,.98]
SEEDS=[20261106,20261107,20261108]
def wall_limit(signum,frame):raise BudgetFinished('Confirmed 24-hour formal training/evaluation budget reached')
def validate_resume():
 """Resume only at a verified seed boundary, without resetting the budget."""
 start=load(OUT/'formal_start.json');assert sha(OUT/'protocol.json')==start['protocol_sha256'],'Protocol changed'
 allowed={'reports.py','run_experiment.py','launch.py'}
 changed=[name for name,h in start['source_sha256'].items() if sha(OUT/'source'/name)!=h]
 assert set(changed)<=allowed,f'Training/inference source changed: {changed}'
 completed=load(OUT/'completed_seeds.json');assert completed==SEEDS[:len(completed)],'Completion markers are not an ordered seed prefix'
 for seed in completed:
  run=OUT/'runs'/str(seed);evaluation=OUT/'evaluation'/str(seed)
  for group in ['A1','A2','A1_2d','A2_2d']:
   assert load(run/group/'history.json')[-1]['epoch']==40
   assert (run/group/'best.pt').is_file()
  for group in ['A3','A4']:
   assert load(run/group/'B_oracle/history.json')[-1]['epoch']==40
   assert (run/group/'inference_bundle.pt').is_file();load(run/group/'stage_decision.json')
  for split in ['test','scaffold_test']:
   n=load(OUT/'data'/f'{split}_summary.json')['usable_conformers']
   assert len(load(evaluation/f'MMFF_{split}.json'))==48
   for group in ['A0','A1','A2','A1_2d','A2_2d','A3','A4']:
    target=evaluation/group/split
    assert len(load(target/'connection.json'))==len(GRID)
    assert len(load(target/'port_cases.json'))==n*len(GRID)
    assert len(load(target/'searched_graphs.json'))==96
    if group in ['A3','A4']:
     assert len(load(target/'geometry_oracle.json'))==n*len(GRID)
     assert len(load(target/'geometry_predicted.json'))==96
     assert (target/'geometry_oracle_coordinates.npz').is_file()
     load(target/'conformer_diversity.json');load(target/'torsions.json')
 for seed in SEEDS[len(completed):]:
  assert not (OUT/'runs'/str(seed)).exists(),f'Partial seed {seed} requires checkpoint-level recovery'
  assert not (OUT/'evaluation'/str(seed)).exists(),f'Partial evaluation {seed}'
 origin=datetime.fromisoformat(start['UTC']);elapsed=(datetime.now(timezone.utc)-origin).total_seconds()
 assert elapsed>=0,'Clock moved before original start'
 if elapsed>=24*3600:raise BudgetFinished('Original 24-hour experiment deadline has passed; no automatic budget reset')
 return dict(completed_seeds=completed,elapsed_seconds=elapsed,remaining_seconds=24*3600-elapsed,deadline_UTC=(origin+timedelta(hours=24)).isoformat(),changed_supervisor_sources=changed)

def refresh_reports():
 """Record report failures without losing scientific work or the first error."""
 try:
  from reports import main as report
  report();return True
 except Exception as exc:
  failure=dict(UTC=datetime.now(timezone.utc).isoformat(),error=repr(exc),traceback=traceback.format_exc());dump(OUT/'report_failures'/f'{time.time_ns()}.json',failure);print('REPORT FAILED',json.dumps(failure),flush=True)
  status=load(OUT/'status.json') if (OUT/'status.json').exists() else {}
  if status.get('training_complete'):status.update(stage='report_failed',training_complete=False,training_and_evaluation_complete=True)
  status['last_report_error']=repr(exc);dump(OUT/'status.json',status);return False
def save_eval(seed,group,split,model,rows,diffuser,mode='3d'):
 target=OUT/'evaluation'/str(seed)/group/split;target.mkdir(parents=True,exist_ok=True);summary,cases=connection_evaluation(model,rows,diffuser,seed+950000,GRID,mode,details=True);dump(target/'connection.json',summary);dump(target/'port_cases.json',cases);return summary
def meanfield(rows,part,key):return float(np.mean([r[part][key] for r in rows if r.get('graph_available') and r[part][key] is not None]))
def main(resume=False):
 assert load(OUT/'preflight.json')['passed'] and load(OUT/'prelaunch_verification.json')['passed'];assert load(OUT/'dataset_verification.json')['passed'];assert load(OUT/'short_workflow_verification.json')['passed'];assert load(OUT/'data/complete.json')['complete']
 assert all(sha(p)==v['sha256'] for p,v in load(OUT/'audit_initial.json')['source_fingerprints'].items())
 if resume:
  recovery=validate_resume();completed=recovery['completed_seeds'];elapsed=recovery['elapsed_seconds'];dump(OUT/'resume_events'/f'{time.time_ns()}.json',dict(UTC=datetime.now(timezone.utc).isoformat(),pid=os.getpid(),**recovery,source_sha256={p.name:sha(p) for p in (OUT/'source').glob('*.py')}));print('RESUME',json.dumps(recovery),flush=True)
 else:
  assert not (OUT/'formal_start.json').exists(),'Formal run already exists; use --resume at a verified seed boundary'
  completed=[];elapsed=0.;dump(OUT/'formal_start.json',dict(UTC=datetime.now(timezone.utc).isoformat(),pid=os.getpid(),physical_GPU=os.environ.get('CUDA_VISIBLE_DEVICES'),wall_budget_hours=24,protocol_sha256=sha(OUT/'protocol.json'),source_sha256={p.name:sha(p) for p in (OUT/'source').glob('*.py')}))
 torch.set_num_threads(2);torch.cuda.set_device(0);trainer.START=time.monotonic()-elapsed;signal.signal(signal.SIGALRM,wall_limit);signal.alarm(max(1,math.ceil(24*3600-elapsed)));dump(OUT/'status.json',dict(stage='initializing_resumed_training' if resume else 'initializing_formal_training',training_complete=False,completed_seeds=completed));device='cuda:0';data={s:torch.load(OUT/'data'/f'{s}.pt',weights_only=False) for s in ['train','val','test','scaffold_test']};diff=make_official_diffuser(OUT/'cache/noise',device);finalbundle=OUT/'runs'/str(completed[-1])/'A4/inference_bundle.pt' if completed else None
 for seed in SEEDS:
  if seed in completed:print('SKIP VERIFIED COMPLETED SEED',seed,flush=True);continue
  budget();models={'A0':ConnectionHead(False).to(device)};histories={}
  baseline_val,_=connection_evaluation(models['A0'],choose_views(data['val'],256,1),diff,seed+810000,[.02,.3,.9]);baseline_score=float(np.mean([m['f1'] for m in baseline_val]));dump(OUT/'runs'/str(seed)/'A0_validation.json',baseline_val)
  for group in ['A1','A2','A1_2d','A2_2d']:
   models[group],histories[group]=train_classifier(group,seed,data['train'],data['val'],diff,device)
  # Collection remains finite and uses no true graph to repair predicted graphs.
  traincases,trainreport=collect_graphs(models['A2'],data['train'],diff,seed+920000,'train',128)
  valcases,valreport=collect_graphs(models['A2'],data['val'],diff,seed+920000,'val',32)
  torch.save(traincases,OUT/'runs'/str(seed)/'predicted_train_views.pt');torch.save(valcases,OUT/'runs'/str(seed)/'predicted_val_views.pt')
  for group in ['A3','A4']:
   geo,history=train_oracle_geometry(group,seed,data['train'],data['val'],diff,device);classifier=copy.deepcopy(models['A2']);valrows=choose_views(data['val'],256,1);oracle,_=geometry_evaluation(geo,valrows,diff,seed+810000,[.02,.3,.9]);dump(OUT/'runs'/str(seed)/group/'B_oracle/validation.json',oracle)
   before=meanfield(oracle,'before','bond_MAE_A');after=meanfield(oracle,'after','bond_MAE_A');a2_score=max(x.get('selection_score',-1) for x in histories['A2']);gateA=a2_score>baseline_score;gateB=after<before*.9
   if group=='A4':gateB=gateB and meanfield(oracle,'after','angle_MAE_deg')<meanfield(oracle,'before','angle_MAE_deg')
   # B must demonstrate independent recovery before adaptation to predicted graphs.
   c=adapt_predicted(group,seed,geo,classifier,traincases,valcases,device) if gateB else dict(skipped=True,reason='Oracle geometry failed predefined validation improvement gate')
   gateC=not c.get('skipped') and c.get('improved',False);d=adapt_predicted(group,seed,geo,classifier,traincases,valcases,device,joint=True) if gateA and gateB and gateC else dict(skipped=True,reason='A/B/C independent validation gates not all passed',gates=dict(A=gateA,B=gateB,C=gateC))
   bundle=OUT/'runs'/str(seed)/group/'inference_bundle.pt';torch.save(dict(classifier=classifier.state_dict(),geometry=geo.state_dict(),local_features=True,group=group,seed=seed,protocol_sha256=sha(OUT/'protocol.json'),baseline_sha256=sha(BASELINE)),bundle);finalbundle=bundle;models[group]=classifier
   dump(OUT/'runs'/str(seed)/group/'stage_decision.json',dict(A_validation_pass=gateA,B_validation_pass=gateB,C=c,D=d,bundle_sha256=sha(bundle)))
   for split in ['test','scaffold_test']:
    budget();target=OUT/'evaluation'/str(seed)/group/split;target.mkdir(parents=True,exist_ok=True);metrics,arrays=geometry_evaluation(geo,data[split],diff,seed+950000,GRID,save_coordinates=True);dump(target/'geometry_oracle.json',metrics);np.savez_compressed(target/'geometry_oracle_coordinates.npz',**arrays);cm,tm=conformers(data[split],arrays,group,seed,split);dump(target/'conformer_diversity.json',cm);dump(target/'torsions.json',tm)
  for group,model in models.items():
   for split in ['test','scaffold_test']:
    budget();mode='2d' if group.endswith('_2d') else '3d';save_eval(seed,group,split,model,data[split],diff,mode)
    cases,report=collect_graphs(model,data[split],diff,seed+950000,f'{group}_{split}',32,mode);target=OUT/'evaluation'/str(seed)/group/split;dump(target/'searched_graphs.json',report)
    if group in ['A3','A4']:
     weights=torch.load(OUT/'runs'/str(seed)/group/'inference_bundle.pt',map_location=device,weights_only=False);geo=GeometryHead().to(device);geo.load_state_dict(weights['geometry']);graphs={(r['id'],r['t']):r['graph'] for r in report if r.get('graph') is not None};metrics,_=geometry_evaluation(geo,choose_views(data[split],32,1),diff,seed+950000,[.02,.3,.9],graphs=graphs);dump(target/'geometry_predicted.json',metrics)
  for split in ['test','scaffold_test']:budget();dump(OUT/'evaluation'/str(seed)/f'MMFF_{split}.json',mmff(data[split],diff,seed+950000,split,device))
  completed.append(seed);dump(OUT/'completed_seeds.json',completed);print('SEED COMPLETE',seed,flush=True)
  # Refresh reports after each full seed, while retaining incomplete status for later seeds.
  refresh_reports();models.clear();torch.cuda.empty_cache()
 assert finalbundle is not None;integration=audit_native(finalbundle);assert integration['passed'];assert all(sha(p)==v['sha256'] for p,v in load(OUT/'audit_initial.json')['source_fingerprints'].items())
 dump(OUT/'status.json',dict(stage='complete',training_complete=True,completed_seeds=completed,elapsed_hours=(time.monotonic()-trainer.START)/3600,official_weights_unchanged=True,full_GEOM_training_started=False));print('FORMAL EXPERIMENT COMPLETE',flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--resume',action='store_true');parser.add_argument('--check-resume',action='store_true');args=parser.parse_args()
 if args.check_resume:print(json.dumps(validate_resume()),flush=True);sys.exit(0)
 try:main(resume=args.resume)
 except BudgetFinished as exc:dump(OUT/'status.json',dict(stage='budget_limit',training_complete=False,reason=str(exc),completed_seeds=load(OUT/'completed_seeds.json') if (OUT/'completed_seeds.json').exists() else []));print(str(exc),flush=True)
 except Exception as exc:dump(OUT/'status.json',dict(stage='failed',training_complete=False,error=repr(exc),traceback=traceback.format_exc()));raise
 finally:
  signal.alarm(0)
  refresh_reports()
