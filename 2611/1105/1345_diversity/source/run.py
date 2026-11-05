"""Frozen target-free strategies, current connection head, official SE(3) steps.

One attempt = one pocket/seed/strategy, with no resampling after any failure.
All complete attempts have five identity queries, five bounded graph searches,
and 24 official force evaluations. Non-updating identity queries are observers.
"""
from common import *
import os,copy,time,argparse,subprocess,traceback,importlib.util
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'matplotlib'))
PROJECT=ROOT.parent/'oracle_rebuild_20261031';sys.path.insert(0,str(PROJECT))
import bootstrap,torch
from rdkit import Chem,RDLogger
from scipy.spatial.transform import Rotation
from vocab_chemistry import VocabLibrary,ChemistryRejected
from trajectory_connections import ConnectionRuntime,load_connection
from hierarchical_topology import hierarchical_search,SearchConfig
from joint_sampler import chemical_graph,NativeState,topology_key
from run_validation import load_official,seed_all
from validation import SamplingAudit
from oracle_graph import assemble
from sigmadock.chem.parsing import get_protein,structure_to_rdkit
from geometry_constraints_20261103.repair import repair
from geometry_constraints_20261103.stereo import derive
sp=importlib.util.spec_from_file_location('frozen_local_optimizer',OLD/'source/optimizer.py');opt=importlib.util.module_from_spec(sp);sp.loader.exec_module(opt)
STRATEGIES={
 'mask_top1':dict(initial='mask',draw='top1',temperature=1.,iterative=False),
 'mask_categorical':dict(initial='mask',draw='categorical',temperature=1.,iterative=False),
 'mask_temp15':dict(initial='mask',draw='categorical',temperature=1.5,iterative=False),
 'mask_topk50':dict(initial='mask',draw='topk',top_k=50,temperature=1.,iterative=False),
 'mask_topp095':dict(initial='mask',draw='topp',top_p=.95,temperature=1.,iterative=False),
 'joint_prior_iterative':dict(initial='frequency_prior',draw='top1',temperature=1.,iterative=True,threshold=.3),
 'alternating_categorical':dict(initial='frequency_prior',draw='categorical',temperature=1.,iterative=True,threshold=0.),
}
def softmax(x):
 z=x-x.max(axis=-1,keepdims=True);p=np.exp(z);return p/p.sum(-1,keepdims=True)
def draw_distribution(logits,policy):
 p=softmax(logits/policy['temperature']);q=p.copy()
 if policy['draw']=='top1':q[:]=0;q[np.arange(len(q)),p.argmax(-1)]=1.
 elif policy['draw']=='topk':
  ids=np.argsort(-p,axis=-1,kind='stable')[:,policy['top_k']:];np.put_along_axis(q,ids,0,axis=-1);q/=q.sum(-1,keepdims=True)
 elif policy['draw']=='topp':
  ids=np.argsort(-p,axis=-1,kind='stable');rankp=np.take_along_axis(p,ids,-1);keep=rankp.cumsum(-1)-rankp<policy['top_p'];ordered=rankp*keep;q[:]=0;np.put_along_axis(q,ids,ordered,-1);q/=q.sum(-1,keepdims=True)
 return p,q
def entropy(p):return -(p*np.log(np.maximum(p,1e-300))).sum(-1)
def weight_paths():
 return dict(identity=ROOT.parent/'core497_joint_experiments_20261019/training3684/runs/geom_ipa3684_100ep_20261029_cuda01234_bs372_fast/best_loss.pt',connection=PROJECT/'results/connection_robust_20261101/training_v2/best.pt',pose=bootstrap.CHECKPOINT)
def make_plan():
 paths=weight_paths()
 plan=dict(version=1,training=False,gpu_count=1,visible_gpu=os.environ.get('CUDA_VISIBLE_DEVICES'),cohort=['12AV','7ILZ','9RT3','9SYL','9T4A','9T51'],cohort_role='Previously examined diagnostic pockets only; no generalization claim',pilot_pockets=['12AV'],pilot_seeds=[40001,40002],seeds=list(range(41001,41013)),strategies=STRATEGIES,N_prior=dict(values=[3,4,5,6],probabilities=[.25]*4,source='Predefined uniform diagnostic prior; not fit to target or training size',seed_rule='numpy.default_rng(seed+912711).choice(values)'),initial_pose=dict(center='Centroid of supplied model protein heavy atoms; no ligand COM read',sigma_A=2.7,rotations='scipy Rotation.random(seed)'),pose_steps=24,grid='Official 25-point cubic t=1 to .005, noise_scale=0',revision_steps=[5,10,15,20],search_config=SearchConfig().__dict__,search_k=1,connection_policy='wide',connection_head='full TimeConnectionPredictor; never legacy_staged',budget='5 identity queries, 5 connection forward/search calls per full attempt; 24 pose steps; no retries',identity='All-MASK current relative poses, no receptor, no connection-prior channel; source geom',joint_adaptation='Reuse frequency initialization, .3 confidence updates, legal fallback and shape transfer from sample_joint. Replace early Port head with CURRENT TimeConnectionPredictor + identical hierarchical search in all arms. Fixed current random-frame initialization shared by all arms.',topology_policy='Each arm proposes a graph at identical update steps. One-shot arms freeze identity but can update graph. Identity updates and connections accepted together only after legal graph+reference+shape validation.',repair='Frozen rigid and local-flex B1 independently from native_final; no water',QC_config_sha256=sha256(QC),PB_config_sha256=sha256(OUT/'posebusters_no_water.json'),weights={k:dict(path=str(v),sha256=sha256(v)) for k,v in paths.items()},history_baseline_path=str(ROOT/'2611/1105/1039+no_water_test/free_probe'),selection='All settings a priori; no sweep, deduplication, quality filtering, test-selected parameters, or training',geometry_transfer_max_fit_A=.5)
 path=OUT/'frozen_plan.json'
 if path.exists():assert load(path)==plan,'Frozen plan changed'
 else:write_json(path,plan)
 # Protect all previous task files, not just the old summary.
 old=ROOT/'2611/1105/1039+no_water_test';finger={str(p):sha256(p) for p in old.rglob('*') if p.is_file()}
 fp=OUT/'historical_fingerprints.json'
 if fp.exists():assert load(fp)==finger
 else:write_json(fp,finger)
 return plan
class Client:
 def __init__(self,worker='pilot'):
  self.log=(OUT/f'identity_server_{worker}.log').open('a');self.proc=subprocess.Popen([sys.executable,'-B',str(OUT/'source/identity_server.py'),'--assets',str(OUT/'identity_assets')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,bufsize=1)
  line=self.proc.stdout.readline();assert line and json.loads(line).get('ready'),line
 def query(self,T,R,t,seed,folder,step,policy):
  self.proc.stdin.write(json.dumps(dict(T=T.tolist(),R=R.tolist(),t=float(t),seed=seed))+'\n');self.proc.stdin.flush();line=self.proc.stdout.readline();assert line,'identity subprocess exited';response=json.loads(line);assert 'error' not in response,response
  logits=np.asarray(response['logits'],dtype=np.float64);raw=softmax(logits);p,q=draw_distribution(logits,policy)
  np.savez_compressed(folder/f'identity_{step:02d}.npz',logits=logits.astype('float32'),model_probabilities=raw.astype('float32'),sampling_probabilities=q.astype('float32'),T=T,R=R)
  rank=np.argsort(-raw,axis=-1,kind='stable')[:,:10]
  stats=dict(step=step,t=float(t),entropy=entropy(raw).tolist(),sampling_entropy=entropy(q).tolist(),top1_probability=raw.max(-1).tolist(),top1_class=raw.argmax(-1).tolist(),benzene_probability=raw[:,3577].tolist(),top10=rank.tolist(),top10_probability=np.take_along_axis(raw,rank,-1).tolist(),protein_conditioned=False,connection_prior_used=False)
  return raw,q,stats
 def close(self):
  if self.proc.poll() is None:self.proc.stdin.write('{"close":true}\n');self.proc.stdin.flush();self.proc.wait(timeout=30)
  self.log.close()
def context_from_pocket(pid):
 path=ROOT/'2611/1105/1039+no_water_test/pockets'/f'{pid}_no_water.pdb'
 pocket=structure_to_rdkit(get_protein(path));xyz=pocket.GetConformer().GetPositions();assert all(a.GetPDBResidueInfo().GetResidueName().strip() not in WATER for a in pocket.GetAtoms())
 return dict(pocket=pocket,pocket_com=xyz.mean(0).astype('float32'),pdb_path=str(path),ligand_path=None)
def save_stage(folder,stage,mol,spec,rec,context):
 d=folder/stage;sdf(d/'raw.sdf',mol);write_json(d/'original_spec.json',spec)
 from denovo_quality.io import load_molecule
 q,_,_=load_molecule(dict(record_index=0,sdf_path=str(d/'raw.sdf'),chemistry_spec_path=str(d/'original_spec.json')))
 q,ds,st=derive(q,spec);sdf(d/'derived.sdf',q);write_json(d/'spec.json',ds);write_json(d/'stereo_decisions.json',st)
 rec['manifest'].append(dict(molecule_id=f'{rec["attempt_id"]}__{stage}',group=stage,receptor_id=rec['receptor_id'],strategy=rec['strategy'],seed=rec['seed'],phase=rec['phase'],sdf_path=str(d/'derived.sdf'),chemistry_spec_path=str(d/'spec.json'),protein_path=context['pdb_path'],record_index=0))
def frames_for(lib,labels,centers,R):
 templates=[lib.template(v)[1] for v in labels];T=np.array([centers[i]-x.mean(0)@R[i].T for i,x in enumerate(templates)]);return T,[x@R[i].T+T[i] for i,x in enumerate(templates)]
def attempt(den,lib,head,client,context,pid,seed,strategy,phase,plan):
 policy=STRATEGIES[strategy];aid=f'{phase}_{pid}_{seed}_{strategy}';folder=OUT/phase/'structures'/aid;folder.mkdir(parents=True,exist_ok=True)
 record_path=folder/'result.json'
 if record_path.exists():return load(record_path)
 started=time.perf_counter();n=int(np.random.default_rng(seed+912711).choice(plan['N_prior']['values']));rng=np.random.default_rng(seed);rel_centers=rng.normal(0,2.7,(n,3));centers=rel_centers+context['pocket_com'];R=Rotation.random(n,random_state=seed).as_matrix();id_rng=np.random.default_rng(seed+17331)
 rec=dict(attempt_id=aid,phase=phase,receptor_id=pid,strategy=strategy,seed=seed,N=n,status='failed',failure_stage='identity',manifest=[],identity_queries=0,connection_queries=0,pose_NFE=0,events=[],identity_statistics=[],repair_errors={})
 try:
  raw,q,stats=client.query(rel_centers,R,1.,seed,folder,0,policy);rec['identity_queries']+=1;rec['identity_statistics'].append(stats)
  if policy['initial']=='frequency_prior':
   prior=np.array([p['train_occurrences']*(p['n_heavy']>=2 and bool(p['candidate_sites'])) for p in lib.properties],dtype=float);prior/=prior.sum();labels=id_rng.choice(len(prior),n,p=prior).tolist()
  else:labels=[int(id_rng.choice(len(p),p=p)) for p in q]
  rec.update(initial_labels=labels.copy(),labels=labels.copy(),failure_stage='connection_decode');T,desired=frames_for(lib,labels,centers,R)
  np.savez_compressed(folder/'raw_cloud.npz',xyz=np.concatenate(desired),T=T,R=R,z=[a.GetAtomicNum() for v in labels for a in lib.template(v)[0].GetAtoms()])
  def search_for(types,tt,rr,t,step):
   rt=ConnectionRuntime(lib,head,types,[None]*n,'cuda:0',policy='wide',head='full');pred=rt.predict(tt,rr,float(t),policy='wide',head='full');rec['connection_queries']+=1
   search=hierarchical_search(rt,pred,k=1,config=SearchConfig());assert not search['truth_used'];search['raw_head_complete']=pred['assembled'] is not None;search['raw_complete_chemical_graph']=pred['raw_complete_chemical_graph'];write_json(folder/f'search_{step:02d}.json',search)
   assert search['candidates'],'no legal complete graph within frozen search budget'
   return search['candidates'][0]['connections']
  connections=search_for(labels,T,R,1.,0);rec['failure_stage']='graph_build';graph,groups,spec=chemical_graph(lib,labels,connections,context,[None]*n,seed=seed)
  mol=assemble(spec)
  for g,coords in zip(groups,desired):
   for a,x in zip(g,coords):mol.GetConformer().SetAtomPosition(a,x.tolist())
  save_stage(folder,'decoded_initial',mol,spec,rec,context)
  seed_all(seed);state=NativeState(den,graph,groups,'cuda:0',desired=desired,preferred=R)
  rec['initial_transfer_fit_rmsd_A']=state.transfer_errors
  rec['failure_stage']='native_sampling';grid=(1.+torch.arange(25,device='cuda:0')/24*(.005**(1/3)-1.))**3
  with torch.inference_mode(),SamplingAudit(den) as audit:
   for step,t in enumerate(grid[:-1]):
    if step in plan['revision_steps']:
     tt,rr,fit=state.slot_frames(lib,labels,[None]*n);raw,q,stats=client.query(tt-context['pocket_com'],rr,float(t),seed,folder,step,policy);rec['identity_queries']+=1;rec['identity_statistics'].append(stats)
     proposed=labels.copy()
     if policy['iterative']:
      guesses=[int(id_rng.choice(len(p),p=p)) for p in q]
      for i,g in enumerate(guesses):
       if raw[i].max()>=policy['threshold']:proposed[i]=g
     event=dict(step=step,t=float(t),old_labels=labels.copy(),proposed_labels=proposed,accepted=False)
     try:
      nc=search_for(proposed,tt,rr,t,step)
      if topology_key(proposed,nc)==topology_key(labels,connections):event['reason']='same_chemistry'
      else:
       ng,gg,ss=chemical_graph(lib,proposed,nc,context,[None]*n,seed=seed);previous=state.physical();wanted=[]
       for i,label in enumerate(proposed):
        old_actual=previous[state.pn+np.asarray(state.groups[i])]
        if label==labels[i]:wanted.append(old_actual.copy())
        else:
         _,x=lib.template(label);wanted.append((x-x.mean(0))@rr[i].T+old_actual.mean(0))
       replacement=NativeState(den,ng,gg,'cuda:0',desired=wanted,preferred=rr);event['transfer_fit_rmsd_A']=replacement.transfer_errors
       if max(replacement.transfer_errors,default=0)>.5:raise ChemistryRejected('Reference transfer fit exceeds frozen 0.5 A bound')
       state,graph,groups,spec=replacement,ng,gg,ss;labels,connections=proposed,nc;audit.static_hash=None;event['accepted']=True
     except (ChemistryRejected,ValueError,RuntimeError,AssertionError) as exc:
      if isinstance(exc,torch.OutOfMemoryError):raise
      event.update(reason=repr(exc),fallback='retain_last_legal_graph')
     rec['events'].append(event)
    state.step(t,grid[step]-grid[step+1]);rec['pose_NFE']+=1
   state.batch=den._update_batch(state.batch,state.pos0,state.pos)
  write_json(folder/'sampling_audit.json',audit.report());xyz=state.physical()[state.pn:state.pn+sum(map(len,state.groups))];mol=assemble(spec)
  for a,x in enumerate(xyz):mol.GetConformer().SetAtomPosition(a,x.tolist())
  rec.update(labels=labels,connections=connections,groups=groups,neural_success=True,failure_stage=None)
  save_stage(folder,'native_final',mol,spec,rec,context)
  np.savez_compressed(folder/'model_pocket.npz',xyz=graph.ref_pos[:state.pn].cpu().numpy(),center=graph.pocket_com.cpu().numpy())
  fids=np.empty(len(xyz),dtype=int)
  for i,g in enumerate(groups):fids[g]=i
  for stage,fun in [('rigid',lambda:repair(mol,fids,context['pdb_path'])),('local_flex',lambda:opt.repair(mol,fids,context['pdb_path'],local=load(OLD/'fixed_parameters.json')['local']))]:
   try:
    fixed,meta=fun();write_json(folder/stage/'repair.json',meta);save_stage(folder,stage,fixed,spec,rec,context)
   except Exception as exc:rec['repair_errors'][stage]=dict(error=repr(exc),traceback=traceback.format_exc())
  rec['status']='ok' if not rec['repair_errors'] else 'repair_incomplete'
 except Exception as exc:rec.update(error=repr(exc),traceback=traceback.format_exc())
 rec['seconds']=time.perf_counter()-started;write_json(record_path,rec);return rec
def main():
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=['pilot','main'],required=True);p.add_argument('--shard',type=int,default=0);p.add_argument('--world',type=int,default=1);a=p.parse_args();assert 0<=a.shard<a.world;torch.set_num_threads(1);torch.cuda.set_device(0);torch.use_deterministic_algorithms(True);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False;RDLogger.DisableLog('rdApp.*')
 plan=make_plan();phase=OUT/a.phase;phase.mkdir(exist_ok=True);lib=VocabLibrary();den=load_official(bootstrap.CHECKPOINT,OUT/('official_cache' if a.world==1 else f'official_cache_{a.shard}'),'cuda:0');head=load_connection(weight_paths()['connection'],'cuda:0');assert type(head['runtime_model']).__name__=='TimeConnectionPredictor'
 suffix='' if a.world==1 else f'_shard_{a.shard}';client=Client(a.phase+suffix);results=[];case_index=0
 try:
  for pid in plan['pilot_pockets'] if a.phase=='pilot' else plan['cohort']:
   context=context_from_pocket(pid)
   for seed in plan['pilot_seeds'] if a.phase=='pilot' else plan['seeds']:
    selected=case_index%a.world==a.shard;case_index+=1
    if not selected:continue
    for strategy in STRATEGIES:
     r=attempt(den,lib,head,client,context,pid,seed,strategy,a.phase,plan);results.append(r);write_json(phase/f'attempts{suffix}.json',results);write_json(phase/f'manifest{suffix}.json',[m for row in results for m in row['manifest']]);print(a.phase,len(results),pid,seed,strategy,r['N'],r['status'],r.get('error',''),round(r['seconds'],2),flush=True)
 finally:client.close()
 write_json(phase/f'weight_integrity{suffix}.json',{k:sha256(v)==plan['weights'][k]['sha256'] for k,v in weight_paths().items()})
if __name__=='__main__':main()
