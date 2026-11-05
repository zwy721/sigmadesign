"""Small explicitly new target-free adapter using frozen existing components.
Not a claim that the historical independent runner supported identity generation.
"""
from common import *
import os,subprocess,copy,time,traceback
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
PROJECT=ROOT.parent/'oracle_rebuild_20261031';sys.path.insert(0,str(PROJECT))
import bootstrap,torch
from rdkit import Chem
from scipy.spatial.transform import Rotation
from vocab_chemistry import VocabLibrary
from trajectory_connections import ConnectionRuntime,load_connection
from hierarchical_topology import hierarchical_search,SearchConfig
from joint_sampler import chemical_graph,NativeState
from run_validation import load_official,seed_all
from oracle_graph import assemble
from geometry_constraints_20261103.repair import repair
from geometry_constraints_20261103.stereo import derive
import importlib.util
module_spec=importlib.util.spec_from_file_location('old_local_optimizer',OLD/'source/optimizer.py');opt=importlib.util.module_from_spec(module_spec);sys.path.insert(0,str(OLD/'source'));module_spec.loader.exec_module(opt)
# common.py above is already imported; optimizer only uses ROOT and denovo imports.

def main():
 out=OUT/'free_probe';out.mkdir(exist_ok=True);torch.set_num_threads(1);device=torch.device('cuda:0');torch.cuda.set_device(device);torch.use_deterministic_algorithms(True);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
 weights=load(ROOT/'multi_candidates_20261103/plan.json')['weight_sha256'];identity_ck=ROOT.parent/'core497_joint_experiments_20261019/training3684/runs/geom_ipa3684_100ep_20261029_cuda01234_bs372_fast/best_loss.pt';cw=next(Path(p) for p in weights if 'training_v2' in p)
 plan=dict(cohort=sorted({r['receptor_id'] for r in load(OUT/'raw_cloud_manifest.json')}),seeds=[1801,1802,1803],attempts=18,slots=5,slot_policy='Fixed five slots for every receptor, never read true count.',identity='All MASK, GEOM matched existing checkpoint, top1 once at t=1; no target data or target lookup.',identity_checkpoint_sha256=sha256(identity_ck),connection_checkpoint_sha256=sha256(cw),official_checkpoint_sha256=sha256(bootstrap.CHECKPOINT),decoder='Existing hierarchical_search SearchConfig(), k=1, wide/full; no retries or quality selection.',pose='Sampled T normal sigma2.7 A around supplied pocket COM, seeded uniform R; official NativeState 25-point original cubic grid, 24 steps, noise0.',repair='Original rigid and frozen B1 each independently from native final, no water.',limitations='New minimal component adapter, not the original frozen native runner. GEOM identity does not read receptor. Pocket region and COM remain ligand-defined historical site. No variable fragment count, binding or prospective success claim.',training=False)
 if not (out/'plan.json').exists():write_json(out/'plan.json',plan)
 else:assert load(out/'plan.json')==plan
 den=load_official(bootstrap.CHECKPOINT,out/'official_cache',device);lib=VocabLibrary();head=load_connection(cw,device)
 log=(out/'identity_server.log').open('w');server=subprocess.Popen([str(ROOT.parent/'.venv-cutoff-train/bin/python'),'-B',str(OUT/'source/free_identity_server.py'),'--assets',str(out/'identity_assets'),'--device','cuda:1'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=log,text=True,bufsize=1)
 ready=server.stdout.readline();assert ready and json.loads(ready).get('ready'),ready
 results=[];manifest=[]
 try:
  for pid in plan['cohort']:
   # Load only context, then discard fixture labels/target; inference APIs never receive these.
   context=copy.deepcopy(torch.load(ROOT/f'fbdd_efficiency_generalization_20261103/independent_inputs/{pid}/fixture.pt',weights_only=False,map_location='cpu')['context']);context['pdb_path']=str(OUT/'pockets'/f'{pid}_no_water.pdb');context['ligand_path']=None
   for seed in plan['seeds']:
    folder=out/'structures'/f'{pid}_{seed}';folder.mkdir(parents=True,exist_ok=True);rec=dict(receptor_id=pid,seed=seed,candidate_rank=1,status='failed',failure_stage='identity',manifest=[])
    t0=time.perf_counter()
    try:
     rng=np.random.default_rng(seed);n=5;centers=rng.normal(0,2.7,(n,3))+np.asarray(context['pocket_com']).reshape(3);R=Rotation.random(n,random_state=seed).as_matrix()
     request=dict(T=centers.tolist(),R=R.tolist(),t=1.,seed=seed);server.stdin.write(json.dumps(request)+'\n');server.stdin.flush();response=json.loads(server.stdout.readline());assert 'error' not in response,response;labels=response['types'];assert len(labels)==5 and all(0<=v<3684 for v in labels);write_json(folder/'identity.json',response)
     templates=[lib.template(v)[1] for v in labels];T=np.array([centers[i]-x.mean(0)@R[i].T for i,x in enumerate(templates)]);cloud=np.concatenate([x@R[i].T+T[i] for i,x in enumerate(templates)]);z=[a.GetAtomicNum() for v in labels for a in lib.template(v)[0].GetAtoms()];np.savez_compressed(folder/'raw_cloud.npz',xyz=cloud,z=z,T=T,R=R)
     rec.update(labels=labels,identity_confidence=response['confidence'],failure_stage='connection_decode')
     rt=ConnectionRuntime(lib,head,labels,[None]*n,device,policy='wide',head='full');pred=rt.predict(T,R,1.,policy='wide',head='full');search=hierarchical_search(rt,pred,k=1,config=SearchConfig());assert not search['truth_used'];write_json(folder/'search.json',search)
     assert search['candidates'],'no legal graph';candidate=search['candidates'][0];rec['failure_stage']='graph_build';graph,groups,spec=chemical_graph(lib,labels,candidate['connections'],context,[None]*n,seed=seed)
     seed_all(seed);state=NativeState(den,graph,groups,device,desired=[x@R[i].T+T[i] for i,x in enumerate(templates)],preferred=R);flat=[a for g in groups for a in g];mol=assemble(spec);initial=np.empty_like(cloud);initial[flat]=cloud
     def save_stage(stage,m,sp):
      d=folder/stage;sdf(d/'raw.sdf',m);q,_,_=__import__('denovo_quality.io',fromlist=['load_molecule']).load_molecule(dict(record_index=0,sdf_path=str(d/'raw.sdf'),chemistry_spec_path=str(folder/'original_spec.json')));q,ds,st=derive(q,sp);sdf(d/'derived.sdf',q);write_json(d/'spec.json',ds);write_json(d/'stereo_decisions.json',st)
      mr=dict(molecule_id=f'free_{stage}__{pid}_{seed}',group=stage,receptor_id=pid,seed=seed,candidate_rank=1,sdf_path=str(d/'derived.sdf'),chemistry_spec_path=str(d/'spec.json'),protein_path=context['pdb_path'],record_index=0);manifest.append(mr);rec['manifest'].append(mr)
     write_json(folder/'original_spec.json',spec)
     for a,x in enumerate(initial):mol.GetConformer().SetAtomPosition(a,x.tolist())
     save_stage('decoded_initial',mol,spec);rec['failure_stage']='native_sampling'
     grid=(1.+torch.arange(25,device=device)/24*(.005**(1/3)-1.))**3
     with torch.inference_mode():
      for step,t in enumerate(grid[:-1]):state.step(t,grid[step]-grid[step+1])
     xyz=state.physical()[state.pn:state.pn+sum(map(len,groups))]
     for a,x in enumerate(xyz):mol.GetConformer().SetAtomPosition(a,x.tolist())
     save_stage('native_final',mol,spec);rec['neural_NFE']=24;rec['failure_stage']='repair'
     fids=np.empty(len(xyz),dtype=int)
     for i,g in enumerate(groups):fids[g]=i
     for stage,fun in [('rigid',lambda:repair(mol,fids,context['pdb_path'])),('local_flex',lambda:opt.repair(mol,fids,context['pdb_path'],local=load(OLD/'fixed_parameters.json')['local']))]:
      q,meta=fun();write_json(folder/stage/'repair.json',meta);save_stage(stage,q,spec)
     rec.update(status='ok',failure_stage=None,groups=groups,source='Unconditional identities; fixed user slot count; site-conditioned neural pose; searched graph.')
    except Exception as exc:rec.update(error=repr(exc),traceback=traceback.format_exc())
    rec['seconds']=time.perf_counter()-t0;write_json(folder/'result.json',rec);results.append(rec);write_json(out/'attempts.json',results);write_json(out/'manifest.json',manifest);print('FREE',pid,seed,rec['status'],rec.get('failure_stage'),rec.get('error',''),flush=True)
 finally:
  if server.poll() is None:server.stdin.write('{"close":true}\n');server.stdin.flush();server.wait(timeout=30)
  log.close()
 write_json(out/'weight_integrity.json',dict(unchanged=sha256(identity_ck)==plan['identity_checkpoint_sha256'] and sha256(cw)==plan['connection_checkpoint_sha256'] and sha256(bootstrap.CHECKPOINT)==plan['official_checkpoint_sha256']))
if __name__=='__main__':main()
