"""Existing GEOM identity weights; inference entry made solely from supplied poses.
No dataset, ligand record, target identity, adjacency, or target coordinate lookup.
The feature branch intentionally matches existing identity_server geom_matched.
"""
import os,sys,json,contextlib,argparse,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];PROJECT=ROOT.parent/'oracle_rebuild_20261031';LEGACY=ROOT.parent/'core497_joint_experiments_20261019/training3684'
def main():
 p=argparse.ArgumentParser();p.add_argument('--assets',type=Path,required=True);p.add_argument('--device',default='cuda:1');a=p.parse_args()
 import torch
 torch.set_num_threads(1)
 ck=LEGACY/'runs/geom_ipa3684_100ep_20261029_cuda01234_bs372_fast/best_loss.pt';saved=torch.load(ck,map_location='cpu',weights_only=False);cfg=saved['config'];cfg['joint_options']['connection_graph']='none';cfg['train']['packed_geom']=False;cfg['train']['prepared_geom']=False
 a.assets.mkdir(parents=True,exist_ok=True);cp=a.assets/'config.json';cp.write_text(json.dumps(cfg));os.environ['CUTOFF_TRAIN_CONFIG']=str(cp)
 sys.path=[str(LEGACY)]+[x for x in sys.path if Path(x or '.').resolve() not in {PROJECT,Path(__file__).resolve().parent}]
 with contextlib.redirect_stdout(sys.stderr):
  import runtime
  for nm in ['mapping.json','templates.json','properties.json']:
   if not (a.assets/nm).exists():shutil.copyfile(runtime.ASSETS/nm,a.assets/nm)
  if not (a.assets/'diffusion_cache').exists():shutil.copytree(runtime.ASSETS/'diffusion_cache',a.assets/'diffusion_cache')
  runtime.ASSETS=a.assets
  from data import input_graph
  from geom_data import initialize
  from model3684 import ConnectionIPA3684
  model=ConnectionIPA3684(a.device);model.load_state_dict(saved['model'],strict=True);model.eval();model.requires_grad_(False);pg=initialize()['pg']
 print(json.dumps(dict(ready=True,K=cfg['K'],checkpoint=str(ck))),flush=True)
 for line in sys.stdin:
  try:
   r=json.loads(line)
   if r.get('close'):break
   assert set(r)<= {'T','R','t','seed','close'}
   with contextlib.redirect_stdout(sys.stderr):
    T=torch.tensor(r['T'],dtype=torch.float32);T-=T.mean(0);R=torch.tensor(r['R'],dtype=torch.float32);n=len(T);tokens=torch.full((n,),cfg['MASK'],dtype=torch.long)
    e=dict(id='target_free',source='geom',pg=pg,origin=torch.zeros(3),tokens=tokens,labels=tokens.clone(),Tt=T,Rt=R,t=float(r['t']),seed=r['seed'],mode='all_mask',oov_templates=[None]*n,pose_supervision=False,geometry_training=False,edges=torch.zeros(n,n,dtype=torch.long),anchors=torch.full((n,n,2),-100,dtype=torch.long),stereo=torch.full((n,n),-100,dtype=torch.long))
    e['graph']=input_graph(e)
    assert not e['graph']['case']['native_anchor_pairs'] and not len(pg['pos'])
    with torch.inference_mode():prob=model._run([e])[0]['logits'].softmax(-1).cpu();values,indices=prob.topk(5,-1)
   print(json.dumps(dict(types=indices[:,0].tolist(),confidence=values[:,0].tolist(),top5=indices.tolist(),K=cfg['K'],protein_conditioned_identity=False,target_information_used=False)),flush=True)
  except Exception as exc:
   import traceback;traceback.print_exc(file=sys.stderr);print(json.dumps(dict(error=repr(exc))),flush=True)
if __name__=='__main__':main()
