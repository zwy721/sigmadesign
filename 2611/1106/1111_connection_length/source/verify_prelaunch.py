"""Independent validation and deployment checks before formal model fitting."""
from common import *
import torch
from models import *
from metrics import *
from train import class_loss
from robust_connection_data import make_official_diffuser
from inference import audit_native

def main():
 torch.set_num_threads(2);assert load(OUT/'preflight.json')['passed'];data=load(OUT/'data/complete.json');assert all(sha(OUT/'data'/f'{s}.pt')==h for s,h in data['split_pt_sha256'].items())
 selection=load(OUT/'data/selected_molecules.json');sets=[{r['canonical'] for r in v} for v in selection.values()];assert all(not a&b for a,b in __import__('itertools').combinations(sets,2));assert not {r['scaffold'] for r in selection['scaffold_test']}&{r['scaffold'] for k,v in selection.items() if k!='scaffold_test' for r in v}
 rows=torch.load(OUT/'data_smoke/val.pt',weights_only=False);rows=choose_views(rows,4,1);b=pack(rows,'cuda:0');diff=make_official_diffuser(OUT/'cache/noise','cuda:0');times=torch.full((len(rows),),.3,device='cuda:0');T,R=noisy_frames(b,diff,times,seed=7171);saved=torch.load(OUT/'preflight_smoke.pt',map_location='cpu',weights_only=False);classifier=ConnectionHead(True).cuda();classifier.load_state_dict(saved['classifier']);geo=GeometryHead().cuda();geo.load_state_dict(saved['geometry']);checks={}
 f1=features(b,T,R,times,True,'2d');f2=features(b,T*7+3,R,times,True,'2d')
 with torch.no_grad():
  p,_=classifier(b,f1);q,_=classifier(b,f2);p=p[torch.argsort(f1['indices'])];q=q[torch.argsort(f2['indices'])];assert torch.allclose(p,q,atol=1e-5,rtol=1e-5)
 checks['2D_ablation_removes_dynamic_and_template_coordinate_channels']=True
 edges=edges_for(b,[r['connections'] for r in rows]);outT,outR=geo(b,T,R,times,edges);loss=class_loss(classifier,b,features(b,outT,outR,times,True));loss.backward();norm=sum(float(p.grad.norm()) for p in geo.parameters() if p.grad is not None);assert norm>1e-8;checks['joint_CE_gradient_can_reach_geometry_without_search_gradient']=norm
 with torch.no_grad():
  t1,r1=geo(b,T,R,times,edges);b['noise']=dict(b['noise'],T0=torch.full_like(b['noise']['T0'],173),R0=torch.zeros_like(b['noise']['R0']));t2,r2=geo(b,T,R,times,edges);error=max(float((t1-t2).abs().max()),float((r1-r2).abs().max()));assert torch.allclose(t1,t2,atol=1e-5,rtol=1e-5) and torch.allclose(r1,r2,atol=1e-5,rtol=1e-5),error
 checks['geometry_runtime_does_not_read_clean_reference_states']=dict(passed=True,max_GPU_roundoff=error,tolerance=1e-5)
 fulltrain=torch.load(OUT/'data/train.pt',weights_only=False);double=next(r for r in fulltrain if any(e['order']==2 and e['stereo']>0 for e in r['connections']));bd=pack([double],'cuda:0');td=torch.tensor([.15],device='cuda:0');dt,dr=noisy_frames(bd,diff,td,seed=7174);gd=GeometryHead().cuda();ot,orr=gd(bd,dt,dr,td,edges_for(bd,[double['connections']]));ld,parts=geometry_loss(bd,ot,orr,'A4');ld.backward();assert len(bd['torsion_index']) and torch.isfinite(ld) and any(p.grad is not None and torch.isfinite(p.grad).all() and p.grad.abs().sum()>0 for p in gd.parameters());checks['specified_double_bond_stereo_and_nonrotatable_torsion_gradient']=dict(passed=True,torsion_constraints=len(bd['torsion_index']),loss=float(ld));del fulltrain
 metric,_=connection_evaluation(classifier,rows,diff,7172,[.02,.3,.9],details=True);checks['independent_validation_metrics_computed']=metric
 one=rows[0];bb=pack([one],'cuda:0');tt,rr=noisy_frames(bb,diff,torch.tensor([.3],device='cuda:0'),seed=7173);graph,info=searched_graph(classifier,one,tt,rr,.3);assert not info['truth_used'];checks['bounded_search_does_not_consult_truth']=info
 bundle=OUT/'preflight_inference_bundle.pt';torch.save(dict(classifier=saved['classifier'],geometry=saved['geometry'],local_features=True),bundle)
 checks['native_interface']=audit_native(bundle);assert checks['native_interface']['passed'];dump(OUT/'prelaunch_verification.json',dict(passed=True,checks=checks));print('PRELAUNCH VERIFICATION PASSED',json.dumps(checks),flush=True)
if __name__=='__main__':main()
