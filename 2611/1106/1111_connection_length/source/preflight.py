"""Executable gradient, rigidity, equivariance, leakage, and tiny-overfit checks."""
from common import *
import time,traceback
import numpy as np
import torch
from torch.nn import functional as F
from scipy.spatial.transform import Rotation
from models import *
from robust_connection_data import make_official_diffuser

def f1(logits,y):
 p=logits.argmax(-1);tp=int(((p==y)&(y>0)).sum());return 2*tp/max(1,int((p>0).sum())+int((y>0).sum()))

def main():
 torch.set_num_threads(2);torch.manual_seed(20261106);began=time.time();device='cuda:0'
 rows=torch.load(OUT/'data_smoke/train.pt',weights_only=False);unique={}
 for r in rows:unique.setdefault(r['molecule_id'],r)
 rows=list(unique.values());b=pack(rows,device);diff=make_official_diffuser(OUT/'cache/noise',device)
 times=torch.linspace(.02,.6,len(rows),device=device);T,R=noisy_frames(b,diff,times,seed=8711);f=features(b,T,R,times,local=True);y=b['y'][f['indices']].long();st=b['stereo'][f['indices']].long();assert int((y>0).sum())==sum(r['target']['true_bonds'] for r in rows)
 checks={};checks['all_original_positive_ports_retained']=dict(passed=True,positive_count=int((y>0).sum()),all_pairs=len(y))
 original_y=b['y'].clone();b['y'].zero_();f_no_target=features(b,T,R,times,local=True);assert torch.equal(f['x'],f_no_target['x']) and torch.equal(f['extra'],f_no_target['extra']);b['y']=original_y;checks['classifier_features_independent_of_targets']=dict(passed=True)
 classifier=ConnectionHead(local=True).to(device);opt=torch.optim.AdamW([p for p in classifier.parameters() if p.requires_grad],lr=.002,weight_decay=.0001);counts=torch.bincount(y,minlength=3).float();w=(counts.max()/counts.clamp_min(1)).sqrt().clamp(max=30)
 initial_logits,_=classifier(b,f);initial_f1=f1(initial_logits,y)
 # Fixed noisy poses are intentional for the memorization diagnostic only.
 for step in range(1200):
  opt.zero_grad(set_to_none=True);logits,s=classifier(b,f);loss=F.cross_entropy(logits,y,weight=w)
  double=y==2
  if double.any():loss=loss+.2*F.cross_entropy(s[double],st[double])
  loss.backward();torch.nn.utils.clip_grad_norm_(classifier.parameters(),5.);opt.step()
  if (step+1)%100==0:print('connection-overfit',step+1,float(loss),f1(logits,y),flush=True)
  if step>=199 and f1(logits,y)>=.98:break
 with torch.no_grad():logits,_=classifier(b,f);final_f1=f1(logits,y)
 checks['tiny_connection_overfit']=dict(passed=final_f1>=.9,initial_f1=initial_f1,final_f1=final_f1,steps=step+1,molecules=len(rows),noise_times=times.cpu().tolist())
 assert final_f1>=.9,checks['tiny_connection_overfit']
 # Geometry checks use eight distinct training molecules, with true graphs explicitly Oracle.
 b=pack(rows[:8],device);times=torch.full((8,),.15,device=device);T,R=noisy_frames(b,diff,times,seed=8712);edges=edges_for(b,[r['connections'] for r in b['rows']]);geo=GeometryHead().to(device)
 outT,outR=geo(b,T,R,times,edges);loss,parts=geometry_loss(b,outT,outR,'A3');loss.backward();norm=sum(float(p.grad.norm()) for p in geo.parameters() if p.grad is not None);checks['bond_loss_to_geometry_parameters_nonzero']=dict(passed=norm>1e-6,gradient_norm=norm);assert norm>1e-6
 initial_bond=float(parts['bond']);initial_full=float(geometry_loss(b,T,R,'A4')[0]);opt=torch.optim.AdamW(geo.parameters(),lr=.0005,weight_decay=.0001)
 for step in range(1000):
  opt.zero_grad(set_to_none=True);outT,outR=geo(b,T,R,times,edges);loss,parts=geometry_loss(b,outT,outR,'A4');loss.backward();torch.nn.utils.clip_grad_norm_(geo.parameters(),5.);opt.step()
  if (step+1)%100==0:print('geometry-overfit',step+1,float(loss),{k:float(v) for k,v in parts.items()},flush=True)
  if step>=299 and float(parts['bond'])<.01 and float(loss)<.05:break
 with torch.no_grad():
  outT,outR=geo(b,T,R,times,edges);full,parts=geometry_loss(b,outT,outR,'A4');checks['tiny_geometry_overfit']=dict(passed=float(parts['bond'])<initial_bond*.2 and float(full)<initial_full*.2,initial_bond_mse=initial_bond,final_bond_mse=float(parts['bond']),initial_full_loss=initial_full,final_full_loss=float(full),steps=step+1)
  before=cloud(b,T,R).double();after=cloud(b,outT,outR).double();valid=b['atom_mask'][:,:,None]&b['atom_mask'][:,None,:];rigid_error=float((torch.cdist(before,before,compute_mode='donot_use_mm_for_euclid_dist')-torch.cdist(after,after,compute_mode='donot_use_mm_for_euclid_dist')).abs()[valid].max());checks['rigid_internal_geometry_preserved']=dict(passed=rigid_error<1e-4,max_distance_change_A=rigid_error,distance_evaluation='FP64 direct pair differences; no diagonal cancellation from Gram matrix')
  change=float((geometry_values(b,T,R)[1]-geometry_values(b,outT,outR)[1]).abs().max());checks['network_actually_changes_bond_lengths']=dict(passed=change>1e-3,max_change_A=change)
  G=torch.tensor(Rotation.from_rotvec([.3,-.2,.7]).as_matrix(),device=device,dtype=T.dtype);shift=T.new_tensor([3.,-2.,1.]);gt,gr=geo(b,T@G.T+shift,G[None]@R,times,edges);errorT=float((gt-(outT@G.T+shift)).abs().max());errorR=float((gr-G[None]@outR).abs().max());checks['global_SE3_equivariance']=dict(passed=errorT<2e-4 and errorR<2e-4,max_T_error_A=errorT,max_R_error=errorR)
 assert all(v['passed'] for v in checks.values()),checks
 frozen=load(OUT/'audit_initial.json')['source_fingerprints'];assert all(sha(p)==h['sha256'] for p,h in frozen.items());checks['upstream_sources_and_weights_unchanged']=dict(passed=True,files=len(frozen))
 torch.save(dict(classifier=classifier.state_dict(),geometry=geo.state_dict()),OUT/'preflight_smoke.pt');report=dict(passed=True,checks=checks,elapsed_seconds=time.time()-began,training_scope='Fixed tiny train-only batches; no validation/test use; not formal performance')
 dump(OUT/'preflight.json',report);print('PREFLIGHT PASSED',json.dumps(report),flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:dump(OUT/'preflight_failure.json',dict(error=repr(exc),traceback=traceback.format_exc()));raise
