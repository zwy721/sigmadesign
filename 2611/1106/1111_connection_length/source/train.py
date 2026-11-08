"""Controlled A0/A1/A2/2D ablations and Oracle -> predicted-graph geometry."""
from common import *
import copy,time,random
import numpy as np
import torch
from torch.nn import functional as F
from models import *
from metrics import *

class BudgetFinished(Exception):pass
START=None
def budget():
 if time.monotonic()-START>24*3600:raise BudgetFinished('Confirmed 24-hour formal training budget reached')
def checkpoint(path,model,epoch,score,history,optimizer=None,**extra):
 payload=dict(model=model.state_dict(),epoch=epoch,selection_score=score,history=history,protocol_sha256=sha(OUT/'protocol.json'),source_sha256={p.name:sha(p) for p in (OUT/'source').glob('*.py')},baseline_sha256=sha(BASELINE),**extra)
 if optimizer is not None:payload['optimizer']=optimizer.state_dict()
 path.parent.mkdir(parents=True,exist_ok=True);torch.save(payload,path)
def epoch_rows(rows,epoch,seed):
 groups={}
 for r in rows:groups.setdefault(r['molecule_id'],[]).append(r)
 ids=list(groups);np.random.default_rng(seed+epoch).shuffle(ids)
 return [groups[k][epoch%len(groups[k])] for k in ids]
def noise_times(n,device,seed):
 gen=torch.Generator(device=device).manual_seed(seed);u=torch.rand(n,device=device,generator=gen);v=torch.rand(n,device=device,generator=gen)
 return torch.where(u<.25,.0001+v*.0499,torch.where(u<.60,.05+v*.35,.4+v*.60))
def class_loss(model,b,f):
 logits,stereo=model(b,f);y=b['y'][f['indices']].long();counts=torch.bincount(y,minlength=3).float();weights=(counts.max()/counts.clamp_min(1)).sqrt().clamp(max=30)
 mol=b['chem']['pair_mol'][f['indices']];w=torch.bincount(mol,minlength=len(b['rows'])).float().clamp_min(1).reciprocal()[mol];w=w/w.mean();loss=(F.cross_entropy(logits,y,weight=weights,reduction='none')*w).sum()/(weights[y]*w).sum().clamp_min(1e-8)
 double=y==2
 if double.any():loss=loss+.2*F.cross_entropy(stereo[double],b['stereo'][f['indices']][double].long())
 return loss

def train_classifier(group,seed,train,val,diffuser,device):
 folder=OUT/'runs'/str(seed)/group;local=group.startswith('A2');mode='2d' if group.endswith('_2d') else '3d';torch.manual_seed(seed);model=ConnectionHead(local).to(device);opt=torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],lr=.0003,weight_decay=.0001);history=[];best=-1;val=choose_views(val,256,1)
 for epoch in range(40):
  budget();model.train();rows=epoch_rows(train,epoch,seed);losses=[];began=time.time()
  for begin in range(0,len(rows),16):
   budget();batch=rows[begin:begin+16];b=pack(batch,device);ns=seed+epoch*100000+begin;times=noise_times(len(batch),device,ns);T,R=noisy_frames(b,diffuser,times,seed=ns+193);f=features(b,T,R,times,local,mode);opt.zero_grad(set_to_none=True);loss=class_loss(model,b,f);assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.,error_if_nonfinite=True);opt.step();losses.append(float(loss.detach()))
  rec=dict(epoch=epoch+1,train_loss=float(np.mean(losses)),seconds=time.time()-began,molecules=len(rows),optimizer_steps=len(losses));score=None
  if (epoch+1)%5==0:
   metrics,_=connection_evaluation(model,val,diffuser,seed+810000,[.02,.3,.9],mode);score=float(np.mean([x['f1'] for x in metrics]));rec.update(validation=metrics,selection_score=score)
  history.append(rec);checkpoint(folder/'latest.pt',model,epoch+1,score,history,opt,group=group,local_features=local,mode=mode)
  if score is not None and score>best:best=score;checkpoint(folder/'best.pt',model,epoch+1,score,history,group=group,local_features=local,mode=mode)
  dump(folder/'history.json',history);dump(OUT/'status.json',dict(stage='connection_training',seed=seed,group=group,epoch=epoch+1,target_epochs=40,training_complete=False));print('CONNECTION',seed,group,json.dumps(rec),flush=True)
 model.load_state_dict(torch.load(folder/'best.pt',map_location=device,weights_only=False)['model']);return model,history

def geo_selection(model,rows,diffuser,seed,group):
 model.eval();losses=[]
 with torch.no_grad():
  rows=choose_views(rows,256,1)
  for ti,t in enumerate([.02,.3,.9]):
   for begin in range(0,len(rows),16):
    batch=rows[begin:begin+16];b=pack(batch,next(model.parameters()).device);times=torch.full((len(batch),),t,device=b['local'].device);T,R=evaluation_frames(b,diffuser,times,seed+810000);oT,oR=model(b,T,R,times,edges_for(b,[r['connections'] for r in batch]));loss,_=geometry_loss(b,oT,oR,group);losses.append(float(loss))
 return float(np.mean(losses))

def train_oracle_geometry(group,seed,train,val,diffuser,device):
 folder=OUT/'runs'/str(seed)/group/'B_oracle';torch.manual_seed(seed+93);model=GeometryHead().to(device);opt=torch.optim.AdamW(model.parameters(),lr=.0003,weight_decay=.0001);history=[];best=float('inf')
 for epoch in range(40):
  budget();model.train();rows=epoch_rows(train,epoch,seed);losses=[];partsums={};began=time.time()
  for begin in range(0,len(rows),16):
   budget();batch=rows[begin:begin+16];b=pack(batch,device);ns=seed+epoch*100000+begin;times=noise_times(len(batch),device,ns);T,R=noisy_frames(b,diffuser,times,seed=ns+193);opt.zero_grad(set_to_none=True);oT,oR=model(b,T,R,times,edges_for(b,[r['connections'] for r in batch]));loss,parts=geometry_loss(b,oT,oR,group);assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.,error_if_nonfinite=True);opt.step();losses.append(float(loss.detach()))
   for k,v in parts.items():partsums[k]=partsums.get(k,0)+float(v.detach())
  rec=dict(epoch=epoch+1,train_loss=float(np.mean(losses)),parts={k:v/len(losses) for k,v in partsums.items()},seconds=time.time()-began,molecules=len(rows));score=None
  if (epoch+1)%5==0:score=geo_selection(model,val,diffuser,seed,group);rec['selection_score']=score
  history.append(rec);checkpoint(folder/'latest.pt',model,epoch+1,score,history,opt,group=group,oracle_connections=True)
  if score is not None and score<best:best=score;checkpoint(folder/'best.pt',model,epoch+1,score,history,group=group,oracle_connections=True)
  dump(folder/'history.json',history);dump(OUT/'status.json',dict(stage='oracle_geometry_training',seed=seed,group=group,epoch=epoch+1,target_epochs=40,training_complete=False));print('GEOMETRY_B',seed,group,json.dumps(rec),flush=True)
 model.load_state_dict(torch.load(folder/'best.pt',map_location=device,weights_only=False)['model']);return model,history

def collect_graphs(classifier,rows,diffuser,seed,split,limit,mode='3d'):
 device=next(classifier.parameters()).device;cases=[];report=[];rows=choose_views(rows,limit,1)
 for ti,t in enumerate([.02,.3,.9]):
  for k,row in enumerate(rows):
   budget();b=pack([row],device);times=torch.tensor([t],device=device);T,R=evaluation_frames(b,diffuser,times,seed)
   try:graph,info=searched_graph(classifier,row,T,R,t,mode)
   except Exception as exc:graph=None;info=dict(error=repr(exc),legal_complete=False,exact_graph=False)
   report.append(dict(id=row['id'],t=t,split=split,graph=graph,**info))
   if graph is not None:
    bm,am=supervision_masks(row,graph)
    if bm.any():cases.append(dict(row=row,T=T.detach().cpu(),R=R.detach().cpu(),t=t,graph=graph,bond_mask=bm,angle_mask=am))
   if (k+1)%16==0:print('PREDICTED_GRAPH',seed,split,t,k+1,'/',len(rows),flush=True)
 dump(OUT/'runs'/str(seed)/f'predicted_graphs_{split}.json',report);return cases,report

def predicted_batch(cases,device):
 b=pack([v['row'] for v in cases],device);T=torch.cat([v['T'] for v in cases]).to(device);R=torch.cat([v['R'] for v in cases]).to(device);times=torch.tensor([v['t'] for v in cases],device=device);edges=edges_for(b,[v['graph'] for v in cases]);bm=torch.cat([v['bond_mask'] for v in cases]).to(device);am=torch.cat([v['angle_mask'] for v in cases]).to(device);return b,T,R,times,edges,bm,am

def adapt_predicted(group,seed,geo,classifier,cases,valcases,device,joint=False):
 stage='D_joint' if joint else 'C_predicted';folder=OUT/'runs'/str(seed)/group/stage
 if not cases or not valcases:return dict(skipped=True,reason='No legal predicted graph containing a correct supervised edge in train or validation')
 parameters=list(geo.parameters())+([p for p in classifier.parameters() if p.requires_grad] if joint else []);opt=torch.optim.AdamW(parameters,lr=.00005,weight_decay=.0001);history=[];best_state=copy.deepcopy(geo.state_dict());best_class=copy.deepcopy(classifier.state_dict());epochs=5 if joint else 10
 def validation():
  values=[]
  with torch.no_grad():
   for start in range(0,len(valcases),8):
    b,T,R,times,edges,bm,am=predicted_batch(valcases[start:start+8],device);oT,oR=geo(b,T,R,times,edges);loss,_=geometry_loss(b,oT,oR,group,bm,am)
    if joint:loss=loss+class_loss(classifier,b,features(b,oT,oR,times,classifier.local))
    values.append(float(loss))
  return float(np.mean(values))
 best=validation();initial_score=best;history.append(dict(epoch=0,validation_selection_loss=best));checkpoint(folder/'best.pt',geo,0,best,history,classifier_state=classifier.state_dict(),group=group)
 for epoch in range(epochs):
  budget();order=np.random.default_rng(seed+epoch).permutation(len(cases));losses=[]
  for start in range(0,len(order),8):
   batch=[cases[k] for k in order[start:start+8]];b,T,R,times,edges,bm,am=predicted_batch(batch,device);opt.zero_grad(set_to_none=True);oT,oR=geo(b,T,R,times,edges);loss,_=geometry_loss(b,oT,oR,group,bm,am)
   if joint:loss=loss+class_loss(classifier,b,features(b,oT,oR,times,classifier.local))
   assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(parameters,5.,error_if_nonfinite=True);opt.step();losses.append(float(loss.detach()))
  score=validation();rec=dict(epoch=epoch+1,train_loss=float(np.mean(losses)),validation_selection_loss=score,train_legal_views=len(cases),validation_legal_views=len(valcases),joint_gradients='CE -> classifier and differentiable corrected coordinates -> geometry; geometry targets -> geometry only; discrete search fixed' if joint else 'Geometry -> geometry head only; classifier frozen');history.append(rec)
  if score<best:best=score;best_state=copy.deepcopy(geo.state_dict());best_class=copy.deepcopy(classifier.state_dict());checkpoint(folder/'best.pt',geo,epoch+1,score,history,classifier_state=classifier.state_dict(),group=group)
  checkpoint(folder/'latest.pt',geo,epoch+1,score,history,opt,classifier_state=classifier.state_dict(),group=group);dump(folder/'history.json',history);print(stage,seed,group,json.dumps(rec),flush=True)
 geo.load_state_dict(best_state);classifier.load_state_dict(best_class);return dict(skipped=False,history=history,initial_validation_loss=initial_score,best_validation_loss=best,improved=best<initial_score*.99)
