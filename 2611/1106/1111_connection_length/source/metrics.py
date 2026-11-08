"""Raw prediction, bounded legal search, and directly reconstructed geometry."""
from common import *
import time,math
import numpy as np
import torch
from rdkit import Chem
from sklearn.metrics import average_precision_score
from models import *
from trajectory_connections import ConnectionRuntime,load_connection
from hierarchical_topology import hierarchical_search,SearchConfig
from vocab_chemistry import VocabLibrary,assemble_slots

def key(e,stereo=False):return tuple(int(e[k]) for k in ['i','j','a','b','order']+(['stereo'] if stereo else []))
def choose_views(rows,limit=None,conformers=1):
 groups={}
 for r in rows:groups.setdefault(r['molecule_id'],[]).append(r)
 ids=sorted(groups,key=lambda s:stable('independent-evaluation',s))
 if limit is not None:ids=ids[:limit]
 return [r for sid in ids for r in sorted(groups[sid],key=lambda r:r['conformer_index'])[:conformers]]

_EVAL_NOISE={}
def evaluation_frames(b,diffuser,times,seed):
 """Common random numbers indexed by conformer/time, independent of batch/subset."""
 frames=[];offset=0;device=b['chem']['labels'].device
 for m,row in enumerate(b['rows']):
  n=len(row['chem']['labels']);t=round(float(times[m]),6);cachekey=(id(diffuser),seed,row['id'],t,str(device))
  if cachekey not in _EVAL_NOISE:
   single=dict(chem=dict(labels=b['chem']['labels'][offset:offset+n],slot_mol=torch.zeros(n,device=device,dtype=torch.long)),noise={k:v[offset:offset+n] for k,v in b['noise'].items()})
   _EVAL_NOISE[cachekey]=noisy_frames(single,diffuser,times[m:m+1],seed=stable('fixed-evaluation',seed,row['id'],t)%(2**31))
  frames.append(_EVAL_NOISE[cachekey]);offset+=n
 return torch.cat([f[0] for f in frames]),torch.cat([f[1] for f in frames])

def rates(p,y):
 pred=p.argmax(-1);active=pred>0;true=y>0;typed=active&true&(pred==y);binary=active&true
 tp=int(typed.sum());n=int(active.sum());total=int(true.sum());bt=int(binary.sum())
 return dict(precision=tp/max(1,n),recall=tp/max(1,total),f1=2*tp/max(1,n+total),binary_precision=bt/max(1,n),binary_recall=bt/max(1,total),AP=float(average_precision_score(true,p[:,1:].sum(-1))) if total else None,typed_true_positive=tp,predicted=n,true_bonds=total,port_pairs=len(y),bond_order_accuracy_on_true_ports=float((pred[true]==y[true]).mean()) if total else None)

@torch.no_grad()
def connection_evaluation(model,rows,diffuser,seed,grid,mode='3d',batch_size=8,details=False):
 global _LIB
 if _LIB is None:_LIB=VocabLibrary()
 device=next(model.parameters()).device;model.eval();aggregate=[];cases=[]
 for ti,t in enumerate(grid):
  probabilities=[];ys=[];rawexact=rawlegal=0;stereo_correct=stereo_total=0;top={k:0 for k in [1,3,5]};totalpositive=0;start=time.perf_counter()
  for begin in range(0,len(rows),batch_size):
   batch=rows[begin:begin+batch_size];b=pack(batch,device);times=torch.full((len(batch),),t,device=device);T,R=evaluation_frames(b,diffuser,times,seed);f=features(b,T,R,times,model.local,mode);logits,sl=model(b,f);p=logits.softmax(-1).cpu().numpy();sp=sl.softmax(-1).cpu().numpy();index=f['indices'].cpu().numpy();y=b['y'][f['indices']].cpu().numpy();probabilities.append(p);ys.append(y)
   true_stereo=b['stereo'][f['indices']].cpu().numpy();sm=(y==2)&(true_stereo>0);stereo_total+=int(sm.sum());stereo_correct+=int((sp.argmax(-1)[sm]==true_stereo[sm]).sum())
   c=b['chem'];a,bb=c['pairs'][f['indices']].T;meta=torch.stack([c['slot'][a],c['slot'][bb],c['atom'][a],c['atom'][bb]],-1).cpu().numpy();molecule=c['pair_mol'][f['indices']].cpu().numpy();offset=0
   for m,row in enumerate(batch):
    pick=molecule==m;pp=p[pick];yy=y[pick];mm=meta[pick].copy();mm[:,:2]-=offset;orders=pp.argmax(-1);pred=[dict(i=int(z[0]),j=int(z[1]),a=int(z[2]),b=int(z[3]),order=int(o),stereo=int(s.argmax()) if o==2 else 0) for z,o,s in zip(mm,orders,sp[pick]) if o>0];truth={key(e) for e in row['connections']};exact={key(e) for e in pred}==truth;rawexact+=exact
    legal=False
    if len(pred)==len(row['chem']['labels'])-1:
     try:assemble_slots(_LIB,row['chem']['labels'].tolist(),pred);legal=True
     except Exception:pass
    rawlegal+=legal;ranks=[]
    for v in row['connections']:
     group=(mm[:,0]==v['i'])&(mm[:,1]==v['j']);where=np.flatnonzero(group);score=pp[where,v['order']];target=(mm[where,2]==v['a'])&(mm[where,3]==v['b']);rank=1+int((score>score[target][0]).sum());ranks.append(rank)
     for k in top:top[k]+=int(rank<=k)
    totalpositive+=len(ranks)
    if details:cases.append(dict(id=row['id'],molecule_id=row['molecule_id'],conformer_index=row['conformer_index'],t=t,N=len(row['chem']['labels']),raw_exact_graph=exact,raw_legal_complete_graph=legal,port_ranks=ranks,**rates(pp,yy)))
    offset+=len(row['chem']['labels'])
  result=rates(np.concatenate(probabilities),np.concatenate(ys));result.update(t=t,conformers=len(rows),molecules=len({r['molecule_id'] for r in rows}),raw_exact_graph_count=rawexact,raw_legal_complete_graph_count=rawlegal,raw_exact_graph_rate=rawexact/max(1,len(rows)),raw_legal_complete_graph_rate=rawlegal/max(1,len(rows)),topk={str(k):v/max(1,totalpositive) for k,v in top.items()},stereo_accuracy_specified_true_ports=stereo_correct/max(1,stereo_total) if stereo_total else None,stereo_denominator=stereo_total,classifier_forward_batches=math.ceil(len(rows)/batch_size),classifier_molecule_forwards=len(rows),seconds=time.perf_counter()-start);aggregate.append(result)
 return aggregate,cases

_LIB=None;_BASE=None
@torch.no_grad()
def searched_graph(model,row,T,R,t,mode='3d'):
 global _LIB,_BASE
 if _LIB is None:_LIB=VocabLibrary()
 if _BASE is None:_BASE=load_connection(BASELINE,'cpu')
 began=time.perf_counter();device=next(model.parameters()).device;b=pack([row],device);times=torch.tensor([t],device=device);f=features(b,T,R,times,model.local,mode);logits,sl=model(b,f);p=logits.softmax(-1).cpu().numpy();sp=sl.softmax(-1).cpu().numpy();classifier_seconds=time.perf_counter()-began;c=b['chem'];a,bb=c['pairs'][f['indices']].T;meta=torch.stack([c['slot'][a],c['slot'][bb],c['atom'][a],c['atom'][bb]],-1).cpu().numpy();rt=ConnectionRuntime(_LIB,_BASE,row['chem']['labels'].tolist(),[None]*len(row['chem']['labels']),'cpu',policy='wide',head='full');pool=[]
 for k,(i,j,a,bb) in enumerate(meta):
  for order in rt.allowed.get((int(i),int(a)),set())&rt.allowed.get((int(j),int(bb)),set()):
   confidence=float(p[k,order]);pool.append(dict(i=int(i),j=int(j),a=int(a),b=int(bb),order=order,stereo=int(sp[k].argmax()) if order==2 else 0,confidence=confidence,log_odds=float(np.log(max(confidence,1e-8))-np.log(max(float(p[k,0]),1e-8)))))
 started=time.perf_counter();search=hierarchical_search(rt,dict(meta=meta,p=p,sp=sp,pool=pool),k=1,config=SearchConfig());seconds=time.perf_counter()-started
 graph=search['candidates'][0]['connections'] if search['candidates'] else None
 return graph,dict(seconds=seconds,classification_and_features_seconds=classifier_seconds,classifier_forwards=1,legal_complete=graph is not None,exact_graph=graph is not None and {key(e,True) for e in graph}=={key(e,True) for e in row['connections']},search_pool=len(search['pool']),truth_used=search['truth_used'])

def supervision_masks(row,graph):
 selected={key(e,True) for e in graph};correct=[key(e,True) in selected for e in row['connections']];bm=torch.tensor(correct,dtype=torch.bool);am=[];slot=row['atom_slot'].tolist();bi={tuple(sorted(x.tolist())):k for k,x in enumerate(row['bond_index'])}
 for a,j,c in row['angle_index'].tolist():
  required=[bi[tuple(sorted((i,k)))] for i,k in [(a,j),(j,c)] if slot[i]!=slot[k]];am.append(all(correct[k] for k in required))
 return bm,torch.tensor(am,dtype=torch.bool)

def geometry_measure(b,T,R,bond_mask=None,angle_mask=None):
 x,d,cos=geometry_values(b,T,R);m=b['bond_mask'] if bond_mask is None else b['bond_mask']&bond_mask;am=torch.ones(len(cos),device=T.device,dtype=torch.bool) if angle_mask is None else angle_mask
 be=(d-b['bond_target']).abs();ae=(cos.acos()-b['angle_target']).abs();ni=b['nonbond_index'];rd=b['radii'][b['atom_mask']]
 severe=torch.linalg.vector_norm(x[ni[:,0]]-x[ni[:,1]],dim=-1)<.6*(rd[ni[:,0]]+rd[ni[:,1]]) if len(ni) else torch.zeros(0,dtype=torch.bool,device=T.device)
 trueT,trueR=b['noise']['T0'],b['noise']['R0'];slot=torch.repeat_interleave(torch.arange(len(T),device=T.device),b['atom_mask'].sum(1));i,j=slot[b['bond_index']].T
 relativeT=torch.einsum('nji,nj->ni',R[i],T[j]-T[i]);refT=torch.einsum('nji,nj->ni',trueR[i],trueT[j]-trueT[i]);relativeR=R[i].transpose(-1,-2)@R[j];refR=trueR[i].transpose(-1,-2)@trueR[j];errorR=relativeR.transpose(-1,-2)@refR;angleR=((errorR.diagonal(dim1=-2,dim2=-1).sum(-1)-1)/2).clamp(-1,1).acos()
 return dict(bond_MAE_A=float(be[m].mean()) if m.any() else None,bond_count=int(m.sum()),angle_MAE_deg=float(torch.rad2deg(ae[am]).mean()) if am.any() else None,angle_count=int(am.sum()),GEOM_reference_bond_MAE_A=float((d-b['geom_bond']).abs().mean()),relative_translation_MAE_A=float(torch.linalg.vector_norm(relativeT[m]-refT[m],dim=-1).mean()) if m.any() else None,relative_rotation_MAE_deg=float(torch.rad2deg(angleR[m]).mean()) if m.any() else None,severe_interfragment_clash=bool(severe.any()),severe_clash_pairs=int(severe.sum()),bond_pass=bool((be[m]<=.15).all()) if m.any() else False,angle_pass=bool((torch.rad2deg(ae[am])<=15).all()) if am.any() else False)

@torch.no_grad()
def geometry_evaluation(model,rows,diffuser,seed,grid,graphs=None,save_coordinates=False):
 device=next(model.parameters()).device;model.eval();results=[];coordinates_saved={}
 for ti,t in enumerate(grid):
  for k,row in enumerate(rows):
   b=pack([row],device);times=torch.tensor([t],device=device);T,R=evaluation_frames(b,diffuser,times,seed);graph=row['connections'] if graphs is None else graphs.get((row['id'],t))
   if graph is None:results.append(dict(id=row['id'],molecule_id=row['molecule_id'],t=t,N=len(row['chem']['labels']),graph_available=False));continue
   mask,angles=supervision_masks(row,graph);mask=mask.to(device);angles=angles.to(device);outT,outR=model(b,T,R,times,edges_for(b,[graph]));before=geometry_measure(b,T,R,mask,angles);after=geometry_measure(b,outT,outR,mask,angles);x=coordinates(b,T,R);out=coordinates(b,outT,outR);correct={key(e,True) for e in graph}=={key(e,True) for e in row['connections']}
   record=dict(id=row['id'],molecule_id=row['molecule_id'],t=t,N=len(row['chem']['labels']),graph_available=True,graph_exact=correct,correct_edges=int(mask.sum()),predicted_edges=len(graph),true_edges=len(row['connections']),before=before,after=after,geometry_head_calls=1,geometry_network_rounds=model.iterations,coordinate_change_RMS_A=float((out-x).square().sum(-1).mean().sqrt()),joint_success=correct and after['bond_pass'] and after['angle_pass'] and not after['severe_interfragment_clash'])
   record['clean']=geometry_measure(b,b['noise']['T0'],b['noise']['R0'],mask,angles);results.append(record)
   if save_coordinates:coordinates_saved[f'{row["id"]}__{t}']=out.cpu().numpy()
 return results,coordinates_saved
