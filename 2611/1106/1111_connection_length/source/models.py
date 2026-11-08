"""Frozen-compatible classification plus differentiable rigid SE(3) correction."""
from common import *
import torch
from rdkit import Chem
from torch import nn
from torch.nn import functional as F
from robust_connection_model import TimeConnectionPredictor
from robust_connection_data import pack_records,candidate_features,noisy_frames

def pack(rows,device):
 b=pack_records(rows,device);loc=[];z=[];radius=[];mask=[];atoms=[];bond=[];angle=[];nonbond=[];torsions=[];torsion_bonds=[];offset=0;bond_offset=0
 for row in rows:
  start=0
  for xyz in row['local']:
   n=len(xyz);loc.append(F.pad(xyz,(0,0,0,44-n)));z.append(F.pad(row['atom_z'][start:start+n],(0,44-n)));radius.append(F.pad(row['atom_radii'][start:start+n],(0,44-n)));mask.append(torch.arange(44)<n);start+=n
  if row.get('molblock'):
   mol=Chem.MolFromMolBlock(row['molblock'],removeHs=False);parent=row['parent_order'];inv={p:i for i,p in enumerate(parent)};slots=row['atom_slot']
   for ei,e in enumerate(row['connections']):
    if e['order']!=2:continue
    a,bb=row['bond_index'][ei].tolist();left=[inv[v.GetIdx()] for v in mol.GetAtomWithIdx(parent[a]).GetNeighbors() if slots[inv[v.GetIdx()]]==slots[a]];right=[inv[v.GetIdx()] for v in mol.GetAtomWithIdx(parent[bb]).GetNeighbors() if slots[inv[v.GetIdx()]]==slots[bb]]
    if left and right:torsions.append([min(left)+offset,a+offset,bb+offset,min(right)+offset]);torsion_bonds.append(bond_offset+ei)
  bond.append(row['bond_index']+offset);angle.append(row['angle_index']+offset);nonbond.append(row['nonbond_index']+offset);offset+=len(row['atom_z']);bond_offset+=len(row['bond_index'])
 b.update(local=torch.stack(loc).to(device),atom_z=torch.stack(z).to(device),radii=torch.stack(radius).to(device),atom_mask=torch.stack(mask).to(device),bond_index=torch.cat(bond).to(device),bond_target=torch.cat([r['bond_target'] for r in rows]).to(device),bond_mask=torch.cat([r['bond_mask'] for r in rows]).to(device),geom_bond=torch.cat([r['geom_bond'] for r in rows]).to(device),angle_index=torch.cat(angle).to(device),angle_target=torch.cat([r['angle_target'] for r in rows]).to(device),nonbond_index=torch.cat(nonbond).to(device))
 b.update(torsion_index=torch.tensor(torsions,device=device,dtype=torch.long).reshape(-1,4),torsion_bonds=torch.tensor(torsion_bonds,device=device,dtype=torch.long));return b

def cloud(b,T,R):return torch.einsum('nij,naj->nai',R,b['local'])+T[:,None]
def coordinates(b,T,R):return cloud(b,T,R)[b['atom_mask']]

def extra_features(b,T,R,index):
 c=b['chem'];pair=c['pairs'][index];a,bb=pair.T;i,j=c['slot'][a],c['slot'][bb]
 delta=T[j]-T[i];left=torch.einsum('nji,nj->ni',R[i],delta)/10;right=torch.einsum('nji,nj->ni',R[j],-delta)/10;relative=(R[i].transpose(-1,-2)@R[j]).flatten(1)
 world=cloud(b,T,R);values=[]
 for start in range(0,len(index),256):
  si,sj=i[start:start+256],j[start:start+256];ds=torch.cdist(world[si],world[sj]);valid=b['atom_mask'][si,:,None]&b['atom_mask'][sj,None,:];radii=b['radii'][si,:,None]+b['radii'][sj,None,:]
  minimum=ds.masked_fill(~valid,1000).flatten(1).min(-1).values/10;overlap=(.75*radii-ds).clamp_min(0)*valid
  near_i=torch.linalg.vector_norm(b['local'][si]-c['xyz'][a[start:start+256],None],dim=-1);near_j=torch.linalg.vector_norm(b['local'][sj]-c['xyz'][bb[start:start+256],None],dim=-1)
  weight=torch.exp(-near_i.square()/4)[:,:,None]*torch.exp(-near_j.square()/4)[:,None,:]*valid
  values.append(torch.stack([minimum,overlap.sum((1,2))/valid.sum((1,2)).clamp_min(1),(overlap*weight).sum((1,2))/weight.sum((1,2)).clamp_min(1e-6)],-1))
 def env(port,slot):
  d=torch.linalg.vector_norm(b['local'][slot]-c['xyz'][port,None],dim=-1);scale=d.new_tensor([1.,2.,4.]);return (torch.exp(-d.square()[...,None]/scale.square())*b['atom_mask'][slot,:,None]).mean(1)
 return torch.cat([left,right,relative,torch.cat(values),env(a,i),env(bb,j)],-1)

def features(b,T,R,t,local=False,mode='3d'):
 # Enumerate all identity-permitted pairs. No distance pruning or true-label insertion.
 f=candidate_features(b['chem'],T,R,t,topk=1000000000,cutoff=1e9,half=False)
 assert len(f['indices'])==len(b['y'])
 if mode=='2d':
  f['x']=torch.cat([f['x'][:,:14],torch.zeros_like(f['x'][:,14:])],-1)
 if local:f['extra']=extra_features(b,T,R,f['indices']) if mode=='3d' else T.new_zeros((len(f['indices']),24))
 f['mode']=mode
 return f

class ConnectionHead(nn.Module):
 def __init__(self,local=False,baseline=None):
  super().__init__();saved=baseline or torch.load(BASELINE,map_location='cpu',weights_only=False);self.base=TimeConnectionPredictor(saved['legacy_state']);self.base.load_state_dict(saved['model'],strict=True);self.local=local
  if local:
   self.geometry=nn.Sequential(nn.LayerNorm(24),nn.Linear(24,96),nn.SiLU(),nn.Linear(96,64),nn.SiLU());self.bond=nn.Linear(64,3);self.stereo=nn.Linear(64,3)
   nn.init.zeros_(self.bond.weight);nn.init.zeros_(self.bond.bias);nn.init.zeros_(self.stereo.weight);nn.init.zeros_(self.stereo.bias)
 def forward(self,b,f):
  c=b['chem']
  if f.get('mode')=='2d':c=dict(c,xyz=torch.zeros_like(c['xyz']))
  selection=torch.arange(len(f['indices']),device=f['indices'].device);p,s=self.base(c,f,self.base.context(c),selection)
  if self.local:
   h=self.geometry(f['extra']);p=p+self.bond(h);s=s+self.stereo(h)
  return p,s

def skew(w):
 z=torch.zeros_like(w[:,0]);x,y,q=w.unbind(-1);return torch.stack([z,-q,y,q,z,-x,-y,x,z],-1).reshape(-1,3,3)

def edges_for(b,graphs):
 c=b['chem'];ports={(int(s),int(a)):k for k,(s,a) in enumerate(zip(c['slot'].cpu(),c['atom'].cpu()))};edges=[];orders=[];stereo=[];molecules=[];offset=0
 for m,(row,graph) in enumerate(zip(b['rows'],graphs)):
  for e in graph:
   edges.append([ports[offset+e['i'],e['a']],ports[offset+e['j'],e['b']]]);orders.append(e['order']);stereo.append(e.get('stereo',0));molecules.append(m)
  offset+=len(row['chem']['labels'])
 device=c['labels'].device
 return dict(ports=torch.tensor(edges,device=device,dtype=torch.long).reshape(-1,2),order=torch.tensor(orders,device=device),stereo=torch.tensor(stereo,device=device),molecule=torch.tensor(molecules,device=device))

class GeometryHead(nn.Module):
 def __init__(self,iterations=4):
  super().__init__();self.iterations=iterations;self.identity=nn.Embedding(3684,16);self.element=nn.Embedding(128,8);self.atom=nn.Sequential(nn.Linear(11,48),nn.SiLU(),nn.Linear(48,32));self.bag=nn.Linear(16,32);self.net=nn.Sequential(nn.LayerNorm(165),nn.Linear(165,256),nn.SiLU(),nn.Linear(256,128),nn.SiLU(),nn.Linear(128,6));nn.init.zeros_(self.net[-1].weight);nn.init.zeros_(self.net[-1].bias)
 def forward(self,b,T,R,t,edges):
  c=b['chem'];p=edges['ports']
  if not len(p):return T,R
  a,bb=p.T;a,bb=torch.cat([a,bb]),torch.cat([bb,a]);i,j=c['slot'][a],c['slot'][bb];mol=torch.cat([edges['molecule']]*2);order=torch.cat([edges['order']]*2);stereo=torch.cat([edges['stereo']]*2);labels=c['labels'];identity=self.identity(labels)
  atom=self.atom(torch.cat([self.element(b['atom_z'].clamp_max(127)),b['local']/5],-1));atom=(atom*b['atom_mask'][...,None]).sum(1)/b['atom_mask'].sum(1).clamp_min(1)[:,None]
  bag=identity.new_zeros((len(t),16)).index_add(0,c['slot_mol'],identity)/c['n'][:,None];bag=self.bag(bag)
  tm=t[mol,None];time=tm*tm.new_tensor([1.,2.,4.,8.,16.,32.]);temb=torch.cat([time.sin(),time.cos(),tm,tm.square()],-1)
  for _ in range(self.iterations):
   local_delta=torch.einsum('nji,nj->ni',R[i],T[j]-T[i])/10;relative=(R[i].transpose(-1,-2)@R[j]).flatten(1);other=torch.einsum('nji,nj->ni',R[i],torch.einsum('nij,nj->ni',R[j],c['xyz'][bb])+T[j]-T[i])/5
   x=torch.cat([identity[i],identity[j],local_delta,relative,c['xyz'][a]/5,other,F.one_hot(order-1,2).to(T.dtype),F.one_hot(stereo,3).to(T.dtype),temb,atom[i],atom[j],bag[mol]],-1)
   twist=self.net(x);translation=12*torch.tanh(twist[:,:3]);omega=.6*torch.tanh(twist[:,3:]);degree=torch.bincount(i,minlength=len(T)).clamp_min(1).to(T.dtype)
   dT=T.new_zeros(T.shape).index_add(0,i,torch.einsum('nij,nj->ni',R[i],translation))/degree[:,None];dW=T.new_zeros(T.shape).index_add(0,i,omega)/degree[:,None]
   # Remove a pure molecular translation, without constraining any relative pose.
   mean=dT.new_zeros((len(t),3)).index_add(0,c['slot_mol'],dT)/c['n'][:,None];T=T+dT-mean[c['slot_mol']];R=R@torch.matrix_exp(skew(dW))
  return T,R

def geometry_values(b,T,R):
 x=coordinates(b,T,R);bi=b['bond_index'];distance=torch.linalg.vector_norm(x[bi[:,0]]-x[bi[:,1]],dim=-1)
 ai=b['angle_index'];v=x[ai[:,0]]-x[ai[:,1]];w=x[ai[:,2]]-x[ai[:,1]];cos=F.cosine_similarity(v,w,dim=-1).clamp(-1+1e-7,1-1e-7)
 return x,distance,cos

def geometry_loss(b,T,R,group,bond_mask=None,angle_mask=None):
 x,d,cos=geometry_values(b,T,R);mask=b['bond_mask'] if bond_mask is None else b['bond_mask']&bond_mask
 bond=(d[mask]-b['bond_target'][mask]).square().mean() if mask.any() else d.sum()*0
 am=torch.ones(len(cos),device=T.device,dtype=torch.bool) if angle_mask is None else angle_mask
 angle=(cos[am]-b['angle_target'][am].cos()).square().mean() if am.any() else d.sum()*0
 # Relative bond axes in each incident fragment frame; invariant to free axial torsion.
 trueT=b['noise']['T0'];trueR=b['noise']['R0'];truth=coordinates(b,trueT,trueR);bi=b['bond_index'];vec=x[bi[:,1]]-x[bi[:,0]];ref=truth[bi[:,1]]-truth[bi[:,0]]
 slot=torch.repeat_interleave(torch.arange(len(T),device=T.device),b['atom_mask'].sum(1))
 bm=torch.ones(len(bi),device=T.device,dtype=torch.bool) if bond_mask is None else bond_mask
 pose=d.sum()*0
 if bm.any():
  for end in [0,1]:
   si=slot[bi[:,end]];pred=F.normalize(torch.einsum('nji,nj->ni',R[si],vec),dim=-1);target=F.normalize(torch.einsum('nji,nj->ni',trueR[si],ref),dim=-1);pose=pose+(pred[bm]-target[bm]).square().mean()/2
 ti=b['torsion_index'];tm=bm[b['torsion_bonds']]
 if len(ti) and tm.any():
  def circular(xyz):
   q=xyz[ti];u=q[:,0]-q[:,1];v=F.normalize(q[:,2]-q[:,1],dim=-1);w=q[:,3]-q[:,2];u=F.normalize(u-(u*v).sum(-1)[:,None]*v,dim=-1);w=F.normalize(w-(w*v).sum(-1)[:,None]*v,dim=-1);return torch.stack([(u*w).sum(-1),(torch.cross(u,w,dim=-1)*v).sum(-1)],-1)
  pose=pose+(circular(x)[tm]-circular(truth)[tm]).square().mean()
 local=x.sum()*0;ni=b['nonbond_index']
 if len(ni):
  rd=b['radii'][b['atom_mask']];nd=torch.linalg.vector_norm(x[ni[:,0]]-x[ni[:,1]],dim=-1);local=(.6*(rd[ni[:,0]]+rd[ni[:,1]])-nd).clamp_min(0).square().mean()
 loss=bond if group=='A3' else bond+2*angle+.5*pose
 return loss,dict(bond=bond,angle=angle,relative_axis=pose,clash=local)
