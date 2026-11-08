"""Independent MMFF baseline and multi-conformer/torsion diagnostics."""
from common import *
import time,itertools
from collections import defaultdict
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolTransforms
from models import *
from metrics import choose_views,geometry_measure,evaluation_frames

def aligned_rmsd(a,b):
 a=a-a.mean(0);b=b-b.mean(0);u,_,v=np.linalg.svd(a.T@b);d=np.eye(3);d[-1,-1]=np.linalg.det(u@v);return float(np.sqrt(np.mean(np.sum((a@u@d@v-b)**2,axis=-1))))

def conformers(rows,arrays,group,seed,split):
 grouped=defaultdict(list);lookup={r['id']:r for r in rows};records=[];torsions=[]
 for name,xyz in arrays.items():
  sid,ts=name.rsplit('__',1);r=lookup[sid];grouped[r['molecule_id'],float(ts)].append((r,xyz))
  parent=Chem.MolFromMolBlock(r['molblock'],removeHs=False);cf=Chem.Conformer(parent.GetNumAtoms());inverse=np.argsort(r['parent_order']);canonical_xyz=xyz[inverse]
  for i,p in enumerate(canonical_xyz):cf.SetAtomPosition(i,tuple(map(float,p)))
  for e in r['connections']:
   # Parent atoms follow the exact saved fragment-to-parent map.
   start=sum(len(v) for v in r['local'][:e['i']]);end=sum(len(v) for v in r['local'][:e['j']]);a=r['parent_order'][start+e['a']];b=r['parent_order'][end+e['b']];left=[x.GetIdx() for x in parent.GetAtomWithIdx(a).GetNeighbors() if x.GetIdx()!=b];right=[x.GetIdx() for x in parent.GetAtomWithIdx(b).GetNeighbors() if x.GetIdx()!=a]
   if left and right:
    ids=[min(left),a,b,min(right)];value=rdMolTransforms.GetDihedralDeg(cf,*ids);reference=rdMolTransforms.GetDihedralDeg(parent.GetConformer(),*ids);error=abs((value-reference+180)%360-180);torsions.append(dict(group=group,seed=seed,split=split,id=sid,molecule_id=r['molecule_id'],t=float(ts),order=e['order'],dihedral_deg=value,reference_dihedral_deg=reference,circular_error_deg=error,free_single_bond_not_supervised=e['order']==1))
 refs=defaultdict(list)
 for r in rows:
  T=r['noise']['T0'];R=r['noise']['R0'];refs[r['molecule_id']].append(np.concatenate([x.numpy()@R[i].numpy().T+T[i].numpy() for i,x in enumerate(r['local'])]))
 for (m,t),values in grouped.items():
  predictions=[x for _,x in values];reference=refs[m];dist=np.array([[aligned_rmsd(x,y) for y in reference] for x in predictions]);pairs=[aligned_rmsd(a,b) for a,b in itertools.combinations(predictions,2)]
  refpairs=[aligned_rmsd(a,b) for a,b in itertools.combinations(reference,2)];records.append(dict(group=group,seed=seed,split=split,molecule_id=m,t=t,predicted_conformers=len(predictions),reference_conformers=len(reference),pair_RMSD_mean_A=float(np.mean(pairs)) if pairs else None,reference_pair_RMSD_mean_A=float(np.mean(refpairs)) if refpairs else None,reference_coverage_at_1A=float((dist.min(0)<=1).mean()),nearest_reference_RMSD_mean_A=float(dist.min(1).mean()),symmetry_correction=False))
 return records,torsions

def mmff(rows,diffuser,seed,split,device):
 result=[]
 for ti,t in enumerate([.02,.3,.9]):
  for k,row in enumerate(choose_views(rows,16,1)):
   b=pack([row],device);times=torch.tensor([t],device=device);T,R=evaluation_frames(b,diffuser,times,seed);coords=coordinates(b,T,R).cpu().numpy();mol=Chem.MolFromMolBlock(row['molblock'],removeHs=False);cf=Chem.Conformer(mol.GetNumAtoms());inverse=np.argsort(row['parent_order'])
   for i,p in enumerate(coords[inverse]):cf.SetAtomPosition(i,tuple(map(float,p)))
   mol.RemoveAllConformers();mol.AddConformer(cf);mol=Chem.AddHs(mol,addCoords=True);began=time.time()
   try:
    status=AllChem.MMFFOptimizeMolecule(mol,mmffVariant='MMFF94s',maxIters=200);assert status>=0,'MMFF unavailable';xyz=mol.GetConformer().GetPositions()[:len(inverse)][row['parent_order']];x=torch.tensor(xyz,device=device,dtype=T.dtype)
    bi=b['bond_index'];d=torch.linalg.vector_norm(x[bi[:,0]]-x[bi[:,1]],dim=-1);ai=b['angle_index'];cos=F.cosine_similarity(x[ai[:,0]]-x[ai[:,1]],x[ai[:,2]]-x[ai[:,1]],dim=-1).clamp(-1+1e-7,1-1e-7);m=b['bond_mask'];ni=b['nonbond_index'];rd=b['radii'][b['atom_mask']];clash=bool((torch.linalg.vector_norm(x[ni[:,0]]-x[ni[:,1]],dim=-1)<.6*(rd[ni[:,0]]+rd[ni[:,1]])).any()) if len(ni) else False
    rec=dict(status='evaluated',optimizer_status=status,converged=status==0,bond_MAE_A=float((d[m]-b['bond_target'][m]).abs().mean()),angle_MAE_deg=float(torch.rad2deg((cos.acos()-b['angle_target']).abs()).mean()),severe_interfragment_clash=clash)
   except Exception as exc:rec=dict(status='failed',error=repr(exc))
   result.append(dict(seed=seed,split=split,id=row['id'],t=t,topology='Oracle true chemical graph',method='Independent full-atom MMFF94s from noisy input; maxIters=200; no neural attribution',seconds=time.time()-began,**rec))
 return result
