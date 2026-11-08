"""Canonical graph splits and explicitly audited original GEOM conformers."""
from common import *
import argparse,sqlite3,zlib,pickle,itertools,traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import torch
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.spatial.transform import Rotation

def database(name):return sqlite3.connect(f'file:{DATA}/geom_database/{name}?mode=ro&immutable=1',uri=True)
def scaffold(s):
 mol=Chem.MolFromSmiles(s);Chem.RemoveStereochemistry(mol)
 return MurckoScaffold.MurckoScaffoldSmiles(mol=mol,includeChirality=False)
def choose(counts):
 c=database('samples.sqlite');g=database('graphs.sqlite')
 info={gid:smiles for gid,smiles in g.execute('SELECT graph_id,canonical_smiles FROM graphs')}
 groups={s:{} for s in ['train','dev','test']};excluded=Counter();rawcounts=Counter()
 for sid,gid,split,eligible in c.execute('SELECT sample_id,graph_id,split,identity_connection_eligible FROM samples'):
  rawcounts[split]+=1
  if not eligible:excluded[split+':not_in_frozen_vocabulary']+=1;continue
  if gid in groups[split]:excluded[split+':duplicate_canonical_graph']+=1;continue
  groups[split][gid]=sid
 def pool(split):
  for gid in sorted(groups[split],key=lambda x:stable('graph-split-20261106',info[x])):
   graph=json.loads(zlib.decompress(g.execute('SELECT graph_json FROM graphs WHERE graph_id=?',(gid,)).fetchone()[0]))
   if graph['n_fragments']<2:excluded[split+':single_fragment_no_connection_target']+=1;continue
   yield dict(id=groups[split][gid],graph_id=gid,canonical=info[gid],scaffold=scaffold(info[gid]),source_split=split)
 selected={};selected['scaffold_test']=list(itertools.islice(pool('test'),counts['scaffold_test']))
 heldout={r['scaffold'] for r in selected['scaffold_test']}
 for name,split in [('train','train'),('val','dev'),('test','test')]:
  rows=[]
  for row in pool(split):
   if row['scaffold'] in heldout:excluded[split+':reserved_test_scaffold']+=1;continue
   rows.append(row)
   if len(rows)==counts[name]:break
  selected[name]=rows
 assert all(len(selected[k])==v for k,v in counts.items()),{k:len(v) for k,v in selected.items()}
 sets=[{r['canonical'] for r in v} for v in selected.values()]
 assert all(not (a&b) for a,b in itertools.combinations(sets,2))
 assert not heldout&{r['scaffold'] for k,v in selected.items() if k!='scaffold_test' for r in v}
 return selected,dict(source_records=dict(rawcounts),selection_exclusions=dict(excluded),canonical_cross_split=0,scaffold_test_overlap=0,selection='Deterministic hash ordering, preserve original train/dev/test; reserve scaffold-test scaffold families before selecting all other sets; no outcome-based selection')

def fit_frame(template,xyz):
 x=template-template.mean(0);center=xyz.mean(0);y=xyz-center
 u,_,vh=np.linalg.svd(x.T@y);q=u@vh
 if np.linalg.det(q)<0:u[:,-1]*=-1;q=u@vh
 return center,q.T,(xyz-center)@q.T,float(np.sqrt(np.mean(np.sum((x@q-y)**2,axis=-1))))

_LIB=None;_CODEC=None
def process(item):
 global _LIB,_CODEC
 torch.set_num_threads(1);RDLogger.DisableLog('rdApp.*')
 if _LIB is None:
  from vocab_chemistry import VocabLibrary
  from cutoff_data.codec import get_codec
  _LIB=VocabLibrary();_CODEC=get_codec()
 from robust_connection_data import compile_chemistry
 c=database('samples.sqlite');g=database('graphs.sqlite');s=json.loads(zlib.decompress(c.execute('SELECT sample_json FROM samples WHERE sample_id=?',(item['id'],)).fetchone()[0]));graph=json.loads(zlib.decompress(g.execute('SELECT graph_json FROM graphs WHERE graph_id=?',(item['graph_id'],)).fetchone()[0]));c.close();g.close()
 result=dict(selection=item,source_sha256=s['source_sha256'],conformers=[],failures=[])
 try:
  assert sha(s['pickle_path'])==s['source_sha256'],'raw pickle SHA mismatch'
  with open(s['pickle_path'],'rb') as f:raw=pickle.load(f)
  count=len(raw['conformers']);indices=[0]+sorted(range(1,count),key=lambda x:stable('conformer-20261106',s['sample_id'],x))[:3]
  result['available_conformers']=count;result['selected_indices']=indices
  labels=[f['class_id'] for f in graph['fragments']];assert all(0<=x<3684 for x in labels)
  parent_order=[a for f in graph['fragments'] for a in f['core_to_parent']];assert sorted(parent_order)==list(range(graph['n_heavy']))
  inverse=np.argsort(parent_order);connections=[]
  for v in graph['connections']:
   a,b=v['endpoints'];st={0:0,int(Chem.BondStereo.STEREOE):1,int(Chem.BondStereo.STEREOZ):2}[v['stereo']]
   connections.append(dict(i=a['fragment'],j=b['fragment'],a=a['anchor_core'],b=b['anchor_core'],order=int(v['bond_order']),stereo=st))
  for ci in indices:
   try:
    rawmol=raw['conformers'][ci]['rd_mol'];canonical,source_map,parent_xyz,stereo=_CODEC.prepare_molecule(rawmol)
    assert canonical==item['canonical'],'conformer canonical graph/identity mismatch'
    assert stereo['status']=='specified_stereo_consistent',f'3D stereo mismatch: {stereo}'
    cleanmol=Chem.MolFromSmiles(canonical);cf=Chem.Conformer(cleanmol.GetNumAtoms())
    for ai,xyz in enumerate(parent_xyz):cf.SetAtomPosition(ai,xyz)
    cleanmol.AddConformer(cf);molh=Chem.AddHs(cleanmol,addCoords=True);props=AllChem.MMFFGetMoleculeProperties(molh,mmffVariant='MMFF94s')
    chem=compile_chemistry(_LIB,labels);T=[];R=[];cloud=[];atomic=[];radii=[];fit=[];atom_slots=[];parent_to_slot={}
    for si,f in enumerate(graph['fragments']):
     tmpl,ref=_LIB.template(labels[si]);xyz=parent_xyz[f['core_to_parent']];center,rot,local,error=fit_frame(ref,xyz)
     assert np.max(np.abs(local@rot.T+center-xyz))<1e-6
     T.append(center);R.append(rot);cloud.append(local);fit.append(error);atomic.extend(a.GetAtomicNum() for a in tmpl.GetAtoms());radii.extend(Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in tmpl.GetAtoms());atom_slots.extend([si]*len(local))
     for a in f['core_to_parent']:parent_to_slot[a]=si
     for pi in torch.where(chem['slot']==si)[0].tolist():
      ai=int(chem['atom'][pi]);nb=[x.GetIdx() for x in tmpl.GetAtomWithIdx(ai).GetNeighbors()];vec=local[ai]-local[nb] if nb else np.zeros((1,3));out=(vec/np.linalg.norm(vec,axis=-1,keepdims=True).clip(1e-8)).sum(0);out/=max(np.linalg.norm(out),1e-8);chem['xyz'][pi]=torch.tensor(local[ai],dtype=torch.float32);chem['direction'][pi]=torch.tensor(out,dtype=torch.float32)
    T=np.asarray(T);R=np.asarray(R);T-=T.mean(0)
    n=len(labels);edges=torch.zeros(n,n,dtype=torch.long);anchors=torch.full((n,n,2),-100,dtype=torch.long);st=torch.zeros(n,n,dtype=torch.long)
    bidx=[];btarget=[];bmask=[];bgeom=[]
    for v in connections:
     i,j,a,b,o=(v[k] for k in ['i','j','a','b','order']);edges[i,j]=edges[j,i]=o;anchors[i,j]=torch.tensor([a,b]);anchors[j,i]=torch.tensor([b,a]);st[i,j]=st[j,i]=v['stereo']
     pa=graph['fragments'][i]['core_to_parent'][a];pb=graph['fragments'][j]['core_to_parent'][b];param=props.GetMMFFBondStretchParams(molh,pa,pb) if props else None
     bidx.append([int(inverse[pa]),int(inverse[pb])]);btarget.append(float(param[2]) if param else 0.);bmask.append(param is not None);bgeom.append(float(np.linalg.norm(parent_xyz[pa]-parent_xyz[pb])))
    ia,ib=chem['pairs'].long().T;si,sj=chem['slot'][ia],chem['slot'][ib];match=(anchors[si,sj,0]==chem['atom'][ia])&(anchors[si,sj,1]==chem['atom'][ib]);y=edges[si,sj].to(torch.int8);y[~match]=0
    assert int((y>0).sum())==len(connections),'true port not enumerated'
    angles=[];atarget=[]
    for a in cleanmol.GetAtoms():
     j=a.GetIdx()
     for i,k in itertools.combinations([x.GetIdx() for x in a.GetNeighbors()],2):
      if parent_to_slot[i]==parent_to_slot[j]==parent_to_slot[k]:continue
      v,w=parent_xyz[i]-parent_xyz[j],parent_xyz[k]-parent_xyz[j];angle=np.arccos(np.clip(np.dot(v,w)/(np.linalg.norm(v)*np.linalg.norm(w)),-1,1));angles.append([int(inverse[i]),int(inverse[j]),int(inverse[k])]);atarget.append(float(angle))
    dist=Chem.GetDistanceMatrix(cleanmol);nonbond=[]
    for a in range(len(parent_order)):
     for b in range(a+1,len(parent_order)):
      if atom_slots[a]!=atom_slots[b] and dist[parent_order[a],parent_order[b]]>2:nonbond.append([a,b])
    row=dict(id=f"{item['id']}:{ci}",molecule_id=item['id'],canonical=canonical,scaffold=item['scaffold'],source_split=item['source_split'],conformer_index=ci,canonical_to_source_atom=source_map,chem=chem,target=dict(y=y,stereo=st[si,sj].to(torch.int8),true_bonds=len(connections),enumerated_positive=int((y>0).sum())),connections=connections,noise=dict(T0=torch.tensor(T,dtype=torch.float32),R0=torch.tensor(R,dtype=torch.float32),C0=torch.tensor(T/2.7,dtype=torch.float32),Q0=torch.tensor(R,dtype=torch.float32),refR=torch.eye(3).repeat(n,1,1),offset=torch.zeros(n,3)),local=[torch.tensor(x,dtype=torch.float32) for x in cloud],atom_z=torch.tensor(atomic),atom_radii=torch.tensor(radii),atom_slot=torch.tensor(atom_slots),bond_index=torch.tensor(bidx).reshape(-1,2),bond_target=torch.tensor(btarget),bond_mask=torch.tensor(bmask),geom_bond=torch.tensor(bgeom),angle_index=torch.tensor(angles).reshape(-1,3),angle_target=torch.tensor(atarget),nonbond_index=torch.tensor(nonbond).reshape(-1,2),template_fit_rmsd_A=fit,parent_order=parent_order,molblock=Chem.MolToMolBlock(cleanmol),MMFF94s_parameterized=props is not None)
    result['conformers'].append(row)
   except Exception as exc:result['failures'].append(dict(conformer_index=ci,error=repr(exc),traceback=traceback.format_exc()))
 except Exception as exc:result['failures'].append(dict(error=repr(exc),traceback=traceback.format_exc()))
 return result

def main():
 p=argparse.ArgumentParser();p.add_argument('--smoke',action='store_true');p.add_argument('--workers',type=int,default=8);a=p.parse_args()
 counts=dict(train=16,val=8,test=8,scaffold_test=8) if a.smoke else dict(train=2000,val=256,test=256,scaffold_test=256)
 output=OUT/('data_smoke' if a.smoke else 'data');output.mkdir(exist_ok=True)
 selected,audit=choose(counts);dump(output/'selected_molecules.json',selected);dump(output/'selection_audit.json',audit)
 for split,items in selected.items():
  rows=[];allresults=[];failed=[]
  with ProcessPoolExecutor(max_workers=a.workers) as executor:
   for k,res in enumerate(executor.map(process,items,chunksize=1)):
    rows.extend(res.pop('conformers'));allresults.append(res)
    for f in res['failures']:failed.append(dict(molecule_id=res['selection']['id'],**f))
    if (k+1)%50==0:print(split,k+1,'molecules;',len(rows),'conformers;',len(failed),'failures',flush=True)
  torch.save(rows,output/f'{split}.pt');dump(output/f'{split}_sources.json',allresults);dump(output/f'{split}_failures.json',failed)
  dump(output/f'{split}_summary.json',dict(raw_selected_molecules=len(items),usable_molecules=len({r['molecule_id'] for r in rows}),usable_conformers=len(rows),failures=len(failed),N_distribution=dict(Counter(str(len(r['chem']['labels'])) for r in rows)),true_bond_orders=dict(Counter(str(c['order']) for r in rows for c in r['connections'])),mmff_bonds_available=sum(int(r['bond_mask'].sum()) for r in rows),mmff_bonds_total=sum(len(r['bond_mask']) for r in rows),pt_sha256=sha(output/f'{split}.pt')))
  print(split,load(output/f'{split}_summary.json'),flush=True)
 dump(output/'complete.json',dict(complete=True,selected=counts,split_pt_sha256={s:sha(output/f'{s}.pt') for s in counts},source_fingerprints=load(OUT/'audit_initial.json')['source_fingerprints']))
if __name__=='__main__':main()
