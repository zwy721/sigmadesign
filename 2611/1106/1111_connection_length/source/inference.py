"""Deployable input contract: types, optional local atoms/ports, T/R, t only."""
from common import *
import numpy as np
import torch
from rdkit import Chem
from models import *
from metrics import searched_graph
from vocab_chemistry import VocabLibrary
from robust_connection_data import compile_chemistry

def input_row(types,local_atoms=None,candidate_ports=None):
 lib=VocabLibrary();types=list(map(int,types));assert all(0<=x<3684 for x in types),'Expected known predicted fragment classes, never MASK'
 chem=compile_chemistry(lib,types);local=[];zs=[];radii=[];slots=[]
 for i,typ in enumerate(types):
  mol,xyz=lib.template(typ);xyz=np.asarray(local_atoms[i]) if local_atoms is not None else xyz
  assert xyz.shape==(mol.GetNumAtoms(),3) and np.isfinite(xyz).all()
  local.append(torch.tensor(xyz,dtype=torch.float32));zs.extend(a.GetAtomicNum() for a in mol.GetAtoms());radii.extend(Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms());slots.extend([i]*len(xyz))
  for p in torch.where(chem['slot']==i)[0].tolist():
   ai=int(chem['atom'][p]);nb=[a.GetIdx() for a in mol.GetAtomWithIdx(ai).GetNeighbors()];v=xyz[ai]-xyz[nb] if nb else np.zeros((1,3));direction=(v/np.linalg.norm(v,axis=-1,keepdims=True).clip(1e-8)).sum(0);direction/=max(np.linalg.norm(direction),1e-8);chem['xyz'][p]=torch.tensor(xyz[ai],dtype=torch.float32);chem['direction'][p]=torch.tensor(direction,dtype=torch.float32)
 if candidate_ports is not None:
  allowed={tuple(map(int,v)) for v in candidate_ports};ports={(int(s),int(a)):p for p,(s,a) in enumerate(zip(chem['slot'],chem['atom']))};assert allowed<=set(ports),'Invalid supplied candidate port';allowed_ids={ports[v] for v in allowed};chem['pairs']=torch.tensor([p.tolist() for p in chem['pairs'] if int(p[0]) in allowed_ids and int(p[1]) in allowed_ids],dtype=torch.int32).reshape(-1,2)
 n=len(types);zero=torch.zeros(len(chem['pairs']),dtype=torch.int8)
 return dict(id='runtime_input',molecule_id='runtime_input',chem=chem,local=local,atom_z=torch.tensor(zs),atom_radii=torch.tensor(radii),atom_slot=torch.tensor(slots),target=dict(y=zero,stereo=zero.clone(),true_bonds=0),connections=[],noise=dict(T0=torch.zeros(n,3),R0=torch.eye(3).repeat(n,1,1),C0=torch.zeros(n,3),Q0=torch.eye(3).repeat(n,1,1),refR=torch.eye(3).repeat(n,1,1),offset=torch.zeros(n,3)),bond_index=torch.zeros(0,2,dtype=torch.long),bond_target=torch.zeros(0),bond_mask=torch.zeros(0,dtype=torch.bool),geom_bond=torch.zeros(0),angle_index=torch.zeros(0,3,dtype=torch.long),angle_target=torch.zeros(0),nonbond_index=torch.zeros(0,2,dtype=torch.long))

class Runtime:
 def __init__(self,bundle,device='cuda:0'):
  saved=torch.load(bundle,map_location='cpu',weights_only=False);self.classifier=ConnectionHead(local=saved['local_features']).to(device);self.classifier.load_state_dict(saved['classifier']);self.geometry=GeometryHead().to(device);self.geometry.load_state_dict(saved['geometry']);self.classifier.eval();self.geometry.eval();self.device=device
 @torch.no_grad()
 def predict(self,types,T,R,t,*,candidate_ports=None,local_atoms=None,legal_connections=None):
  row=input_row(types,local_atoms,candidate_ports);T=torch.as_tensor(T,device=self.device,dtype=torch.float32);R=torch.as_tensor(R,device=self.device,dtype=torch.float32);times=torch.tensor([t],device=self.device,dtype=torch.float32);b=pack([row],self.device)
  if legal_connections is None:graph,search=searched_graph(self.classifier,row,T,R,t)
  else:
   from vocab_chemistry import assemble_slots
   assemble_slots(VocabLibrary(),types,legal_connections);graph=legal_connections;search=dict(supplied_legal_graph=True)
  if graph is None:return dict(status='no_legal_graph',T=T.cpu().numpy(),R=R.cpu().numpy(),connections=None,search=search)
  newT,newR=self.geometry(b,T,R,times,edges_for(b,[graph]));return dict(status='complete',T=newT.cpu().numpy(),R=newR.cpu().numpy(),connections=graph,search=search,geometry_iterations=self.geometry.iterations)

def audit_native(bundle):
 from run_validation import load_official
 from joint_sampler import chemical_graph,NativeState
 from sigmadock.chem.parsing import get_protein,structure_to_rdkit
 from validation import SamplingAudit
 lib=VocabLibrary();pocket_path=ROOT/'2611/1105/1039+no_water_test/pockets/12AV_no_water.pdb';pocket=structure_to_rdkit(get_protein(pocket_path));context=dict(pocket=pocket,pocket_com=pocket.GetConformer().GetPositions().mean(0).astype('float32'),pdb_path=str(pocket_path),ligand_path=None)
 runtime=Runtime(bundle);types=[3577,3577,1214];random=np.random.default_rng(20261106);T=random.normal(0,2.7,(3,3))+context['pocket_com'];from scipy.spatial.transform import Rotation
 R=Rotation.random(3,random_state=20261106).as_matrix();pred=runtime.predict(types,T,R,.3);assert pred['status']=='complete',pred['status']
 graph,groups,spec=chemical_graph(lib,types,pred['connections'],context,[None]*3,seed=20261106);den=load_official(bootstrap.CHECKPOINT,OUT/'cache/native','cuda:0')
 with torch.inference_mode(),SamplingAudit(den) as sampling:
  state=NativeState(den,graph,groups,'cuda:0');state.step(torch.tensor(.3,device='cuda:0'),torch.tensor(.01,device='cuda:0'));tt,rr,fit=state.slot_frames(lib,types,[None]*3);xyz=state.physical();local=[(xyz[state.pn+np.asarray(g)]-tt[i])@rr[i] for i,g in enumerate(groups)];p=runtime.predict(types,tt,rr,.29,local_atoms=local,legal_connections=pred['connections']);desired=[np.asarray(v)@p['R'][i].T+p['T'][i] for i,v in enumerate(local)];replacement=NativeState(den,graph,groups,'cuda:0',desired=desired,preferred=p['R']);assert max(replacement.transfer_errors)<.01;replacement.step(torch.tensor(.29,device='cuda:0'),torch.tensor(.01,device='cuda:0'));assert np.isfinite(replacement.physical()).all()
 report=dict(passed=True,scope='API compatibility smoke with fixed supplied predicted-vocabulary types, generated legal graph and random pocket poses; not docking efficacy',types=types,ground_truth_ligand_or_GEOM_used_at_runtime=False,geometry_transfer_errors_A=replacement.transfer_errors,sampling_audit=sampling.report(),official_checkpoint_sha256=sha(bootstrap.CHECKPOINT),upstream_weights_unchanged=all(sha(p)==v['sha256'] for p,v in load(OUT/'audit_initial.json')['source_fingerprints'].items()))
 dump(OUT/'native_integration.json',report);return report
