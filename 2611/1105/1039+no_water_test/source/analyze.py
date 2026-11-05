from common import *
from collections import Counter,defaultdict
from itertools import combinations
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator,rdMolAlign
from denovo_quality.io import load_molecule,heavy_copy
from denovo_quality.protein import read_protein
from denovo_quality.chemistry import standardized_smiles
from denovo_quality.diversity import distribution,training_reference
from geometry_constraints_20261103.stereo import derive
STAGES=['decoded_initial','native_final','rigid','local_flex']

def cloud_metrics(xyz,z,protein):
 atoms=[a for a in read_protein(protein)['atoms'] if a['kind']=='protein'];coords=np.array([a['xyz'] for a in atoms]);pt=Chem.GetPeriodicTable();r=np.array([pt.GetRvdw(int(v)) for v in z]);pr=np.array([a['vdw_radius_A'] for a in atoms]);d=np.linalg.norm(xyz[:,None]-coords[None,:],axis=-1);ratio=d/(r[:,None]+pr[None,:]);bad=ratio<load(QC)['protein']['severe_vdw_ratio'];pairs=np.argwhere(bad)
 return dict(protein_pass=not bad.any(),protein_clash_pairs=int(bad.sum()),protein_clashing_atoms=int(bad.any(1).sum()),protein_min_vdw_ratio=float(ratio.min()),protein_min_distance_A=float(d.min())),[dict(ligand_atom=int(i),protein_atom=int(j),protein_residue=atoms[j]['residue'],protein_atom_name=atoms[j]['name'],distance_A=float(d[i,j]),vdw_ratio=float(ratio[i,j])) for i,j in pairs]

def rate(rows,k,n):return dict(pass_count=sum(r.get(k) is True for r in rows),fail_count=sum(r.get(k) is False for r in rows),unavailable_count=n-sum(isinstance(r.get(k),(bool,np.bool_)) for r in rows),denominator=n)
def stereo_change(before,after):
 a={(d['type'],d['center']):d for d in before['decisions']};b={(d['type'],d['center']):d for d in after['decisions']};out=dict(flips=0,lost=0,regained=0)
 for k in a.keys()|b.keys():
  x,y=a.get(k,{}),b.get(k,{})
  if x.get('reliable') and y.get('reliable'):out['flips']+=x['geometry_CIP']!=y['geometry_CIP']
  out['lost']+=bool(x.get('reliable') and not y.get('reliable'));out['regained']+=bool(y.get('reliable') and not x.get('reliable'))
 return out

def diversity(records,mols,n):
 gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);valid=[r for r in records if r['checks'].get('chemically_valid_connected') is True];chem=[r['modules']['chemistry'] for r in valid];cs=[c['constitution_smiles'] for c in chem];iso=[c['canonical_smiles'] for c in chem];fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in cs];sim=[DataStructs.TanimotoSimilarity(a,b) for a,b in combinations(fps,2)];same=defaultdict(list);conformers=[]
 for r in valid:same[r['metadata']['receptor_id'],r['modules']['chemistry']['canonical_smiles']].append(r)
 for (pid,smi),rr in same.items():
  ar=[];world=[]
  for a,b in combinations(rr,2):
   ma,mb=mols[a['molecule_id']],mols[b['molecule_id']];x=ma.GetConformer().GetPositions();y=mb.GetConformer().GetPositions();matches=ma.GetSubstructMatches(mb,useChirality=True,uniquify=False,maxMatches=4096)
   if matches:world.append(min(float(np.sqrt(np.square(x[list(m)]-y).sum(1).mean())) for m in matches));ar.append(float(rdMolAlign.GetBestRMS(Chem.Mol(ma),Chem.Mol(mb),maxMatches=4096)))
  if ar:conformers.append(dict(receptor_id=pid,smiles=smi,n=len(rr),aligned_RMSD_A=distribution(ar),receptor_frame_RMSD_A=distribution(world),distinct_aligned_1A_fraction=float(np.mean(np.array(ar)>=1))))
 return dict(denominator=n,valid=len(valid),unique_constitutions=len(set(cs)),unique_isomeric=len(set(iso)),unique_constitution_fraction_all=len(set(cs))/n,internal_diversity=1-float(np.mean(sim)) if sim else None,pair_tanimoto=distribution(sim),conformer_groups=conformers)

def cohort(label,records,raw_clouds,denominator,attempts=None):
 target=OUT if label=='paired_oracle' else OUT/'free_probe';target.mkdir(exist_ok=True)
 per=[];bonds=[];angles=[];clashes=[];molmap={};coord={};stereos={};specs={};expected=[]
 for r in records:
  m=r['metadata'];stage=m['group'];ident=dict(stage=stage,receptor_id=m['receptor_id'],seed=int(m['seed']),candidate_rank=int(m['candidate_rank']));mol,info,fids=load_molecule(m);heavy,ff,_=heavy_copy(mol,fids);xyz=heavy.GetConformer().GetPositions();z=[a.GetAtomicNum() for a in heavy.GetAtoms()];protein,cp=cloud_metrics(xyz,z,m['protein_path']);molmap[r['molecule_id']]=heavy;coord[stage,key(r)]=xyz
  assert protein['protein_pass']==r['checks'].get('no_severe_protein_clash')
  stpath=Path(m['sdf_path']).with_name('stereo_decisions.json')
  if stpath.exists():st=load(stpath)
  else:_,_,st=derive(mol,load(m['chemistry_spec_path']))
  stereos[stage,key(r)]=st;specs[stage,key(r)]=load(m['chemistry_spec_path']);g=r['modules'].get('geometry',{});pb=r['modules'].get('posebusters',{})
  assert not WATER_TESTS&set(pb.get('tests',{}))
  row=dict(**ident,HQ_no_water=r['checks'].get('high_quality'),**protein,chemical_valid=r['checks'].get('chemical_valid'),fully_connected=r['checks'].get('fully_connected'),**{'QC_'+k:v for k,v in r['checks'].items()},abnormal_bonds=g.get('abnormal_bond_count'),abnormal_angles=g.get('abnormal_angle_count'),internal_severe_pairs=g.get('severe_internal_clash_pairs'),stereo_reliable=st['reliable'],stereo_unavailable=st['unavailable'],PB_failed=[k for k,v in pb.get('tests',{}).items() if v is False],PB_unavailable=[k for k,v in pb.get('tests',{}).items() if v is None],PB_status=pb.get('status'),PB_errors=pb.get('errors'),failed_checks=r.get('failed_checks'),unavailable_checks=r.get('unavailable_checks'),errors=r.get('errors'),sdf_path=m['sdf_path'])
  for relation in ['cross_fragment','within_fragment']:
   bb=[b for b in g.get('bonds',[]) if b['fragment_relation']==relation];row[relation+'_abnormal_bonds']=sum(b['abnormal'] for b in bb);row[relation+'_bond_MAE_A']=float(np.mean([abs(b['deviation_A']) for b in bb])) if bb else None
  repair_path=Path(m['sdf_path']).with_name('repair.json')
  if repair_path.exists():
   o=load(repair_path);row.update(optimizer_converged=o['success'],optimizer_iterations=o['iterations'],max_local_displacement_A=o.get('max_local_displacement_A',0),optimizer_termination=o['termination'])
  per.append(row)
  for b in g.get('bonds',[]):bonds.append(dict(**ident,**b))
  for a in g.get('angles',[]):angles.append(dict(**ident,fragment_relation='within_fragment' if len({ff[i] for i in a['atoms']})==1 else 'cross_fragment',**a))
  clashes.extend(dict(**ident,**p) for p in cp)
 for cloud in raw_clouds:
  ident={k:cloud[k] for k in ['receptor_id','seed','candidate_rank']};protein,cp=cloud_metrics(np.asarray(cloud['xyz']),cloud['z'],cloud['protein_path']);per.append(dict(stage='raw_predecode',**ident,HQ_no_water=None,chemical_valid=None,fully_connected=None,**protein,shared_raw_cloud_id=cloud.get('shared_raw_cloud_id'),graph_metrics_status='not_evaluable_no_connected_graph'));clashes.extend(dict(stage='raw_predecode',**ident,**p) for p in cp)
 for row in per:
  if row['stage']=='raw_predecode':continue
  k=key(row);x=coord[row['stage'],k]
  for ref in ['decoded_initial','native_final','rigid']:
   if (ref,k) not in coord or coord[ref,k].shape!=x.shape:continue
   dx=x-coord[ref,k];row['RMSD_from_'+ref+'_A']=float(np.sqrt(np.square(dx).sum(1).mean()));row['max_atom_from_'+ref+'_A']=float(np.linalg.norm(dx,axis=1).max());row['centroid_from_'+ref+'_A']=float(np.linalg.norm(dx.mean(0)))
   row['stereo_from_'+ref]=stereo_change(stereos[ref,k],stereos[row['stage'],k])
 summary={};by_receptor=[]
 for stage in ['raw_predecode']+STAGES:
  rr=[r for r in per if r['stage']==stage];s=dict(denominator=denominator,records=len(rr),missing=denominator-len(rr),rates={k:rate(rr,k,denominator) for k in ['protein_pass','HQ_no_water','chemical_valid','fully_connected']+[f'QC_{c}' for c in load(QC)['quality']['required_checks']]},abnormal_bonds=sum(r.get('abnormal_bonds') or 0 for r in rr),abnormal_angles=sum(r.get('abnormal_angles') or 0 for r in rr),protein_clash_pairs=sum(r.get('protein_clash_pairs') or 0 for r in rr),internal_severe_pairs=sum(r.get('internal_severe_pairs') or 0 for r in rr),stereo_reliable=sum(r.get('stereo_reliable') or 0 for r in rr),stereo_unavailable=sum(r.get('stereo_unavailable') or 0 for r in rr),optimizer_converged=sum(r.get('optimizer_converged') is True for r in rr),PB_failures=Counter(x for r in rr for x in r.get('PB_failed',[])))
  for k in ['RMSD_from_native_final_A','max_atom_from_native_final_A','RMSD_from_rigid_A','RMSD_from_decoded_initial_A','centroid_from_native_final_A']:s[k]=distribution([r.get(k) for r in rr])
  if stage!='raw_predecode':s['diversity']=diversity([r for r in records if r['metadata']['group']==stage],molmap,denominator)
  if stage=='raw_predecode':
   for field in ['abnormal_bonds','abnormal_angles','internal_severe_pairs','stereo_reliable','stereo_unavailable','optimizer_converged']:s[field]=None
  summary[stage]=s
  for pid in sorted({r['receptor_id'] for r in per}):
   pp=[r for r in rr if r['receptor_id']==pid];n=24 if label=='paired_oracle' else 3;by_receptor.append(dict(stage=stage,receptor_id=pid,denominator=n,records=len(pp),protein_pass=sum(r['protein_pass'] is True for r in pp),HQ_no_water=None if stage=='raw_predecode' else sum(r.get('HQ_no_water') is True for r in pp),chemical_valid=None if stage=='raw_predecode' else sum(r.get('chemical_valid') is True for r in pp)))
 pairs=[];index={(r['stage'],key(r)):r for r in per}
 for a,b in [('raw_predecode','decoded_initial'),('decoded_initial','native_final'),('native_final','rigid'),('native_final','local_flex'),('rigid','local_flex')]:
  for metric in ['protein_pass','HQ_no_water']:
   values=[]
   keys={key(r) for r in per}
   for k in keys:
    x=index.get((a,k),{}).get(metric);y=index.get((b,k),{}).get(metric);values.append((x,y))
   pairs.append(dict(before=a,after=b,metric=metric,denominator=denominator,both_pass=sum(x is True and y is True for x,y in values),gain=sum(x is False and y is True for x,y in values),loss=sum(x is True and y is False for x,y in values),both_fail=sum(x is False and y is False for x,y in values),unavailable=denominator-sum(isinstance(x,(bool,np.bool_)) and isinstance(y,(bool,np.bool_)) for x,y in values)))
 write_json(target/'summary.json',summary);savecsv(target/'per_candidate_results.csv',per);savecsv(target/'by_receptor.csv',by_receptor);savecsv(target/'paired_transitions.csv',pairs);savecsv(target/'bond_details.csv',bonds);savecsv(target/'angle_details.csv',angles);savecsv(target/'protein_clashes.csv',clashes);write_json(target/'per_candidate_results.json',per);write_json(target/'paired_transitions.json',pairs)
 return summary

def main():
 records=jsonl(OUT/'quality/molecule_metrics.jsonl');clouds=[]
 missing=[]
 for r in records:
  pb=r['modules'].get('posebusters',{});names=[k for k,v in pb.get('tests',{}).items() if v is None]
  if names:missing.append(dict(stage=r['metadata']['group'],id=r['molecule_id'],missing=names,errors=pb.get('errors'),HQ=r['checks'].get('high_quality')))
 write_json(OUT/'PB_unavailable_records.json',missing)
 native={key(r):r for r in load(ROOT/'fbdd_efficiency_generalization_20261103/independent/native/raw_manifest.json')}
 for r in load(OUT/'raw_cloud_manifest.json'):
  xyz=np.load(r['prefix_path'])['xyz'];mol,_,_=load_molecule(native[key(r)]);z=[mol.GetAtomWithIdx(a).GetAtomicNum() for g in r['groups'] for a in g];clouds.append(dict(r,xyz=xyz,z=z))
 s=cohort('paired_oracle',records,clouds,144);print('MAIN', {k:(v['rates']['protein_pass']['pass_count'],v['rates']['HQ_no_water']['pass_count']) for k,v in s.items()},flush=True)
 if (OUT/'free_probe/quality/molecule_metrics.jsonl').exists():
  attempts=load(OUT/'free_probe/attempts.json');clouds=[]
  for r in attempts:
   p=OUT/'free_probe/structures'/f"{r['receptor_id']}_{r['seed']}"/'raw_cloud.npz'
   if p.exists():
    a=np.load(p);clouds.append(dict(receptor_id=r['receptor_id'],seed=r['seed'],candidate_rank=1,xyz=a['xyz'],z=a['z'],protein_path=str(OUT/'pockets'/f"{r['receptor_id']}_no_water.pdb")))
  f=cohort('free_probe',jsonl(OUT/'free_probe/quality/molecule_metrics.jsonl'),clouds,18,attempts);print('FREE',{k:(v['rates']['protein_pass']['pass_count'],v['rates']['HQ_no_water']['pass_count']) for k,v in f.items()},flush=True)
if __name__=='__main__':main()
