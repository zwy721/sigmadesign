"""Attempt-denominated statistics. No output deduplication or quality selection."""
from common import *
from collections import Counter,defaultdict
from itertools import combinations
import argparse
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import rdFingerprintGenerator,rdMolAlign
from rdkit.Chem.Scaffolds import MurckoScaffold
from denovo_quality.diversity import training_reference,structural_fingerprint,distribution
RDLogger.DisableLog('rdApp.*')
STAGES=['decoded_initial','native_final','rigid','local_flex']
def avg(v):
 v=[x for x in v if x is not None];return float(np.mean(v)) if v else None
def H(v):
 p=np.array(list(Counter(v).values()),float);p/=p.sum() if len(p) else 1
 return float(-(p*np.log(p)).sum()) if len(p) else None
def tree_key(r):
 adj=defaultdict(list)
 for e in r.get('connections',[]):adj[e['i']].append(e['j']);adj[e['j']].append(e['i'])
 def rooted(i,p):return '('+''.join(sorted(rooted(j,i) for j in adj[i] if j!=p))+')'
 return min(rooted(i,-1) for i in range(r['N'])) if len(r.get('connections',[]))==r['N']-1 else None
def merge(phase):
 rows=[load(p) for p in sorted((OUT/phase/'structures').glob('*/result.json'))]
 write_json(OUT/phase/'attempts.json',rows);write_json(OUT/phase/'manifest.json',[m for r in rows for m in r['manifest']]);return rows
def main():
 p=argparse.ArgumentParser();p.add_argument('--phase',default='main');p.add_argument('--merge-only',action='store_true');a=p.parse_args();rows=merge(a.phase)
 if a.merge_only:print(len(rows),'attempts merged');return
 plan=load(OUT/'frozen_plan.json');qc=jsonl(OUT/a.phase/'quality/molecule_metrics.jsonl');by={r['molecule_id']:r for r in qc};generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
 exact,training,training_meta=training_reference(ROOT/'de_novo_quality_20261103/training_index.json');write_json(OUT/'novelty_reference.json',training_meta)
 pergen=[];molcache={};checks=[];bonds=[];angles=[];idstats=[];events=[]
 for r in rows:
  row={k:r.get(k) for k in ['attempt_id','receptor_id','seed','strategy','N','status','failure_stage','error','seconds','identity_queries','connection_queries','pose_NFE','initial_labels','labels','repair_errors']}
  events.extend(dict(attempt_id=r['attempt_id'],strategy=r['strategy'],receptor_id=r['receptor_id'],seed=r['seed'],**e) for e in r['events'])
  row.update(initial_fragment_combination=sorted(r.get('initial_labels',[])),terminal_fragment_combination=sorted(r.get('labels',[])),identity_changes=sum(e['old_labels']!=e['proposed_labels'] for e in r['events']),accepted_revisions=sum(e['accepted'] for e in r['events']),rejected_revisions=sum(not e['accepted'] and e.get('reason')!='same_chemistry' for e in r['events']),unlabelled_tree=tree_key(r))
  search_path=OUT/a.phase/'structures'/r['attempt_id']/'search_00.json'
  if search_path.exists():
   search=load(search_path);row.update(initial_raw_head_complete_tree=search['raw_head_complete'],initial_raw_head_complete_chemical_graph=search['raw_complete_chemical_graph'],initial_search_legal_candidates=len(search['candidates']),initial_search_pool_size=len(search['pool']))
  for stats in r['identity_statistics']:idstats.append(dict(attempt_id=r['attempt_id'],strategy=r['strategy'],receptor_id=r['receptor_id'],seed=r['seed'],**stats))
  for stage in STAGES:
   mid=r['attempt_id']+'__'+stage;q=by.get(mid);c=q.get('checks',{}) if q else {};chem=q.get('modules',{}).get('chemistry',{}) if q else {};geom=q.get('modules',{}).get('geometry',{}) if q else {};meta=next((m for m in r['manifest'] if m['group']==stage),None)
   for ck in ['chemical_valid','fully_connected','high_quality','no_severe_protein_clash','no_severe_internal_clash','bond_geometry_valid','angle_geometry_valid','stereo_definition_complete','stereo_geometry_consistent','posebusters_all_pass']:
    row[f'{stage}_{ck}']=c.get(ck)
   row[f'{stage}_constitution_smiles']=chem.get('constitution_smiles');row[f'{stage}_canonical_smiles']=chem.get('canonical_smiles');row[f'{stage}_abnormal_cross_bonds']=geom.get('cross_fragment_abnormal_bond_count');row[f'{stage}_abnormal_bonds']=geom.get('abnormal_bond_count');row[f'{stage}_abnormal_angles']=geom.get('abnormal_angle_count')
   if q:
    checks.append(dict(attempt_id=r['attempt_id'],stage=stage,checks=c,errors=q['errors'],unavailable=q['unavailable_checks']))
    for v in geom.get('bonds',[]):
     if v['abnormal']:bonds.append(dict(attempt_id=r['attempt_id'],stage=stage,**v))
    spec=load(meta['chemistry_spec_path']) if meta else {}
    fragment_of={k:f['slot'] for f in spec.get('fragments',[]) for k in f['atom_ids']}
    cross_angles=0
    for v in geom.get('angles',[]):
     cross=len({fragment_of.get(k,-1) for k in v['atoms']})>1
     if v['abnormal']:
      angles.append(dict(attempt_id=r['attempt_id'],stage=stage,cross_fragment=cross,**v));cross_angles+=cross
    row[f'{stage}_abnormal_cross_angles']=cross_angles
   if meta and c.get('chemical_valid') is True:
    m=Chem.SDMolSupplier(meta['sdf_path'],removeHs=True)[0]
    if m is not None:
     molcache[mid]=m
     row[f'{stage}_scaffold']=MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(chem['constitution_smiles']),includeChirality=False)
     row[f'{stage}_novel_vs_partial_training']=chem.get('constitution_smiles') not in training if training else None
   if stage in ['rigid','local_flex'] and meta:
    rm=load(Path(meta['sdf_path']).parent/'repair.json');row[f'{stage}_repair_RMSD_A']=rm.get('heavy_RMSD_no_alignment_A');row[f'{stage}_max_atom_displacement_A']=rm.get('max_atom_displacement_A');row[f'{stage}_optimizer_success']=rm.get('success')
  pergen.append(row)
 savecsv(OUT/'PER_GENERATION_RESULTS.csv',pergen);write_json(OUT/'PER_GENERATION_RESULTS.json',pergen);savecsv(OUT/'ABNORMAL_BONDS.csv',bonds);savecsv(OUT/'ABNORMAL_ANGLES.csv',angles);write_json(OUT/'all_quality_checks.json',checks);write_json(OUT/'identity_statistics.json',idstats);savecsv(OUT/'ALL_UPDATE_EVENTS.csv',events);savecsv(OUT/'REJECTED_UPDATE_EVENTS.csv',[e for e in events if not e['accepted'] and e.get('reason')!='same_chemistry'])
 def summarize(rr,strategy,pocket):
  n=len(rr);initial=[x for r in rr for x in r.get('initial_labels') or []];terminal=[x for r in rr for x in r.get('labels') or []]
  s=dict(strategy=strategy,receptor_id=pocket,attempts=n,neural_outputs=sum(r['native_final_constitution_smiles'] is not None for r in rr),failure_stages=dict(Counter(r['failure_stage'] for r in rr if r['failure_stage'])),N_distribution=dict(Counter(r['N'] for r in rr)),identity_unique_initial=len(set(initial)),identity_unique_terminal=len(set(terminal)),identity_entropy_initial=H(initial),identity_entropy_terminal=H(terminal),vocabulary_coverage_initial=len(set(initial))/3684,vocabulary_coverage_terminal=len(set(terminal))/3684,identity_token_repeat_fraction_initial=1-len(set(initial))/len(initial) if initial else None,benzene_fraction_initial=initial.count(3577)/len(initial) if initial else None,benzene_fraction_terminal=terminal.count(3577)/len(terminal) if terminal else None,fragment_combination_unique_initial=len({tuple(r['initial_fragment_combination']) for r in rr}),fragment_combination_unique_terminal=len({tuple(r['terminal_fragment_combination']) for r in rr}),unlabelled_trees=len({r['unlabelled_tree'] for r in rr if r['unlabelled_tree']}),seconds_total=sum(r['seconds'] for r in rr),seconds_per_attempt=avg([r['seconds'] for r in rr]),pose_NFE=sum(r['pose_NFE'] for r in rr),connection_queries=sum(r['connection_queries'] for r in rr),identity_queries=sum(r['identity_queries'] for r in rr),accepted_revisions=sum(r['accepted_revisions'] for r in rr),rejected_revisions=sum(r['rejected_revisions'] for r in rr))
  s.update(initial_raw_head_complete_tree=sum(r.get('initial_raw_head_complete_tree') is True for r in rr),initial_search_success=sum(bool(r.get('initial_search_legal_candidates')) for r in rr),prior_eligible_vocabulary=load(OUT/'generator_provenance.json')['support'])
  selected_ids={r['attempt_id'] for r in rr};ii=[x for x in idstats if x['attempt_id'] in selected_ids]
  for label,qq in [('initial',[x for x in ii if x['step']==0]),('all_queries',ii)]:
   for field in ['entropy','sampling_entropy','top1_probability','benzene_probability']:s[f'{label}_{field}_mean']=avg([v for x in qq for v in x[field]])
   top=[v for x in qq for v in x['top1_class']];s[f'{label}_benzene_top1_fraction']=top.count(3577)/len(top) if top else None
  s['initial_actual_draw_entropy']=load(OUT/'generator_provenance.json')['entropy_nats'] if plan['strategies'][strategy]['initial']=='frequency_prior' else s['initial_sampling_entropy_mean']
  for stage in STAGES:
   valid=[r for r in rr if r[f'{stage}_chemical_valid'] is True and r[f'{stage}_fully_connected'] is True and r[f'{stage}_constitution_smiles']];const=[r[f'{stage}_constitution_smiles'] for r in valid];isom=[r[f'{stage}_canonical_smiles'] for r in valid]
   fps=[structural_fingerprint(molcache[r['attempt_id']+'__'+stage],generator) for r in valid];d=[1-DataStructs.TanimotoSimilarity(fps[i],fps[j]) for i,j in combinations(range(len(fps)),2)]
   unique_const=sorted(set(const));ufps=[generator.GetFingerprint(Chem.MolFromSmiles(v)) for v in unique_const];ud=[1-DataStructs.TanimotoSimilarity(ufps[i],ufps[j]) for i,j in combinations(range(len(ufps)),2)]
   for ck in ['chemical_valid','fully_connected','high_quality','no_severe_protein_clash','no_severe_internal_clash','stereo_definition_complete','posebusters_all_pass']:
    s[f'{stage}_{ck}_count']=sum(r[f'{stage}_{ck}'] is True for r in rr);s[f'{stage}_{ck}_unavailable']=sum(r[f'{stage}_{ck}'] is None for r in rr)
   s.update({f'{stage}_valid_connected':len(valid),f'{stage}_unique_constitution':len(set(const)),f'{stage}_unique_isomeric':len(set(isom)),f'{stage}_valid_unique_per_attempt':len(set(const))/n if n else None,f'{stage}_repeat_construct_fraction':1-len(set(const))/len(valid) if valid else None,f'{stage}_cross_seed_pair_repeat':avg([const[i]==const[j] for i,j in combinations(range(len(valid)),2) if valid[i]['seed']!=valid[j]['seed']]),f'{stage}_scaffold_unique':len({r.get(f'{stage}_scaffold') for r in valid}),f'{stage}_scaffold_unique_fraction':len({r.get(f'{stage}_scaffold') for r in valid})/len(valid) if valid else None,f'{stage}_internal_diversity':avg(d),f'{stage}_unique_structure_internal_diversity':avg(ud),f'{stage}_HQ_unique':len({r[f'{stage}_constitution_smiles'] for r in rr if r[f'{stage}_high_quality'] is True}),f'{stage}_novel_unique_partial_training':sum(v not in training for v in unique_const) if training else None})
   for field in ['abnormal_cross_bonds','abnormal_cross_angles','abnormal_bonds','abnormal_angles']:s[f'{stage}_{field}_total']=sum(r.get(f'{stage}_{field}') or 0 for r in rr)
   if stage in ['rigid','local_flex']:
    for field in ['repair_RMSD_A','max_atom_displacement_A']:s[f'{stage}_{field}']=distribution([r.get(f'{stage}_{field}') for r in rr])
    s[f'{stage}_optimizer_success_count']=sum(r.get(f'{stage}_optimizer_success') is True for r in rr)
  return s
 pockets=[];summary=[]
 for strategy in plan['strategies']:
  sr=[r for r in pergen if r['strategy']==strategy];summary.append(summarize(sr,strategy,'ALL'))
  for pid in plan['cohort']:pockets.append(summarize([r for r in sr if r['receptor_id']==pid],strategy,pid))
 savecsv(OUT/'PER_POCKET_DIVERSITY.csv',pockets);savecsv(OUT/'STRATEGY_SUMMARY.csv',summary)
 nstrata=[]
 for strategy in plan['strategies']:
  for pocket in ['ALL']+plan['cohort']:
   for count in sorted({r['N'] for r in pergen}):
    rr=[r for r in pergen if r['strategy']==strategy and r['N']==count and (pocket=='ALL' or r['receptor_id']==pocket)]
    if rr:
     ns=summarize(rr,strategy,pocket);ns['N_stratum']=count;nstrata.append(ns)
 savecsv(OUT/'PER_N_DIVERSITY.csv',nstrata);write_json(OUT/'PER_N_DIVERSITY.json',nstrata)
 cross=[]
 for strategy in plan['strategies']:
  sr=[r for r in pergen if r['strategy']==strategy]
  for p1,p2 in combinations(plan['cohort'],2):
   left={r['seed']:r for r in sr if r['receptor_id']==p1};right={r['seed']:r for r in sr if r['receptor_id']==p2};seeds=sorted(set(left)&set(right));validseeds=[s for s in seeds if left[s]['native_final_constitution_smiles'] and right[s]['native_final_constitution_smiles']]
   L={r['native_final_constitution_smiles'] for r in left.values() if r['native_final_constitution_smiles']};R={r['native_final_constitution_smiles'] for r in right.values() if r['native_final_constitution_smiles']};deltas=[]
   for seed in seeds:
    aa=OUT/a.phase/'structures'/left[seed]['attempt_id']/'identity_00.npz';bb=OUT/a.phase/'structures'/right[seed]['attempt_id']/'identity_00.npz';deltas.append(float(np.max(np.abs(np.load(aa)['model_probabilities']-np.load(bb)['model_probabilities']))))
   cross.append(dict(strategy=strategy,pocket_1=p1,pocket_2=p2,paired_seeds=len(seeds),initial_identity_equal_fraction=avg([left[s]['initial_labels']==right[s]['initial_labels'] for s in seeds]),initial_probability_max_abs_delta=max(deltas,default=None),terminal_identity_equal_fraction=avg([left[s]['labels']==right[s]['labels'] for s in seeds]),valid_seed_pairs=len(validseeds),same_seed_constitution_equal_fraction=avg([left[s]['native_final_constitution_smiles']==right[s]['native_final_constitution_smiles'] for s in validseeds]),constitution_Jaccard=len(L&R)/len(L|R) if L|R else None))
 savecsv(OUT/'CROSS_POCKET_COMPARISON.csv',cross)
 conformers=[]
 for strategy in plan['strategies']:
  for pid in plan['cohort']:
   for stage in ['native_final','rigid','local_flex']:
    groups=defaultdict(list)
    for r in pergen:
     if r['strategy']==strategy and r['receptor_id']==pid and r.get(f'{stage}_canonical_smiles') and r['attempt_id']+'__'+stage in molcache:groups[r[f'{stage}_canonical_smiles']].append(r)
    for smiles,rr in groups.items():
     aligned=[];world=[];errors=[]
     for x,y in combinations(rr,2):
      try:
       m=Chem.Mol(molcache[x['attempt_id']+'__'+stage]);n=Chem.Mol(molcache[y['attempt_id']+'__'+stage]);matches=m.GetSubstructMatches(n,useChirality=True,uniquify=False,maxMatches=4096);xx=m.GetConformer().GetPositions();yy=n.GetConformer().GetPositions()
       if not matches:raise ValueError('No stereo-preserving atom mapping')
       world.append(min(float(np.sqrt(np.mean(np.sum((xx[list(match)]-yy)**2,-1)))) for match in matches));aligned.append(float(rdMolAlign.GetBestRMS(m,n,maxMatches=4096)))
      except Exception as exc:errors.append(repr(exc))
     if len(rr)>1:conformers.append(dict(strategy=strategy,receptor_id=pid,stage=stage,canonical_smiles=smiles,outputs=len(rr),aligned_RMSD_A=distribution(aligned),receptor_frame_RMSD_A=distribution(world),errors=errors,symmetry_match_cap=4096))
 write_json(OUT/'CONFORMER_DIVERSITY.json',conformers);savecsv(OUT/'CONFORMER_DIVERSITY.csv',conformers)
 data=dict(attempts=len(rows),summary=summary,pockets=pockets,cross_pocket=cross,training_reference=training_meta,conformer_groups=len(conformers),rejection_reasons=dict(Counter(e.get('reason') for e in events if not e['accepted'] and e.get('reason')!='same_chemistry')));write_json(OUT/'summary.json',data)
 print('analysis complete',len(rows),flush=True)
if __name__=='__main__':main()

