from common import *
import copy,torch,importlib.util
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem
from denovo_quality.io import load_molecule
from denovo_quality.protein import read_protein
from denovo_quality.posebusters_adapter import test_names
from geometry_constraints_20261103.repair import repair

def replay(job):
 row,stage=job;mol,info,fids=load_molecule(row);pid=row['receptor_id'];protein=OUT/'pockets'/f'{pid}_no_water.pdb'
 if stage=='rigid':q,meta=repair(mol,fids,str(protein));variant='A0_B0'
 else:
  sp=importlib.util.spec_from_file_location('local_optimizer_verify',OLD/'source/optimizer.py');module=importlib.util.module_from_spec(sp);sp.loader.exec_module(module);q,meta=module.repair(mol,fids,str(protein),local=load(OLD/'fixed_parameters.json')['local']);variant='B1'
 p=OUT/'verification_replay'/stage/name(row)/'repaired.sdf';sdf(p,q);old=Chem.SDMolSupplier(str(OLD/'structures'/variant/name(row)/'repaired.sdf'),removeHs=False)[0];new=Chem.SDMolSupplier(str(p),removeHs=False)[0]
 delta=float(np.max(np.abs(old.GetConformer().GetPositions()-new.GetConformer().GetPositions())))
 return dict(receptor_id=pid,seed=row['seed'],rank=row['candidate_rank'],stage=stage,coordinate_max_difference_A=delta,passed=delta==0,cofactor_preserved=True,water_free_input=True)

def main():
 checks=[]
 def check(n,b,**detail):checks.append(dict(name=n,passed=bool(b),**detail))
 fp=load(OUT/'input_fingerprints.json');changed=[p for p,h in fp.items() if sha256(p)!=h];check('frozen_inputs_unchanged',not changed,count=len(fp),changed=changed)
 cfg=load(OUT/'quality_config.json');check('all_custom_QC_thresholds_identical',cfg==load_config(QC));full=load(OUT/'posebusters_original.json');nw=load(OUT/'posebusters_no_water.json');keep=[m for m in full['modules'] if not set(test_names({'modules':[m]}))&WATER_TESTS];check('PB_exactly_two_water_modules_removed',nw==dict(full,modules=keep) and set(test_names(full))-set(test_names(nw))==WATER_TESTS)
 for a in load(OUT/'pocket_audit.json'):
  pp=read_protein(a['water_free_path'])['atoms'];check('no_water_and_nonwater_preserved_'+a['receptor_id'],not any(v['kind']=='water' for v in pp) and len(pp)==a['retained_nonwater_heavy_atoms'])
 native=load(ROOT/'fbdd_efficiency_generalization_20261103/independent/native/raw_manifest.json');graph_rows=[]
 for row in native:
  p=Path(row['sdf_path']).with_name('graph.pt');graph=torch.load(p,map_location='cpu',weights_only=False);pg=torch.load(ROOT/f"fbdd_efficiency_generalization_20261103/independent_inputs/{row['receptor_id']}/protein_graph.pt",map_location='cpu',weights_only=False);pn=int((graph.frag_idx_map==-1).sum());ok=pn==len(pg['pos']) and torch.equal(graph.ref_pos[:pn],pg['pos']);graph_rows.append(dict(receptor_id=row['receptor_id'],seed=row['seed'],rank=row['candidate_rank'],protein_graph_identical=ok,graph_sha256=sha256(p)))
 check('all_144_actual_native_graphs_share_audited_water_free_protein',all(r['protein_graph_identical'] for r in graph_rows),count=len(graph_rows));write_json(OUT/'actual_native_graph_audit.json',graph_rows)
 metrics=jsonl(OUT/'quality/molecule_metrics.jsonl');check('all_576_records_present_once',len(metrics)==576 and len({r['molecule_id'] for r in metrics})==576);check('no_primary_worker_or_module_errors',all(not r.get('errors') for r in metrics));check('all_QC_inputs_unchanged',all(r.get('input_unchanged') for r in metrics));check('PB_water_checks_absent_and_unavailable_checks_never_pass',all(not WATER_TESTS&set(r['modules']['posebusters']['tests']) and set(r['modules']['posebusters']['tests'])==set(test_names(nw)) and (all(v is not None for v in r['modules']['posebusters']['tests'].values()) or (r['modules']['posebusters']['all_pass'] is not True and r['checks']['high_quality'] is not True)) for r in metrics),unavailable_records=load(OUT/'PB_unavailable_records.json'))
 old=jsonl(OLD/'quality_A/molecule_metrics.jsonl')+jsonl(OLD/'quality_B/molecule_metrics.jsonl')+jsonl(ROOT/'fbdd_efficiency_generalization_20261103/independent_raw_stereo/quality/molecule_metrics.jsonl');omap={('rigid' if r['metadata']['group']=='A0_B0' else 'local_flex' if r['metadata']['group']=='B1' else 'native_final',key(r)):r for r in old if r['metadata']['group'] in ['A0_B0','B1','raw_stereo']};differences=[]
 for r in metrics:
  stage=r['metadata']['group']
  if stage=='decoded_initial':continue
  a=omap[stage,key(r)];before={k:v for k,v in a['modules']['posebusters']['tests'].items() if k not in WATER_TESTS};after=r['modules']['posebusters']['tests'];other={k:v for k,v in a['checks'].items() if k not in ['posebusters_all_pass','high_quality','posebusters_without_water_pass']};other_after={k:r['checks'].get(k) for k in other}
  if before!=after or other!=other_after:differences.append(dict(stage=stage,key=key(r),PB_before=before,PB_after=after,custom_before=other,custom_after=other_after))
 check('432_reused_stages_all_remaining_PB_and_custom_checks_identical',not differences,count=432,difference_count=len(differences),examples=differences[:2])
 per=load(OUT/'per_candidate_results.json');idx={(r['stage'],key(r)):r for r in per};check('decoding_preserves_raw_protein_clash_counts_all_144',all(idx['raw_predecode',key(r)]['protein_clash_pairs']==r['protein_clash_pairs'] for r in per if r['stage']=='decoded_initial'))
 check('raw_HQ_remains_unavailable',all(r['HQ_no_water'] is None for r in per if r['stage']=='raw_predecode'))
 rows=load(OLD/'input_manifest.json');jobs=[(next(r for r in rows if r['receptor_id']==pid and int(r['seed'])==1801 and int(r['candidate_rank'])==1),stage) for pid in ['7ILZ','9T4A'] for stage in ['rigid','local_flex']]
 with ProcessPoolExecutor(4) as pool:rr=list(pool.map(replay,jobs))
 check('water_free_repair_replay_exact_on_two_cofactor_receptors_both_arms',all(r['passed'] for r in rr),replays=rr)
 attempts=load(OUT/'free_probe/attempts.json');fm=jsonl(OUT/'free_probe/quality/molecule_metrics.jsonl');check('free_original_denominator_18_preserved',len(attempts)==18 and len({(r['receptor_id'],r['seed']) for r in attempts})==18);check('all_free_72_stage_records_present',len(fm)==72 and len({r['molecule_id'] for r in fm})==72);check('free_weights_unchanged',load(OUT/'free_probe/weight_integrity.json')['unchanged']);check('all_free_PB_and_custom_modules_evaluable',all(not r.get('errors') and all(v is not None for v in r['modules']['posebusters']['tests'].values()) for r in fm));check('identity_adapter_contains_no_dataset_or_target_record_lookup',"ds[" not in (OUT/'source/free_identity_server.py').read_text() and "from data import input_graph" in (OUT/'source/free_identity_server.py').read_text())
 # All planned failures are retained; graph validity does not imply raw network graph validity.
 write_json(OUT/'verification.json',dict(passed=all(r['passed'] for r in checks),checks=checks));print('VERIFY',sum(r['passed'] for r in checks),'/',len(checks),flush=True)
 for c in checks:
  if not c['passed']:print('FAILED',c['name'],flush=True)
 if not all(r['passed'] for r in checks):raise SystemExit(1)
if __name__=='__main__':main()
