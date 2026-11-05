"""Independent integrity checks, attempt accounting and sampler boundary tests."""
from common import *
import itertools,ast
from collections import Counter
from denovo_quality.common import conjunction
def main():
 results=[];failures=[]
 def check(name,fn):
  try:
   detail=fn();results.append(dict(name=name,passed=True,detail=detail))
  except Exception as exc:
   failures.append(name);results.append(dict(name=name,passed=False,error=repr(exc)))
 plan=load(OUT/'frozen_plan.json');rows=load(OUT/'main/attempts.json');pilot=load(OUT/'pilot/attempts.json');manifest=load(OUT/'main/manifest.json');qc=jsonl(OUT/'main/quality/molecule_metrics.jsonl')
 def assert_true(v,msg=''):
  assert v,msg;return True
 def cohort():
  expected=set(itertools.product(plan['cohort'],plan['seeds'],plan['strategies']));actual=[(r['receptor_id'],r['seed'],r['strategy']) for r in rows];assert len(actual)==len(set(actual))==len(expected) and set(actual)==expected;return len(actual)
 check('exact_frozen_attempt_grid_no_duplicates',cohort)
 check('pilot_complete_separate_denominator',lambda:assert_true(len(pilot)==14 and all(r['phase']=='pilot' for r in pilot) and all(r['phase']=='main' for r in rows)))
 check('N_is_predefined_prior_only',lambda:assert_true(all(r['N']==int(np.random.default_rng(r['seed']+912711).choice([3,4,5,6])) for r in rows)))
 def budgets():
  complete=[r for r in rows if r.get('neural_success')];assert all((r['identity_queries'],r['connection_queries'],r['pose_NFE'])==(5,5,24) for r in complete);return dict(complete=len(complete),early_failed=len(rows)-len(complete))
 check('equal_complete_attempt_neural_and_search_budgets',budgets)
 def coverage():
  expected={m['molecule_id'] for m in manifest};actual={r['molecule_id'] for r in qc};missing=[m for m in manifest if m['molecule_id'] not in actual]
  assert len(manifest)==len(expected) and len(qc)==len(actual) and actual<=expected
  assert all(m['group']=='decoded_initial' for m in missing)
  pending=list(csv.DictReader((OUT/'PENDING_QUALITY_RECORDS.csv').open()))
  assert {m['molecule_id'] for m in pending}=={m['molecule_id'] for m in missing}
  assert sum(r['metadata']['group'] in ['native_final','rigid','local_flex'] for r in qc)==len(rows)*3
  return dict(required_terminal_QC=len(rows)*3,completed_QC=len(qc),optional_initial_QC_pending=len(missing))
 check('required_terminal_QC_complete_optional_initial_pending_explicit',coverage)
 check('all_QC_inputs_unchanged',lambda:assert_true(all(r.get('input_unchanged') is True for r in qc)))
 check('HQ_conjunction_no_missing_as_pass',lambda:assert_true(all(r['checks']['high_quality']==conjunction(r['checks'].get(k) for k in load_config(QC)['quality']['required_checks']) for r in qc)))
 def protocol():
  from denovo_quality import posebusters_adapter as pb
  install_pb();current=pb.effective_config('dock',1);assert current==load(OUT/'posebusters_no_water.json')
  old=ROOT/'2611/1105/1039+no_water_test';assert sha256(QC)==plan['QC_config_sha256']==sha256(old/'quality_config.json')
  assert sha256(OUT/'posebusters_no_water.json')==plan['PB_config_sha256']==sha256(old/'posebusters_no_water.json')
  expected=set(pb.test_names(current));assert len(expected)==20
  for r in qc:
   tests=r['modules'].get('posebusters',{}).get('tests',{})
   assert not (set(tests)&WATER_TESTS)
  return sorted(expected)
 check('frozen_HQ_no_water_protocol_exact',protocol)
 check('weights_unchanged',lambda:assert_true(all(sha256(v['path'])==v['sha256'] for v in plan['weights'].values())))
 check('all_historical_files_unchanged',lambda:assert_true(all(sha256(p)==v for p,v in load(OUT/'historical_fingerprints.json').items())))
 check('generator_sources_templates_unchanged',lambda:assert_true(all(sha256(p)==v for p,v in load(OUT/'generator_provenance.json')['files'].items())))
 def artifacts():
  complete=0;rejected=0;probs=0;searches=0
  for phase,records in [('pilot',pilot),('main',rows)]:
   for r in records:
    folder=OUT/phase/'structures'/r['attempt_id']
    if r.get('neural_success'):
     audit=load(folder/'sampling_audit.json');assert len(audit['forwards'])==24 and audit['all_dynamic_recomputations_equal'] and audit['static_graph_unchanged'];complete+=1
    if r['status']=='failed':assert r.get('error') and r.get('traceback') and r.get('failure_stage')
    for p in folder.glob('identity_*.npz'):
     d=np.load(p);raw=d['model_probabilities'];q=d['sampling_probabilities'];assert np.isfinite(raw).all() and (raw>=0).all() and np.allclose(raw.sum(-1),1,atol=2e-6);assert np.isfinite(q).all() and (q>=0).all() and np.allclose(q.sum(-1),1,atol=2e-6);assert len(raw)==r['N'] and raw.shape[1]==3684;probs+=1
    for p in folder.glob('search_*.json'):
     s=load(p);assert not s['truth_used'];assert s['config']==plan['search_config'];searches+=1
    for e in r['events']:
     if not e['accepted'] and e.get('reason')!='same_chemistry':assert e.get('fallback')=='retain_last_legal_graph';rejected+=1
  return dict(complete_official_audits=complete,probability_arrays=probs,searches=searches,rejected_updates_retained=rejected)
 check('all_probabilities_searches_failures_and_official_audits',artifacts)
 def model_pockets():
  from denovo_quality.protein import read_protein
  for r in rows:
   if not r.get('neural_success'):continue
   p=OUT/'main/structures'/r['attempt_id']/'model_pocket.npz';d=np.load(p);atoms=[a for a in read_protein(str(ROOT/'2611/1105/1039+no_water_test/pockets'/f'{r["receptor_id"]}_no_water.pdb'))['atoms'] if a['kind']=='protein'];xyz=np.array([a['xyz'] for a in atoms]);assert d['xyz'].shape==xyz.shape and np.allclose(d['xyz'],xyz,atol=1e-5);assert np.allclose(d['center'].reshape(3),xyz.mean(0),atol=1e-5)
  return True
 check('actual_model_pockets_water_free_and_centers_protein_only',model_pockets)
 def boundary():
  server=(OUT/'source/identity_server.py').read_text();runner=(OUT/'source/run.py').read_text();tree=ast.parse(server)
  assert 'PoseDataset(' not in server and "assert not len(pg['pos'])" in server and "set(r)<=" in server
  assert "fixture.pt" not in runner and "shared_meta.json" not in runner and "ligand_path=None" in runner
  heads=[keyword.value.value for call in ast.walk(ast.parse(runner)) if isinstance(call,ast.Call) for keyword in call.keywords if keyword.arg=='head' and isinstance(keyword.value,ast.Constant)]
  assert heads and set(heads)=={'full'}
  return 'Static input contract plus runtime empty-protein assertions; no target dataset construction.'
 check('target_free_boundary_and_current_head',boundary)
 def sampling():
  import run
  logits=np.array([[3.,1.,-2.,-5.]])
  raw=run.softmax(logits);p,q=run.draw_distribution(logits,dict(draw='top1',temperature=1.));assert q.argmax()==raw.argmax() and q.max()==1
  p,q=run.draw_distribution(logits,dict(draw='topk',top_k=2,temperature=1.));assert (q>0).sum()==2 and np.allclose(q.sum(-1),1)
  p,q=run.draw_distribution(logits,dict(draw='topp',top_p=.95,temperature=1.));idx=np.argsort(-p[0]);keep=q[0,idx]>0;assert p[0,idx][keep].sum()>=.95 and p[0,idx][keep][:-1].sum()<.95
  p2,q=run.draw_distribution(logits,dict(draw='categorical',temperature=1.5));assert run.entropy(p2)>run.entropy(raw)
  # Proposals have positive support for every model class before truncation.
  assert (raw>0).all()
  return 'argmax, normalized top-k, inclusive top-p threshold, temperature entropy tests passed'
 check('probability_sampling_boundary_tests',sampling)
 def replay_identity():
  prior=np.load(OUT/'frequency_prior.npz')['probabilities'];paired={}
  for r in rows:
   raw=np.load(OUT/'main/structures'/r['attempt_id']/'identity_00.npz');rng=np.random.default_rng(r['seed']+17331)
   if plan['strategies'][r['strategy']]['initial']=='frequency_prior':expected=rng.choice(len(prior),r['N'],p=prior).tolist()
   else:expected=[int(rng.choice(len(q),p=q.astype(float)/q.sum(dtype=float))) for q in raw['sampling_probabilities']]
   assert expected==r['initial_labels']
   key=r['seed'];prob=raw['model_probabilities']
   if key in paired:assert np.allclose(paired[key],prob,atol=1e-6,rtol=1e-5)
   else:paired[key]=prob
  return 'Initial draws reproduced from saved probabilities/seed; raw probabilities agree across all strategies and pockets within atol=1e-6, rtol=1e-5. GPU roundoff is reported.'
 check('independent_initial_draw_replay_and_common_distribution',replay_identity)
 check('CSV_preserves_every_attempt',lambda:assert_true(len(list(csv.DictReader((OUT/'PER_GENERATION_RESULTS.csv').open())))==len(rows)))
 check('per_pocket_rows_complete',lambda:assert_true(len(list(csv.DictReader((OUT/'PER_POCKET_DIVERSITY.csv').open())))==len(plan['cohort'])*len(plan['strategies'])))
 check('each_initial_distribution_agrees_within_float_tolerance_across_pockets',lambda:assert_true(all(r['initial_probability_max_abs_delta']<=1e-6 and r['initial_identity_equal_fraction']==1 for r in load(OUT/'summary.json')['cross_pocket'])))
 write_json(OUT/'verification.json',dict(passed=len(results)-len(failures),total=len(results),failed=failures,optional_initial_QC_pending=len(manifest)-len(qc),checks=results))
 print('VERIFICATION',len(results)-len(failures),'/',len(results),failures,flush=True)
 if failures:raise SystemExit(1)
if __name__=='__main__':main()

