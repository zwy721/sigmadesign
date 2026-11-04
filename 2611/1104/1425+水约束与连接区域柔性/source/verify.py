"""Full-cohort regression without sampling, reference alignment or retuning."""
from study_common import *
from rdkit import Chem, RDLogger
from denovo_quality.io import load_molecule
from geometry_constraints_20261103.stereo import derive
RDLogger.DisableLog('rdApp.*')


def graph(m):
    return ([tuple((a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsotope(),a.GetNumRadicalElectrons())) for a in m.GetAtoms()], sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),str(b.GetBondType())) for b in m.GetBonds()))


def main():
    checks={}
    def check(label,value):
        checks[label]=bool(value)
        assert value,label
    fingerprints=load(OUT/'input_fingerprints.json')
    check('all_historical_input_hashes_unchanged',all(Path(p).is_file() and sha256(p)==s for p,s in fingerprints.items()))
    lock=load(OUT/'execution_lock.json')
    check('fixed_parameters_unchanged',sha256(CONFIG)==lock['fixed_parameters_sha256'])
    check('numerical_sources_unchanged',all(sha256(ROOT/p)==s for p,s in lock['numerical_sources'].items()))
    check('QC_file_byte_identical',sha256(QC)==sha256(OUT/'quality_config.json')==load(CONFIG)['quality_config_sha256'])
    for exp in ['A','B']:
        check(f'QC_{exp}_configuration_identical',load(OUT/f'quality_{exp}/config_used.json')==load_config(QC))
        check(f'PB_{exp}_configuration_identical',load(OUT/f'quality_{exp}/posebusters_effective_config.json')==load(OLD/'quality/posebusters_effective_config.json'))
    rows=load(OUT/'input_manifest.json');raw={key(r):r for r in rows}
    metrics=jsonl(OUT/'quality_A/molecule_metrics.jsonl')+jsonl(OUT/'quality_B/molecule_metrics.jsonl')
    processing=load(OUT/'processing_A.json')+load(OUT/'processing_B.json')
    check('144_original_candidates',len(rows)==144 and len(raw)==144)
    check('432_processing_records_including_failures',len(processing)==432)
    check('no_new_neural_NFE',load(CONFIG)['additional_neural_NFE']==0)
    check('all_optimizers_produced_structures',all(r['status']=='ok' for r in processing))
    check('all_432_QC_records_retained',len(metrics)==432)
    check('all_144_pairs_for_each_variant',all({key(r) for r in metrics if r['metadata']['group']==v}==set(raw) for v in ['A0_B0','A1','B1']))
    previous={key(r):r for r in jsonl(OLD/'quality/molecule_metrics.jsonl') if r['metadata']['group']=='stereo'}
    controls=[r for r in metrics if r['metadata']['group']=='A0_B0']
    check('baseline_strict_7_of_144',sum(r['checks']['high_quality'] is True for r in controls)==7)
    check('baseline_auxiliary_103_of_144',sum(aux_hq(r) for r in controls)==103)
    check('baseline_water_only_96',sum(aux_hq(r) and r['checks']['high_quality'] is not True for r in controls)==96)
    req=load(QC)['quality']['required_checks']
    check('baseline_required_QC_exact_per_candidate',all(all(r['checks'].get(k)==previous[key(r)]['checks'].get(k) for k in req+['high_quality']) for r in controls))
    check('baseline_PB_exact_per_candidate',all(r['modules']['posebusters']['tests']==previous[key(r)]['modules']['posebusters']['tests'] for r in controls))
    pp=[r for r in processing if r['variant']=='A0_B0']
    check('baseline_coordinates_serialized_exact',all(r['baseline_serialized_coordinate_max_difference_A']==0. for r in pp))
    check('baseline_stereo_spec_exact',all(r['baseline_stereo_spec_equal'] for r in pp))
    graph_ok=spec_ok=coords_ok=stereo_ok=separate=True
    max_local=0.;max_rigid=0.;max_fixed=0.
    expected_waters=load(OUT/'fixed_waters.json')
    for r in metrics:
        rr=raw[key(r)];variant=r['metadata']['group'];folder=Path(r['metadata']['sdf_path']).parent
        m,info,fids=load_molecule(r['metadata']);x,_,xfids=load_molecule(rr)
        graph_ok &= graph(m)==graph(x) and fids==xfids
        original=OUT/'structures/original'/name(rr)
        spec_ok &= sha256(original/'raw.sdf')==sha256(rr['sdf_path']) and sha256(original/'spec.json')==sha256(rr['chemistry_spec_path'])
        repaired,_,_=load_molecule(dict(rr,sdf_path=str(folder/'repaired.sdf')))
        q,ds,decisions=derive(repaired,load(rr['chemistry_spec_path']))
        stereo_ok &= ds==load(folder/'derived_spec.json') and decisions==load(folder/'stereo_decisions.json')
        coords_ok &= np.array_equal(repaired.GetConformer().GetPositions(),m.GetConformer().GetPositions())
        opt=load(folder/'repair.json')
        check_label='parameters_match'
        separate &= opt['parameters']==load(CONFIG)['baseline_parameters']
        if variant=='A1':
            separate &= opt['water_enabled'] and not opt['local_enabled']
            separate &= opt['water_ids']==[w['water_id'] for w in expected_waters if w['receptor_id']==rr['receptor_id']]
            max_rigid=max(max_rigid,opt['max_fragment_pair_distance_change_A'])
        elif variant=='B1':
            separate &= opt['local_enabled'] and not opt['water_enabled']
            max_local=max(max_local,opt['max_local_displacement_A'])
            max_fixed=max(max_fixed,opt['fixed_region_pair_distance_error_A'])
        else:
            max_rigid=max(max_rigid,opt['max_fragment_pair_distance_change_A'])
    check('connectivity_bond_orders_atom_numbering_and_fragment_partitions_unchanged',graph_ok)
    check('all_raw_SDF_and_spec_copies_byte_exact',spec_ok)
    check('same_stereo_assignment_reproducible_all_432',stereo_ok)
    check('stereo_step_never_moves_coordinates_all_432',coords_ok)
    check('single_factors_separate_and_all_waters_retained',separate)
    check('A0_A1_fragments_rigid_to_1e_minus_8_A',max_rigid<1e-8)
    check('B1_nonlocal_region_rigid_to_1e_minus_8_A',max_fixed<1e-8)
    check('B1_local_displacements_at_most_0p10_A',max_local<=.10+1e-12)
    check('no_primary_QC_unavailable',all(all(r['checks'].get(k) is not None for k in req) for r in metrics))
    result=dict(passed=all(checks.values()),checks=checks,check_count=len(checks),historical_fingerprints=len(fingerprints),candidates=144,QC_records=len(metrics),max_local_displacement_A=max_local,max_rigid_pair_distance_error_A=max_rigid,max_fixed_region_pair_distance_error_A=max_fixed,additional_neural_NFE=0,unconverged_retained=sum(not load(Path(r['metadata']['sdf_path']).with_name('repair.json'))['success'] for r in metrics))
    write_json(OUT/'verification.json',result)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
