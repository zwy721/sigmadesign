"""Read immutable inputs and freeze one descriptive protocol before optimization."""
from study_common import *
from collections import Counter
from datetime import datetime, timezone
from rdkit import Chem
from denovo_quality.io import load_molecule
from denovo_quality.common import versions
from denovo_quality.posebusters_adapter import effective_config
from geometry_constraints_20261103.repair import PARAMETERS
from geometry_constraints_20261103.stereo import PARAMETERS as STEREO
from posebusters.tools.protein import get_atom_type_mask


def main():
    assert not CONFIG.exists(), 'Protocol already frozen; do not overwrite this study.'
    rows = load(OLD / 'native/raw_manifest.json')
    previous = [r for r in jsonl(OLD / 'quality/molecule_metrics.jsonl') if r['metadata']['group'] == 'stereo']
    idx = {key(r): r for r in previous}
    assert len(rows) == len(idx) == 144
    assert len({key(r) for r in rows}) == 144
    hq = sum(r['checks']['high_quality'] is True for r in previous)
    aux = sum(aux_hq(r) for r in previous)
    water_only = sum(aux_hq(r) and r['checks']['high_quality'] is not True for r in previous)
    assert (hq, aux, water_only) == (7, 103, 96)
    assert load(OLD / 'config.json')['repair'] == PARAMETERS
    assert load(OLD / 'config.json')['stereo'] == STEREO
    assert sha256(QC) == load(OLD / 'config.json')['quality_config_sha256']
    fingerprints = {}
    audit = []
    for r in rows:
        paths = [r['sdf_path'], r['chemistry_spec_path'], r['protein_path']]
        old = Path(idx[key(r)]['metadata']['sdf_path']).parent
        paths += [str(old / f) for f in ['raw.sdf', 'repaired.sdf', 'derived.sdf', 'derived_spec.json', 'repair.json', 'stereo_decisions.json']]
        for p in paths:
            if not Path(p).is_file():
                raise FileNotFoundError(p)
            fingerprints[p] = sha256(p)
        mol, info, fids = load_molecule(r)
        assert info['sanitize_success'] and not info.get('declared_chemistry_error')
        assert all(a.GetAtomicNum() > 1 for a in mol.GetAtoms())
        spec = load(r['chemistry_spec_path'])
        assert len(fids) == mol.GetNumAtoms()
        audit.append(dict(receptor_id=r['receptor_id'], seed=r['seed'], rank=r['candidate_rank'], heavy_atoms=mol.GetNumAtoms(), fragments=len(set(fids)), cross_bonds=len(spec['connections']), raw_sdf=r['sdf_path'], original_spec=r['chemistry_spec_path'], protein=r['protein_path'], baseline_derived=str(old/'derived.sdf')))
    pb = effective_config('dock', 1)
    assert pb == load(OLD/'quality/posebusters_effective_config.json')
    dp = next(m['parameters'] for m in pb['modules'] if m['name'] == 'Distance to waters')
    catalog = list(csv.DictReader((DIAG/'water_catalog.csv').open()))
    water_rows = []
    for pid in sorted({r['receptor_id'] for r in rows}):
        r = next(r for r in rows if r['receptor_id'] == pid)
        p = Chem.MolFromPDBFile(r['protein_path'], removeHs=False, sanitize=False, proximityBonding=False)
        mask = get_atom_type_mask(p, dp['ignore_types'])
        ids = [i for i, yes in enumerate(mask) if yes]
        cat = [w for w in catalog if w['receptor_id'] == pid]
        assert ids == [int(w['prepared_atom_index']) for w in cat]
        for w in cat:
            i = int(w['prepared_atom_index'])
            xyz = p.GetConformer().GetPositions()[i]
            assert p.GetAtomWithIdx(i).GetAtomicNum() == 8
            np.testing.assert_array_equal(xyz, json.loads(w['xyz']))
            water_rows.append(dict(receptor_id=pid, water_id=w['water_id'], source_water_id=w['source_water_id'], prepared_atom_index=i, prepared_serial=int(w['prepared_serial']), xyz=xyz.tolist(), radius_A=float(w['radius_A'])))
    assert len(water_rows) == 79
    sources = [QC, OLD/'native/raw_manifest.json', OLD/'quality/molecule_metrics.jsonl', OLD/'quality/posebusters_effective_config.json', OLD/'config.json', DIAG/'water_catalog.csv']
    sources += list((ROOT/'denovo_quality').glob('*.py'))
    sources += [ROOT/'geometry_constraints_20261103/repair.py', ROOT/'geometry_constraints_20261103/stereo.py']
    for p in sources:
        fingerprints[str(p)] = sha256(p)
    protocol = dict(frozen_at_UTC=datetime.now(timezone.utc).isoformat(), purpose='Single-setting descriptive paired numerical repair on previously inspected test structures', additional_neural_NFE=0, training=False, reference_conformers_used=False, new_candidates=0,
        development_decision='Only one existing development receptor and previously inspected labels; insufficient independent development support. No sweep, no selection, no subsequent retuning.',
        baseline_parameters=PARAMETERS, stereo_parameters=STEREO, quality_config_sha256=sha256(QC),
        water=dict(weight_kcal_mol_A2=200., radius_ratio=.80, epsilon_A=1e-6, radius_definition='RDKit van der Waals radii, exactly the original PB water catalog; water O radius 1.55 A', potential='E=k*sum(max(ratio*(r_i+r_w)-sqrt(||x_i-w||^2+epsilon^2),0)^2)', gradient='dE/dx_i=-2*k*sum(overlap*(x_i-w)/sqrt(||x_i-w||^2+epsilon^2))', smoothness='C1 squared hinge with regularized distance; all catalog waters retained without outcome filtering', rationale='Same k and ratio as existing non-water repulsion; one prior setting, 0.05 ratio margin above unchanged PB 0.75 gate'),
        local=dict(max_norm_A=.10, component_bound_A=.10/np.sqrt(3), position_k_kcal_mol_A2=200., atom_selection='Cross-fragment connection endpoints and their first graph neighbors', coordinates='delta in original fragment local Cartesian frame; x=(local0+delta) R.T + center + translation', bound='Each delta component in [-0.10/sqrt(3),+0.10/sqrt(3)] A, guaranteeing norm<=0.10 A; this is a cube, not a spherical constraint', penalty='200*sum(||delta||^2); gradient=400*delta', nonlocal_atoms='Exactly original fragment SE(3), no independent displacement', rationale='One predefined 0.10 A diagnostic; not selected from test outcomes'),
        optimizer=dict(method='L-BFGS-B', maxiter=800, ftol=1e-10, gtol=1e-4, maxls=40, hydrogen_iterations=200, raw_initial_coordinates='identical for A0, A1 and B1', initial_rigid_parameters='zero', initial_local_parameters='zero', stopping='Original numerical convergence or budget; no outcome-dependent stopping'),
        variants=dict(A0_B0='Unmodified original repair function, no water/local flexibility; one shared control', A1='Only add fixed-water energy and analytic gradient', B1='Only add local port/neighbor coordinates and their position regularizer; no water energy'),
        excluded_optional_B2='Not run; no extra all-atom comparison requested as mandatory', groups=Counter(r['receptor_id'] for r in rows), versions=versions(), reporting='Strict HQ and auxiliary HQ ignoring exactly two water checks always separate. All are postprocessed saved candidates, not raw network output.')
    write_json(OUT/'input_manifest.json', rows)
    write_json(OUT/'input_fingerprints.json', fingerprints)
    write_json(OUT/'fixed_waters.json', water_rows)
    savecsv(OUT/'input_audit.csv', audit)
    write_json(OUT/'baseline_audit.json', dict(n=144, receptors=6, seeds_per_receptor=3, candidates_per_seed=8, fixed_water_atoms=79, strict_HQ=hq, auxiliary_no_water_HQ=aux, water_only=water_only, water_distance_fail=sum(r['modules']['posebusters']['tests']['minimum_distance_to_waters'] is False for r in previous), water_volume_fail=sum(r['modules']['posebusters']['tests']['volume_overlap_with_waters'] is False for r in previous)))
    write_json(OUT/'quality_config.json', load(QC))
    write_json(OUT/'posebusters_effective_config.json', pb)
    write_json(CONFIG, protocol)
    print(json.dumps(load(OUT/'baseline_audit.json'), ensure_ascii=False), flush=True)

if __name__ == '__main__':
    main()
