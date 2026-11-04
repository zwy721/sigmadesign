"""CPU only: A0/A1 paired replay, original QC, then B1 with original QC."""
import argparse
import copy
import os
import shutil
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from study_common import *
from rdkit import Chem, RDLogger
from denovo_quality.io import load_molecule
from denovo_quality.common import object_hash, versions
from denovo_quality.pipeline import task as quality_task, COMPUTE_MODULES
from geometry_constraints_20261103.repair import repair as original_repair
from geometry_constraints_20261103.stereo import derive
from optimizer import repair as experimental_repair
RDLogger.DisableLog('rdApp.*')


def process(args):
    row, variant = args
    folder = OUT/'structures'/variant/name(row)
    marker = folder/'result.json'
    if marker.exists():
        return load(marker)
    folder.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    cpu_start = time.process_time()
    result = dict(variant=variant, receptor_id=key(row)[0], seed=key(row)[1], rank=key(row)[2], status='failed')
    try:
        mol, info, fids = load_molecule(row)
        assert info['sanitize_success'] and not info.get('declared_chemistry_error')
        cfg = load(CONFIG)
        spec = load(row['chemistry_spec_path'])
        original = OUT/'structures/original'/name(row)
        original.mkdir(parents=True, exist_ok=True)
        # One shared byte-exact copy for all arms, with original atom numbering.
        for src, dest in [(row['sdf_path'], original/'raw.sdf'), (row['chemistry_spec_path'], original/'spec.json')]:
            if not dest.exists():shutil.copyfile(src,dest)
            assert sha256(src) == sha256(dest)
        if variant == 'A0_B0':
            q, meta = original_repair(mol, fids, row['protein_path'])
        elif variant == 'A1':
            w = dict(cfg['water'], atoms=waters_for(row))
            q, meta = experimental_repair(mol, fids, row['protein_path'], water=w)
        elif variant == 'B1':
            q, meta = experimental_repair(mol, fids, row['protein_path'], local=cfg['local'])
        else:
            raise ValueError(variant)
        repair_cpu = time.process_time()-cpu_start
        # Match historical SDF round-trip before applying the exact same stereo function.
        save_sdf(folder/'repaired.sdf', q)
        q, info, _ = load_molecule(dict(row, sdf_path=str(folder/'repaired.sdf')))
        assert info['sanitize_success']
        stereo, ds, decisions = derive(q, spec)
        save_sdf(folder/'derived.sdf', stereo)
        write_json(folder/'derived_spec.json', ds)
        write_json(folder/'stereo_decisions.json', decisions)
        write_json(folder/'repair.json', meta)
        if variant == 'A0_B0':
            old = OLD/'cases'/f"{row['receptor_id']}_{row['seed']}"/f"candidate_{int(row['candidate_rank']):02d}"
            original_derived = Chem.SDMolSupplier(str(old/'derived.sdf'), removeHs=False)[0]
            new_derived = Chem.SDMolSupplier(str(folder/'derived.sdf'), removeHs=False)[0]
            error = float(np.max(np.abs(original_derived.GetConformer().GetPositions()-new_derived.GetConformer().GetPositions())))
            # Keep a record even if platform/version differences arise; strict check at final verification.
            result.update(baseline_serialized_coordinate_max_difference_A=error, baseline_stereo_spec_equal=ds == load(old/'derived_spec.json'))
        result.update(status='ok', repair_cpu_seconds=repair_cpu, repair_wall_seconds=meta['seconds'], postprocess_cpu_seconds=time.process_time()-cpu_start, postprocess_wall_seconds=time.perf_counter()-started, raw_sha256=sha256(row['sdf_path']), spec_sha256=sha256(row['chemistry_spec_path']), derived_sdf=str(folder/'derived.sdf'), derived_spec=str(folder/'derived_spec.json'))
    except Exception as e:
        result.update(error=f'{type(e).__name__}: {e}', traceback=traceback.format_exc(), postprocess_cpu_seconds=time.process_time()-cpu_start, postprocess_wall_seconds=time.perf_counter()-started)
    write_json(marker, result)
    return result


def evaluate(args):
    arguments = args
    t = time.process_time()
    result = quality_task(arguments)
    result['quality_cpu_seconds'] = time.process_time()-t
    # Separate timing file from immutable package cache; main result JSONL preserves CPU measurements.
    return result


def run(variants, workers):
    cfg = load(CONFIG)
    assert sha256(QC) == cfg['quality_config_sha256']
    fingerprints = load(OUT/'input_fingerprints.json')
    assert all(sha256(p) == s for p,s in fingerprints.items())
    rows = load(OUT/'input_manifest.json')
    start = time.perf_counter()
    jobs = [(row, variant) for row in rows for variant in variants]
    results = []
    with ProcessPoolExecutor(workers) as pool:
        futures = [pool.submit(process, job) for job in jobs]
        for future in as_completed(futures):
            r = future.result()
            results.append(r)
            if len(results)%24 == 0 or r['status'] != 'ok':
                print('repair',variants,len(results),'/',len(jobs),r['status'],flush=True)
    results.sort(key=lambda r: (r['variant'],key(r)))
    label = 'A' if variants == ['A0_B0','A1'] else 'B'
    write_json(OUT/f'processing_{label}.json', results)
    manifest = []
    failures = []
    for r in results:
        row = next(x for x in rows if key(x) == key(r))
        if r['status'] != 'ok':
            failures.append(r)
            continue
        q = copy.deepcopy(row)
        q.update(molecule_id=f"{r['variant']}__{name(row)}", group=r['variant'], sdf_path=r['derived_sdf'], chemistry_spec_path=r['derived_spec'])
        q.pop('reference_path', None)
        manifest.append(q)
    write_json(OUT/f'quality_manifest_{label}.json', manifest)
    write_json(OUT/f'processing_failures_{label}.json', failures)
    quality_out = OUT/f'quality_{label}'
    quality_out.mkdir(exist_ok=True)
    config = load_config(QC)
    package = ROOT/'denovo_quality'
    code = {p.name:sha256(p) for p in package.glob('*.py')}
    compute = object_hash({n:code[n] for n in COMPUTE_MODULES})
    deps = versions()
    write_json(quality_out/'protocol.json',dict(config=config,config_file_sha256=sha256(QC),source_sha256=code,computation_sha256=compute,versions=deps,manifest_sha256=sha256(OUT/f'quality_manifest_{label}.json'),input_count=len(manifest),reference_conformer_supplied=False,PB_energy_ensemble_conformers=50,additional_neural_NFE=0))
    write_json(quality_out/'config_used.json',config)
    from denovo_quality.posebusters_adapter import effective_config
    write_json(quality_out/'posebusters_effective_config.json',effective_config('dock',1))
    metrics = []
    qstart = time.perf_counter()
    with ProcessPoolExecutor(workers) as pool:
        futures = [pool.submit(evaluate,(r,config,str(quality_out),compute,deps,False)) for r in manifest]
        for future in as_completed(futures):
            r = future.result()
            metrics.append(r)
            if len(metrics)%24 == 0:
                print('QC',label,len(metrics),'/',len(manifest),flush=True)
    idx = {r['molecule_id']:r for r in metrics}
    metrics = [idx[r['molecule_id']] for r in manifest]
    with (quality_out/'molecule_metrics.jsonl').open('w') as f:
        for r in metrics:f.write(json.dumps(clean(r),ensure_ascii=False,allow_nan=False)+'\n')
    write_json(OUT/f'timing_{label}.json',dict(wall_seconds=time.perf_counter()-start,quality_wall_seconds=time.perf_counter()-qstart,repair_cpu_seconds=sum(r.get('repair_cpu_seconds',0) for r in results),postprocess_cpu_seconds=sum(r['postprocess_cpu_seconds'] for r in results),quality_cpu_seconds=sum(r['quality_cpu_seconds'] for r in metrics),workers=workers,additional_neural_NFE=0, note='Per-process CPU sums; wall time separately measured. Cache CPU is measured on cache retrieval if a run is resumed.'))
    print(label,'complete; strict counts', {v:sum(r['checks'].get('high_quality') is True for r in metrics if r['metadata']['group']==v) for v in variants},flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--experiment', choices=['A','B'],required=True)
    p.add_argument('--workers',type=int,default=8)
    a = p.parse_args()
    run(['A0_B0','A1'] if a.experiment=='A' else ['B1'],a.workers)
