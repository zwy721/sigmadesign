from common import *
import copy,torch
from collections import Counter
from rdkit import Chem,RDLogger
from denovo_quality.io import load_molecule
from denovo_quality.protein import read_protein
from geometry_constraints_20261103.stereo import derive
RDLogger.DisableLog('rdApp.*')

def main():
 rows=load(OLD/'input_manifest.json');assert len(rows)==144 and len({key(r) for r in rows})==144
 cases=load(ROOT/'fbdd_efficiency_generalization_20261103/test_cases.json')
 audit=[];proteins={};fingerprints={}
 def track(p):fingerprints[str(Path(p).resolve())]=sha256(p)
 for pid in sorted({r['receptor_id'] for r in rows}):
  row=next(r for r in rows if r['receptor_id']==pid);path=Path(row['protein_path']);folder=path.parent
  lines=path.read_text().splitlines(keepends=True);new=[];removed=[]
  for line in lines:
   if line.startswith(('ATOM  ','HETATM','ANISOU')) and line[17:20].strip() in WATER:removed.append(line)
   else:new.append(line)
  dest=OUT/'pockets'/f'{pid}_no_water.pdb';dest.parent.mkdir(exist_ok=True);dest.write_text(''.join(new));proteins[pid]=str(dest)
  aa=read_protein(str(path))['atoms'];bb=read_protein(str(dest))['atoms'];assert [a for a in aa if a['kind']!='water']==bb
  f=torch.load(folder/'fixture.pt',map_location='cpu',weights_only=False);pg=torch.load(folder/'protein_graph.pt',map_location='cpu',weights_only=False);pocket=f['context']['pocket']
  residues=[a.GetPDBResidueInfo().GetResidueName().strip() for a in pocket.GetAtoms()];assert not set(residues)&WATER
  coords=pocket.GetConformer().GetPositions();assert np.max(np.abs(coords-pg['pos'].numpy()))<1e-4
  nearest=np.linalg.norm(coords[:,None]-np.array([a['xyz'] for a in aa if a['kind']=='protein'])[None,:],axis=-1).min(1);assert nearest.max()<.002
  graph=torch.load(folder/'official_true_graph.pt',map_location='cpu',weights_only=False);pn=int((graph.frag_idx_map==-1).sum());assert pn==len(coords) and torch.equal(graph.ref_pos[:pn],pg['pos'])
  audit.append(dict(receptor_id=pid,original_path=str(path),water_free_path=str(dest),original_heavy_atoms=len(aa),water_heavy_atoms=sum(a['kind']=='water' for a in aa),nonwater_cofactors=Counter(a['kind'] for a in aa if a['kind'] not in ['water','protein']),retained_nonwater_heavy_atoms=len(bb),model_pocket_heavy_atoms=len(coords),model_water_atoms=0,model_nonwater_cofactor_atoms=sum(r not in {a['residue'].split(':')[-1].rstrip('0123456789') for a in aa if a['kind']=='protein'} for r in []),serialized_model_residues=sorted(set(residues)),protein_graph_equal_to_fixture=True,official_graph_protein_equal=True,original_and_model_cofactor_policy='Original PDB and new no-water QC retain cofactors; historical official model removes all heterogens. No claim of cofactor-conditioned inference.'))
  for p in [path,dest,folder/'fixture.pt',folder/'protein_graph.pt',folder/'official_true_graph.pt']:track(p)
 track(QC)
 from denovo_quality.posebusters_adapter import effective_config,test_names
 full=copy.deepcopy(effective_config('dock',1));install_pb();from denovo_quality.posebusters_adapter import effective_config
 nw=effective_config('dock',1);removed_modules=[m for m in full['modules'] if set(test_names({'modules':[m]}))&WATER_TESTS];assert nw==dict(full,modules=[m for m in full['modules'] if m not in removed_modules])
 write_json(OUT/'posebusters_original.json',full);write_json(OUT/'posebusters_no_water.json',nw)
 cfg=load_config(QC);write_json(OUT/'quality_config.json',cfg)
 raw_stereo={key(r):r for r in load(ROOT/'fbdd_efficiency_generalization_20261103/independent_raw_stereo/manifest.json')}
 manifest=[];cloud=[]
 for r in rows:
  k=key(r);nm=name(r);pid,seed,rank=k;source=Path(r['sdf_path']).parent
  # Original raw manifest uses local copies: resolve actual saved neural trajectory by canonical manifest.
  native=next(x for x in load(ROOT/'fbdd_efficiency_generalization_20261103/independent/native/raw_manifest.json') if key(x)==k);source=Path(native['sdf_path']).parent
  short=np.load(source/'short.npz');initial=short['trajectory'][0];groups=load(source/'build.json')['groups'];flat=[a for g in groups for a in g]
  prefix=ROOT/'fbdd_efficiency_generalization_20261103/independent/native/prefix'/f'{pid}_{seed}'/'t0.6.npz'
  pre=np.load(prefix)['xyz'];assert np.max(np.abs(initial-pre))<1e-4
  mol,info,fids=load_molecule(native);assert info['sanitize_success'];x=np.empty_like(initial);x[flat]=pre
  for a,xyz in enumerate(x):mol.GetConformer().SetAtomPosition(a,xyz.tolist())
  folder=OUT/'structures/decoded_initial'/nm;sdf(folder/'raw.sdf',mol)
  # Same SDF rounding and geometry-derived stereo for every complete graph stage.
  mol,_,_=load_molecule(dict(native,sdf_path=str(folder/'raw.sdf')));derived,spec,st=derive(mol,load(native['chemistry_spec_path']));sdf(folder/'derived.sdf',derived);write_json(folder/'spec.json',spec);write_json(folder/'stereo_decisions.json',st)
  stage_rows=[('decoded_initial',dict(native,sdf_path=str(folder/'derived.sdf'),chemistry_spec_path=str(folder/'spec.json'))),('native_final',raw_stereo[k])]
  for stage,variant in [('rigid','A0_B0'),('local_flex','B1')]:
   d=OLD/'structures'/variant/nm;stage_rows.append((stage,dict(native,sdf_path=str(d/'derived.sdf'),chemistry_spec_path=str(d/'derived_spec.json'))))
  for stage,base in stage_rows:
   q=dict(base,molecule_id=f'{stage}__{nm}',group=stage,protein_path=proteins[pid]);q.pop('reference_path',None);manifest.append(q)
   track(q['sdf_path']);track(q['chemistry_spec_path'])
  cloud.append(dict(receptor_id=pid,seed=seed,candidate_rank=rank,prefix_path=str(prefix),spec_path=native['chemistry_spec_path'],groups=groups,protein_path=proteins[pid],coordinates='unconnected prefix heavy atoms; flattened fragment order',shared_raw_cloud_id=f'{pid}_{seed}',raw_to_decoded_coordinate_max_A=float(np.max(np.abs(initial-pre)))))
  for p in [prefix,source/'short.npz',source/'final.npz',source/'build.json',native['chemistry_spec_path']]:track(p)
 for p in Path(__file__).parent.glob('*.py'):track(p)
 for p in [ROOT/'fbdd_efficiency_generalization_20261103/prepare_fixtures.py',ROOT/'fbdd_efficiency_generalization_20261103/legacy_prefix.py',ROOT/'fbdd_efficiency_generalization_20261103/native_generate.py',ROOT.parent/'oracle_rebuild_20261031/joint_sampler.py',ROOT.parent/'oracle_rebuild_20261031/identity_server.py']:track(p)
 weights=load(ROOT/'multi_candidates_20261103/plan.json')['weight_sha256'];assert all(sha256(p)==h for p,h in weights.items())
 write_json(OUT/'pocket_audit.json',audit);savecsv(OUT/'pocket_audit.csv',audit);write_json(OUT/'manifest.json',manifest);write_json(OUT/'raw_cloud_manifest.json',cloud);write_json(OUT/'input_fingerprints.json',fingerprints)
 protocol=dict(name='SigmaDock-FBDD no-crystal-water v1 20261105',historical_candidate_reuse=True,reason='Actual serialized model pocket and shared protein graph are water free for all six receptors.',primary='HQ_no_water',PB_removed_exactly=sorted(WATER_TESTS),PB_remaining_config_unchanged=True,custom_thresholds_unchanged=True,water_repulsion=False,original_denominator=144,independent_receptors=6,seeds=[1801,1802,1803],K=8,unique_predecode_clouds=18,stages=['raw_predecode_t0.6','decoded_initial_t0.6','native_final_t0.02','rigid','local_flex'],raw_graph_metrics='Not evaluable before chemical graph; no fabricated HQ.',conditioning='Historical pool: true fragment identities and true fragment count, ligand-defined pocket region. Not free de novo.',decoder='One shared raw cloud per receptor/seed; eight graph candidates. Decoder changes chemical graph, not raw coordinates. Native continuation is a separate neural stage.',repair='Both start native_final independently. Existing MMFF94s rigid protocol and B1 0.1A local bound; no water energy. Same geometry-derived stereo.',additional_neural_NFE=0,weights=weights,source_config_sha256=sha256(QC),free_generation_audit='Separate capability audit; any independent free-generation cohort reported separately.',historical_results_unchanged=True)
 write_json(OUT/'protocol.json',protocol);print('prepared',len(manifest),'QC records; actual model water atoms',sum(a['model_water_atoms'] for a in audit),flush=True)
if __name__=='__main__':main()
