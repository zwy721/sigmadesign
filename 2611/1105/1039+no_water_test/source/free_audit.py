from common import *
import torch,copy
PROJECT=ROOT.parent/'oracle_rebuild_20261031';sys.path.insert(0,str(PROJECT))
import bootstrap
from oracle_graph import rebuild
from pipeline_utils import rigid_fit
from collections import Counter

def main():
 attempts=load(OUT/'free_probe/attempts.json');rows=[]
 for r in attempts:
  folder=OUT/'free_probe/structures'/f"{r['receptor_id']}_{r['seed']}";spec=load(folder/'original_spec.json');context=copy.deepcopy(torch.load(ROOT/f"fbdd_efficiency_generalization_20261103/independent_inputs/{r['receptor_id']}/fixture.pt",weights_only=False,map_location='cpu')['context']);context['pdb_path']=str(OUT/'pockets'/f"{r['receptor_id']}_no_water.pdb");context['ligand_path']=None;g=rebuild(spec,**context);pn=int((g.frag_idx_map==-1).sum());xyz=np.load(folder/'raw_cloud.npz')['xyz'];start=0;fits=[];fitted=[];groups=[f['atom_ids'] for f in spec['fragments']]
  for ids in groups:
   n=len(ids);src=g.ref_pos[pn+np.asarray(ids)].numpy();desired=xyz[start:start+n];T,R,error=rigid_fit(src,desired);fits.append(float(error));fitted.append(src@R.T+T);start+=n
  diff=np.concatenate(fitted)-xyz;rows.append(dict(receptor_id=r['receptor_id'],seed=r['seed'],fragment_fit_RMSD_A=fits,max_fragment_fit_RMSD_A=max(fits),native_start_from_saved_decoded_RMSD_A=float(np.sqrt(np.square(diff).sum(1).mean())),native_start_max_atom_difference_A=float(np.linalg.norm(diff,axis=1).max())))
 ids=[label for r in attempts for label in r.get('labels',[])];confidence=[c for r in attempts for c in r.get('identity_confidence',[])];summary=dict(attempts=len(attempts),failed=sum(r['status']!='ok' for r in attempts),identity_class_counts=Counter(ids),unique_identity_classes=len(set(ids)),identity_top1_confidence_min=min(confidence),identity_top1_confidence_max=max(confidence),all_slots=90,site_identity_uses_receptor=False,identity_generation_mode='one-shot all-MASK at maximum noise t=1, diagnostic only; does not benchmark all possible iterative samplers',maximum_native_start_projection_RMSD_A=max(r['native_start_from_saved_decoded_RMSD_A'] for r in rows),maximum_native_start_projection_atom_displacement_A=max(r['native_start_max_atom_difference_A'] for r in rows),additional_official_pose_NFE=sum(r.get('neural_NFE',0) for r in attempts),identity_forwards=18,connection_forwards=18)
 write_json(OUT/'free_probe/capability_diagnostic.json',summary);savecsv(OUT/'free_probe/native_handoff.csv',rows);print(summary,flush=True)
if __name__=='__main__':main()
