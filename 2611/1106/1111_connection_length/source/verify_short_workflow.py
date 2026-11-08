"""Exercise predicted-graph optimization, evaluation, MMFF and reporting paths."""
from common import *
import time,torch
import train as trainer
from models import *
from metrics import *
from train import adapt_predicted
from benchmarks import mmff,conformers
from robust_connection_data import make_official_diffuser

def main():
 torch.set_num_threads(2);trainer.START=time.monotonic();rows=choose_views(torch.load(OUT/'data_smoke/train.pt',weights_only=False),2,2);saved=torch.load(OUT/'preflight_smoke.pt',map_location='cpu',weights_only=False);classifier=ConnectionHead(True).cuda();classifier.load_state_dict(saved['classifier']);geo=GeometryHead().cuda();geo.load_state_dict(saved['geometry']);diff=make_official_diffuser(OUT/'cache/noise','cuda:0');row=rows[0];b=pack([row],'cuda:0');t=torch.tensor([.3],device='cuda:0');T,R=noisy_frames(b,diff,t,seed=7881);bm,am=supervision_masks(row,row['connections']);case=dict(row=row,T=T.cpu(),R=R.cpu(),t=.3,graph=row['connections'],bond_mask=bm,angle_mask=am)
 # The same sample must retain the same noise across evaluation batch sizes and subsets.
 bb=pack(rows,'cuda:0');times=torch.full((len(rows),),.3,device='cuda:0');whole=evaluation_frames(bb,diff,times,7877);first=pack([rows[0]],'cuda:0');one=evaluation_frames(first,diff,times[:1],7877);n=len(rows[0]['chem']['labels']);assert torch.equal(whole[0][:n],one[0]) and torch.equal(whole[1][:n],one[1])
 import metrics as metrics_module
 metrics_module._EVAL_NOISE.clear();fresh=evaluation_frames(first,diff,times[:1],7877);assert torch.equal(fresh[0],one[0]) and torch.equal(fresh[1],one[1])
 c=adapt_predicted('preflight_C',7880,geo,classifier,[case],[case],'cuda:0');assert not c['skipped'];d=adapt_predicted('preflight_D',7880,geo,classifier,[case],[case],'cuda:0',joint=True);assert not d['skipped']
 metrics,arrays=geometry_evaluation(geo,rows,diff,7882,[.02,.9],save_coordinates=True);assert len(metrics)==len(rows)*2;div,tors=conformers(rows,arrays,'preflight',-1,'train_only');assert div
 repaired=mmff(rows[:1],diff,7883,'train_only','cuda:0');assert len(repaired)==3 and any(r.get('status')=='evaluated' for r in repaired)
 dump(OUT/'short_workflow_verification.json',dict(passed=True,evaluation_noise_batch_subset_invariance=True,evaluation_noise_cache_clear_reproducible=True,scope='Tiny train-only software smoke, not independent model performance',predicted_optimizer_steps=len(c['history'])-1,joint_optimizer_steps=len(d['history'])-1,geometry_views=len(metrics),conformer_rows=len(div),MMFF_views=repaired))
 print('SHORT WORKFLOW PASSED',flush=True)
if __name__=='__main__':main()
