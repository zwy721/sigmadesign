"""CPU preparation workers share only CUDA_VISIBLE_DEVICES=0.
Precompute the last three pockets before the original four workers reach them.
The original workers resume existing result.json without rerunning any attempt.
"""
from common import *
import argparse
import run
def main():
 p=argparse.ArgumentParser();p.add_argument('--pocket',choices=['9SYL','9T4A','9T51'],required=True);p.add_argument('--shard',type=int,required=True);a=p.parse_args()
 run.torch.set_num_threads(1);run.torch.cuda.set_device(0);run.torch.use_deterministic_algorithms(True);run.torch.backends.cuda.matmul.allow_tf32=False;run.torch.backends.cudnn.allow_tf32=False;run.RDLogger.DisableLog('rdApp.*')
 plan=load(OUT/'frozen_plan.json');tag=f'extra_{a.pocket}_{a.shard}';phase=OUT/'main'
 # Model, algorithm and seeds are identical; only CPU scheduling differs.
 lib=run.VocabLibrary();den=run.load_official(run.bootstrap.CHECKPOINT,OUT/f'official_cache_{tag}','cuda:0');head=run.load_connection(run.weight_paths()['connection'],'cuda:0');client=run.Client(tag);results=[]
 try:
  context=run.context_from_pocket(a.pocket)
  for index,seed in enumerate(plan['seeds']):
   if index%4!=a.shard:continue
   for strategy in plan['strategies']:
    r=run.attempt(den,lib,head,client,context,a.pocket,seed,strategy,'main',plan);results.append(r);write_json(phase/f'attempts_{tag}.json',results);print(tag,len(results),seed,strategy,r['status'],r.get('error',''),round(r['seconds'],2),flush=True)
 finally:client.close()
 write_json(phase/f'weight_integrity_{tag}.json',{k:sha256(v)==plan['weights'][k]['sha256'] for k,v in run.weight_paths().items()})
if __name__=='__main__':main()

