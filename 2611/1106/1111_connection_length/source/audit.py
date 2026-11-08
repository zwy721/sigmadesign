"""Read-only source/data audit. This entry point performs no training."""
from common import *
import sqlite3,zlib,pickle,inspect
from collections import Counter
import numpy as np
import torch
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
from cutoff_data.codec import get_codec

def main():
 torch.set_num_threads(2)
 db=DATA/'geom_database'
 con=sqlite3.connect(f'file:{db}/samples.sqlite?mode=ro&immutable=1',uri=True)
 counts=con.execute('SELECT split,count(*),count(distinct graph_id),sum(identity_connection_eligible) FROM samples GROUP BY split').fetchall()
 crossing=con.execute('SELECT count(*) FROM (SELECT graph_id FROM samples GROUP BY graph_id HAVING count(distinct split)>1)').fetchone()[0]
 sample=json.loads(zlib.decompress(con.execute('SELECT sample_json FROM samples ORDER BY sample_id LIMIT 1').fetchone()[0]))
 assert sha(sample['pickle_path'])==sample['source_sha256']
 with open(sample['pickle_path'],'rb') as f:raw=pickle.load(f)
 codec=get_codec();conf=raw['conformers'][0];mol=conf['rd_mol']
 canonical,original,xyz,stereo=codec.prepare_molecule(mol)
 assert canonical==sample['canonical_smiles']
 parent=Chem.MolFromSmiles(canonical);cf=Chem.Conformer(parent.GetNumAtoms())
 for i,p in enumerate(xyz):cf.SetAtomPosition(i,p)
 parent.AddConformer(cf);parent=Chem.AddHs(parent,addCoords=True)
 props=AllChem.MMFFGetMoleculeProperties(parent,mmffVariant='MMFF94s')
 bond=next(iter(parent.GetBonds()));a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
 params=props.GetMMFFBondStretchParams(parent,a,b) if props else None
 files=[BASELINE,bootstrap.CHECKPOINT,ASSETS/'mapping.json',ASSETS/'templates.json',ASSETS/'properties.json',db/'database_manifest.json',BUNDLE/'src/cutoff_data/codec.py']
 files += [PROJECT/f for f in ['trajectory_connections.py','robust_connection_model.py','robust_connection_data.py','train_robust_connections.py','hierarchical_topology.py','vocab_chemistry.py','joint_sampler.py']]
 files += [Path(codec.__file__)]
 files += [Path(load(ROOT/'2611/1105/1345_diversity/frozen_plan.json')['weights']['identity']['path'])]
 baseline=torch.load(BASELINE,map_location='cpu',weights_only=False)
 report=dict(training_started=False,rdkit_version=rdBase.rdkitVersion,torch_version=torch.__version__,source_database=str(db),source_counts=[dict(split=s,records=n,canonical_graphs=g,eligible=e) for s,n,g,e in counts],canonical_graph_cross_split_count=crossing,baseline_checkpoint=str(BASELINE),baseline_epoch=baseline.get('epoch'),baseline_architecture=baseline.get('architecture'),source_fingerprints={str(p):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files},raw_multiconformer_probe=dict(sample_id=sample['sample_id'],source_sha256=sample['source_sha256'],conformers_available=len(raw['conformers']),first_conformer_identity_exact=True,canonical_to_source_atom=original,stereo=stereo),MMFF94s_probe=dict(parameterized=props is not None,bond_stretch_parameters=params,parameter_order='bondType, forceConstant, equilibriumLength_A'),baseline_audit=dict(full_atom_dynamic_geometry=False,features='Port endpoint chemistry, candidate port distance/RBF, port direction cosines, endpoint local coordinates, identities, global chemistry bag and time. No whole-fragment neighborhood/collision encoder.',geometry_correction_output=False,classification_outputs=['no_bond','single','double'],stereo_outputs=['unspecified','E','Z'],previous_training_source='PDBBind record cache, not GEOM',search_is_differentiable_training_target=False),native_interface='NativeState.slot_frames -> correction(T,R) -> rigid coordinate reconstruction/transfer -> NativeState.step; official denoiser frozen')
 dump(OUT/'audit_initial.json',report)
 (OUT/'GEOM_DATASET_AUDIT.md').write_text(f'''# GEOM 数据与现有模型初审\n\n状态：只读审计完成，尚未开始训练。\n\n现有数据库包含240,686条来源记录、240,495种canonical分子图；同一分子图跨原始split为{crossing}。已有词表3684类，仅支持单/双跨片段键。原始多构象pickle可读取且SHA与来源清单一致；探针分子有{len(raw['conformers'])}个构象，原数据库只冻结index 0。正式队列将按canonical完整分子去重并隔离全部构象；另设与训练/验证骨架隔离的测试集。未映射、立体或化学检查失败的分子/构象均记录，不静默替换标签。\n\n当前连接头有端口距离、方向余弦、端口局部化学、类型、全局片段组成和噪声时间，没有全片段原子邻域或动态局部碰撞编码，也不输出可微位姿修正。原带噪训练入口读取PDBBind缓存，本轮建立独立GEOM数据入口。\n\nMMFF94s完整分子键参数读取探针通过：{params}，分别是键类型、力常数、平衡键长（Å）。正式数据将按完整化学环境逐键取参数；不可参数化项单列。\n\n原始文件及权重指纹见audit_initial.json。现有连接头、身份模型与官方SigmaDock权重均作为只读依赖。\n''')
 print(json.dumps({k:v for k,v in report.items() if k not in ['source_fingerprints','raw_multiconformer_probe']},ensure_ascii=False),flush=True)
if __name__=='__main__':main()
