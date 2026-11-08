"""Completion-aware reports; empty metrics remain unavailable, never zero-filled."""
from common import *
import csv,itertools
from collections import defaultdict,Counter
import numpy as np

def csvfile(name,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status']
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def average(values):
 v=[x for x in values if x is not None];return float(np.mean(v)) if v else None
def fmt(value):return 'NA' if value is None else f'{value:.4f}'
def table(rows,fields):return '| '+' | '.join(fields)+' |\n| '+' | '.join(['---']*len(fields))+' |\n'+''.join('| '+' | '.join(str(r.get(k,'NA')) for k in fields)+' |\n' for r in rows)

def main():
 status=load(OUT/'status.json') if (OUT/'status.json').exists() else dict(stage='preparing',training_complete=False);complete=status.get('training_complete',False);state='全部完成' if complete else f'尚未完成；当前阶段：{status.get("stage")}';summaries={s:load(OUT/'data'/f'{s}_summary.json') for s in ['train','val','test','scaffold_test'] if (OUT/'data'/f'{s}_summary.json').exists()};data_rows=[dict(split=s,selected_molecules=x['raw_selected_molecules'],usable_molecules=x['usable_molecules'],usable_conformers=x['usable_conformers'],failed_conformers=x['failures'],MMFF_bonds=x['mmff_bonds_available']) for s,x in summaries.items()]
 conn=[];ports=[];geo=[];diversity=[];torsion=[];search=[];repair=[]
 for p in sorted((OUT/'evaluation').glob('*/*/*/connection.json')):
  seed,group,split=p.relative_to(OUT/'evaluation').parts[:3];conn.extend(dict(seed=int(seed),group=group,split=split,**r) for r in load(p))
  ports.extend(dict(seed=int(seed),group=group,split=split,**r) for r in load(p.with_name('port_cases.json')))
 for p in sorted((OUT/'evaluation').glob('*/*/*/geometry_*.json')):
  seed,group,split=p.relative_to(OUT/'evaluation').parts[:3];condition='Oracle_true_graph' if p.stem=='geometry_oracle' else 'predicted_legal_graph'
  for r in load(p):
   record=dict(seed=int(seed),group=group,split=split,condition=condition,**{k:v for k,v in r.items() if k not in ['clean','before','after']})
   for part in ['clean','before','after']:record.update({part+'_'+k:v for k,v in r.get(part,{}).items()})
   geo.append(record)
 for p in sorted((OUT/'evaluation').glob('*/*/*/searched_graphs.json')):
  seed,group,split=p.relative_to(OUT/'evaluation').parts[:3];search.extend({**r,'seed':int(seed),'group':group,'split':split} for r in load(p))
 for p in (OUT/'evaluation').glob('*/*/*/conformer_diversity.json'):diversity.extend(load(p))
 for p in (OUT/'evaluation').glob('*/*/*/torsions.json'):torsion.extend(load(p))
 for p in (OUT/'evaluation').glob('*/MMFF_*.json'):repair.extend(load(p))
 csvfile('NOISE_LEVEL_COMPARISON.csv',conn);csvfile('PORT_PREDICTION_RESULTS.csv',ports);csvfile('GEOMETRY_METRICS.csv',geo);csvfile('SEARCH_RESULTS.csv',search);csvfile('CONFORMER_DIVERSITY.csv',diversity);csvfile('TORSION_DISTRIBUTIONS.csv',torsion);csvfile('MMFF_BASELINE.csv',repair)
 per_n=[];by_n=defaultdict(list)
 for r in ports:by_n[r['seed'],r['group'],r['split'],r['t'],r['N']].append(r)
 for (seed,group,split,t,n),rows in sorted(by_n.items()):
  tp=sum(r['typed_true_positive'] for r in rows);pred=sum(r['predicted'] for r in rows);truth=sum(r['true_bonds'] for r in rows)
  per_n.append(dict(seed=seed,group=group,split=split,t=t,N=n,views=len(rows),molecules=len({r['molecule_id'] for r in rows}),true_bonds=truth,predicted=pred,typed_true_positive=tp,precision=tp/max(1,pred),recall=tp/max(1,truth),F1=2*tp/max(1,pred+truth),macro_AP=average([r['AP'] for r in rows]),raw_exact_graph_rate=average([r['raw_exact_graph'] for r in rows]),raw_legal_complete_graph_rate=average([r['raw_legal_complete_graph'] for r in rows])))
 csvfile('PER_N_METRICS.csv',per_n)
 cr=[]
 for split,group in sorted({(r['split'],r['group']) for r in conn}):
  rows=[r for r in conn if r['split']==split and r['group']==group];cr.append(dict(split=split,group=group,seeds=len({r['seed'] for r in rows}),noise_strata=len({r['t'] for r in rows}),mean_F1=fmt(average([r['f1'] for r in rows])),mean_AP=fmt(average([r['AP'] for r in rows])),mean_raw_exact_graph_rate=fmt(average([r['raw_exact_graph_rate'] for r in rows])),mean_raw_legal_graph_rate=fmt(average([r['raw_legal_complete_graph_rate'] for r in rows]))))
 gr=[]
 for split,group,condition in sorted({(r['split'],r['group'],r['condition']) for r in geo}):
  rows=[r for r in geo if (r['split'],r['group'],r['condition'])==(split,group,condition)];available=[r for r in rows if r.get('graph_available')];gr.append(dict(split=split,group=group,condition=condition,total_views=len(rows),available_views=len(available),before_bond_MAE_A=fmt(average([r.get('before_bond_MAE_A') for r in available])),after_bond_MAE_A=fmt(average([r.get('after_bond_MAE_A') for r in available])),before_angle_MAE_deg=fmt(average([r.get('before_angle_MAE_deg') for r in available])),after_angle_MAE_deg=fmt(average([r.get('after_angle_MAE_deg') for r in available])),clash_before=fmt(average([int(r['before_severe_interfragment_clash']) for r in available])),clash_after=fmt(average([int(r['after_severe_interfragment_clash']) for r in available])),joint_success_per_all_view=fmt(sum(bool(r.get('joint_success')) for r in rows)/len(rows))))
 # Paired graph-level uncertainty, after averaging all conformers/noises/seeds per graph.
 paired=[]
 for split in ['test','scaffold_test']:
  per=defaultdict(list)
  for r in ports:
   if r['split']==split:per[r['group'],r['molecule_id']].append(r['f1'])
  for left,right in [('A2','A1'),('A2','A2_2d'),('A1','A1_2d')]:
   ids=sorted({m for g,m in per if g==left}&{m for g,m in per if g==right})
   if ids:
    delta=np.array([np.mean(per[left,m])-np.mean(per[right,m]) for m in ids]);rng=np.random.default_rng(20261106);bs=np.mean(delta[rng.integers(0,len(delta),(2000,len(delta)))],axis=1);lo,hi=np.quantile(bs,[.025,.975]);paired.append(dict(split=split,comparison=left+' minus '+right,molecules=len(ids),mean_paired_case_F1_difference=float(delta.mean()),CI95_low=float(lo),CI95_high=float(hi),positive_interval=bool(lo>0),scope='Molecule bootstrap after averaging conformers, noise strata and available seeds; interim until all 3 seeds complete'))
 dump(OUT/'paired_bootstrap.json',paired)
 intro=f'状态：**{state}**。只使用GPU 0，首轮正式训练/评估总墙钟预算24小时；不会自动扩大全量GEOM训练。报告在每个完整种子及最终退出时刷新，未完成项不补零。\n\n'
 note='所有组获得给定的真实片段身份和N；本轮是分子内部连接/几何训练，不代表自由生成或蛋白条件docking。Oracle几何组使用真实连接；预测图组只监督预测正确的连接，错误图不强制拉向真实拓扑。\n\n'
 (OUT/'GEOM_DATASET_AUDIT.md').write_text('# GEOM 数据集审计\n\n'+intro+table(data_rows,list(data_rows[0]) if data_rows else ['split'])+'\n来源为现有GEOM原始多构象文件；原SQLite仅冻结第0构象。本轮按完整canonical/isomeric分子图去重，保留原train/dev/test边界，先预留scaffold测试骨架，再排除其他集合中的这些骨架。每分子最多4构象，索引0加3个固定哈希抽样索引，不按能量/质量筛选。\n\n2768个分子均保留；构象失败及原始索引逐项保存在data/*_failures.json、*_sources.json。失败构象不参加拟合，各组使用完全相同可用数据。切分规则为冻结历史切分器；完整图重建、原子映射、键级及立体一致性核验见dataset_verification.json。MMFF94s按完整分子化学环境逐键取平衡键长，GEOM原构象键长独立保存；参考角度为对应GEOM构象。\n\n分子/骨架隔离针对本轮GEOM训练及验证；没有证明这些骨架未出现在既有PDBBind基线预训练中，不作全历史训练宇宙未见声明。\n')
 (OUT/'CONNECTION_TRAINING_RESULTS.md').write_text('# 连接训练结果\n\n'+intro+note+(table(cr,list(cr[0])) if cr else '尚无正式heldout结果。\n')+'\nF1按正确端口和键级计真阳性；AP按是否成键计。统计所有合法端口候选，不按带噪距离删除正例；困难负样本包含同一分子内全部其他端口对。原始图准确率和原始合法完整率不包含搜索修补。分层搜索固定128骨架beam、24骨架、64端口beam，搜索结果另列SEARCH_RESULTS.csv；每集合固定32分子×第一个可用构象×3噪声，仅作为有明确分母的有限搜索对照。\n\n表中为噪声层与已完成种子的均值，不是把不同层混成一个分母；逐层记录、分子N与端口rank见CSV。检查点每5epoch按验证集256分子和3噪声选择，共40epoch；测试集不参与调参或选择。\n')
 (OUT/'GEOMETRY_RECOVERY_RESULTS.md').write_text('# 可微局部几何恢复\n\n'+intro+note+(table(gr,list(gr[0])) if gr else '尚无正式heldout几何结果。\n')+'\n几何头输出片段局部坐标系下的平移和轴角修正，通过矩阵指数构成刚体旋转，四轮更新；不改变片段内部坐标。每轮分量平移限12Å、旋转限0.6rad。类型/原子、全片段局部原子编码、相对位姿、时间、连接键级和立体状态参与网络输入。\n\nA3只用MMFF94s键长MSE。A4增加2倍键角cosine MSE和0.5倍相对连接轴监督，非旋转双键附加圆周扭转监督；可旋转单键不加目标扭转或全局参考构象坐标MSE。没有将离散搜索或MMFF当作可微训练路径。\n\n几何表按构象—噪声样本平均；预测图不可用视图保留在总分母。联合成功＝连接图精确且合法、全部可评价跨片段键误差≤0.15Å、全部跨片段角误差≤15°、无严重跨片段非键碰撞（排除1–2/1–3邻居，距离<0.6×范德华半径和）。这是本轮内部指标，不是历史HQ_no_water。\n\nMMFF_BASELINE.csv单列从同一噪声起点出发的独立Oracle拓扑、全原子MMFF94s修复，最多200迭代，保留未收敛和失败；其效果不计作神经网络效果。CONFORMER_DIVERSITY.csv报告同分子不同构象RMSD、最多4参考构象的1Å覆盖和参考多样性；未做原子对称校正。TORSION_DISTRIBUTIONS.csv保留跨连接扭转分布，单键参考偏差不是训练损失。\n')
 ablation='|组|连接分类|局部3D特征|几何监督|\n|---|---|---|---|\n|A0|冻结原头|原特征|无|\n|A1|微调当前TimeConnectionPredictor|原特征|无|\n|A2|同A1|新增全原子邻域、相对T/R、位阻|无|\n|A3|从同seed的A2开始|同A2|仅键长|\n|A4|从同seed的A2开始|同A2|键长+键角+合理相对位姿|\n|A1_2d/A2_2d|同组参数与训练预算|动态及模板坐标通道归零|无|\n\n'
 (OUT/'ABLATION_RESULTS.md').write_text('# 消融与梯度归属\n\n'+intro+ablation+'A/B/C/D阶段及选择记录见runs/<seed>/<group>/stage_decision.json。B只训练几何头，C冻结分类器并在预测合法图中仅监督正确连接，D仅在A/B/C验证门槛均通过后尝试联合微调。D的CE梯度可通过可微坐标特征进入几何头；几何损失进入几何头；离散选图不接收梯度。\n\n独立预检见preflight.json和prelaunch_verification.json，包含真实非零梯度、刚体内部距离、全局SE(3)等变性、清洁参考不进入运行时、二维消融坐标不变性及官方接口。\n\n配对分子bootstrap结果见paired_bootstrap.json，区间跨0时不作显著改进结论；三随机种子全部完成前只作中间观察。\n')
 integration=load(OUT/'native_integration.json') if (OUT/'native_integration.json').exists() else None
 (OUT/'SIGMADOCK_INTEGRATION_AUDIT.md').write_text('# SigmaDock 接口审计\n\n'+intro+'现有TimeConnectionPredictor仅分类，无几何修正输出；原有端口距离/方向特征不等于全片段原子环境。新增模块独立保存，旧权重和接口文件保持只读。\n\n调用入口：source/inference.py 的 Runtime.predict(types,T,R,t,candidate_ports=None,local_atoms=None,legal_connections=None)。类型可来自身份模型；未给图时先分类和合法搜索，给定合法图时直接几何修正。运行时不读取GEOM原构象、真实配体或参考坐标；真实图只用于明确标注的Oracle训练/评价。\n\nNativeState.slot_frames提取当前状态→保留实际片段内部坐标→神经SE(3)修正→同一合法图的NativeState坐标转移→官方NativeState.step继续去噪。仅做接口兼容，不训练或评价docking效能。\n\n接口运行核验：'+('通过。完整记录见native_integration.json；动态边与几何重新计算已审计。\n' if integration and integration.get('passed') else '尚未完成，见prelaunch.log。\n'))
 decision='尚不能给出最终训练结论。等待三个种子的独立测试和scaffold测试完成，禁止据小样本过拟合结果启动全量训练。'
 if complete:
  positive=[p for p in paired if p['comparison']=='A2 minus A1' and p['positive_interval']];decision=f'本轮小规模训练已结束。局部3D特征的配对F1正区间覆盖{len(positive)}/2种heldout条件；详细效应和区间见paired_bootstrap.json。下一步应结合原始端口rank、键级错误、Oracle/预测图几何差距和构象覆盖决定，不能以搜索合法率或MMFF结果代替模型能力。全量GEOM训练仍需另行确认。'
 (OUT/'NEXT_STEP_DECISION.md').write_text('# 下一步决策\n\n'+intro+decision+'\n\n本轮不重新训练片段身份模型，不修改官方SigmaDock，也没有蛋白结合位姿监督。连接错误归因需比较端口rank/键级准确率及2D/3D消融；几何能力需比较直接输出的键长、角度和碰撞；预测图错误保留为连接失败，不将其当作几何头失败掩盖。\n')
 dump(OUT/'summary.json',dict(status=status,datasets=summaries,connection=cr,geometry=gr,paired_bootstrap=paired,formal_connection_rows=len(conn),formal_geometry_rows=len(geo),formal_search_views=len(search)))
 # File hashes capture completed checkpoints and inputs, avoiding changing logs/status.
 if complete:
  files={str(p.relative_to(OUT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in OUT.rglob('*') if p.is_file() and p.suffix in ['.pt','.json','.csv','.py','.md'] and p.name!='delivery_manifest.json' and 'cache' not in p.parts};dump(OUT/'delivery_manifest.json',dict(complete=True,files=files))
 if conn:
  import matplotlib
  matplotlib.use('Agg');import matplotlib.pyplot as plt
  fig,ax=plt.subplots(1,2,figsize=(11,4))
  for split,a in zip(['test','scaffold_test'],ax):
   for group in ['A0','A1','A2','A3','A4','A1_2d','A2_2d']:
    rows=[r for r in conn if r['split']==split and r['group']==group]
    if rows:
     ts=sorted({r['t'] for r in rows});a.plot(ts,[average([r['f1'] for r in rows if r['t']==t]) for t in ts],marker='.',label=group)
   a.set(title=split,xlabel='Noise t',ylabel='Raw typed-port F1',ylim=(0,1));a.legend(fontsize=7)
  fig.tight_layout();fig.savefig(OUT/'NOISE_PERFORMANCE.png',dpi=160);plt.close(fig)
 print('REPORTS REFRESHED',state,flush=True)
if __name__=='__main__':main()
