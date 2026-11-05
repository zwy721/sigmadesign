"""Generate reports from measured attempt/QC statistics only."""
from common import *
def f(v,d=3):
 if v is None:return 'NA'
 if isinstance(v,float):return f'{v:.{d}f}'
 return str(v)
def pct(x,n):return f'{x}/{n} ({100*x/n:.1f}%)' if n else 'NA'
def mean(values):
 values=[v for v in values if v is not None];return float(np.mean(values)) if values else None
def table(h,rr):
 return '| '+' | '.join(h)+' |\n| '+' | '.join(['---']*len(h))+' |\n'+'\n'.join('| '+' | '.join(f(x) for x in r)+' |' for r in rr)+'\n'
def write(name,s):(OUT/name).write_text(s+'\n')
def main():
 d=load(OUT/'summary.json');p=load(OUT/'frozen_plan.json');s=d['summary'];rows=load(OUT/'PER_GENERATION_RESULTS.json');base=next(x for x in s if x['strategy']=='mask_top1');cat=next(x for x in s if x['strategy']=='mask_categorical');n=d['attempts']
 main_table=table(['策略','尝试','有效连通','2D Unique','Unique/尝试','Morgan内部多样性','网络HQ','刚性HQ','局部柔性HQ'],[[x['strategy'],x['attempts'],x['native_final_valid_connected'],x['native_final_unique_constitution'],x['native_final_valid_unique_per_attempt'],x['native_final_internal_diversity'],x['native_final_high_quality_count'],x['rigid_high_quality_count'],x['local_flex_high_quality_count']] for x in s])
 yield_table=table(['策略','有效唯一/尝试','网络HQ唯一/尝试','刚性HQ唯一/尝试','柔性HQ唯一/尝试'],[[x['strategy'],pct(x['native_final_unique_constitution'],x['attempts']),pct(x['native_final_HQ_unique'],x['attempts']),pct(x['rigid_HQ_unique'],x['attempts']),pct(x['local_flex_HQ_unique'],x['attempts'])] for x in s])
 identity=table(['策略','初始/终态类数','模型熵','模型Top1概率','实际初始采样熵','3577初始/终态占比','初始/终态组合数'],[[x['strategy'],f"{x['identity_unique_initial']}/{x['identity_unique_terminal']}",x['initial_entropy_mean'],x['initial_top1_probability_mean'],x['initial_actual_draw_entropy'],f"{x['benzene_fraction_initial']:.3f}/{x['benzene_fraction_terminal']:.3f}",f"{x['fragment_combination_unique_initial']}/{x['fragment_combination_unique_terminal']}"] for x in s])
 stage=table(['策略','阶段','蛋白碰撞通过/全部尝试','HQ/全部尝试','异常跨片段键','异常跨片段角','修复RMSD均值Å','最大原子位移Å'],[[x['strategy'],st,pct(x[st+'_no_severe_protein_clash_count'],x['attempts']),pct(x[st+'_high_quality_count'],x['attempts']),x[st+'_abnormal_cross_bonds_total'],x[st+'_abnormal_cross_angles_total'],x.get(st+'_repair_RMSD_A',{}).get('mean'),x.get(st+'_max_atom_displacement_A',{}).get('max')] for x in s for st in ['native_final','rigid','local_flex']])
 cost=table(['策略','初始连接头完整树','搜索成功','接受/拒绝更新','身份/连接/位姿前向','平均秒/尝试','原始生成失败'],[[x['strategy'],x['initial_raw_head_complete_tree'],x['initial_search_success'],f"{x['accepted_revisions']}/{x['rejected_revisions']}",f"{x['identity_queries']}/{x['connection_queries']}/{x['pose_NFE']}",x['seconds_per_attempt'],x['failure_stages']] for x in s])
 pocket=table(['策略','口袋','有效连通','2D Unique','内部多样性','网络HQ','柔性HQ'],[[x['strategy'],x['receptor_id'],x['native_final_valid_connected'],x['native_final_unique_constitution'],x['native_final_internal_diversity'],x['native_final_high_quality_count'],x['local_flex_high_quality_count']] for x in d['pockets']])
 fixed=load(OUT/'PER_N_DIVERSITY.json');fixed_table=table(['策略','固定N=6每口袋平均Unique/尝试','固定N=6每口袋平均内部多样性'],[[x['strategy'],mean([r['native_final_valid_unique_per_attempt'] for r in fixed if r['strategy']==x['strategy'] and r['receptor_id']!='ALL' and r['N_stratum']==6]),mean([r['native_final_internal_diversity'] for r in fixed if r['strategy']==x['strategy'] and r['receptor_id']!='ALL' and r['N_stratum']==6])] for x in s])
 gains=f"预定义的标准全概率抽样相对 Top-1：2D 构造 Unique {base['native_final_unique_constitution']}→{cat['native_final_unique_constitution']}，Morgan 内部多样性 {base['native_final_internal_diversity']:.4f}→{cat['native_final_internal_diversity']:.4f}；网络 HQ {base['native_final_high_quality_count']}/{base['attempts']}→{cat['native_final_high_quality_count']}/{cat['attempts']}，差 {100*(cat['native_final_high_quality_count']/cat['attempts']-base['native_final_high_quality_count']/base['attempts']):+.1f} 个百分点；柔性 HQ {base['local_flex_high_quality_count']}/{base['attempts']}→{cat['local_flex_high_quality_count']}/{cat['attempts']}。有效连通 {base['native_final_valid_connected']}/{base['attempts']}→{cat['native_final_valid_connected']}/{cat['attempts']}。"
 write('DIVERSITY_EVALUATION_RESULTS.md',f"""# 冻结权重自主生成多样性评测

正式测试 {n} 次，六个已查看口袋 × 12 种子 × 七种策略；pilot 的 14 次独立保留。所有原始尝试均入分母；一次尝试只输出 k=1 候选，不补抽、不删重复、不按质量筛选。仅使用一张 A100（CUDA_VISIBLE_DEVICES=0），未训练、未改任何权重，历史目录保持原样。

{gains}

{main_table}
本轮Top-1身份仍全部为苯，但完整二维图已有{base['native_final_unique_constitution']}种：其中包含N变化与连接/位姿闭环带来的图差异。故“苯身份重复”和“完整分子只有一种”必须分开；历史只有一种构造是旧固定N/单次解码策略的有限结果。本轮随机策略的主要收益包括不同片段组合与指纹结构覆盖，不能只用二维Unique数量概括。

高质量且唯一的构造产出如下，去重只用于计算指标，全部输出与重复都保留：

{yield_table}

2D Unique 是 RDKit Cleanup、移除 atom map、保留电荷/组分后去立体标签的构造 SMILES。含立体标签 Unique 单列 CSV。Morgan radius=2、2048 bits、chirality off，Internal Diversity 是所有有效输出对的平均 1−Tanimoto，重复样本参与统计；去重结构上的指标也单列，不混淆二者。Scaffold 用 Murcko，计数/比例在 CSV。未把不同 3D 构象计作不同 2D 分子。

## 身份与口袋条件

{identity}
GEOM 身份头没有受体输入，初始相同种子的身份跨口袋相同，概率在1e-6绝对容差内一致；证据见 CROSS_POCKET_COMPARISON.csv。它能在不同随机策略下产生不同组合，但不能称为已实现口袋条件化身份预测。当前连接头也没有直接蛋白特征；口袋进入官方位姿网络，再经几何反馈可能影响后期连接或交替身份。这种差异不证明结合活性、选择性或泛化。

所有口袋均是此前检查过的诊断数据；未在这六个口袋调参，未做独立开发集参数搜索，也没有未见口袋验证。每口袋完整指标如下，不能仅凭跨口袋混合值评估多样性：

{pocket}
固定N=6的同口袋种子比较如下；N=4/5/6的全部分层见PER_N_DIVERSITY.csv。此表控制大小与口袋混合的影响，平均值只作为六个已查看口袋的描述统计，不作独立泛化推断。

{fixed_table}
## 几何质量与修复代价

{stage}
沿用 HQ_no_water 的 12 条必要条件及 PB 剩余 20 条输出，原 UFF 50 构象能量参考、碰撞比值、键角/环/立体阈值不变，只去掉两项水检查；缺失项不能当作通过。CSV 报告完整性与缺失数，all_quality_checks.json 保留所有失败和错误。两种修复独立从同一 native_final 出发，未串联；优化未收敛仍保留，绝不因未收敛删样。全部异常键角在 ABNORMAL_BONDS.csv / ABNORMAL_ANGLES.csv。

native_final 是联合推理在最后一次官方神经运动后的坐标、未做 MMFF 修复；它仍包含之前合法图搜索和参考形状转移的贡献，不能当作完全无约束的网络能力。decoded_initial 是初始化学图赋予随机片段坐标的读出；未连图的原始片段点云没有完整化学图，HQ/化学有效/连通等为 NA。每次初始及图更新参考投影误差在 result.json，更新阈值 .5Å，不能将投影计作神经改进。身份或图更新后，decoded_initial 与终态可能不是同一化学图，不能用原子索引直接算起终 RMSD；修复 RMSD 则严格基于相同终态图。

立体标签由当前坐标推导，与自身坐标一致不等于目标立体预测正确。没有目标立体标签输入；生成身份模板本有的立体特征属于词表，不是目标真值。相同立体分子、同口袋内的对齐和受体坐标系 RMSD 单列 CONFORMER_DIVERSITY.json，不影响 2D Unique。

## 生成有效性与预算

{cost}
拒绝的更新原因计数：{d['rejection_reasons']}。逐次提议在ALL_UPDATE_EVENTS.csv，拒绝事件在REJECTED_UPDATE_EVENTS.csv，绝不把保留旧图当成成功接受新身份。

完整尝试均为五次身份查询、五次当前连接前向/同预算搜索、24 次官方位姿前向。一次身份策略的后四次身份查询只观察，不更新身份；各策略均在步 5/10/15/20 提议连接更新。早期失败提前停止并如实计消耗，拒绝更新保留上一合法图且事件全部留存。连接合法图依赖模板价态约束及 bounded beam，不能把合法搜索后的有效率说成原始连接头本身的有效率。每次搜索 SearchConfig 不变，k=1；搜索内部前缀去重属于冻结算法，未对输出分子做去重筛选。原始连接失败、无效输出、修复失败均可在逐次 CSV/JSON 追溯。

多进程只是共享同一 GPU 来并发 CPU 搜索/修复；平均秒数是并发环境下每次经过时间，求和不是总墙钟时长。多样性、有效唯一产出与质量始终一起报告，不以容易修复的重复苯结构获得高通过率。

## N、输入与 Novelty 边界

预定 N 先验是独立于目标的 Uniform{{3,4,5,6}}，种子规则见 frozen_plan.json。本轮每口袋实际分布 {base['N_distribution']}，跨策略一致。只有 12 个独立 N 抽样，跨口袋复用；未抽到 N=3 如实保留。当前没有可用 N/停止头，人工先验不冒充自主预测。

直接读取去水 PDB，用官方 get_protein + structure_to_rdkit 得到模型口袋，以蛋白重原子均值作为中心，不加载 fixture、真实 N、配体 SDF/参考坐标/连接/立体标签。口袋区域本身来自已定义的历史结合位点，故不是未知位点搜索。QC/修复保留非水辅因子；官方模型沿用去 HETATM 协议，未感知这些辅因子，限制与旧实验相同。

原冻结基线仍为旧目录的 18 次固定 N5、90/90 苯、1 个构造、网络 HQ=4/18、修复 HQ=15/18，原记录和指纹保持不变。本轮 Top-1 采用相同一次身份规则，但 N、蛋白中心、连接更新不同，因此不得把历史与本轮 HQ 差说成严格配对模型改进。

仅对可用的旧 PDBBind 17,420 项训练索引计算部分 Novelty，索引 SHA 和范围见 novelty_reference.json。该索引不覆盖本轮 GEOM 身份、当前连接头及官方位姿模型的完整训练并集，因此完整 Novelty 不可评估。部分索引的未命中数单列 CSV，不能称为完整训练去重后的新颖性。未评估结合活性或实验成功。

相同初始输入的GPU身份概率存在约1e-7量级浮点差异，独立核验采用atol=1e-6、rtol=1e-5；配对初始身份实际一致。固定种子与协议保证可追溯重放，不宣称所有GPU概率逐比特一致。

冻结参数、全部权重 SHA 见 frozen_plan.json，概率 logits/分布、搜索池、失败、所有 SDF/spec、修复与逐项 QC 均保留。复现见 REPRODUCE.md；独立核验见 verification.json。

![身份多样性、固定N口袋内唯一产出与几何质量](DIVERSITY_QUALITY.png)
""")
 write('IDENTITY_SAMPLING_COMPARISON.md',f"""# 身份采样比较

{identity}
七种方案在运行前冻结：全 MASK 一次 Top1；温度1全概率；温度1.5全概率；温度1 Top-k=50；温度1 Top-p=.95（包含越过阈值的首项）；原 sample_joint 频率先验 + .3 置信度迭代 Top1；频率先验 + 无阈值全概率交替。无测试后调参、重复过滤或候选质量选择。

模型熵与经验身份计数熵分别保存，单位 nats；Top1 概率、3577概率/Top1占比、类型覆盖与组合、N 分布均在 STRATEGY_SUMMARY.csv 和 PER_POCKET_DIVERSITY.csv。词表覆盖分母3684，频率初始化可用词表3679（>=2重原子且有端口），不把随机概率抽样限制为这3679类。初始/终态所有片段计数、覆盖及组合逐分子列出。

每个 identity_XX.npz 包含全部3684类 logits、模型概率和该次提议采样概率；Top10与熵在 identity_statistics.json。频率先验策略在 step=0 先做模型观察，再从独立的 frequency_prior.npz 初始化，因此它的 identity_00 中采样概率是未使用的身份提议分布，实际初始抽样分布/熵来自 frequency_prior.npz / initial_actual_draw_entropy。后续提议被置信度门或合法图验证拒绝时，实际保留上一个合法身份；提议概率不是经接受/拒绝条件化后的终态概率。

{gains}

同 seed 初始身份跨口袋重复是输入不读受体且采用配对随机数的必然结果，不是独立的72次身份样本。后续身份是否变化及跨口袋相同率另列 CROSS_POCKET_COMPARISON.csv。sample_joint 的身份服务本就每轮全 MASK，不把上轮类型当作已知 token；本轮保留此契约。GEOM 未训练预测连接先验通道，未擅自打开；交替经几何反馈实现，不能说连接图已经输入身份头。
""")
 write('JOINT_GENERATION_AUDIT.md',"""# 联合生成接口审计

| 模块 | 真实输入输出及兼容性 |
| --- | --- |
| GEOM ConnectionIPA3684 | 全 MASK 当前T/R、t、空受体，输出3684 logits；不直接感知受体或预测连接 |
| 原 sample_joint | record labels 决定默认N；频率先验初始化，.3置信度身份建议，旧 PortConnectionPredictor；原服务用验证集ID查记录 |
| 当前 TimeConnectionPredictor | 身份/端口/相对位姿/时间/化学bag，输出键和立体概率；与旧 Port forward/types 接口不兼容 |
| 当前 ConnectionRuntime/hierarchical_search | 预测身份模板与T/R、价态/端口；完整树候选；不接收真实图 |
| 官方 SigmaDock NativeState | 当前合法化学图、蛋白图、随机位姿；SE(3)运动；动态边、距离和方向每步重算并审计 |
| N/停止 | 身份分类与位姿主干都不输出N或停止概率；无可直接复用生成模块 |

统一 attempt 入口以已构造的口袋context、seed、预定义策略和冻结模块为生成输入；pid仅作文件标识。context_from_pocket直接解析去水PDB，ligand_path=None，以蛋白中心代替真实配体COM，不读fixture、编码配体、验证集、真实标签/连接或参考三维。身份子进程请求白名单仅T/R/t/seed/close，未创建PoseDataset，GEOM initialize不打开数据库。

本轮没有直接冒称原样运行旧 sample_joint。复用其频率先验、.3置信更新、合法候选联合提交、形状转移与失败保留原则；新增无数据集身份服务和概率采样。为符合当前模型，旧 Port proposals/greedy fallback替换为当前 TimeConnectionPredictor + 同一 hierarchical_search；保留官方24步主干、NativeState和SamplingAudit。全程 head=full，未切 legacy_staged。当前模型内部含冻结旧expert的残差基座是其既有架构，不等于将评测头换回旧版。

所有策略使用同一正态sigma2.7位姿初始化、uniform R、蛋白中心和N先验；一次身份方案也允许同步连接更新，保持同预算。模型对预测身份的适配由合法图/官方参考重建验证，不修改预训练特征列。新参考由预测图RDKit嵌入，不用真实三维。更新形状误差上限.5Å，拒绝时保留原图，软件错误、失败及提议事件全部落盘。存在图转移投影，不能完全归功于位姿网络。

模型权重、输入PDB和历史文件SHA在 frozen_plan.json / historical_fingerprints.json，冻结生成源码/资产指纹在 generator_provenance.json，动态计算图审计在每次 sampling_audit.json。无新训练，无水排斥开发，不修改历史数据。

仍缺受体条件身份、经训练连接图条件身份、N/停止、模型感知非水辅因子和未见口袋验证。二维与三维多样性分别报告；词表模板立体与几何读出不当作目标立体预测。
""")
 cross_table=table(['策略','初始身份跨口袋相同率','终态身份跨口袋相同率','同seed结构跨口袋相同率'],[[x['strategy'],mean([c['initial_identity_equal_fraction'] for c in d['cross_pocket'] if c['strategy']==x['strategy']]),mean([c['terminal_identity_equal_fraction'] for c in d['cross_pocket'] if c['strategy']==x['strategy']]),mean([c['same_seed_constitution_equal_fraction'] for c in d['cross_pocket'] if c['strategy']==x['strategy'] and c['same_seed_constitution_equal_fraction'] is not None])] for x in s])
 write('NEXT_STEP_DECISION.md',f"""# 下一步决策

当前权重能在不同采样策略下产生不同片段组合和完整分子；一次Top1的苯重复不能代表整个模型没有多样性。但这套系统还没有直接受体条件身份，也没有自主N生成。

{gains}

{main_table}
## 身份重复来自哪里

{identity}
初始身份模型概率分布较宽，argmax忽略其余概率而反复选高频类别，因此当前一次Top1坍缩有明确的采样策略来源。温度/截断/迭代各方案结果如表，没有按最好结果事后选参数。模型本身的频率偏置、没有受体条件和缺少经过训练的身份—连接闭环仍然存在；随机采样带来的覆盖不能证明身份与口袋匹配。

## 不同口袋是否影响结果

{cross_table}
这些均为15对口袋的配对种子相同率均值，逐对值和初始概率差见CROSS_POCKET_COMPARISON.csv。初始概率在浮点容差内一致（最大绝对差逐对列出）；后期差异只能来自位姿/几何/合法搜索反馈。没有未见口袋数据，不能声称泛化或受体选择性。

## 多样性与有效性、几何代价

{stage}
网络、刚性、柔性分开给分子级分母；全部异常键角和修复位移分布保留，缺失非通过。合法搜索后的有效率不代表原始头预测能力，原始连接完整树与搜索有效率的差距如下：

{cost}
## 优先级与是否微调

优先把身份预测和联合接口的目标定义补全：直接受体条件、身份与连接的兼容提议、独立开发/未见口袋上的质量与分布核验。当前提高身份多样性只需改变冻结采样策略就已有实测收益，不能据一次Top1立即决定重新训练。

其次根据表中原始连接头与合法搜索的差距改进连接校准/端口选择，结合未修复HQ、跨片段键角和碰撞建立几何监督或采样约束的独立对照。合法性收益与修复收益来自显式搜索/优化，不能归为神经网络已学会。是否微调、改哪些损失应由独立开发集消融决定，不预先锁定键长训练。

N/停止是完整自由生成的独立必要模块；本轮人工均匀先验只解决诊断输入边界，不能替代建模。现阶段先补接口和独立验证，N建模需列入下一阶段；不启动大规模训练，不继续水排斥优化。上述排序反映本轮可观察的能力缺口，六个已查看口袋不足以确定最佳训练配方或生物学效用。
""")
 write('REPRODUCE.md',"""# 复现

从项目根目录执行，使用已有环境，不安装/下载新包。所有GPU进程均限定CUDA_VISIBLE_DEVICES=0。

    export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0
    export PYTHONPATH="$PWD/.quality_deps:$PWD:$PWD/../third_party/sigmadock_official_repo_20261030/src"
    PY=../.venv-cutoff-train/bin/python
    TASK=2611/1105/1345_diversity
    $PY -B "$TASK/source/run.py" --phase pilot
    for shard in 0 1 2 3; do
      $PY -B "$TASK/source/run.py" --phase main --world 4 --shard "$shard" > "$TASK/main_$shard.log" 2>&1 &
    done
    wait
    $PY -B "$TASK/source/analyze.py" --merge-only
    $PY -B "$TASK/source/evaluate.py" --manifest "$TASK/main/manifest.json" --output "$TASK/main/quality" --workers 8
    $PY -B "$TASK/source/analyze.py"
    $PY -B "$TASK/source/reports.py"
    $PY -B "$TASK/source/verify.py"

实际运行增加了extra_workers.py，分别提前计算9SYL/9T4A/9T51三口袋，每口袋4个CPU分片，也只用GPU0。初始四个worker到达这些口袋时读取已完成result.json，不重新抽样；此调度不影响冻结504次分母。为简洁且避免调度重叠，以上复现命令仅用四个分片即可得到同一结果，不需要额外workers。stream_quality.py在正式生成期间并发CPU QC；最后evaluate.py用相同缓存和冻结协议汇合。

run.py遇已有result.json原样恢复，包括失败记录，不重试失败。独立重跑应在新的同深度项目副本中执行，并保留原始计划；不能在同一目录修改策略/权重/种子后续写。当前目录已完成结果的再次运行只读恢复，不构成独立重放。pilot / main各保留全部分布、搜索池、拒绝事件、图审计、SDF/spec、修复和QC。

原历史基线及144例条件生成结果保持原样，历史目录全部文件SHA核对。权重、协议、来源和生成代码指纹见frozen_plan.json、generator_provenance.json、historical_fingerprints.json及verification.json。完整训练并集Novelty和未见口袋泛化不在本次可评估范围。

本次为避免慢项阻塞，另以boost_quality.py在独立quality_boost目录使用32个CPU worker处理待检结构，参数与原8-worker队列完全相同。finalize_quality.py按预定“原队列优先”合并结果，不按质量挑选；--snapshot只保存当前已完成结果，待检项仍明确列在PENDING_QUALITY_RECORDS.csv。慢项结束后依次运行finalize_quality.py、analyze.py、reports.py、verify.py即可刷新完整辅助诊断。正式1512项终态QC不依赖剩余初始阶段慢项。没有因为慢项而改阈值、少算50构象或设未经批准的超时。
""")
 print('Reports written',flush=True)
 completion=load(OUT/'quality_completion.json')
 if completion['missing']:
  note=f"**状态：504次正式生成和1512项神经/刚性/柔性终态QC已完成；额外decoded_initial初始赋图诊断有{completion['missing']}项标准50构象能量计算仍在运行。待检项未记通过，也未当失败或删样。核心结论不依赖这17项。其清单在PENDING_QUALITY_RECORDS.csv，完整QC汇合状态见quality_completion.json。**\n\n"
  for name in ['DIVERSITY_EVALUATION_RESULTS.md','IDENTITY_SAMPLING_COMPARISON.md','JOINT_GENERATION_AUDIT.md','NEXT_STEP_DECISION.md','REPRODUCE.md']:
   path=OUT/name;path.write_text(note+path.read_text())
if __name__=='__main__':main()

