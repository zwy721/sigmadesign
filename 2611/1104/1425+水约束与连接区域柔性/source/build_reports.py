"""Build reviewable Chinese reports from final metrics, with separate definitions."""
from study_common import *


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(v) for v in r)+' |' for r in rows])


def f(x):
    return f'{x:.4f}' if isinstance(x,(int,float)) else str(x)


def main():
    s=load(OUT/'summary.json');allrows=load(OUT/'per_candidate_results.json')
    ss={r['variant']:r for r in s['overall']};pairs={r['variant']:r for r in s['paired_changes']}
    a,c,b=ss['A1'],ss['A0_B0'],ss['B1']
    va=load(OUT/'verification.json')
    ta,tb=load(OUT/'timing_A.json'),load(OUT/'timing_B.json')
    common='''本轮仅对已有 144 个候选做 CPU 数值修复，新增神经网络 NFE=0。6 个受体 × 3 个原种子 × 8 个候选；连接图、键级、原子编号、片段划分、固定环境及原 QC 全部保留。所有版本都从相同原始坐标开始，以相同原始立体标签准备 MMFF，再经过相同 SDF 回读和既有 16 次扰动立体赋值，最后执行原完整 QC，包括标准 PB 的 50 构象能量检查。无参考构象进入优化或选择。

严格 HQ 是原 required_checks 的联合通过；辅助 no-water HQ 只移除 PB 中固定水距离与固定水体积两个检查，仍保留其他全部门槛。后者仅用于诊断。所有数字均为本轮后处理结构质量，不能表述为模型原始输出通过率。A0/B0 是同一份控制，只运行及统计一次。

独立开发支持不足：现有仅一个开发受体，并且本批标签此前已被查看。本轮在运行前锁定单组参数，无参数扫描、无结果驱动调参；结论为这 6 个受体上的配对描述，不声称独立泛化。未运行可选 B2，也未组合 A1 与 B1。
'''
    def quality(vs):
        return table(['版本','严格 HQ /144','辅助 no-water HQ /144','PB /144','水距离失败','水体积失败','键长通过','键角通过','内部碰撞通过','环境碰撞通过','立体定义完整','立体几何一致','优化收敛 /144'],[[v,ss[v]['strict_HQ'],ss[v]['auxiliary_HQ_no_water'],ss[v]['PB_dock'],ss[v]['water_distance_fail'],ss[v]['water_volume_fail'],ss[v]['QC_bond_geometry_valid_pass'],ss[v]['QC_angle_geometry_valid_pass'],ss[v]['QC_no_severe_internal_clash_pass'],ss[v]['QC_no_severe_protein_clash_pass'],ss[v]['QC_stereo_definition_complete_pass'],ss[v]['QC_stereo_geometry_consistent_pass'],ss[v]['optimizer_converged']] for v in vs])
    def grouped(vs):
        pp=sorted({r['receptor_id'] for r in s['by_receptor']})
        ix={(r['variant'],r['receptor_id']):r for r in s['by_receptor']}
        return table(['受体','版本','严格 HQ /24','辅助 no-water HQ /24','水距离失败','异常键数(跨/内)','异常角数(跨/内)'],[[p,v,ix[v,p]['strict_HQ'],ix[v,p]['auxiliary_HQ_no_water'],ix[v,p]['water_distance_fail'],f"{ix[v,p]['cross_fragment_abnormal_bonds']}/{ix[v,p]['within_fragment_abnormal_bonds']}",f"{ix[v,p]['cross_fragment_abnormal_angles']}/{ix[v,p]['within_fragment_abnormal_angles']}"] for p in pp for v in vs])
    def movements(vs):
        return table(['版本','相对原始 RMSD 均值/最大 Å','相对控制 RMSD 均值/最大 Å','相对控制最大单原子位移 Å','相对控制质心移动均值/最大 Å','修复 CPU 总秒','修复 CPU 中位秒','QC CPU 总秒'],[[v,f"{f(ss[v]['RMSD_from_raw_no_alignment_A_mean'])}/{f(ss[v]['RMSD_from_raw_no_alignment_A_max'])}",f"{f(ss[v]['RMSD_from_control_no_alignment_A_mean'])}/{f(ss[v]['RMSD_from_control_no_alignment_A_max'])}",f(ss[v]['max_atom_displacement_from_control_A_max']),f"{f(ss[v]['centroid_displacement_from_control_A_mean'])}/{f(ss[v]['centroid_displacement_from_control_A_max'])}",f(ss[v]['repair_cpu_seconds']),f(ss[v]['repair_cpu_seconds_median']),f(ss[v]['QC_cpu_seconds'])] for v in vs])
    def stereo(v):
        r=ss[v]
        return f"相对 A0/B0，原本与修复后都可靠但 CIP/EZ 编码改变的中心总数 {r['stereo_reliable_flips_from_control']}，可靠性丢失 {r['stereo_reliability_lost_from_control']}，恢复可靠 {r['stereo_reliability_regained_from_control']}。相对相同原始坐标的独立立体赋值，可靠编码改变 {r['stereo_reliable_flips_from_raw']}、丢失 {r['stereo_reliability_lost_from_raw']}、恢复 {r['stereo_reliability_regained_from_raw']}。这些是几何稳定性诊断，不代表对外部参考立体构型的正确率。"
    def gains(v):
        p=pairs[v]
        return f"严格 HQ 新通过 {p['strict_HQ_gained']}、原通过丢失 {p['strict_HQ_lost']}；辅助 no-water HQ 新通过 {p['auxiliary_HQ_no_water_gained']}、丢失 {p['auxiliary_HQ_no_water_lost']}。"
    a_report=f'''# 固定水软排斥单因素实验 A

{common}
输入只读审计复核：严格 HQ=7/144，辅助 no-water HQ=103/144，仅水阻断=96/144；水距离失败=133/144，水体积失败=0/144。控制回放的全部序列化坐标、立体 sidecar、每候选 required QC 和 PB 二值结果均与历史完全相同。

## 固定方法

A0/B0 调用未修改的既有 `geometry_constraints_20261103/repair.py`。A1 保留相同 MMFF94s、SE(3) 自由度、零初始刚体参数、加氢准备、非水项及 800 次 L-BFGS-B 预算，只增加固定水氧对全部分子重原子的平滑软排斥：

`d_ij = sqrt(||x_i-w_j||² + ε²)`，`L_ij = 0.80(r_i+r_j)`，`E_water = k Σ max(L_ij-d_ij,0)²`。

`∂E/∂x_i = -2k Σ max(L_ij-d_ij,0) (x_i-w_j)/d_ij`。

预锁定 k=200 kcal mol⁻¹ Å⁻²、ε=10⁻⁶ Å；水 O 的 RDKit vdW 半径 1.55 Å，其余按原 PB 使用同一半径表。平方铰链为 C1 连续，正则距离避免零距离梯度奇点。权重和半径比沿用原非水软排斥设置，未做测试集搜索。优化作用距离按原子半径和 0.80 比例确定，全部 79 个原水氧坐标和编号保留；没有按距离或通过结果删除水。QC 水距离门槛仍为 0.75 半径和，体积门槛仍为 0.075（vdW_scale=0.5）。每次迭代仅在作用半径内有非零惩罚，作用半径外能量和梯度为零。

L-BFGS-B 的 maxiter=800、ftol=1e-10、gtol=1e-4、maxls=40，氢准备最多 200 次，均保持原样。退出原因、能量、迭代次数、初始/最终水能量及耗时见每候选 `repair.json` 和 CSV。未收敛候选仍保留最后坐标并执行全部 QC；优化收敛与质量通过分别统计。

## 结果

{quality(['A0_B0','A1'])}

{gains('A1')}

{movements(['A0_B0','A1'])}

A1 的相对控制无对齐 RMSD 最大 {f(a['RMSD_from_control_no_alignment_A_max'])} Å，最大单原子变化 {f(a['max_atom_displacement_from_control_A_max'])} Å。原有全局软位移项仍相对原始坐标起作用；未把对齐后误差用于掩盖整体移动。

{stereo('A1')}

A1 的水距离改善伴随明确的坐标漂移和非水几何代价：辅助 no-water HQ 丢失 19 个、蛋白碰撞检查新失败 20 个，键角和酰胺项也出现退化。本轮据此停止增强水排斥权重，不将 65/144 表述成无损修复。

{grouped(['A0_B0','A1'])}

受体是来源分组单位；种子分组见 `by_seed.csv`。6 受体整组重抽样 5000 次得到严格通过率配对差的描述性区间 {pairs['A1']['descriptive_cluster_bootstrap_HQ_rate_delta_95_interval']}，受体少且标签已知，不能解释为独立泛化检验。

## 文件与限制

逐候选 `per_candidate_results.csv`，逐水最短距离和碰撞 `water_contacts.csv`，接触对距离/半径/能量/梯度 `water_pair_details.csv`（6 Å 内所有原子对及每个水的最近和最小相对距离对，完整水目录不变），全部失败及未收敛 `failure_cases.csv`，相对控制退化 `regressions.csv`；所有水坐标和编号在 `fixed_waters.json`。原始 SDF 在 `structures/original/`，控制和 A1 的修复前立体赋值 SDF、赋值后 SDF、sidecar、立体中心决定和优化退出信息在 `structures/A0_B0/` 与 `structures/A1/`。半径表达式、梯度、权重和停止条件在 `fixed_parameters.json`，全部原 QC 单项在 `quality_A/molecule_metrics.jsonl`。

本轮只执行该单组水权重；没有因仍失败而增强排斥。CPU 为进程时间之和，修复与完整 QC 分开；8 个进程的实验 A 实测总墙钟 {f(ta['wall_seconds'])} 秒，不能将 CPU 和或历史 GPU 时间当作端到端延迟。

复现：`bash reproduce.sh /tmp/water_port_reproduction`，从本目录执行；详细命令见 `RUN_COMMANDS.md`。
'''
    b_report=f'''# 连接区域有限柔性实验 B

{common}
## 固定方法

B0 为相同 A0/B0。B1 从与控制相同的原始坐标、相同氢准备和 SE(3) 零参数开始，开放跨片段键两端及各自第一图邻居的局部位移。优化坐标为 `(原始片段局部坐标+δ) Rᵀ+片段中心+平移`；其余原子仅随各自片段刚性运动。所有非水目标项、MMFF94s、原姿态与方向正则、优化器及 800 次预算保持原样，没有水排斥项。

单组预先锁定的 δ 每分量限制为 ±0.10/√3 Å，因此每原子局部位移范数保证≤0.10 Å；可行域是立方体，并非球形约束。额外位置正则为 `200 Σ||δ||²` kcal/mol，梯度为 `400δ`，用于协调局部改变。0.10 Å 仅为本轮预定诊断设置，不是物理真值，也未在 0.05/0.10/0.15 Å 中事后选优。限制的是相对自身片段刚体运动的额外位移，整体坐标仍可按原 SE(3) 改变；不得将二者混同。

B2 全原子强约束诊断未运行。B1 与 A1 分开，不报告混合方法。

## 结果

{quality(['A0_B0','B1'])}

{gains('B1')}

{table(['版本','跨片段异常键数','片段内异常键数','跨片段异常角数','片段内异常角数','异常酰胺扭转数','严重内部碰撞对','严重蛋白碰撞对','跨片段键 MAE 均值 Å','片段内键 MAE 均值 Å'],[[v,ss[v]['cross_fragment_abnormal_bonds'],ss[v]['within_fragment_abnormal_bonds'],ss[v]['cross_fragment_abnormal_angles'],ss[v]['within_fragment_abnormal_angles'],ss[v]['abnormal_amides'],ss[v]['severe_internal_clashes'],ss[v]['severe_protein_clashes'],f(ss[v]['cross_fragment_bond_MAE_A_mean']),f(ss[v]['within_fragment_bond_MAE_A_mean'])] for v in ['A0_B0','B1']])}

键偏差是原 RDKit 距离几何界上下限中点的绝对差，不是相对晶体构象的误差。键角门槛仍采用原端点距离界，明细另报角度数值。`torsion_details.csv` 分别保留酰胺平面性和 ETKDG 匹配扭转诊断；只有原 QC 定义的酰胺阈值参与相应检查，未为一般扭转另造通用阈值。

{movements(['A0_B0','B1'])}

B1 额外局部位移最大 {f(va['max_local_displacement_A'])} Å，未开放区域的片段内成对距离最大误差 {va['max_fixed_region_pair_distance_error_A']:.3g} Å（优化内存中计算）。SDF 仅保留 4 位小数，文件回读后约 10⁻⁴ Å 的差异属于序列化舍入，局部上限验证采用未舍入优化参数，同时 QC 始终使用实际 SDF 坐标。

{stereo('B1')}

B1 主要修复内部碰撞和部分跨片段角度。异常跨片段键总数仍为 12；片段内异常键/角仍为 0，但片段内键偏差均值从 0.0297 增至 0.0315 Å，存在轻微代价。立体定义完整从 140/144 增至 144/144，但有 1 个原本与控制均可靠的编码改变，不能声称立体构型完全保持。原本可靠的中心没有变为不可靠。

{grouped(['A0_B0','B1'])}

6 受体整组重抽样的严格通过率配对差描述性区间 {pairs['B1']['descriptive_cluster_bootstrap_HQ_rate_delta_95_interval']}；全部种子和候选保留，不能将 144 个候选当作 144 个独立来源。

## 交付与回归

原始和各版本 SDF 在 `structures/`；逐键、角、扭转、碰撞、立体和原子位移分别在 `bond_details.csv`、`angle_details.csv`、`torsion_details.csv`、`clash_details.csv`、`stereo_centers.csv`、`atom_displacements.csv`。`per_candidate_results.csv` 含每项原 QC、局部/整体位移、优化终止原因和 CPU 时间；`failure_cases.csv` 保留质量失败与全部未收敛样本。`fixed_parameters.json` 为运行前锁定的全部参数。

新增回归覆盖水软排斥解析梯度、含 SE(3)/MMFF/氢的完整目标梯度、关闭新增项与原修复完全相同、局部位移上限和图保持。全队列核验还检查原控制的每候选坐标/立体/QC 复现、原输入哈希、原子编号、片段划分、原水全保留、立体赋值重现及单因素分离。新增与既有相关单测共 11 项通过，全队列 30 项核验通过；结果见 `regression_tests.log`、`verification.json`。

实验 B（含完整 QC）实测墙钟 {f(tb['wall_seconds'])} 秒；与实验 A 顺序执行，CPU 总和与墙钟分列。复现：`bash reproduce.sh /tmp/water_port_reproduction`（在本目录）；执行命令见 `RUN_COMMANDS.md`。
'''
    (OUT/'WATER_AWARE_REPAIR_RESULTS.md').write_text(a_report)
    (OUT/'LOCAL_PORT_FLEXIBILITY_RESULTS.md').write_text(b_report)
    # Decision prose is generated from measured totals, with explicit regressions and evidence limits.
    summary=f'''# 下一步决策汇总

本轮已完成输入审计 → A0/A1 配对 CPU 修复及完整 QC → B0/B1 单因素比较及完整 QC → 回归核验。新增神经网络 NFE=0；原 144 个候选全部保留，原始和历史结果未覆盖。

{quality(['A0_B0','A1','B1'])}

严格 HQ 与辅助 no-water HQ 定义始终分开；以上均为保存候选的后处理通过率，不是神经网络原始输出质量。

## 1. 水软排斥的收益和代价

严格通过 {c['strict_HQ']}→{a['strict_HQ']}/144，水距离失败 {c['water_distance_fail']}→{a['water_distance_fail']}/144。{gains('A1')} 相对控制无对齐 RMSD 均值/最大 {f(a['RMSD_from_control_no_alignment_A_mean'])}/{f(a['RMSD_from_control_no_alignment_A_max'])} Å，最大原子位移 {f(a['max_atom_displacement_from_control_A_max'])} Å；相对原始坐标 RMSD 均值 {f(a['RMSD_from_raw_no_alignment_A_mean'])} Å。异常键总数 {c['abnormal_bonds']}→{a['abnormal_bonds']}，异常角 {c['abnormal_angles']}→{a['abnormal_angles']}，严重内部碰撞对 {c['severe_internal_clashes']}→{a['severe_internal_clashes']}，严重蛋白碰撞对 {c['severe_protein_clashes']}→{a['severe_protein_clashes']}。各受体差异见详细报告与 `by_receptor.csv`。已观察到明确的漂移与非水几何损失，因此只保留预锁定权重结果，停止增强水排斥。

## 2. 连接局部柔性的收益和代价

B1 严格通过 {c['strict_HQ']}→{b['strict_HQ']}/144，辅助 no-water HQ {c['auxiliary_HQ_no_water']}→{b['auxiliary_HQ_no_water']}/144。异常键 {c['abnormal_bonds']}→{b['abnormal_bonds']}、异常角 {c['abnormal_angles']}→{b['abnormal_angles']}；片段内异常键 {c['within_fragment_abnormal_bonds']}→{b['within_fragment_abnormal_bonds']}，片段内异常角 {c['within_fragment_abnormal_angles']}→{b['within_fragment_abnormal_angles']}，严重内部碰撞对 {c['severe_internal_clashes']}→{b['severe_internal_clashes']}。{gains('B1')} 片段内键偏差均值 {f(c['within_fragment_bond_MAE_A_mean'])}→{f(b['within_fragment_bond_MAE_A_mean'])} Å，未触发新增片段内键/角异常。最大额外局部位移 {f(va['max_local_displacement_A'])} Å；最大未开放区域距离误差 {va['max_fixed_region_pair_distance_error_A']:.3g} Å。

{stereo('B1')}

## 3. 质量、稳定性与 CPU 成本

{movements(['A0_B0','A1','B1'])}

严格质量净增：A1 为 {a['strict_HQ']-c['strict_HQ']:+d}，B1 为 {b['strict_HQ']-c['strict_HQ']:+d}。完整 QC 成本单列，优化收敛数分别为 {c['optimizer_converged']}、{a['optimizer_converged']}、{b['optimizer_converged']}/144；未收敛并未被删除，质量通过也不等同数值收敛。具体退化候选和退出原因见 `regressions.csv`、`failure_cases.csv`。

以本轮严格通过数量为目标，A1 的计算收益更高：修复 CPU 为控制的 {a['repair_cpu_seconds']/c['repair_cpu_seconds']:.3f} 倍，仅增加 {a['repair_cpu_seconds']-c['repair_cpu_seconds']:.2f} 个 CPU 秒，共新增 58 个严格通过；B1 为 {b['repair_cpu_seconds']/c['repair_cpu_seconds']:.3f} 倍，但严格通过数未增加。若优先考虑坐标稳定和非水几何，B1 的平均相对控制位移更小，内部碰撞下降明显，适合作为后续局部几何诊断方向。

因此，下一轮优先在独立开发样本上验证 A1 的水距离收益与漂移/碰撞约束，保持本轮权重不再增强；B1 的少量非水改善不足以成为当前严格质量瓶颈的主解法。A1 存在 19 个辅助通过丢失、20 个蛋白碰撞新失败及 1 个可靠编码变化，B1 也有 1 个辅助通过丢失和 1 个可靠编码变化；两者都需要独立证据再决定是否替代控制。CPU 为单次运行测量，不声称微小时间差具有统计显著性。

## 4. 证据边界与后续要求

本批 6 个受体、18 个种子组的标签先前已知，现有开发支持不足以选择参数或证明泛化。每个受体的 24 个候选有共同来源，只做配对描述和受体整组重抽样；受体少时区间不稳定，不能因总计 144 而认定证据充分。

本轮足以检验两个单因素在这些保存数据上的直接数值作用，也识别了新增几何和立体风险；不足以直接支持将更复杂的联合优化作为默认方案。后续具体缺少：独立开发集上预定参数/预算/停止规则，未查看标签的新受体测试，多来源/连接类型的分层覆盖，位移与可靠立体保持的预定验收标准，以及优化收敛和重复执行稳定性的独立验证。后续若研究 A+B 组合，须在参数冻结后另作对照，不能把本轮两个独立效果相加。

原输入及数值代码指纹在 `input_fingerprints.json` 与 `execution_lock.json`，全量核验 {va['check_count']} 项通过，明细见 `verification.json`。复现命令和交付导航见 `RUN_COMMANDS.md`、`README.md`。
'''
    (OUT/'NEXT_DECISION_SUMMARY.md').write_text(summary)
    print('Three reports generated.',flush=True)

if __name__=='__main__':main()
