# 任务三最终报告：独立口袋泛化测试

报告日期：2026-11-04。固定模型和阈值后，完成 **6 个测试受体×3 seeds×Top-8=144 个候选**；K1/4/8 使用同池嵌套前缀。开发集仅 1 个可用受体、3 seeds、24 个候选。本次整理没有重跑神经轨迹。

**几何修复在外部小样本中出现同方向改善，但 CASF 的高通过率没有复制。** 公平地给予相同三维立体赋值后，Top-8 HQ 从 2/144（1.39%）增至 7/144（4.86%），不同构造从 1 增至 6，有产出的 pocket/seed 从 2/18 增至 4/18。最终 HQ 集中在 3/6 个受体。该绝对通过率低于 CASF 70/288（24.31%）；两队列不是配对实验，不能把差异全归因于某一模型或蛋白性质。

## 独立性审计与实际覆盖

排除完整 PDBBind v2020 的 19,443 个成员，而非仅本地 17,420 个训练成员；排除官方 PoseBusters/Astex 数据包全部 513 个成员（比论文 PB v2 308＋Astex 85 的口径更宽）。同时排除冻结 GEOM 所有 split 及本地已建立 PDBBind 配体索引的相同二维构造。蛋白所有链在 RCSB 30% 序列簇均无上述排除集成员，直接序列搜索 identity≥30%、E-value≤1 也无命中。开发、测试及各受体之间没有共享 30% 簇。角色在生成前依身份 hash 固定，不按生成质量分组。

这些检查支持“对已审计模型数据清单和序列规则独立”；没有做完整局部口袋三维相似性比对，不能保证所有远缘结构相似均已排除。近期公开时间只作辅助，不代替成员和序列审计。官方训练/验证来源见 [SigmaDock 论文 §3.1](https://arxiv.org/html/2511.04854v2#S3.SS1)，结构、CCD 和链序列来源见 [RCSB Data API](https://data.rcsb.org/)。原查询、响应、缓存和指纹保留于 `data_audit/`、`independent_inputs/`。

| PDB | 用途 | 蛋白 | 30% 簇 | OOV 片段 | 最大模板拟合误差 Å |
|---|---|---|---|---|---|
| 9YSP | development | dCTP pyrophosphatase 1 | 11038 | 0 | 0.426 |
| 9RT3 | test | Choline O-acetyltransferase | 3916 | 0 | 0.059 |
| 9T4A | test | Diphosphoinositol polyphosphate phosphohydrolase 1 | 3510 | 2 | 1.627 |
| 7ILZ | test | Nucleoprotein | 1277 | 0 | 0.010 |
| 9T51 | test | StAR-related lipid transfer protein 4 | 9710 | 1 | 0.457 |
| 9SYL | test | tRNA (cytosine(38)-C(5))-methyltransferase | 7921 | 0 | 0.163 |
| 12AV | test | Ornithine aminotransferase, mitochondrial | 1408 | 2 | 0.103 |

三个 OOV 测试受体均使用显式片段化学结构和独立未结合模板；9T4A 最大模板拟合误差 1.627 Å 仍保留，未按拟合质量筛除。拟合误差只用于输入表示审计，不是排序或质量优化目标。真实身份、N、OOV 化学及配体邻域口袋仍为 Oracle，因此这不是自由 de novo；也没有验证真实实验结合、活性或药效。

原计划开发 3、测试 12 个输入；最终完整可用为开发 1、测试 6。未进入采样的输入和原因如下，均在生成质量评价之前处理：

| 输入 | 阶段/角色 | 未进入采样的原因 |
|---|---|---|
| 9QJF | development | 配体重原子缺失 |
| 9SOD | development | 配体重原子缺失 |
| 9TAR | test | 共价配体–受体连接不在冻结非共价管线范围 |
| 34JY | test | 配体重原子缺失 |
| 39MR | 源输入范围 | Free amino acid rather than drug-like bound ligand |
| 9YEB | 源输入范围 | Metal-dependent cofactor requires chemistry outside this noncovalent pipeline |
| 9RMV | 源输入范围 | Polyethylene glycol crystallization additive |
| 28LK | 源输入范围 | Membrane lipid rather than drug-like bound ligand |

输入失败不假装是“实际采样后的有效分子”；主质量表限定于完整可用的 18 个已锁定尝试，分母始终为 18K。因缺失/不支持而未进入运行的源输入另表公开，不能据此宣称覆盖所有原计划蛋白。完整证据见 [independent_input_coverage.csv](independent_input_coverage.csv)、[independent_fixture_preparation.json](independent_fixture_preparation.json)、[input_scope_exclusions.json](input_scope_exclusions.json)、[cohort_lock.json](cohort_lock.json)。

## 四版本与 K1/4/8 结果

raw 是网络生成坐标和构图输出的原显式立体标签；raw_stereo 只做可靠三维立体赋值；repaired 只做刚性修复并保留原标签；stereo 为修复＋可靠赋值。**隔离修复作用的主比较为 raw_stereo→stereo**；仅 repaired 缺少同等立体赋值，不能作为公平的最终质量对照。

| K | 版本 | 合法连通 | 无严重内碰撞 | 键长 | 键角 | 扭转/酰胺 | PB mol | PB dock | HQ/计划分母（率） | 不同 HQ 构造 | 有 HQ/18 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | raw | 18 | 10 | 12 | 10 | 16 | 8 | 0 | 0/18（0.0%） | 0 | 0 |
| 1 | raw_stereo | 18 | 10 | 12 | 10 | 16 | 8 | 0 | 0/18（0.0%） | 0 | 0 |
| 1 | repaired | 18 | 17 | 18 | 15 | 18 | 15 | 1 | 0/18（0.0%） | 0 | 0 |
| 1 | stereo | 18 | 17 | 18 | 15 | 18 | 15 | 1 | 1/18（5.6%） | 1 | 1 |
| 4 | raw | 72 | 42 | 39 | 41 | 67 | 28 | 0 | 0/72（0.0%） | 0 | 0 |
| 4 | raw_stereo | 72 | 42 | 39 | 41 | 67 | 28 | 0 | 0/72（0.0%） | 0 | 0 |
| 4 | repaired | 72 | 66 | 66 | 59 | 72 | 56 | 2 | 0/72（0.0%） | 0 | 0 |
| 4 | stereo | 72 | 66 | 66 | 59 | 72 | 56 | 2 | 2/72（2.8%） | 2 | 2 |
| 8 | raw | 144 | 81 | 74 | 82 | 129 | 56 | 2 | 0/144（0.0%） | 0 | 0 |
| 8 | raw_stereo | 144 | 81 | 74 | 82 | 129 | 56 | 2 | 2/144（1.4%） | 1 | 2 |
| 8 | repaired | 144 | 130 | 132 | 118 | 144 | 110 | 7 | 0/144（0.0%） | 0 | 0 |
| 8 | stereo | 144 | 130 | 132 | 118 | 144 | 110 | 7 | 7/144（4.9%） | 6 | 4 |

每个几何列是该检查独立通过的数量，不能相加推导联合通过数；PB mol 是分子内部检查，PB dock 是保留水检查的标准完整 dock 配置。当前 144 个原始候选全部有效且连通，无生成缺额；若发生缺额或失败，仍保留计划分母。网络直接 raw HQ 为 0/144，raw_stereo 的 2/144 已含后处理，最终 7/144 更不能表述为网络直接质量。仅 repaired 各 K 的 HQ 均为 0，主要提示完整立体定义也是联合门槛的一部分。

## 配对提升、逐受体产量和残余问题

| K | HQ 率配对差（百分点） | 受体聚类 95% 描述性区间（百分点） |
|---|---|---|
| 1 | 5.556 | 0.000–16.667 |
| 4 | 2.778 | 0.000–5.556 |
| 8 | 3.472 | 0.694–6.250 |

以 6 个受体的完整 RCSB30 簇为单位，保留该受体所有 seeds 和候选，配对 bootstrap 20,000 次。K8 差值为 +3.472 个百分点，区间约 [0.694, 6.250]。这些是小样本、选择后的受控队列的描述性区间，不能当作总体泛化保证；K1/4 区间包括 0。

| 受体 | raw_stereo HQ/24 | 修复＋赋值 HQ/24 | 不同构造：前→后 | 有 HQ 尝试：前→后/3 | 最终构造对 Tanimoto |
|---|---|---|---|---|---|
| 12AV | 0 | 0 | 0→0 | 0→0 | — |
| 7ILZ | 0 | 0 | 0→0 | 0→0 | — |
| 9RT3 | 2 | 4 | 1→3 | 2→2 | 0.345 |
| 9SYL | 0 | 1 | 0→1 | 0→1 | — |
| 9T4A | 0 | 0 | 0→0 | 0→0 | — |
| 9T51 | 0 | 2 | 0→2 | 0→1 | 0.443 |

最终不同构造平均为 1.0/受体；K1/4 分别为 1/6、2/6。两个相似性数值仅来自至少两个 HQ 构造的受体，其他空白表示不可计算。逐候选/口袋完整表见 [independent_candidate_results.csv](independent_candidate_results.csv)、[independent_pocket_yields.csv](independent_pocket_yields.csv)。

| 冻结检查 | raw_stereo 未通过/144 | stereo 未通过/144 |
|---|---|---|
| 无严重内部碰撞 | 63 | 14 |
| 键长 | 70 | 12 |
| 键角 | 62 | 26 |
| 扭转/酰胺 | 15 | 0 |
| 无严重蛋白碰撞 | 67 | 5 |
| 立体定义完整 | 15 | 4 |
| PB mol | 88 | 34 |
| PB dock（含水） | 142 | 137 |

各失败项可重叠。修复明显降低分子内部几何及蛋白严重碰撞失败，但最终仍只有 7/144 通过标准 PB dock。保持原主指标不变时，“去掉水相关检查”的辅助 PB 指标为 106/144，其中 99 个候选只在放宽水检查后满足该辅助 PB 定义；这不是综合 HQ，也不证明它们具有可接受的真实结合。该差异提示水相关接触判定是当前严格 PB 结果的重要限制，后续应独立审计输入水与接触规则，不能用事后删水替换本次主结果。

## 化学覆盖：源配体分层与实际生成分开

9RT3 的源配体包含跨片段双键与 E/Z，其 24 个候选单独归入源化学分层；其余 5 个测试受体源连接为非环单键。**全部 144 个实际候选均未提出跨片段双键或 E/Z 连接**，包括 9RT3 的全部 24 个候选。因此“输入含 E/Z”不等于“生成的 E/Z 连接几何已经验证”，该受体的合格替代构造不能用作双键/EZ 泛化成功证据。

冻结分割没有提供跨片段环连接样本，当前搜索要求片段图为树，不支持跨片段成环拓扑；片段内部自带的环不构成跨片段环覆盖。当前证据仅支持实际生成域内的结论；非单键、E/Z、跨片段环仍有明确覆盖缺口。源分层见 [independent_chemistry_strata.csv](independent_chemistry_strata.csv)，实际候选连接类型见 [generated_connection_chemistry.csv](generated_connection_chemistry.csv)。不把不支持的类型静默排除后声称推广成功。

## 成本与端到端计时

| K | 姿态 NFE | 连接头前向 | 组件秒数之和 | NFE/不同 HQ 构造 | 组件 s/不同 HQ 构造 | 实际研究墙钟 s |
|---|---|---|---|---|---|---|
| 1 | 1800 | 18 | 983.285 | 1,800.0 | 983.285 | — |
| 4 | 6120 | 18 | 2,052.275 | 3,060.0 | 1,026.137 | — |
| 8 | 11880 | 18 | 3,256.655 | 1,980.0 | 542.776 | 2,030.324 |

| K | 前缀 s | 完整候选搜索 s | 构图 s | 短程 s | 长程 s | 修复＋赋值 CPU s | 最终版 QC CPU s |
|---|---|---|---|---|---|---|---|
| 1 | 108.337 | 519.057 | 9.915 | 52.365 | 182.644 | 11.168 | 99.800 |
| 4 | 108.337 | 519.057 | 40.837 | 223.431 | 789.825 | 40.331 | 330.458 |
| 8 | 108.337 | 519.057 | 84.741 | 396.345 | 1,382.432 | 77.396 | 688.347 |

姿态 NFE 为 18×(20＋80K)，连接头各 18 次前向另计。K1/4 是同一 K8 真实轨迹的嵌套成本核算：保留完整搜索成本，只累加选中前缀候选的构图/短程/长程/后处理/最终版 QC，没有重复计算神经轨迹。没有实际单独测量 K1/4 部署墙钟。

K8 实际研究墙钟 **2,030.324 秒（33.84 分钟）**，包含模型载入、完整生成、I/O、修复和原三版本质量评测。该运行含多 GPU/多进程协作及等待，共享 A100 负载；完整并行记录见 [parallel_execution.json](independent/parallel_execution.json)。组件和混合 CPU 累计时间与分支计时，既不等于这一墙钟，也不等于串行部署延迟。raw_stereo 作为后补对照另用 137.938 秒、新增神经 NFE=0，不含于原 K8 墙钟。开发集完整研究墙钟为 385.884 秒、姿态 NFE=1,980；测试为 11,880。本次任务一至三总新增姿态 NFE=13,860，不重复累计反事实 K1/4。

计时证据：[independent_K_results.csv](independent_K_results.csv)、[independent/end_to_end.json](independent/end_to_end.json)、[independent_raw_stereo/timing.json](independent_raw_stereo/timing.json)、[tasks123_compute_ledger.csv](tasks123_compute_ledger.csv)。完整 16＋64 步续采样与模型未变化证据见 [independent/native/protocol.json](independent/native/protocol.json)。

任务三结论：冻结末端修复带来的局部改善在这组外部受体上重现；高绝对 HQ 产量及完整化学类型的泛化尚未建立。应保留末端修复、严格水检查和覆盖披露，扩大独立样本与输入表示审计后再判断普适性。

只重建报告、图和核验：`bash fbdd_efficiency_generalization_20261103/reproduce_tasks123_reports.sh`。所有统计和冻结指纹由 [tasks123_verification.json](tasks123_verification.json) 对照原始产物核验。
