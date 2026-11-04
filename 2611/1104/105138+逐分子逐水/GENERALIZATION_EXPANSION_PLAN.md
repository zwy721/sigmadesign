# 扩大独立评测的方案与下一轮单项对照

日期：2026-11-04。本文件给出方案；本轮没有启动新的大规模采样、训练或搜索，不修改当前冻结主评测条件。

## 当前覆盖为何只有 1 个开发受体和 6 个测试受体

| 筛选阶段 | 单位 | 数量 |
|---|---|---|
| RCSB query total matches | entry | 5438 |
| Retrieved first page only | entry | 2000 |
| Initial no-overlap 30pct clusters | entry | 917 |
| Fetched entry metadata cache | entry | 300 |
| Chemical inventory after initial entry/ligand filters | ligand_entity | 135 |
| After entry/constitution/cofactor/cluster exclusions | ligand_entity | 57 |
| After direct sequence exclusion | ligand_entity | 38 |
| Greedy cluster and constitution deduplication | ligand_entity | 19 |
| Identity-hash ordered cap | ligand_entity | 15 |
| After input scope exclusion | ligand_entity | 11 |
| Complete fixture development | receptor | 1 |
| Complete fixture test | receptor | 6 |

各阶段单位不同，不能用 300−135 当作“165 个蛋白被化学排除”：一个 entry 可以有多个配体实体，另有蛋白实体数、核酸、元素、重原子数等先验条件。300 个缓存 entry 可证实实际元数据抓取上限；完整流水线当前没有为每个返回空列表的 entry 保存详细排除原因，这是审计缺口。

| 初始 135 个配体实体的首次排除原因 | 数量 |
|---|---|
| PDBBind_index_constitution_overlap | 35 |
| 30pct_cluster_overlap | 9 |
| strict_eligible_before_dedup_cap | 38 |
| excluded_cofactor | 32 |
| direct_sequence_overlap | 19 |
| GEOM_constitution_overlap | 2 |

此表采用首次排除的互斥顺序计数；全部可重叠排除标记见 `inventory_exclusion_reasons.csv`。57 个实体进入直接序列审计，19 个有排除集序列命中，38 个留下；按受体30%簇与二维化学贪心去重后为19个，再受15个身份 hash 截断。这个截断和随后的输入准备损失都发生在生成质量评测前。

| 输入 | 角色 | 状态 | 原因 |
|---|---|---|---|
| 9QJF | development | failed | ValueError('Incomplete observed ligand; missing heavy atoms C13,C14,C15,C17,C19,C20,C9,N1,N18,N2,O1,O16,O2,O3') |
| 9SOD | development | failed | ValueError('Incomplete observed ligand; missing heavy atoms C15,C16,C17,C18,C21,C22,C24,O23,O25') |
| 9TAR | test | failed | ValueError('Covalent ligand-receptor system outside frozen noncovalent pipeline: [{"conn_type_id": "covale", "ptnr1_label_asym_id": "A", "ptnr2_label_asym_id": "D", "ptnr1_label_comp_id": "CYS", "ptnr2_label_comp_id": "A1JUW"}]') |
| 34JY | test | failed | ValueError('Incomplete observed ligand; missing heavy atoms C8,C9,CL1,O1') |
| 39MR | test | scope_excluded | Free amino acid rather than drug-like bound ligand |
| 9YEB | test | scope_excluded | Metal-dependent cofactor requires chemistry outside this noncovalent pipeline |
| 9RMV | test | scope_excluded | Polyethylene glycol crystallization additive |
| 28LK | test | scope_excluded | Membrane lipid rather than drug-like bound ligand |

15 个计划输入中，4 个被归入源化学范围排除，3 个缺失配体重原子，1 个共价体系不支持；剩 7 个完整受体（开发1、测试6）。完整/不完整配体副本与单一 alt conformer 的检查另保存于 `ligand_copy_completeness.csv`；三个重原子缺失案例中，可找到完整单一构象副本的输入为 无。该检查未补原子、未改变当前输入，也未将失败者补回测试。

实际缓存查询从 **2025-04-01** 开始、分辨率≤1.8 Å；旧选择描述称“post-March-2026”，并非缓存查询本身的条件。已完成7个受体的实际公开日期均在2026-03之后，所以这不推翻现有成员/序列审计，但暴露了文字描述与查询 provenance 的不一致。下一轮必须从机器读取的查询生成队列说明。5438 个匹配只下载了第一页2000条；后续又只抓917个初筛合格entry中的300个。因此当前小队列不能解释为已穷尽可用独立结构。

## 扩大评测：先锁输入，再固定角色

建议目标为 **20 个开发受体＋60 个新测试受体**，每受体3 seeds、统一K8实际运行，再对K1/4作嵌套核算。旧9YSP及本轮6个测试受体只作历史审计/回归，不再当作未见测试集；不能把本次已观察的结果用于下一轮测试集选参。

1. 完整分页获取预定时间窗内的全部搜索结果，缓存查询、分页起止、返回总数、抓取失败及每个 entry/配体的每条排除原因。先沿用≤1.8 Å及当前非共价范围；若不足，只能在新测试角色锁定前预先定义额外分辨率层或扩大时间窗，单独报告，不逐个挑选高质量结果。
2. 排除三个冻结模型涉及的完整训练、验证、选模清单：PDBBind全成员、官方PB/Astex包、GEOM所有split及已索引PDBBind相同配体构造。记录清单完整性和来源；不能只排本地训练子集或依赖发布日期。
3. 对全部蛋白链做30%序列簇排除及相同直接序列命中规则；开发、测试之间及新旧队列之间按簇隔离。增加局部口袋结构相似性审计，算法、覆盖度/相似度门槛先在开发侧固定；未完成结构比对者不能标作“完全无同源口袋”。不能按“无序列簇标签”自动判断独立。
4. 配体按完整单一构象/副本选择，保留原B因子、占有率、altloc、链和残基映射；不混合两个alt构象来拼齐分子。先计算完整性、可编码性、OOV槽、片段模板误差、固定水和辅因子情况，再分配新角色。缺失原子修补若另立流程，作为预先定义的输入准备对照，不能混进当前“未修补观察坐标”主队列。
5. 按受体簇分层后采用冻结身份hash分配20/60；禁止按最终HQ、水失败数或模板拟合好坏调换角色。模板误差分布仅作分层与披露，不据测试结果决定排除阈值。锁定病例、权重、配置和源码指纹之后才运行生成。

## 化学类型与统计设计

60个测试受体建议预设40个源跨片段非环单键、20个源跨片段双键/EZ的挑战层，并覆盖词表内/OOV、酰胺、多片段和不同水数量。固定模型能否实际提出非单键/EZ候选应独立计数；源类型不能替代实际生成覆盖。若找不到足量严格独立输入，报告不足并保留分层分母，不用单键样本填满后仍宣称E/Z得到验证。

另建立至少10个源跨片段成环的输入覆盖清单。当前片段树搜索不支持此类拓扑，明确记为不支持；它们属于已规划输入范围的覆盖损失，不能声称已有成功生成。只有在另一个预先定义的方法版本支持成环后，才进行该方法与当前冻结基线的单独比较。

真实片段身份、N、显式OOV化学和已知配体邻域口袋继续作为Oracle条件，清楚记录；解除这些条件应另立分阶段实验。本方案不宣称自由de novo。

主终点仍为原标准PB与原HQ、受体内二维构造产量和有产出的pocket/seed比例；去水敏感性、标准之外的辅助HQ、原始网络与后处理质量分别报告。以受体/同源簇为重采样单位，保留所有seeds；分层报告，预先固定分析而非从多个终点挑赢家。

60个新测试×3 seeds×8候选=1440候选，姿态NFE为60×3×660=118,800；20开发为39,600，合计158,400，连接头240次前向另计。不将反事实K1/4重复计费。不按本轮共享GPU的墙钟简单外推总时间；先在新开发集测模型加载、搜索、采样、修复、QC及调度墙钟，给出独占/共享设备和失败重试预算。以上只是计划预算，尚未消耗。

## 下一轮最值得做的一个配对对照

**在保存的同一批原始候选上，对比“当前末端修复”与“额外考虑固定水排斥的末端修复”，评价始终使用完全相同的标准PB dock和原HQ。** 保持神经模型、候选图、片段内部刚性、已有主阈值以及固定水坐标相同，只在实验分支中加入水接触项；当前末端修复保留为冻结对照，不能把实验分支结果覆盖回现有基线。

优先使用全部保留水作条件，避免依据已看到的冲突删除“碍事的水”。具体权重与停止条件先在新的独立开发输入固定，再在已保存候选做描述性配对诊断，最后用新锁定测试队列确认。六个晶体参考也沿用同条件作输入相容性控制。

同时记录标准HQ、PB逐项、跨片段几何、位姿漂移、立体稳定性、CPU耗时和水相关失败改善/恶化，按受体成组比较。该对照不需要重复神经采样，额外姿态NFE可以为0；仍会产生末端优化和QC成本。若只能通过大幅漂移或破坏内部几何来避免水，便不支持采用该约束；若严格HQ改善且其他门槛和漂移可控，才支持把水条件纳入后续模型/修复研究。

选择该项的依据是96/144个候选仅被水阻断，而六个晶体参考均相容；它能在原评测规则下直接检验“水未进入修复目标”这一具体解释。它不能证明水置换的真实热力学，也不是将去水辅助质量改名为HQ。本轮只提出此对照，没有启动新优化或大型生成实验。

## 产物与复现

原主质量阈值和 PB 配置原样复制至 `configs/`，原输入及轨迹仅作读取；所有 CSV 与源文件 SHA-256 保存在本目录。`diagnostic_config.json` 明确辅助指标和水注释规则；它们仅用于归因，没有用作新筛选器或测试集调参。

```bash
bash independent_diagnostics_20261104/reproduce.sh
```

该命令重建离线分析、四份报告和核验。六个晶体对照有带输入/配置指纹的 CPU QC 缓存；不运行神经采样。本轮新神经 NFE=0。核验见 `verification.json`，导航与综合判断见 `README.md`。
