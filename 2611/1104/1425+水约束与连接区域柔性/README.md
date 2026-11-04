# 固定水与连接区域有限柔性：保存结构的离线修复对照

本目录完成 `prompt.md` 要求的 144 个保存候选的 CPU 配对实验。新增神经 NFE=0，未新增候选、训练或结构检索。所有历史文件只读，新产物集中于本目录。

先读 [下一步决策](NEXT_DECISION_SUMMARY.md)、[水软排斥结果](WATER_AWARE_REPAIR_RESULTS.md)、[连接区域有限柔性结果](LOCAL_PORT_FLEXIBILITY_RESULTS.md)。严格 HQ 与仅忽略两项水检查的辅助 HQ 始终分开；任何修复后通过率都不代表模型原始输出通过率。

- `fixed_parameters.json`、`quality_config.json`、`posebusters_effective_config.json`：预锁定方法参数及原完整 QC。
- `input_audit.csv`、`baseline_audit.json`、`input_manifest.json`：输入核对及 7/144、103/144、96/144 基线复核。
- `fixed_waters.json`：全部 79 个原水氧的坐标、原编号和半径。
- `structures/original/`：144 份字节一致的原始 SDF 和连接 sidecar。
- `structures/A0_B0/`、`structures/A1/`、`structures/B1/`：各 144 份修复 SDF、同流程立体赋值 SDF、sidecar、立体决定、优化退出与耗时。A0/B0 共用一份冻结控制。
- `per_candidate_results.csv`：全部候选、各 QC、位移、收敛和 CPU 时间；`failure_cases.csv` 保留质量失败和未收敛样本，`regressions.csv` 列出相对控制的质量退化。
- `bond_details.csv`、`angle_details.csv`、`torsion_details.csv`、`clash_details.csv`、`stereo_centers.csv`、`atom_displacements.csv`、`water_contacts.csv`、`water_pair_details.csv`：逐项明细，含跨片段/片段内归属或固定水编号。
- `summary.csv`、`by_receptor.csv`、`by_seed.csv`、`paired_changes.csv`：按受体及原种子分组的配对统计。6 个受体的相关候选不作为 144 个独立来源。
- `quality_A/`、`quality_B/`：原完整 QC 的逐候选 JSONL、缓存、评测日志与诊断结构；此处 QC 输出的优化诊断结构不作为实验候选。
- `processing_A.json`、`processing_B.json`、`timing_A.json`、`timing_B.json`：处理状态、CPU 与墙钟时间；优化未收敛与质量失败分开统计。
- `raw_geometry_stereo.json`：相同原始坐标的独立几何立体决定，用于逐中心比较。
- `input_fingerprints.json`、`execution_lock.json`、`verification.json`、`regression_tests.log`：不可变性、原控制复现、同 QC/立体流程及数值梯度回归证据。
- `source/`、`tests/`、`reproduce.sh`、[运行命令](RUN_COMMANDS.md)：代码、测试及完整 CPU 复现。

B1 为原片段 SE(3) 加连接端点及第一邻居的额外局部坐标；预锁定位移上限是每分量 ±0.10/√3 Å，保证范数≤0.10 Å，另加位置正则。A1 仅加水项，保持片段刚性。两种单因素分别比较，未跑可选 B2 或 A+B 联合方法。现有开发支持不足，仅做单参数描述性诊断，无测试集调参或泛化声明。
