# SigmaDock 几何约束与 Top-K 评测

本目录为独立实验，不覆盖原始 SDF、连接图、模型权重和状态。任务原文为 `prompt_used.md`。

使用父目录 `.venv-cutoff-train/bin/python`（工作区上一层），依赖沿用 `.quality_deps` 与 `../generation_experiment/vina_dependencies`。完整复现入口：

```bash
bash geometry_constraints_20261103/run.sh
```

脚本默认使用 GPU 0–3 分片；GPU 分配可按空闲资源改动。所有模型冻结，复用原候选池、受体、刚性片段划分、t=0.6 状态和 t=0.4 原生检查点。存在结果时直接复用。恢复新 raw 分支严格读取原 C/Q/pos/pos0/T0/R0 与 RNG；中途干预读取同一个检查点 RNG，在修复后的物理坐标重建官方刚体状态，并记录物理转移误差。

`fixed_protocol.json` 在评测前冻结：单组先验参数，无 CASF 超参数搜索。`quality_config.json` 保留旧工具所有质量阈值，只关闭不参与质量门槛的可选 Vina 诊断，examples 数设为 0；完整标准 PB dock 和 50 构象能量检查均保留。

输入与输出：

- `pool_completion/`：184 个新增完整后续分支；原来 104 个只读复用。
- `short/`：288 个短程 SDF，显式图由原 spec 保留。
- `post_repair/`：288 个完整候选的最终修复。
- `mid_repair/`、`mid_continuation/`：36 个单候选 t=0.4 修复及其官方后续扩散。
- `true_graph_post/`：36 个已知真图对照的最终修复。
- `stereo/`：独立几何立体图、spec 和逐中心可靠性决定。
- `short_quality/`：288 条完整标准质量记录。
- `final_quality/`：1260 条记录（raw/post 各 288 及各自 stereo 派生版、mid 36 及派生版、真图修复 36）。
- `metric_audit.csv`、`bond_trajectories.csv`、`angle_final.csv`、`clash_final.csv`、`amide_final.csv`：逐分子、逐键、逐阶段归因。
- `paired_quality_displacement.csv`、`repair_summary.csv`、`paired_cluster_bootstrap.json`：配对修复和受体分组统计。
- `topK_pool_quality.csv`、`selection_policies.csv`、`topK_pocket_yields.csv`：产量、策略和计算成本。
- `short_predictivity.csv`、`short_to_final.csv`：全部短程候选对最终质量的预测能力。
- `verification.json`：输入、模型哈希、原始质量不变、阈值不变、图与坐标不变、记录完整性核验。

根目录交付报告：`METRIC_AUDIT.md`、`GEOMETRY_ERROR_ANALYSIS.md`、`GEOMETRY_REPAIR_RESULTS.md`、`MULTICANDIDATE_DENOVO_QUALITY.md`。原先严格 RMSD 统计已纠正，旧版保存在 `pre_audit_quality/`，84 条原主质量结果不变。

NFE 为姿态主干前向，成功新轨迹 11,776+2,304=14,080 NFE。此次共享 GPU 调度迁移有 8 条中断分支，未记录的额外开销界限为 0–512 NFE；成功轨迹耗时不含这些损耗或未单独计时的连接搜索、立体赋值。部署预算与本次研究补全成本分别报告，不将批处理并行时间之和声称为端到端延迟。
