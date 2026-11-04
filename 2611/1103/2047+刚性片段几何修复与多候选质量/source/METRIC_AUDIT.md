# 评价指标一致性审计

官方真实连接组严格 RMSD < 2 Å 的正确结果为 **25/36**。原来的 0/36 来自立体原子映射处理不一致，既不是几何采样失败，也不是综合质量门槛附加到 RMSD 上。

| 组别 | 构造一致 | 显式立体一致 | 构造 RMSD < 2 Å | 严格 RMSD < 2 Å |
|---|---:|---:|---:|---:|
| A_predicted_t0.6 | 1/36 | 0/36 | 1/36 | 0/36 |
| B_official_true_graph | 36/36 | 36/36 | 25/36 | 25/36 |
| C_crystal | 12/12 | 12/12 | 12/12 | 12/12 |

`standardized_smiles` 经 RDKit Cleanup 清除了参考 SDF 在非立体/对称原子上的冗余手性标签；旧 RMSD 代码直接用未清理标签调用 `GetSubstructMatches(useChirality=True)`。36 个官方真图样本规范化后立体相同，但原始立体映射均为空。修复在重原子副本上统一执行 `AssignStereochemistry(cleanIt=True, force=True)` 后再映射，不从几何或参考补写生成图标签。

RMSD 定义：相同重原子构造；允许分子图对称映射；取受体坐标系下的最小重原子 RMSD；不做刚体对齐。严格版本还要求规范化显式立体一致及保持手性的映射；上限 4096 个映射。不要求 PB、能量、碰撞或综合高质量通过。严格图一致也不表示三维构象的实际手性正确：官方真图组仍有 6/36 的声明标签与最终三维不一致，继续由独立质量项判为失败。

自动化测试覆盖全部 36 个真图样本、镜像标签排除、平移 3 Å 不被对齐消除、旋转雅可比、刚性修复及近共面/E/Z 赋值。原始报告快照保留在 `geometry_constraints_20261103/pre_audit_quality/`；重新生成旧质量报告时所有主质量检查保持不变。

| 组别 | 标准 PB dock | 仅移除水检查的 PB 敏感性 | 自定义无严重内部碰撞 | PB 无内部碰撞 |
|---|---:|---:|---:|---:|
| A_predicted_t0.6 | 5/36 | 9/36 | 19/36 | 22/36 |
| B_official_true_graph | 10/36 | 16/36 | 31/36 | 31/36 |
| C_crystal | 9/12 | 11/12 | 12/12 | 12/12 |

PB 仍运行官方 dock 全部 22 项，阈值不变。水敏感性只从同一次输出中移除 `minimum_distance_to_waters` 和 `volume_overlap_with_waters` 两项，再取其余项目的合取；不是修改阈值或替代标准 PB。晶体对照也只有 9/12 标准 PB dock 通过，水敏感性为 11/12。

键长和键角均使用冻结 RDKit 距离几何界限和 25% 容差；逐分子自定义/PB 结果完全一致。自定义碰撞额外检查 1–3 对（下界的 0.70），1–4 使用距离几何下界的 0.70，远程对使用 vdW 半径和的 0.60；PB 排除 1–2 和 1–3，对剩余非键对用距离几何下界的 0.70。因此预测组有 3 个自定义失败而 PB 通过的样本，属于定义差异。逐对记录保留拓扑及距离，不把两套检查合并成一个分数。

自定义蛋白碰撞不含固定晶体水；标准 PB 分别检查蛋白、有机辅因子、无机辅因子和水。无机辅因子使用共价半径，避免把正常金属配位当成蛋白碰撞。所有标准 PB 检查仍进入综合高质量合取。

逐分子审计：[metric_audit.csv](geometry_constraints_20261103/metric_audit.csv)；指标对照：[PB_custom_comparison.csv](geometry_constraints_20261103/PB_custom_comparison.csv)；回归测试：[test_geometry_constraints.py](tests/test_geometry_constraints.py)。
