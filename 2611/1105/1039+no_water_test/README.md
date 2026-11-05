# 无结晶水协议与生成质量：精简归档

先读 [完整评估报告](NO_WATER_EVALUATION_REPORT.md)、[冻结协议](NO_WATER_PROTOCOL.md)、[自主生成能力审计](CAPABILITY_AUDIT.md)。

| 主队列阶段 | 蛋白碰撞通过 | HQ_no_water |
| --- | --- | --- |
| 神经网络终点 | 77/144 | 24/144 |
| 刚性修复后 | 139/144 | 103/144 |
| 局部柔性修复后 | 139/144 | 106/144 |

主队列仍使用真实片段身份与片段数，大部分HQ提升来自修复。独立无真实身份诊断的18次尝试中，90个槽位均预测为苯，仅产生一种分子构造；不能把修复后15/18通过HQ解释为可靠的自主生成。历史模型口袋确已去水，但模型尚未感知保留在环境中的非水辅因子。

本次重新计算648条完整QC，23/23核验通过。完整本地文件仍保留；当前GitHub仅收录必要内容：

- 报告、原始prompt、[复现说明](REPRODUCE.md)及[source/](source/)。
- [冻结协议配置](protocol.json)、[质量阈值](quality_config.json)、PB原始/无水配置。
- [主队列逐候选CSV](per_candidate_results.csv)、[受体分层](by_receptor.csv)、[配对转移](paired_transitions.csv)、[汇总](summary.json)。
- [口袋审计](pocket_audit.csv)、[144个实际模型图审计](actual_native_graph_audit.json)、[上游失败与排除](upstream_attrition.json)。
- [无身份诊断计划](free_probe/plan.json)、[18次尝试摘要](free_probe/attempt_summary.csv)、[身份坍缩诊断](free_probe/capability_diagnostic.json)、[诊断汇总](free_probe/summary.json)和关键CSV。
- [核验结果](verification.json)、[最终完整性核验](final_integrity.json)、[PB不可评估记录](PB_unavailable_records.json)、[阶段图](stage_quality.png)。
- [MANIFEST.json](MANIFEST.json)：本次上传文件大小、SHA-256、来源及仅为GitHub阅读调整的文档链接。

不收录完整SDF、模型权重、扩散缓存、逐分子QC/日志、逐键/逐角/逐原子大表及重复JSON。原始完整实验完全不改动。文档中完整本地文件已标注；仅克隆本精简仓库不足以独立重跑全部实验。历史任务保持原样。
