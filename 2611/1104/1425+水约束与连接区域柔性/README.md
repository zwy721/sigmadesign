# 固定水与连接区域有限柔性：精简结果归档

本目录保留本轮实验的报告、代码、固定参数、关键结果表和验证摘要。144 个保存候选的 CPU 配对实验已经完成，新增神经 NFE=0；本地完整实验结果保持原样。

先读 [下一步决策](NEXT_DECISION_SUMMARY.md)、[水软排斥结果](WATER_AWARE_REPAIR_RESULTS.md)、[连接区域有限柔性结果](LOCAL_PORT_FLEXIBILITY_RESULTS.md)。严格 HQ 分别为 A0/B0 7/144、A1 65/144、B1 7/144；仅忽略两项水检查的辅助 HQ 分别为 103/144、84/144、106/144。以上均为修复后结构质量。

- `prompt.md`、`PROMPT_STATUS.md`：任务原文及来源。
- `fixed_parameters.json`、`quality_config.json`、`posebusters_effective_config.json`、`fixed_waters.json`：冻结的优化参数、原完整 QC 与固定水坐标。
- `baseline_audit.json`：7/144、103/144、仅水阻断 96/144 的基线复核。
- `per_candidate_results.csv`、`failure_cases.csv`、`regressions.csv`：逐候选综合质量、失败与退化清单。
- `summary.csv`、`by_receptor.csv`、`by_seed.csv`、`paired_changes.csv`：总体及原受体/种子分组的配对统计。
- `verification.json`、`delivery_verification.json`、`execution_lock.json`、`regression_tests.log`：完整本地实验的 30 项核验、11 项单测及冻结证据。
- `timing_A.json`、`timing_B.json`、`experiment_A.log`、`experiment_B.log`：CPU/墙钟时间与任务级运行摘要。
- `source/`、`tests/`、`reproduce.sh`、[运行命令](RUN_COMMANDS.md)：数值代码、测试和原工作区复现入口。
- `MANIFEST.json`：本次精简归档的范围、文件大小和哈希。

完整 SDF、逐候选 QC JSONL、分子级缓存/日志、逐键/逐角/逐原子/逐水明细以及重复的大型 JSON 不再放在当前 GitHub 版本中，均保留于原本地实验目录。报告和代码中提到的这些路径指向完整本地结果，不表示文件仍包含在精简仓库。验证摘要中的结构数量和 QC 数量描述已完成的实验，不是当前仓库的文件数量。

复现需要原项目布局、历史输入和已安装依赖；本目录不提供独立的数据与环境副本。按普通 Git 删除提交进行精简，先前完整上传仍可在历史提交中查看。
