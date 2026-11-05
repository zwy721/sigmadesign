> GitHub 精简归档：保留报告、源码、固定参数、关键汇总与核验。完整结构、QC缓存、逐原子明细及完整指纹留在原本地工作区；复现需要原项目依赖。本文仅调整归档链接与此范围说明，实验结论和数值不变。

# 复现与文件说明

从项目根目录运行，使用已有环境，不下载或安装依赖：

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PWD/.quality_deps:$PWD:$PWD/../third_party/sigmadock_official_repo_20261030/src"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
PY=../.venv-cutoff-train/bin/python
TASK=2611/1105/1039+no_water_test
$PY -B "$TASK/source/prepare.py"
$PY -B "$TASK/source/evaluate.py" --workers 8
```

`prepare.py` 审计并复用已经生成的144个候选及两种既有无水修复；不会重新运行历史神经网络。完整依赖路径和SHA在 input_fingerprints.json 中。6个去水口袋、解码起点SDF和spec保存在本任务目录；其余复用结构通过manifest绝对路径追溯。保持这些历史依赖，才可完整复现。

18次无真实身份诊断需要两张可用GPU。其既有计划不可修改参数后继续沿用同一输出目录；如需独立重跑，先在新的项目副本执行。首次运行：

```bash
$PY -B "$TASK/source/free_probe.py"
$PY -B "$TASK/source/evaluate.py" \
  --manifest "$TASK/free_probe/manifest.json" \
  --output "$TASK/free_probe/quality" --workers 8
$PY -B "$TASK/source/free_audit.py"
$PY -B "$TASK/source/analyze.py"
$PY -B "$TASK/source/verify.py"
$PY -B "$TASK/source/build_reports.py"
```

`free_probe.py` 采用最小新适配器，不是原历史native入口。只验证已有GEOM身份头、当前连接头、官方位姿网络的一种固定诊断策略，不代表全部自由生成方案。启动阶段曾修正独立Python进程的同名common模块导入冲突；修正发生在首个候选尝试前，没有据候选结果修改生成策略。

`PB_unavailable_records.json` 保存标准PB未返回有限值的记录。本次唯一一条为解码起点的internal_energy，不跳过、不补成通过。`verify.py` 在两个含辅因子受体上重放两种优化，检查去水前后的输出逐坐标精确一致。

最终报告为 NO_WATER_EVALUATION_REPORT.md；数值表包括主队列与free_probe下的summary、逐候选、键、角及严重蛋白碰撞对。`quality/molecule_metrics.jsonl` 是详细QC原始结果；cache/details为可再生运行产物。原始模型权重、原始来源和旧报告未更改。

未来上传GitHub宜只保留报告、源码、配置、审计与关键汇总CSV/JSON、图和指纹；完整SDF、模型缓存和逐原子大表留在本地。当前任务未要求上传，未执行GitHub写入。
