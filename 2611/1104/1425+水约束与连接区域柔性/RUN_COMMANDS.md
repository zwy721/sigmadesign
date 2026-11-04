> GitHub 精简归档说明：以下命令用于原始本地工作区，依赖原项目模块、历史输入和既有 Python 环境。完整结构、QC 明细和输入指纹保留在本地；仅克隆本仓库不足以独立重跑实验。

# 执行与复现命令

工作目录为项目根目录 `/home/zhouwanyao/drive_download/releases/casf_paired_eval`。Python 使用现有 `../.venv-cutoff-train/bin/python`，依赖使用项目已有目录；无需安装或网络检索。

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PWD/.quality_deps:$PWD/../generation_experiment/vina_dependencies:$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
study_python="$PWD/../.venv-cutoff-train/bin/python"
study_task="$PWD/2611/1104/1425+水约束与连接区域柔性"
```

本轮实际顺序：

```bash
"$study_python" -B "$study_task/source/audit.py"
"$study_python" -B -m pytest -q -p no:cacheprovider "$study_task/tests/test_repair.py"
"$study_python" -B "$study_task/source/run_experiments.py" --experiment A --workers 8
"$study_python" -B "$study_task/source/run_experiments.py" --experiment B --workers 8
"$study_python" -B "$study_task/source/analyze.py"
"$study_python" -B "$study_task/source/verify.py"
"$study_python" -B "$study_task/source/build_reports.py"
```

参数在 audit 阶段冻结，并在首次 A 优化启动后、读取新 QC 结果前记录到 `execution_lock.json`。固定参数和数值源文件的哈希均进行最终核验。优化前先运行新增的解析梯度/控制回归测试；完成后复跑新增测试，并运行原刚性优化、旋转及立体赋值相关回归。每个阶段日志单独保存。

完整重新运行请指定一个新输出目录，脚本会拒绝覆盖已有冻结实验：

```bash
bash "$study_task/reproduce.sh" /tmp/water_port_reproduction
```

不传参数时写入任务目录下带 UTC 时间的新 `replications/` 子目录。脚本依次完成审计、冻结、测试、A、B、统计、全队列核验及报告生成。数据和固定环境继续从已记录路径读取；不存在的输入直接报出路径，不推测补齐。没有神经网络前向、训练、GPU 作业、外部结构检索或新增候选。

单独重新生成统计和检查时，使用相同 Python、PYTHONPATH，运行 `analyze.py`、`verify.py`、`build_reports.py`。这只重写本轮派生表和报告，不改历史输入或 SDF；若要保留本轮每个派生产物的原样快照，请使用新输出目录完整复现。
