# 复现与进度检查

以下命令在 `/home/zhouwanyao/drive_download/releases/casf_paired_eval` 执行。解释器为 `../.venv-cutoff-train/bin/python`，版本及依赖审计见 `audit_initial.json`。

```bash
# 查看真实阶段与日志；process.json 记录PID、启动命令、日志路径。
cat geom_connection_geometry_training/status.json
cat geom_connection_geometry_training/process.json
tail -n 30 geom_connection_geometry_training/training.log
nvidia-smi --id=0
```

现有数据已准备完毕且全部分子图核验通过。以下是从准备到正式训练的顺序；已有正式运行时不要重复启动。数据、原权重及依赖哈希保存在 `data/complete.json`、`dataset_verification.json`、`audit_initial.json`，正式源码和配置哈希保存在 `formal_start.json`。

```bash
export CUDA_VISIBLE_DEVICES=0
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=2
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/audit.py
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/prepare.py --smoke --workers 8
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/prepare.py --workers 8
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/verify_dataset.py
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/preflight.py
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/verify_prelaunch.py
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/verify_short_workflow.py
# 在新的实验副本中启动；原目录已存在formal_start.json时会拒绝重启。
# launch.py后台脱离终端运行，并记录process.json、training.log和PID。
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/launch.py
# 可单独刷新已完成部分的报告，不会补齐未完成结果。
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/reports.py
```

随机种子为20261106、20261107、20261108。A1/A2及各自二维消融各40epoch；A3/A4真实连接几何各40epoch。每epoch训练分子均为2000个，按epoch轮换可用构象；每5epoch在全部256验证分子、固定三个噪声层选checkpoint。阶段C用固定128训练/32验证分子的预测合法图适配10epoch；阶段D只在A/B/C验证门槛通过后联合微调5epoch。阶段门槛及选择结果保存在 `stage_decision.json`。

评价噪声由种子、构象ID及时间共同固定，与batch大小、数据子集和处理次序无关；各组共享同一坐标起点，MMFF与相应神经评价样本成对。完整测试采用全部可用构象与7个噪声层；离散搜索固定每集合32分子、首个构象、3噪声层；MMFF固定16分子、首个构象、3噪声层，最多200迭代并保留不收敛记录。

相对平移和旋转误差是对具体GEOM参考构象的诊断，不能将可旋转单键的参考扭转偏差视为化学失败。多构象RMSD/覆盖不做原子对称校正。GPU上的index_add等可能产生微小浮点差异；软件不宣称全部训练位级确定性。

## 2026-11-06 报告故障后的续跑

第一个种子全部训练/评估完成后，搜索报告合并重复传入split字段，导致13:28退出。修复仅涉及reports.py、run_experiment.py和launch.py，训练目标、数据、模型、噪声与超参数不变。原始代码、错误状态、日志和首个种子103个文件的哈希保存在 `recovery/report_split_failure_20261106/`，`recovery_verification.json`记录实际报告回归、结果不变和报告错误隔离检查。

已加入种子边界续跑；完成种子的全部训练历史、评价分母与文件必须验证通过，未完成种子目录必须尚不存在，才能继续。`formal_start.json`不变，24小时截止仍为2026-11-07 12:07:17（北京时间），包含此前运行与中断时间。续跑源码指纹与进程分别记录在 `resume_events/`、`process_history/` 和当前 `process.json`。

```bash
# 只检查可否恢复，不启动GPU任务
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/run_experiment.py --check-resume
# 原进程退出后，在已验证完成的种子边界后台恢复；禁止重复启动活跃任务
../.venv-cutoff-train/bin/python -B geom_connection_geometry_training/source/launch.py --resume
```

报告故障单独记录在report_failures/，不再使后续种子训练退出；最终报告失败会明确标为report_failed，不会声称整轮完整交付。当前实现只支持已完成种子之间的恢复，部分种子的checkpoint恢复需另行处理。
