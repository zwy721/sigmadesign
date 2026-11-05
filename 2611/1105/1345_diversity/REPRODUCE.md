# 复现

从项目根目录执行，使用已有环境，不安装/下载新包。所有GPU进程均限定CUDA_VISIBLE_DEVICES=0。

    export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0
    export PYTHONPATH="$PWD/.quality_deps:$PWD:$PWD/../third_party/sigmadock_official_repo_20261030/src"
    PY=../.venv-cutoff-train/bin/python
    TASK=2611/1105/1345_diversity
    $PY -B "$TASK/source/run.py" --phase pilot
    for shard in 0 1 2 3; do
      $PY -B "$TASK/source/run.py" --phase main --world 4 --shard "$shard" > "$TASK/main_$shard.log" 2>&1 &
    done
    wait
    $PY -B "$TASK/source/analyze.py" --merge-only
    $PY -B "$TASK/source/evaluate.py" --manifest "$TASK/main/manifest.json" --output "$TASK/main/quality" --workers 8
    $PY -B "$TASK/source/analyze.py"
    $PY -B "$TASK/source/reports.py"
    $PY -B "$TASK/source/verify.py"

实际运行增加了extra_workers.py，分别提前计算9SYL/9T4A/9T51三口袋，每口袋4个CPU分片，也只用GPU0。初始四个worker到达这些口袋时读取已完成result.json，不重新抽样；此调度不影响冻结504次分母。为简洁且避免调度重叠，以上复现命令仅用四个分片即可得到同一结果，不需要额外workers。stream_quality.py在正式生成期间并发CPU QC；最后evaluate.py用相同缓存和冻结协议汇合。

run.py遇已有result.json原样恢复，包括失败记录，不重试失败。独立重跑应在新的同深度项目副本中执行，并保留原始计划；不能在同一目录修改策略/权重/种子后续写。当前目录已完成结果的再次运行只读恢复，不构成独立重放。pilot / main各保留全部分布、搜索池、拒绝事件、图审计、SDF/spec、修复和QC。

原历史基线及144例条件生成结果保持原样，历史目录全部文件SHA核对。权重、协议、来源和生成代码指纹见frozen_plan.json、generator_provenance.json、historical_fingerprints.json及verification.json。完整训练并集Novelty和未见口袋泛化不在本次可评估范围。

本次为避免慢项阻塞，另以boost_quality.py在独立quality_boost目录使用32个CPU worker处理待检结构，参数与原8-worker队列完全相同。finalize_quality.py按预定“原队列优先”合并结果，不按质量挑选；--snapshot只保存当前已完成结果，待检项仍明确列在PENDING_QUALITY_RECORDS.csv。慢项结束后依次运行finalize_quality.py、analyze.py、reports.py、verify.py即可刷新完整辅助诊断。正式1512项终态QC不依赖剩余初始阶段慢项。没有因为慢项而改阈值、少算50构象或设未经批准的超时。

