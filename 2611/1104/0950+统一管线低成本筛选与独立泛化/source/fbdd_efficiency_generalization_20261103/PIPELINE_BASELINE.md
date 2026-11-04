# 任务一最终报告：统一管线与回归基线

报告日期：2026-11-04。范围为任务一；全部统计来自已完成的运行记录。本次报告整理没有重新采样或训练模型。

统一入口 `pipeline.py` 已完成历史输入回归及新独立输入的实际运行。相同 CASF Top-8 输入得到 **70/288（24.31%）综合高质量候选、41 个按受体去重的二维构造、19/36 个有产出的 pocket/seed**；标准 PoseBusters dock 为 71/288。288 个最终 SDF 的序列化坐标与历史结果完全相同，规范化立体图及质量配置一致。

## 固定流程与输入条件

`t=0.6 无连接前缀 → 多候选合法连接 → 官方完整续采样 → 刚性片段受约束修复 → 可靠三维立体赋值 → 冻结质量评价 → 受体内二维构造去重`。

旧前缀执行 20 步；官方续采样执行 16 步至 t=0.4，再执行 64 步至最终时间。连接候选采用冻结连接头和原有层级合法搜索。K=1/4/8 取同一 Top-8 排序的嵌套前缀；不补复制候选。每个 pocket/seed 的姿态预算为 `20 + 80K` NFE，连接头前向另计。搜索仍采用当前完整候选搜索，不能把减少 K 自动等同于减少搜索费用。

真实片段身份、片段数 N、OOV 片段的显式化学结构，以及已知配体邻域定义的口袋仍是受控输入；不输入真实跨片段连接作为候选答案。这是 Oracle 条件下的 FBDD 评价，尚不能代表自由 de novo 或未知口袋发现。

`generate.py` 用独立解释器分别调用 `legacy_prefix.py` 和 `native_generate.py`，隔离旧模型与官方代码的同名模块；实际导入来源、模型指纹、随机种子和原生断点保存在各运行的 protocol 中。新开发集已完成 24 个候选，独立测试已完成 144 个候选，均沿用本入口和冻结配置。

## 输出版本与质量定义

每个候选保留 `raw.sdf`、`repaired.sdf`、`derived.sdf`，以及原始/派生 spec、修复诊断和立体可靠性决定。raw 通过来源链接保留；既有历史轨迹没有被覆盖。另保存 `raw_stereo` 对照：仅给未修复的原始坐标做相同可靠三维立体赋值，供公平评价修复效果。

HQ 要求化学有效、完全连通、无严重分子内碰撞、键长/键角/环/扭转几何合格、无严重蛋白碰撞、位于目标口袋、标准 PB dock 全检查通过、立体定义完整且与几何一致。每项必须明确为 True；缺失或失败均不算通过。不以参考构造、参考 RMSD、QED、Vina 或应变阈值作为主质量门槛。PB 保留水相关检查；优化器不含水，二者口径不同且均未更改。

不同构造按 `(receptor_id, constitution_smiles)` 去重，跨三个 seed 合并；同构造不同立体异构体或不同构象不重复增加二维构造计数。不同受体的相同构造分别计数，因此 41 是受体内去重后求和，不是全局去重。

## 冻结修复、赋值与指纹

修复沿用 `a_priori_v1`：MMFF94s，位姿/朝向系数 5/2，蛋白排斥系数 200、半径比例 0.8，酰胺取最近 cis/trans ±5°，L-BFGS-B 最多 800 次，氢优化 200 次。各片段只允许刚体移动，片段内重原子距离误差要求 <1e-8 Å。每轴平移范围 ±2 Å、每轴旋转范围 ±π/2；重原子移动超过 2 Å 是软惩罚，不能解释为硬位移上限。

立体赋值要求归一化四面体体积 ≥0.10，在 16 次 σ=0.01 Å 扰动下 CIP 稳定；双键偏平面 ≤15°。不读取参考配体手性。不可靠中心保留未定义，由相同质量检查判定。

| 冻结模型 | SHA-256 |
|---|---|
| best_loss.pt（/data/zhouwanyao/cutoff_geom_pdbbind_joint_200ep_20261026_cuda156） | `4c45e8b310f5fa4cf3dceee737af635ac9723cf2e5795e6c31a09244a47a3ac7` |
| official_sample_checkpoint_0.ckpt（/home/zhouwanyao/drive_download/releases/pretrained/sigmadock_official） | `db15427ca349e6f1e5f894bff841112c7360384886aa472667d8011307cad382` |
| best.pt（/home/zhouwanyao/drive_download/releases/oracle_rebuild_20261031/results/connection_robust_20261101/training_v2） | `b79dd58d80174d1b90105c1d7c4e4f9d30f71962033c8787fd4c9be908cedf93` |

质量配置 SHA-256：`26311524167940621690c83e3c5521dadaf71fab4276989b61689dbab7d6b4ad`。完整参数及脚本指纹见 [baseline_replay/config.json](baseline_replay/config.json)、[cohort_lock.json](cohort_lock.json)；最终交付文件与源码指纹见 [tasks123_verification.json](tasks123_verification.json)。修复和阈值在独立测试前已冻结，本阶段没有新网络训练。

## CASF K=1/4/8 回归

| K | 版本 | 候选 | HQ | HQ 率 | 不同 HQ 构造 | PB dock | 有 HQ 尝试/36 |
|---|---|---|---|---|---|---|---|
| 1 | raw：网络显式图 | 36 | 0 | 0.0% | 0 | 5 | 0 |
| 1 | raw_stereo：仅赋值 | 36 | 3 | 8.3% | 2 | 5 | 3 |
| 1 | repaired：仅修复 | 36 | 3 | 8.3% | 2 | 11 | 3 |
| 1 | stereo：修复＋赋值 | 36 | 10 | 27.8% | 7 | 11 | 10 |
| 4 | raw：网络显式图 | 144 | 1 | 0.7% | 1 | 12 | 1 |
| 4 | raw_stereo：仅赋值 | 144 | 6 | 4.2% | 4 | 12 | 6 |
| 4 | repaired：仅修复 | 144 | 7 | 4.9% | 4 | 38 | 6 |
| 4 | stereo：修复＋赋值 | 144 | 37 | 25.7% | 21 | 38 | 18 |
| 8 | raw：网络显式图 | 288 | 2 | 0.7% | 2 | 25 | 1 |
| 8 | raw_stereo：仅赋值 | 288 | 14 | 4.9% | 10 | 25 | 9 |
| 8 | repaired：仅修复 | 288 | 11 | 3.8% | 7 | 71 | 7 |
| 8 | stereo：修复＋赋值 | 288 | 70 | 24.3% | 41 | 71 | 19 |

网络直接输出 Top-8 HQ 为 **2/288**；加立体赋值而不修复为 **14/288**；最终修复加赋值为 **70/288**。隔离几何修复应比较 raw_stereo 与 stereo（14→70），不能把全部增益都归给网络，也不能将“仅修复”与“修复加赋值”的差异算作几何改善。K1/4 表取历史同池排序前缀；新入口对完整 K8 进行了逐候选精确回归。

回归仅在坐标、化学图和质量配置完全一致后复用历史 QC；没有重复运行神经轨迹，也没有声称重新独立计算了全部 QC。证据见 [baseline_regression.json](baseline_replay/baseline_regression.json)、[pipeline_results.json](baseline_replay/pipeline_results.json)、[quality_manifest.json](baseline_replay/quality_manifest.json)。

## 成本与复现

本次 CASF 回放新增神经 NFE=0；完整 Top-8 的逻辑预算为 23,760 姿态 NFE。重新修复、赋值及 I/O 实测墙钟为 28.405 秒，因神经轨迹及 QC 被复用，不能当作完整推理延迟。历史缺少搜索/模型加载计时；阶段耗时与 QC CPU 累计秒数分开报告。

新独立 K8 从输入至三版本质量筛选的实际研究墙钟为 2,030.324 秒，另补 raw_stereo 诊断 137.938 秒，均见任务三报告。K1/4 由同一已保存 K8 轨迹作成本核算，没有单独部署墙钟。任务一至三合计新增姿态 NFE 为 13,860（开发 1,980＋测试 11,880）；CASF 回放为 0。任务四费用另列于原完整账本，不混入本次范围。

只重新生成本次报告、图和核验记录，不运行模型：

```bash
bash fbdd_efficiency_generalization_20261103/reproduce_tasks123_reports.sh
```

新输入接口为 `pipeline.py --K 1/4/8 --cases cases.json --output NEW_DIRECTORY`；已有原生结果可用 `--raw-manifest` 接入。历史精确回归使用 `--replay-casf`。完整实验复现入口 `run_study.sh` 会运行包括任务四在内的实验，不能作为仅整理报告的命令。本次交付以只读结果重建命令为准。

任务一结论：统一生成、后处理与评价入口已形成，并通过完整历史回归和新输入实跑。该结果固定了后续比较基础；独立性、化学覆盖和生成质量的外推限制见任务三。
