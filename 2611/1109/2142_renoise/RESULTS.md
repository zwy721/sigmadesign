# 两阶段固定图二次去噪：结果与决策

本轮已完成，7组各32例，所有分母保留。**B改善几何质量，但尚不能证明结合位姿更准确。**

A为原始动态24步；B在其终态局部加噪至0.2再去噪24步；三个C分别使用历史初始化、官方原seed先验和官方新seed先验；D在原流程首次t≤0.2时提前冻结；E从A终态加噪至1再去噪，匹配额外24 NFE。开发选择仅使用6例独立开发集，未按32例择优。

B平均内部/蛋白严重碰撞对数0.969/0.844→0.313/0.219；参考键长MAE0.179→0.062 Å、角度MAE11.71→4.49度。B平均未对齐重原子RMSD1.99 Å，接触保持率54.4%。D仅10/32例有可比较原子，38原子，平均每例覆盖率4.29%；不能将其映射RMSD与B/C/E全分子RMSD直接比较。

九项测试通过；224条终态、4192条逐步记录及2790个历史文件哈希核验通过。C_replay全部32×24步与历史逐位一致。原图/坐标推导立体口径分列；所有修正、旧源码及冻结版本均保留。两种非规范片段切分及checkpoint缺少ligand–ligand空间编码器的限制明确记录。历史压力集并非独立泛化数据；未使用MMFF修复这些采样输出。

压缩包包括十项要求交付、完整224条QC、源代码、上游依赖源、冻结配置、逐分子/逐步表格、计算开销和版本修正证据。外部权重、输入结构及完整轨迹仍是复现前提。`run/DELIVERY_VERIFICATION.json`等字段为发布前快照。

---

# Two-stage fixed-graph results

32 historical pressure cases per method; all denominators retained. HQ uses frozen coordinate-derived stereo; original-declared HQ is separate. No training or sampling repair was used. These cases are not an independent generalization set.

Original stereo definition complete: A 12/32, B 12/32, C_replay 12/32, C_official_seed 12/32, C_official_new 12/32, D 14/32, E 12/32.

|Method|HQ /32|Original-stereo HQ /32|PB /32|Internal pairs mean|Protein pairs mean|RMSD to A, Å|Logical NFE|
|---|---:|---:|---:|---:|---:|---:|---:|
|A|2|1|3|0.969|0.844|0.000|24|
|B|9|3|16|0.312|0.219|1.991|48|
|C_replay|11|3|15|0.281|0.219|3.989|48|
|C_official_seed|9|3|11|0.094|0.281|4.788|48|
|C_official_new|8|1|11|0.406|0.312|5.100|48|
|D|8|4|10|0.844|0.500|1.780|24|
|E|8|2|12|0.281|0.438|4.633|48|

D RMSD is available in 10/32 cases only, covering 38 atoms; mean per-molecule mapping coverage=4.29%. Its mapped RMSD must not be compared with full-molecule B/C/E drift.

Frozen development choice: B t_start=0.2, reverse noise_scale=0.0; D freeze threshold=0.2. C_replay is the historical initialization control; C_official_seed/new use IGSO3 official priors.

## Paired HQ changes versus A

|Method|Wins|Losses|Exact two-sided McNemar p|
|---|---:|---:|---:|
|B|8|1|0.0390625|
|C_replay|10|1|0.0117188|
|C_official_seed|8|1|0.0390625|
|C_official_new|7|1|0.0703125|
|D|6|0|0.03125|
|E|7|1|0.0703125|

P values are descriptive, unadjusted across six contrasts and correlated pockets; no independent confirmatory claim. D changes predicted chemistry, so its contrast is not an isolated denoising effect. E holds fixed chemistry and 24 extra NFE; it does not isolate freezing itself.

Geometry MAE in this study uses saved network reference geometry. Frozen HQ geometry instead uses its existing RDKit bounds/amide checks. These targets are explicitly separate. Protein collision counts use protein atoms only; cofactors remain in full PB and QC diagnostic outputs. Contacts are heavy-atom protein pairs within 4.5 Å.

Collision columns and trajectories use the original frozen predicted graph and unrounded network coordinates. HQ_derived_* columns use the serialized coordinate-derived QC graph. Stereo-sensitive RDKit bounds can differ: 9RT3_41005/C_official_new has1 original-graph severe pair versus5 derived-QC severe pairs at identical serialized coordinates. No threshold or endpoint was changed; see QUALITY_GRAPH_SCOPE_AUDIT.json. Original-stereo HQ applies the historical original-declaration gates to the frozen QC, without claiming network stereo prediction.

Full raw quality, errors, unavailable checks, path files and failures are retained under main/.

---

# Next-step decision

1. Fixed-graph local second pass: B HQ=9/32 versus A=2/32, paired wins/losses=8/1, descriptive exact p=0.0390625. Do not extrapolate from this repeatedly inspected pressure set.
2. Collision recovery: B mean internal/protein severe pairs=0.312/0.219, A=0.969/0.844. Inspect both new and removed pair counts, not only HQ.
3. Local versus full resampling: B=9/32; C_replay=11/32, C_official_seed=9/32, C_official_new=8/32; E=8/32. No per-molecule best-of-N selection was used.
4. Early freezing: D=8/32 at24 logical NFE, with changed chemistry in 32/32 cases; B/C/E each cost48. This is an engineering tradeoff, not a fixed-chemistry causal comparison.
5. Pseudoforces and torque: see the paired same-input directional projections, deterministic rigid shadows and actual pair-distance changes. Outward force does not imply outward rigid displacement, and neither guarantees final HQ.
6. Input representation: saved chemistry/dummy references and actual masked auxiliary edges pass their audits. inter_fragments is unavailable in this checkpoint by design. Additive edge encoders would require a separate architecture/training study, not silently enabling missing weights.
7. Pose retention: B mean drift=1.991 Å; C/E and contact changes are reported separately. This study contains no crystal pose accuracy validation.
8. Efficiency: this run records synchronized A100 path wall time and NFE; historic CPU MMFF timing has different scopes. It does not establish that an extra neural pass is cheaper than terminal MMFF.
9. Framework implication: two-stage dynamic topology plus stable-graph recovery is testable, but favorable pressure-test counts alone do not justify a new training run. First confirm on untouched pockets/seeds under the frozen protocol.
10. Failure attribution: graph construction is audited, but out-of-training-distribution graph partitions and treating the .005 terminal state as a clean proxy remain limitations. Absent direct learned ligand-ligand spatial edges is a concrete checkpoint restriction. The present comparisons cannot uniquely separate input distribution mismatch from learned score capacity.

Recommended next action is an independent confirmatory cohort with the already frozen B/D settings and complete original-stereo reporting, before model/weight changes. Keep the historical MMFF route separate as an engineering reference.
