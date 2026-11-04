# Top-K 的真实 de novo 质量与成本

全部 36 个 pocket/seed 的 8 个不同二维构造候选均完成官方后续采样，共 288 个。复用原有 288 条短程与 104 条完整后续结果，新增 184×64=11,776 次官方主干前向；中途修复另用 36×64=2,304 次。全部候选完整验证避免了只评估旧选择器挑中候选的选择偏差。

下表保留完整生成后通过同一质量门槛的所有分子，并按受体和二维构造跨三个 seed 去重。“任一高质量”以 36 个 pocket/seed 为分母；产量以 12 个受体为分母。未使用参考构造或 RMSD 进行选择。

| 输出/K | 有效连通 | 无严重内碰撞 | 键长 | 键角 | PB mol/dock | 立体完整 | 高质量 | 任一高质量池 | 高质量不同结构/口袋 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pool_raw/1 | 36/36 | 19/36 | 23/36 | 18/36 | 12/5 | 8/36 | 0/36 | 0/36 | 0.000 |
| pool_raw/4 | 144/144 | 70/144 | 84/144 | 68/144 | 43/12 | 28/144 | 1/144 | 1/36 | 0.083 |
| pool_raw/8 | 288/288 | 157/288 | 171/288 | 149/288 | 97/25 | 55/288 | 2/288 | 1/36 | 0.167 |
| pool_raw_stereo/1 | 36/36 | 19/36 | 23/36 | 18/36 | 12/5 | 33/36 | 3/36 | 3/36 | 0.167 |
| pool_raw_stereo/4 | 144/144 | 70/144 | 84/144 | 68/144 | 43/12 | 134/144 | 6/144 | 6/36 | 0.333 |
| pool_raw_stereo/8 | 288/288 | 157/288 | 171/288 | 149/288 | 97/25 | 264/288 | 14/288 | 9/36 | 0.833 |
| pool_post/1 | 36/36 | 36/36 | 31/36 | 33/36 | 30/11 | 8/36 | 3/36 | 3/36 | 0.167 |
| pool_post/4 | 144/144 | 140/144 | 119/144 | 122/144 | 112/38 | 28/144 | 7/144 | 6/36 | 0.333 |
| pool_post/8 | 288/288 | 283/288 | 237/288 | 244/288 | 220/71 | 55/288 | 11/288 | 7/36 | 0.583 |
| pool_post_stereo/1 | 36/36 | 36/36 | 31/36 | 33/36 | 30/11 | 35/36 | 10/36 | 10/36 | 0.583 |
| pool_post_stereo/4 | 144/144 | 140/144 | 119/144 | 122/144 | 112/38 | 140/144 | 37/144 | 18/36 | 1.750 |
| pool_post_stereo/8 | 288/288 | 283/288 | 237/288 | 244/288 | 220/71 | 281/288 | 70/288 | 19/36 | 3.417 |

| 输出/K | 高质量不同构造数（按口袋去重） | 不在所供训练清单内 | 与对应参考构造不同 |
|---|---:|---:|---:|
| pool_raw/1 | 0 | 0 | 0 |
| pool_raw/4 | 1 | 1 | 1 |
| pool_raw/8 | 2 | 2 | 2 |
| pool_post_stereo/1 | 7 | 7 | 6 |
| pool_post_stereo/4 | 21 | 21 | 20 |
| pool_post_stereo/8 | 41 | 40 | 39 |

| 输出/K | 部署策略 | 高质量输出 | 不同结构/口袋 | 每个不同高质量结构 NFE | 采样+修复秒/不同结构 |
|---|---|---:|---:|---:|---:|
| pool_raw/1 | generation_rank | 0 | 0.000 | 无产出 | 无产出 |
| pool_raw/1 | initial_heuristic | 0 | 0.000 | 无产出 | 无产出 |
| pool_raw/1 | short_geometry | 0 | 0.000 | 无产出 | 无产出 |
| pool_raw/1 | full_quality_filter | 0 | 0.000 | 无产出 | 无产出 |
| pool_raw/4 | generation_rank | 0 | 0.000 | 无产出 | 无产出 |
| pool_raw/4 | initial_heuristic | 1 | 0.083 | 3600.0 | 411.8 |
| pool_raw/4 | short_geometry | 1 | 0.083 | 5328.0 | 618.6 |
| pool_raw/4 | full_quality_filter | 1 | 0.083 | 12240.0 | 1846.0 |
| pool_raw/8 | generation_rank | 0 | 0.000 | 无产出 | 无产出 |
| pool_raw/8 | initial_heuristic | 1 | 0.083 | 3600.0 | 411.5 |
| pool_raw/8 | short_geometry | 1 | 0.083 | 7632.0 | 889.7 |
| pool_raw/8 | full_quality_filter | 2 | 0.167 | 11880.0 | 1849.6 |
| pool_post_stereo/1 | generation_rank | 10 | 0.583 | 514.3 | 63.3 |
| pool_post_stereo/1 | initial_heuristic | 10 | 0.583 | 514.3 | 63.3 |
| pool_post_stereo/1 | short_geometry | 10 | 0.583 | 514.3 | 63.3 |
| pool_post_stereo/1 | full_quality_filter | 10 | 0.583 | 514.3 | 63.3 |
| pool_post_stereo/4 | generation_rank | 10 | 0.583 | 514.3 | 63.3 |
| pool_post_stereo/4 | initial_heuristic | 10 | 0.750 | 400.0 | 49.1 |
| pool_post_stereo/4 | short_geometry | 12 | 0.833 | 532.8 | 64.7 |
| pool_post_stereo/4 | full_quality_filter | 37 | 1.750 | 582.9 | 93.7 |
| pool_post_stereo/8 | generation_rank | 10 | 0.583 | 514.3 | 63.3 |
| pool_post_stereo/8 | initial_heuristic | 10 | 0.750 | 400.0 | 49.1 |
| pool_post_stereo/8 | short_geometry | 11 | 0.833 | 763.2 | 91.6 |
| pool_post_stereo/8 | full_quality_filter | 70 | 3.417 | 579.5 | 96.2 |

generation_rank 采用原连接候选顺序；initial_heuristic 采用已有 t=0.6 连接/几何评分；short_geometry 采用已有 t=0.4 几何评分；三者各保留一个。full_quality_filter 在全部候选完成采样、所列后处理和质量检查后保留所有合格不同结构，是可部署的最终过滤，不是提前排序，也不是依赖参考的 oracle。候选池中“任一高质量”和通过数是对应输出管线的覆盖上界，单个早期选择器不保证达到。

部署 NFE：按初始评分选择一个为 20+80=100；短程筛选一个为 20+16K+64（K=1/4/8：100/148/212）；保留全部完整候选为 20+80K（100/340/660）。20 步旧模型前缀在每个 pocket/seed 内共享。原图与三维立体后处理不增加神经网络 NFE，修复增加 CPU 时间。NFE 专指姿态主干，不包括连接头前向。实测秒数为保存的前缀、短程、后续采样与修复耗时之和，受共享 GPU 负载影响；未计未单独留时的 Top-K 搜索和立体赋值，不是端到端延迟。质量评测 CPU 累计耗时另存表中，不能与并行墙钟耗时混同。

本次成功完成的新轨迹合计 14,080 NFE；调度迁移还中断过 8 条未完成后续分支，每条最多 64 NFE，额外损耗的精确前向数未保留，范围为 0–512 NFE。这一操作开销不属于上述部署算法预算，也不能忽略为零；实测成功轨迹耗时不含中断重试损耗。

| 短程特征 | 最终标签 | 正例/288 | ROC AUC |
|---|---|---:|---:|
| initial_score | final_PB | 25/288 | 0.865 |
| initial_score | final_high_3D_stereo | 14/288 | 0.879 |
| short_score | final_PB | 25/288 | 0.842 |
| short_score | final_high_3D_stereo | 14/288 | 0.844 |
| short_bad_bonds | final_PB | 25/288 | 0.737 |
| short_bad_bonds | final_high_3D_stereo | 14/288 | 0.771 |
| short_severe_pairs | final_PB | 25/288 | 0.664 |
| short_severe_pairs | final_high_3D_stereo | 14/288 | 0.714 |

这些 AUC 只衡量固定特征对最终原始三维质量的事后预测能力；没有在 CASF 上拟合阈值或训练选择器。完整 288 个短程质量记录和逐候选短程→最终表同时提供。后处理后的结果与 raw 分开，不将其作为网络直接生成通过率。

所有化学不同均以冻结标准化二维构造定义，未把不同 seed 或同构造立体标签当成新的二维分子。相对于实际 PDBBind 训练清单的缺失只覆盖该训练阶段；不宣称覆盖 GEOM 预训练或官方 SigmaDock 完整训练宇宙。参考恢复仅在通用辅助 JSON 中，未进入评分或质量门槛。

数据：[完整 Top-K 质量](geometry_constraints_20261103/topK_pool_quality.csv)、[逐口袋产量](geometry_constraints_20261103/topK_pocket_yields.csv)、[选择策略与成本](geometry_constraints_20261103/selection_policies.csv)、[短程预测能力](geometry_constraints_20261103/short_predictivity.csv)、[短程逐候选对应](geometry_constraints_20261103/short_to_final.csv)、[完整质量 JSONL](geometry_constraints_20261103/final_quality/molecule_metrics.jsonl)。

混合不同口袋的 AUC 可能含分子大小/口袋差异，应同时看同一 pocket/seed 内的排序。对最终三维立体后处理高质量标签，仅 9/36 个候选池同时具有正负样本：

| 特征 | 池内 AUC 均值 | 正负样本对加权 AUC |
|---|---:|---:|
| initial_score | 0.518 | 0.512 |
| short_score | 0.719 | 0.707 |
| short_bad_bonds | 0.583 | 0.616 |
| short_bad_angles | 0.602 | 0.610 |
| short_severe_pairs | 0.726 | 0.732 |

池内排序明细见 `short_predictivity_within_attempt.csv`。高 AUC 不等于每个候选池的第一名都合格，也不改变完整生成过滤的成本。
