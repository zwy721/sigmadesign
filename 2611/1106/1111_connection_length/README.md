# GEOM 连接与局部几何小规模训练（全部完成）

2000/256/256/256 个训练、验证、普通测试与骨架隔离测试分子，最多4构象，3种子，单GPU0；正式训练与评估约4.91小时。完整产物450项及冻结输入16项的SHA-256已复核，评价分母核验通过。

## 结果与限制

普通测试原始typed-port F1：A0=0.4455，A1=0.5194，A2=0.5194；骨架测试：0.4575、0.5321、0.5323。新增局部3D特征A2相对A1未显示显著收益；已有3D通道与二维消融的比较另列。

A4在Oracle真实连接的普通测试中，将键长MAE从2.8473Å降到0.3751Å、角度MAE从32.2712°降到23.6352°，但碰撞比例从58.89%升至70.02%，联合成功率约0.03%。A3仅键长监督使角度与碰撞恶化。不得据此宣称已解决分子几何生成。MMFF独立诊断结果单列，不归于网络收益。

所有组给定真实片段身份与N；Oracle几何使用真实连接。SigmaDock只通过接口兼容性smoke，未验证docking能力。阶段D依验证门槛执行，3/6个组种子组合未通过门槛并按协议跳过；全部决策保存在runs/。未启动全量GEOM训练。

## 文件导航

- [数据审计](GEOM_DATASET_AUDIT.md)、[连接结果](CONNECTION_TRAINING_RESULTS.md)、[几何结果](GEOMETRY_RECOVERY_RESULTS.md)、[消融](ABLATION_RESULTS.md)
- [接口审计](SIGMADOCK_INTEGRATION_AUDIT.md)、[下一步决策](NEXT_STEP_DECISION.md)、[复现](REPRODUCE.md)、[原始任务](prompt.md)
- [连接汇总](CONNECTION_SUMMARY.csv)、[几何汇总](GEOMETRY_SUMMARY.csv)、[噪声分层](NOISE_LEVEL_COMPARISON.csv)、[片段数分层](PER_N_METRICS.csv)
- [全部有限搜索结果](SEARCH_RESULTS.csv)、[MMFF对照](MMFF_BASELINE.csv)、[构象多样性](CONFORMER_DIVERSITY.csv)、[扭转分层汇总](TORSION_SUMMARY.csv)
- [本次完成核验](COMPLETION_VERIFICATION.json)、[归档清单](MANIFEST.json)、[完整本地产物指纹](delivery_manifest.json)、[源码](source/)

## 归档范围

报告、源码、冻结配置、拆分索引、失败记录、训练历史、阶段决策和关键表逐字节保留。新增汇总表从完整结果生成，MANIFEST注明来源。大型PORT_PREDICTION_RESULTS.csv、GEOMETRY_METRICS.csv、TORSION_DISTRIBUTIONS.csv以及权重、张量数据、坐标、缓存和日志保留在本地geom_connection_geometry_training/；完整指纹用于追溯，不代表这些文件均已上传。

此前首种子结束后的报告故障已修复并续跑完成，recovery/保留故障与恢复核验。复现需要原项目布局、外部源码、词表、GEOM数据和已有权重；精简仓库不是独立运行包。历史科学产物未修改。
