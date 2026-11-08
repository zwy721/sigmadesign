# GEOM连接与局部几何训练（小规模三种子已完成）

2000/256/256/256 个训练、验证、普通测试与骨架隔离测试分子，最多4构象，3种子，单GPU0，正式训练与评估约4.91小时。450项本地产物及16项冻结输入的哈希复核通过，测试分母完整。

## 核心结论

普通测试原始typed-port F1：A0=0.4455，A1=0.5194，A2=0.5194；骨架测试：0.4575、0.5321、0.5323。新增局部3D特征相对A1未显示显著收益。

A4在Oracle真实连接的普通测试上，键长MAE从2.8473Å降至0.3751Å、角度MAE从32.2712°降至23.6352°，但碰撞比例从58.89%升至70.02%，联合成功率约0.03%。仅键长监督A3导致角度和碰撞恶化。

本轮给定真实片段身份/N，Oracle几何给定真实连接；MMFF修复结果独立标记。SigmaDock仅验证接口兼容，未验证自由生成或docking能力。阶段D有3/6个组种子组合未过验证门槛，按协议跳过，详见核验文件；未启动全量GEOM训练。

## 报告与结果

- [数据审计](GEOM_DATASET_AUDIT.md)
- [连接训练](CONNECTION_TRAINING_RESULTS.md)
- [几何恢复](GEOMETRY_RECOVERY_RESULTS.md)
- [消融](ABLATION_RESULTS.md)
- [SigmaDock接口审计](SIGMADOCK_INTEGRATION_AUDIT.md)
- [下一步决策](NEXT_STEP_DECISION.md)
- [结果汇总](summary.json)、[噪声分层表](NOISE_LEVEL_COMPARISON.csv)
- [冻结配置](protocol.json)、[完成核验](COMPLETION_VERIFICATION.json)、[原始任务](prompt.md)、[归档清单](MANIFEST.json)

## 归档范围

按用户要求仅保留13个核心文件。完整源码、训练历史、拆分索引、失败明细、搜索/MMFF/构象表、权重、数据与缓存均保留在本地geom_connection_geometry_training/，未修改。本地完整交付清单的哈希记录于MANIFEST；报告中提及的其他文件请在本地查阅。本目录是结果归档，复现需要原项目环境与本地完整产物。
