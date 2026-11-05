# 冻结权重自主生成多样性评测（全部完成）

正式队列：六个诊断口袋 × 12种子 × 七种策略 = 504次尝试，2016/2016条阶段QC已完成；另有14次pilot及56条QC。此前慢能量检查均已返回，没有新增超时、跳过能量或放宽阈值。最终19/19项核验通过，待检为0。只用一张A100，未训练、未修改权重。

正式神经/刚性/局部柔性终点HQ_no_water分别为75/504、406/504、420/504。各策略和全部失败/拒绝记录保留，不能以修复质量替代模型原始能力。一次Top-1仍全部选择苯身份；随机或交替采样能产生不同片段组合。身份头本身不直接读取口袋，不宣称未见口袋泛化。

## 六项指定交付

- [完整多样性评估](DIVERSITY_EVALUATION_RESULTS.md)
- [身份采样比较](IDENTITY_SAMPLING_COMPARISON.md)
- [联合生成接口审计](JOINT_GENERATION_AUDIT.md)
- [逐口袋多样性](PER_POCKET_DIVERSITY.csv)
- [全部504次生成结果](PER_GENERATION_RESULTS.csv)
- [下一步决策](NEXT_STEP_DECISION.md)

## 复现与核验

[原始任务](prompt.md) · [复现说明](REPRODUCE.md) · [冻结计划与权重SHA](frozen_plan.json) · [质量配置](quality_config.json) · [PoseBusters无水配置](posebusters_no_water.json) · [19项核验](verification.json) · [完成状态](RUN_STATE.json) · [QC汇合状态](quality_completion.json) · [源码](source/)

[结果图](DIVERSITY_QUALITY.png) · [SVG](DIVERSITY_QUALITY.svg) · [策略汇总](STRATEGY_SUMMARY.csv) · [按N分层](PER_N_DIVERSITY.csv) · [跨口袋比较](CROSS_POCKET_COMPARISON.csv) · [构象RMSD](CONFORMER_DIVERSITY.csv) · [全部更新事件](ALL_UPDATE_EVENTS.csv) · [166条拒绝更新](REJECTED_UPDATE_EVENTS.csv)

[身份查询摘要](IDENTITY_QUERY_SUMMARY.csv)逐查询保留全部槽位的概率熵、Top1与苯概率；[初始频率先验](frequency_prior.npz)保留精确数值。[pilot汇总](pilot/summary.json)、[pilot全部尝试摘要](pilot/attempt_summary.csv)、[正式QC协议](main/quality/protocol.json)和[pilot QC协议](pilot/quality/protocol.json)单列。完整训练并集Novelty不可评估，部分索引的范围见[说明](novelty_reference.json)。

## 归档范围

本GitHub目录是精简研究归档：科学报告、提示词、源码、冻结配置、关键CSV/JSON、概率查询摘要、全部拒绝更新事件、图与核验记录。所有已有本地报告/配置/源码/关键表均逐字节复制，不改写实验结论。MANIFEST记录大小、SHA、来源与完整本地交付索引的SHA。

完整本地实验约2.9GB，包括全部SDF/spec、全3684类逐槽位logits及概率NPZ、搜索池、逐原子键角明细、QC JSONL/cache/details、日志和完整18,699文件指纹，全部原样保留在原工作区，未上传这些大体积数据。模型权重和训练数据也不上传。核心失败和接受/拒绝信息已由所收录CSV完整保留；全部504次技术生成均完成。

复现仍需要原本地项目目录布局、已有模型、词表和历史去水口袋；仅克隆这个精简仓库不能独立重跑。过去任务目录保持原样。
