> GitHub 精简归档：保留报告、源码、固定参数、关键汇总与核验。完整结构、QC缓存、逐原子明细及完整指纹留在原本地工作区；复现需要原项目依赖。本文仅调整归档链接与此范围说明，实验结论和数值不变。

# 自主生成能力审计

| 路径 | 身份 | 连接 | 位姿 | 槽位/受体限制 |
| --- | --- | --- | --- | --- |
| 历史 native_generate | fixture真实labels/oov | 当前 robust 连接头+合法搜索 | 历史20步前段+官方80步 | 真实N，已知配体定义口袋 |
| 仓库 sample_joint 原型 | joint模式频率先验+可迭代建议 | 早期 PortConnectionPredictor | 官方SE3，可重建图 | 默认真实N；显式slots可避开；旧server要求数据集ID |
| 本次 free_probe 适配器 | 全MASK GEOM头一次top1，不读target record | 当前TimeConnectionPredictor，k1 | 官方24步 | 固定N5；已知口袋；身份头不读受体 |

主队列不得标成自由生成。仓库有自主生成原型，因此“完全没有实现身份生成”也不准确；缺的是当前主评估入口对真实独立输入的完整、经验证集成。

本次适配器创建全新 masked entry，只使用 T/R/t/seed，labels 自身也是 MASK，占位边和立体标签为空，空蛋白 GEOM 分支匹配原 identity_server 的 geom_matched 设置。原 fixture 中只取 pocket/context，未将 labels、真实连接或 target_mol 传给模型。口袋 COM 仍沿用已知结合位点，不构成未知位点搜索。

冻结策略详见 [free_probe/plan.json](free_probe/plan.json)。90 个身份均为苯，说明这个一次 argmax 的高噪声诊断发生坍缩；不能外推为所有其他策略均失败。repair后 HQ 不惩罚重复，也不验证结合，因此必须与1种构造/0指纹多样性一起报告。

未完成能力：受体条件身份模型、可变N生成、当前连接头与迭代身份/位姿的统一闭环、模型感知非水辅因子、显式且经验证的立体生成、冻结参数的多受体盲测和活性/选择性验证。没有以这些缺失为理由虚报完整自主生成成功。
