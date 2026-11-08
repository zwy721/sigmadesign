# GEOM 数据集审计

状态：**全部完成**。只使用GPU 0，首轮正式训练/评估总墙钟预算24小时；不会自动扩大全量GEOM训练。报告在每个完整种子及最终退出时刷新，未完成项不补零。

| split | selected_molecules | usable_molecules | usable_conformers | failed_conformers | MMFF_bonds |
| --- | --- | --- | --- | --- | --- |
| train | 2000 | 2000 | 7845 | 36 | 32532 |
| val | 256 | 256 | 1000 | 6 | 4211 |
| test | 256 | 256 | 1005 | 6 | 4409 |
| scaffold_test | 256 | 256 | 1004 | 4 | 4141 |

来源为现有GEOM原始多构象文件；原SQLite仅冻结第0构象。本轮按完整canonical/isomeric分子图去重，保留原train/dev/test边界，先预留scaffold测试骨架，再排除其他集合中的这些骨架。每分子最多4构象，索引0加3个固定哈希抽样索引，不按能量/质量筛选。

2768个分子均保留；构象失败及原始索引逐项保存在data/*_failures.json、*_sources.json。失败构象不参加拟合，各组使用完全相同可用数据。切分规则为冻结历史切分器；完整图重建、原子映射、键级及立体一致性核验见dataset_verification.json。MMFF94s按完整分子化学环境逐键取平衡键长，GEOM原构象键长独立保存；参考角度为对应GEOM构象。

分子/骨架隔离针对本轮GEOM训练及验证；没有证明这些骨架未出现在既有PDBBind基线预训练中，不作全历史训练宇宙未见声明。
