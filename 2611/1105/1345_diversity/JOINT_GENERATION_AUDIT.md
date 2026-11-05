# 联合生成接口审计

| 模块 | 真实输入输出及兼容性 |
| --- | --- |
| GEOM ConnectionIPA3684 | 全 MASK 当前T/R、t、空受体，输出3684 logits；不直接感知受体或预测连接 |
| 原 sample_joint | record labels 决定默认N；频率先验初始化，.3置信度身份建议，旧 PortConnectionPredictor；原服务用验证集ID查记录 |
| 当前 TimeConnectionPredictor | 身份/端口/相对位姿/时间/化学bag，输出键和立体概率；与旧 Port forward/types 接口不兼容 |
| 当前 ConnectionRuntime/hierarchical_search | 预测身份模板与T/R、价态/端口；完整树候选；不接收真实图 |
| 官方 SigmaDock NativeState | 当前合法化学图、蛋白图、随机位姿；SE(3)运动；动态边、距离和方向每步重算并审计 |
| N/停止 | 身份分类与位姿主干都不输出N或停止概率；无可直接复用生成模块 |

统一 attempt 入口以已构造的口袋context、seed、预定义策略和冻结模块为生成输入；pid仅作文件标识。context_from_pocket直接解析去水PDB，ligand_path=None，以蛋白中心代替真实配体COM，不读fixture、编码配体、验证集、真实标签/连接或参考三维。身份子进程请求白名单仅T/R/t/seed/close，未创建PoseDataset，GEOM initialize不打开数据库。

本轮没有直接冒称原样运行旧 sample_joint。复用其频率先验、.3置信更新、合法候选联合提交、形状转移与失败保留原则；新增无数据集身份服务和概率采样。为符合当前模型，旧 Port proposals/greedy fallback替换为当前 TimeConnectionPredictor + 同一 hierarchical_search；保留官方24步主干、NativeState和SamplingAudit。全程 head=full，未切 legacy_staged。当前模型内部含冻结旧expert的残差基座是其既有架构，不等于将评测头换回旧版。

所有策略使用同一正态sigma2.7位姿初始化、uniform R、蛋白中心和N先验；一次身份方案也允许同步连接更新，保持同预算。模型对预测身份的适配由合法图/官方参考重建验证，不修改预训练特征列。新参考由预测图RDKit嵌入，不用真实三维。更新形状误差上限.5Å，拒绝时保留原图，软件错误、失败及提议事件全部落盘。存在图转移投影，不能完全归功于位姿网络。

模型权重、输入PDB和历史文件SHA在 frozen_plan.json / historical_fingerprints.json，冻结生成源码/资产指纹在 generator_provenance.json，动态计算图审计在每次 sampling_audit.json。无新训练，无水排斥开发，不修改历史数据。

仍缺受体条件身份、经训练连接图条件身份、N/停止、模型感知非水辅因子和未见口袋验证。二维与三维多样性分别报告；词表模板立体与几何读出不当作目标立体预测。

