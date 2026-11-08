# SigmaDock 接口审计

状态：**全部完成**。只使用GPU 0，首轮正式训练/评估总墙钟预算24小时；不会自动扩大全量GEOM训练。报告在每个完整种子及最终退出时刷新，未完成项不补零。

现有TimeConnectionPredictor仅分类，无几何修正输出；原有端口距离/方向特征不等于全片段原子环境。新增模块独立保存，旧权重和接口文件保持只读。

调用入口：source/inference.py 的 Runtime.predict(types,T,R,t,candidate_ports=None,local_atoms=None,legal_connections=None)。类型可来自身份模型；未给图时先分类和合法搜索，给定合法图时直接几何修正。运行时不读取GEOM原构象、真实配体或参考坐标；真实图只用于明确标注的Oracle训练/评价。

NativeState.slot_frames提取当前状态→保留实际片段内部坐标→神经SE(3)修正→同一合法图的NativeState坐标转移→官方NativeState.step继续去噪。仅做接口兼容，不训练或评价docking效能。

接口运行核验：通过。完整记录见native_integration.json；动态边与几何重新计算已审计。
