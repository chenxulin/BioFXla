不建议现在直接大规模多路径生成，主要不是算力问题，而是输入证据还没到可放大的质量。
当前我们有的是 seed 数据底座：ChEMBL 活性、Murcko scaffold、PDB metadata、临床 anchor。它能支持“选候选方向”和“验证 REINVENT4 流程”，但还不能可靠支持百万级生成，因为几个关键约束还没闭合：
专利风险没过滤
asundexian/milvexian 及大量 Bayer/BMS/Janssen 相关 chemotype 可能很拥挤。现在直接大规模 Mol2Mol/LibInvent，容易生成一堆 patent-near neighbors，后面大批报废。

选择性证据没接上
FXIa 必须看 FXa、thrombin、kallikrein、trypsin 等反筛。现在只按 FXIa potency 放大，会偏向泛 serine protease binder，尤其是 amidine/benzamidine 类 S1 binders。

PDB 还只是 metadata
我们已有 PDB ligand ID 和 SMILES，但还没解析 S1/S2/S4 occupancy、关键残基、水介导作用和 pose quality。没有这些约束，LinkInvent/de novo scoring 会比较盲。

scaffold seed 不是最终优先级
Tier A 只说明公开 ChEMBL 中有强活性记录，不代表有 patent white space、可合成、可生成、可选择性优化。直接多路径会把“活性强但拥挤/不可用”的方向一起放大。

REINVENT4 多路径需要 scoring function 先定好
大规模生成真正花时间的是后处理：去重、性质过滤、相似度排除、dock/pharmacophore、选择性反筛。没有先跑 pilot 校准阈值，百万级输出会变成很大的噪声池。

我的建议是：不是慢慢做，而是先用 1-2 天级别的 pilot 把闸门调好，再并行放大。比如：
先选 20-50 个非临床候选代表分子跑 Mol2Mol 小批量。
同时准备 LibInvent scaffold templates 和 clinical/patent similarity exclusion set。
用 pilot 输出校准 MW/TPSA/cLogP/charge/SA/similarity 阈值。
然后再开多路径：Mol2Mol + LibInvent + LinkInvent + REINVENT de novo，每条路径 50k-150k 起步。
这样最终放大的不是“会生成很多分子”的流程，而是“有机会留下可用分子”的流程。