# DIANN_Data_Preprocessing
totally same with the Spectronaut_data_preprocessing, but for DIANN
DIANN 蛋白质组学数据分析流程
这是一个标准化的分析流程，用于处理由 DIANN 软件分析得到的 DIA 标签自由定量 (Label-Free Quantification, LFQ) 蛋白质组学数据。项目从 DIANN 输出的报告文件开始，完成数据质控、预处理、统计分析和核心结果的可视化。

# 项目结构
.
├── data/               # 存放原始数据文件 (如 DIANN 的 .tsv 报告)。
│   └── raw_data/       # DIANN 输出的原始报告文件。
├── figures/            # 存放生成的所有图表 (如PCA图、Upset图、热图)。
├── processed_data/     # 存放处理后的中间数据 (如干净的定量矩阵)。
├── scripts/            # 存放所有的分析代码 (Jupyter Notebook)。
│   └── preprocessing/  # 预处理相关的Notebook。
├── LICENSE             # 项目许可证文件。
└── README.md           # 项目说明文件 (就是你正在看的这个)。

# 环境设置与安装
1. 克隆本仓库
```
git clone https://github.com/Jasper-delong/DIANN_Data_Preprocessing.git
cd [您的仓库目录名]
```

2. 安装所需依赖
本项目依赖于多个 Python 库。建议在虚拟环境 (如 Conda) 中安装。
```
pip install pandas matplotlib seaborn matplotlib-venn upsetplot scikit-learn
```

# 分析流程大纲（未完成）
本流程主要分为两个核心阶段：数据导入与预过滤，以及数据质量控制与重复性评估。
目录 (Table of Contents)
1. 阶段一：数据导入与预处理 (Data Import & Preprocessing)
1.1. 数据导入 (Data Import)
1.2. 数据清洗 (Data Cleaning)
2. 阶段二：质量控制与重复性评估 (Quality Control & Reproducibility Assessment)
2.1. 鉴定数量重复性评估 (Identification Reproducibility)

1. 阶段一：数据导入与预处理 (Data Import & Preprocessing)
本阶段的核心目标是读取 DIANN 的输出结果，并进行初步的清理工作，为后续分析准备好干净的数据矩阵。
1.1. 数据导入 (Data Import)
状态: ✅ 已完成
操作描述:
使用 Python 的 pandas 库读取 DIANN 生成的报告文件。DIANN 会为蛋白质组和肽段分别生成独立的定量矩阵文件。
蛋白质组数据: 读取 report.pg_matrix.tsv 文件。该文件包含了每个蛋白质组在所有样本中的定量信息（通常是 MaxLFQ 强度）。
肽段数据: 读取 report.pr_matrix.tsv 文件。该文件包含了每个肽段母离子在所有样本中的定量信息。
列名重命名: 原始的列名可能包含完整的文件路径，不易阅读。本步骤中包含一个关键的重命名操作，将冗长的路径名映射为简洁的样本名（例如：JDL_1, JDL_2...），方便后续分析和图表展示。
1.2. 数据清洗 (Data Cleaning)
状态: ✅ 已完成 (作为可选项)
操作描述:
此步骤的目的是剔除在样品制备过程中引入的常见污染物。
移除已知的污染物 (Contaminants):
原因: 样品制备过程中几乎不可避免地会引入角蛋白 (Keratin, KRT)、胰蛋白酶 (Trypsin) 等常见污染物。它们并非源自目标生物样本，会干扰后续的定量归一化和差异分析。
操作: 通过关键词匹配（如 KRT, Trypsin）筛选蛋白质的 Genes 或 Protein.Ids 列，移除这些污染物对应的行。在工作流中，此步骤被设计为可选项，允许用户根据需要灵活执行。
关于 Decoy 和 Q-value 过滤的说明:
与某些其他软件不同，DIANN 在生成最终的 report.pg_matrix.tsv 报告时，已经在内部完成了基于 Q-value (FDR) 的过滤（通常为1% FDR）。因此，输出文件中的数据已经是经过筛选的高可信度结果，无需手动进行 Q-value 过滤或移除 Decoy 匹配。
2. 阶段二：质量控制与重复性评估 (Quality Control & Reproducibility Assessment)
本阶段旨在评估实验的技术稳定性，是衡量数据是否可靠的关键步骤。
2.1. 鉴定数量重复性评估 (Identification Reproducibility)
状态: ✅ 已完成
操作描述:
本步骤旨在评估技术重复样本之间鉴定结果的一致性。通过一个高度集成的函数实现，该函数能根据分析需求（组内分析 vs 组间比较）和重复样本数量自动选择最合适的可视化方法。
蛋白质/肽段鉴定数量统计:
目的: 检查在各个样本中鉴定到的蛋白质和肽段数量是否大致相当。
操作: 统计每个样本的有效定量数目（非空值），并生成**条形图 (Bar Plot)**进行可视化比较。
组内重复性重叠分析 (Intra-Group Overlap):
目的: 评估同一组内的技术重复样本鉴定结果的一致性。
操作 (智能选择):
当组内重复数 ≤ 3 时，自动绘制韦恩图 (Venn Diagram)，直观展示少量样本间的重叠情况。
当组内重复数 > 3 时，自动绘制Upset Plot，清晰展示多样本间的复杂交集情况。
组间重叠分析 (Inter-Group Overlap):
目的: 比较不同实验组别之间鉴定的蛋白质有何异同。
操作: 绘制韦恩图。为确保比较的严谨性，引入了阈值逻辑：只有当一个蛋白质在某一组内的至少 N 个重复中被鉴定到时（N可调，例如2次），才认为该蛋白质代表该组参与比较。
如何使用
将您的 DIANN 输出文件（特别是 report.pg_matrix.tsv 和 report.pr_matrix.tsv）放入 /data/raw_data 文件夹。
打开 /scripts/preprocessing 文件夹中的主分析 Notebook。
在 Notebook 的**“参数配置区”**中：
修改 column_mapping 字典以匹配您的原始文件名和期望的短名称。
配置 GROUP_DEFINITIONS 字典以定义您的实验分组。
设置 GROUP_TO_ANALYZE 来选择是进行组内分析还是组间比较。
按顺序运行 Notebook 中的所有代码块，生成的图表会自动保存在 /figures 目录。
作者
[Jasper]
联系邮箱: [您的邮箱地址]
GitHub: @Jasper-delong