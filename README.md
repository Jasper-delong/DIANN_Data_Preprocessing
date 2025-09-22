
# DIANN_Data_Preprocessing
A standardized workflow for processing Data-Independent Acquisition (DIA) proteomics data analyzed by DIANN.

# DIANN 蛋白质组学数据分析流程

这是一个标准化的分析流程，用于处理由 DIANN 软件分析得到的 DIA 标签自由定量 (Label-Free Quantification, LFQ) 蛋白质组学数据。项目从 DIANN 输出的报告文件开始，完成数据质控、预处理、统计分析和核心结果的可视化。

## 项目结构
```
.
├── data/               # 存放原始数据文件 (如 DIANN 的 .tsv 报告)。
│   └── raw_data/       # DIANN 输出的原始报告文件。
├── figures/            # 存放生成的所有图表 (如PCA图、Upset图、热图)。
├── processed_data/     # 存放处理后的中间数据 (如干净的定量矩阵)。
├── scripts/            # 存放所有的分析代码 (Jupyter Notebook)。
│   └── preprocessing/  # 预处理相关的Notebook。
├── LICENSE             # 项目许可证文件。
└── README.md           # 项目说明文件 (就是你正在看的这个)。
```

## 环境设置与安装
1.  克隆本仓库
    ```bash
    git clone https://github.com/Jasper-delong/DIANN_Data_Preprocessing.git
    cd DIANN_Data_Preprocessing
    ```

2.  安装所需依赖
    本项目依赖于多个 Python 库。建议在虚拟环境 (如 Conda) 中安装。
    ```bash
    pip install pandas matplotlib seaborn matplotlib-venn upsetplot scikit-learn pyarrow
    ```

## 分析流程大纲
本流程主要分为两个核心阶段：数据导入与预处理，以及数据质量控制与重复性评估。

### 目录 (Table of Contents)
*   [1. 阶段一：数据导入与预处理 (Data Import & Preprocessing)](#1-阶段一数据导入与预处理-data-import--preprocessing)
    *   [1.1. 数据导入 (Data Import)](#11-数据导入-data-import)
    *   [1.2. 数据清洗 (Data Cleaning)](#12-数据清洗-data-cleaning)
*   [2. 阶段二：质量控制与重复性评估 (Quality Control & Reproducibility Assessment)](#2-阶段二质量控制与重复性评估-quality-control--reproducibility-assessment)
    *   [2.1. 鉴定数量重复性评估 (Identification Reproducibility)](#21-鉴定数量重复性评估-identification-reproducibility)
    *   [2.2. 蛋白质组定量重复性检验 (Quantitative Reproducibility)](#22-蛋白质组定量重复性检验-quantitative-reproducibility)
    *   [2.3. 绝对质量评估 (Absolute Quality Assessment)](#23-绝对质量评估-absolute-quality-assessment)

---

## 1. 阶段一：数据导入与预处理 (Data Import & Preprocessing)
本阶段的核心目标是读取 DIANN 的输出结果，并进行初步的清理工作，为后续分析准备好干净的数据矩阵。

### 1.1. 数据导入 (Data Import)
**状态**: ✅ 已完成

**操作描述**:
使用 Python 的 `pandas` 库读取 DIANN 生成的报告文件。DIANN 会为蛋白质组和肽段分别生成独立的定量矩阵文件。
*   **蛋白质组数据**: 读取 `report.pg_matrix.tsv` 文件。该文件包含了每个蛋白质组在所有样本中的定量信息（通常是 MaxLFQ 强度）。
*   **肽段数据**: 读取 `report.pr_matrix.tsv` 文件。该文件包含了每个肽段母离子在所有样本中的定量信息。
*   **列名重命名**: 原始的列名可能包含完整的文件路径或冗长的文件名，不易阅读。本步骤中包含一个关键的重命名操作，将这些名称映射为简洁的样本名（例如：`JDL_1`, `JDL_2`...），方便后续分析和图表展示。

### 1.2. 数据清洗 (Data Cleaning)
**状态**: ✅ 已完成 (作为可选项)

**操作描述**:
此步骤的目的是剔除在样品制备过程中引入的常见污染物。
*   **移除已知的污染物 (Contaminants)**:
    *   **原因**: 样品制备过程中几乎不可避免地会引入角蛋白 (Keratin, KRT)、胰蛋白酶 (Trypsin) 等常见污染物。它们并非源自目标生物样本，会干扰后续的定量归一化和差异分析。
    *   **操作**: 通过关键词匹配（如 `KRT`, `Trypsin`）筛选蛋白质的 `Genes` 或 `Protein.Ids` 列，移除这些污染物对应的行。在工作流中，此步骤被设计为可选项。

*   **关于 Decoy 和 Q-value 过滤的说明**:
    *   与某些其他软件不同，DIANN 在生成最终的 `report.pg_matrix.tsv` 报告时，**已经在内部完成了基于 Q-value (FDR) 的过滤**（通常为1% FDR）。因此，输出文件中的数据已经是经过筛选的高可信度结果，无需手动进行 Q-value 过滤或移除 Decoy 匹配。

---

## 2. 阶段二：质量控制与重复性评估 (Quality Control & Reproducibility Assessment)
本阶段旨在评估实验的技术稳定性与数据可靠性，是整个分析流程的基石。

### 2.1. 鉴定数量重复性评估 (Identification Reproducibility)
**状态**: ✅ 已完成

**操作描述**:
本步骤旨在评估技术重复样本之间**鉴定结果**的一致性。
*   **蛋白质/肽段鉴定数量统计**:
    *   **目的**: 检查在各个样本中鉴定到的蛋白质和肽段数量是否大致相当。
    *   **操作**: 统计每个样本的有效定量数目（非空值），并生成**条形图 (Bar Plot)**进行可视化比较。
*   **组内重复性重叠分析 (Intra-Group Overlap)**:
    *   **目的**: 评估同一组内的技术重复样本鉴定结果的一致性。
    *   **操作 (智能选择)**:
        *   当组内重复数 **≤ 3** 时，自动绘制**韦恩图 (Venn Diagram)**，直观展示少量样本间的重叠情况。
        *   当组内重复数 **> 3** 时，自动绘制**Upset Plot**，清晰展示多样本间的复杂交集情况。
*   **组间重叠分析 (Inter-Group Overlap)**:
    *   **目的**: 比较不同实验组别之间鉴定的蛋白质有何异同。
    *   **操作**: 绘制**韦恩图**。为确保比较的严谨性，引入了**阈值逻辑**：只有当一个蛋白质在某一组内的至少 N 个重复中被鉴定到时（N可调），才认为该蛋白质代表该组参与比较。

### 2.2. 蛋白质组定量重复性检验 (Quantitative Reproducibility)
**状态**: ✅ 已完成

**操作描述**:
本步骤旨在评估重复样本间**蛋白质丰度测量**的相关性和一致性。分析在对数据进行 Log2 转换后进行，以确保统计的可靠性。
*   **相关性热图 (Correlation Heatmap)**:
    *   **目的**: 计算并展示样本间皮尔逊相关系数，直观地量化样本间的相似度。
    *   **操作**: 生成一个颜色编码的矩阵，数值越接近1（颜色越暖），代表样本间的定量重复性越好。
*   **主成分分析图 (PCA Plot)**:
    *   **目的**: 通过降维将高维的蛋白质组数据在二维空间中展示，检查样本是否按照预设的生物学分组自然聚类。
    *   **操作**: 生成散点图，其中每个点代表一个样本。重复性好的样本会聚集在一起，不同组别的样本应能清晰分开。
*   **定量分布箱线图 (Box Plot)**:
    *   **目的**: 比较各样本定量值的整体分布情况，检查是否存在需要校正的系统性偏差。
    *   **操作**: 为每个样本生成一个箱线图。理想情况下，所有样本的中位数应大致在同一水平线上。

### 2.3. 绝对质量评估 (Absolute Quality Assessment)
**状态**: ✅ 已完成

**操作描述**:
本步骤旨在评估数据本身的绝对质量，确保鉴定结果是可靠的。此分析模块主要基于 DIANN 输出的最全面的主报告文件 **`report.parquet`**。
*   **总体鉴定效率总结**:
    *   **目的**: 从最高层面概览整个实验的成果。
    *   **操作**: 统计在所有（或指定的）样本中，总共鉴定到了多少种不重复的蛋白质组和肽段序列。
*   **每个蛋白质组的肽段支持数**:
    *   **目的**: 评估蛋白质鉴定的证据强度。一个蛋白质被越多的唯一肽段所支持，其鉴定结果就越可靠。
    *   **操作**: 统计每个蛋白质组对应的唯一肽段数量，并绘制**直方图 (Histogram)**，展示其整体分布。通常，大部分蛋白质应由至少2个肽段支持，中位数是一个关键的衡量指标。
*   **酶切效率分析 (Missed Cleavages)**:
    *   **目的**: 评估胰蛋白酶的酶切效率。
    *   **操作**: 智能检查 `report.parquet` 文件中是否存在漏切位点信息。如果存在，则统计唯一肽段中包含 0, 1, 2... 个漏切位点的比例，并生成**条形图**。一个高质量的实验，通常超过70%的肽段应为0个漏切位点。

## 如何使用
1.  将您的 DIANN 输出文件（特别是 `report.pg_matrix.tsv`, `report.pr_matrix.tsv` 和 `report.parquet`）放入 `/data/raw_data` 文件夹。
2.  打开 `/scripts/preprocessing` 文件夹中的主分析 Notebook。
3.  在 Notebook 的**“参数配置区”**中：
    *   修改 `column_mapping` 字典以匹配您的原始文件名和期望的短名称。
    *   配置 `GROUP_DEFINITIONS` 字典以定义您的实验分组。
    *   设置 `GROUP_TO_ANALYZE` 来选择是进行组内分析还是组间比较。
4.  按顺序运行 Notebook 中的所有代码块，生成的图表会自动保存在 `/figures` 目录。

## 作者
*   **[Jasper]**
*   **联系邮箱**: [您的邮箱地址]
*   **GitHub**: @Jasper-delong