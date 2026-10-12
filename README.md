# 16S 扩增子分析流程 README

本文档说明本目录下 `16s.sh` 全流程分析的环境搭建方法、所用软件，以及每一步的
原理、命令参数与结果文件说明。

---

## 目录

1. [分析概览](#1-分析概览)
2. [目录结构](#2-目录结构)
3. [分析环境搭建](#3-分析环境搭建)
4. [所需软件说明](#4-所需软件说明)
5. [每一步的结果说明](#5-每一步的结果说明)
6. [参数调整](#6-参数调整)
7. [常见问题与注意事项](#7-常见问题与注意事项)

---

## 1. 分析概览

测序区域为 **16S rRNA V3–V4**（Illumina 双端），双端读长 271/262 bp，
拼接后扩增子片段约 **429 bp**，**引物已在测序数据中去除**。

| 项目 | 内容 |
|---|---|
| 样本数 | 9 |
| 分组 | `C`（3 例）、`M`（3 例）、`H`（3 例） |
| 输入数据 | `seq/1_1.fastq.gz … seq/9_2.fastq.gz` |
| 分析环境 | conda 环境 `16s`（R 4.4.3） |
| 运行命令 | `bash /home/cz/16s/16s.sh` |

分析流程为：

```
原始双端序列
   ↓ 双端拼接 (vsearch)
   ↓ 质控过滤 / 去冗余 / UNOISE3 去噪 / 去嵌合体
ASV 代表序列 + 特征表
   ↓ RDP 物种注释 + 去除非细菌序列
   ↓ 等量抽平 (10,000 reads)
   ├─→ Alpha 多样性（6 个指数 + 稀释曲线 + 维恩图）
   ├─→ Beta 多样性（Bray-Curtis/Jaccard/Euclidean/UniFrac + PCoA/CPCoA + 热图）
   ├─→ 物种组成（门/纲/目/科/属 堆叠柱状图 + 弦图 + 树图）
   ├─→ GreenGenes 有参比对
   ├─→ FAPROTAX 功能预测（7 张可视化图）
   └─→ 差异分析（edgeR + 火山图 + 热图 + 曼哈顿图 + STAMP + LEfSe）
```

**本次运行关键结果**：拼接 861,957 条 → 质控后 812,251 条 → 去噪 1,599 个
ASV → 去嵌合体后 **298 个 ASV**，分属 5 门 / 8 纲 / 15 目 / 29 科 / 49 属，
FAPROTAX 共注释出 **92 个功能组**。

---

## 2. 目录结构

```
/home/cz/16s/
├── 16s.sh                  主分析脚本（全流程入口）
├── README.md               本文档
├── seq/                    输入数据：<样本>_1.fastq.gz / <样本>_2.fastq.gz
├── result/                 所有分析结果（详见第 5 节）
│   └── metadata.txt        样本分组表（唯一的必需输入）
├── temp/                   中间文件与运行日志
├── script/                 运行期生成的临时 R 脚本（由辅助 shell 脚本写出）
└── db/                     数据库与辅助脚本
    ├── env16s.yml          conda 环境配置文件
    ├── script/             全部 R / shell / python 辅助脚本
    │   ├── FAPROTAX_1.2.7/ FAPROTAX 数据库与 collapse_table.py
    │   ├── FAPROTAX_1.2.6/
    │   ├── sintax_rank_summary.R      ★ 新增：各分类级汇总
    │   ├── function_faprotax_plot.R   ★ 新增：功能预测可视化
    │   ├── alpha_barplot_fix.R        ★ 新增：α 多样性柱状图（修复版）
    │   └── tax_maptree_treemap.R      ★ 新增：树图备用实现
    ├── usearch/            usearch 二进制 + RDP 16S 数据库
    └── gg/                 GreenGenes 97% 参考序列
```

> **metadata.txt 格式要求**：制表符分隔，第一列必须名为 `Sample`，且必须
> 包含 `Group` 列。数值型样本名（`1`、`2`…）可直接使用。
>
> ```text
> Sample	Group	Repetition
> 1	C	1
> 2	C	2
> ...
> ```

---

## 3. 分析环境搭建

### 3.1 环境已就绪（推荐）

本机 conda 环境 `16s` 已按下方步骤构建完成，直接运行脚本即可：

```bash
bash /home/cz/16s/16s.sh
```

脚本内部会自动执行 `conda activate 16s`，**不会使用 base 环境**。

### 3.2 从零重建环境

```bash
# 1) 创建基础环境（工具链 + R + R 包）
conda env create -f /home/cz/16s/db/env16s.yml

# 2) 激活
source /home/cz/miniconda3/etc/profile.d/conda.sh
conda activate 16s
```

```bash
# 3) 补装 amplicon 包
#    amplicon 已从 CRAN 下架、conda-forge 也没有打包，需手动安装。
#    它是纯 R 包（无编译代码），可从主 conda 的 R 库复制：
cp -r /home/cz/miniconda3/lib/R/library/amplicon \
      "$(conda info --base)/envs/16s/lib/R/library/"
```

### 3.3 环境验证

```bash
source /home/cz/miniconda3/etc/profile.d/conda.sh && conda activate 16s

# 检查命令行工具
for t in vsearch mafft FastTree Rscript python seqkit; do
    printf "%-10s " "$t"; command -v $t || echo MISSING
done

# 检查 R 包
Rscript -e 'for (p in c("amplicon","agricolae","vegan","edgeR","ggplot2",
  "circlize","treemap","pheatmap","rstatix","ggsci","ggpubr","patchwork",
  "picante","ape","multcompView","data.tree","VennDiagram"))
  cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly=TRUE)))'
```

预期：6 个命令行工具全部有路径，R 包全部为 `TRUE`。

### 3.4 数据库文件

数据库位于 `db/` 下。**注意：为控制在 GitHub 的仓库体积限制内，两个大体积
参考库未纳入版本控制**，需按 [第 8 节](#8-数据库获取) 自行下载。

| 文件 | 用途 | 大小 | 是否随仓库提供 |
|---|---|---|---|
| `db/usearch/rdp_16s_v18.fa.gz` | RDP 16S 注释 + 去嵌合体参考（压缩版） | 6.5 MB | ✅ 提供，脚本自动解压 |
| `db/usearch/usearch11.0.667_linux32` | usearch 二进制（免费版） | 2.9 MB | ✅ 提供 |
| `db/script/FAPROTAX_1.2.7/` | FAPROTAX 功能数据库 | 1.1 MB | ✅ 提供 |
| `db/gg/97_otus.fa` | GreenGenes 97% 有参比对 | 137 MB | ❌ 需下载 |
| `db/usearch/utax_reference_dataset_all_25.07.2023.fasta.gz` | UTAX 参考库（流程未使用） | 38 MB | ❌ 不提供 |

> 脚本启动时会自动把 `rdp_16s_v18.fa.gz` 解压为 `rdp_16s_v18.fa`，无需手动处理。

---

## 4. 所需软件说明

### 4.1 命令行工具

| 软件 | 版本 | 用途 |
|---|---|---|
| **vsearch** | 2.31.0 | 核心工具：双端拼接、质控、去冗余、UNOISE3 去噪、去嵌合体、OTU 表构建、SINTAX 注释 |
| **mafft** | 7.526 | ASV 代表序列多序列比对 |
| **FastTree** | — | 由比对结果构建系统发育树（近似最大似然） |
| **R** | 4.4.3 | 全部统计分析与绘图 |
| **Python** | 3.10 | 运行 FAPROTAX 的 `collapse_table.py` |
| **seqkit** | — | FASTA/FASTQ 辅助处理 |
| **usearch** | 11.0.667 | **仅**用于 `-fastx_getseqs`、`-alpha_div` 两个子命令 |

> #### ⚠️ 关于 usearch 的重要说明
>
> `db/usearch/` 中的 usearch 是**免费授权版**，它**拒绝**以下子命令：
>
> | 被拒绝的子命令 | 本流程的替代方案 |
> |---|---|
> | `-unoise3`（ASV 去噪） | `vsearch --cluster_unoise` |
> | `-sintax_summary`（分级汇总） | 自写脚本 `sintax_rank_summary.R` |
> | `-otutab_stats`（表统计） | shell/R 直接统计 |
>
> 因此本流程中 usearch 只承担两个免费功能，其余全部使用 vsearch + R 实现。

### 4.2 R 包

**核心分析包**

| R 包 | 用途 |
|---|---|
| `amplicon` | 多样性/物种组成绘图函数库（1.19.0） |
| `vegan` | 稀释抽平、距离矩阵、稀释曲线 |
| `edgeR` | 差异丰度分析（负二项模型） |
| `ape` | 系统发育树读写与加根 |
| `picante` | UniFrac 系统发育距离 |

**绘图与统计包**

`ggplot2`、`reshape2`、`dplyr`、`tidyr`、`patchwork`、`ggrepel`、`scales`、
`RColorBrewer`、`pheatmap`、`circlize`、`treemap`、`rstatix`、`ggsci`、`ggpubr`、
`multcompView`（显著性字母）、`agricolae`（LSD 检验）、`data.tree`、
`VennDiagram`、`optparse`

### 4.3 分析脚本目录中的辅助脚本

以下脚本是为打通流程而**新增或修改**的，均在 `db/script/` 下：

| 脚本 | 说明 |
|---|---|
| `sintax_rank_summary.R` | ★新增。按门/纲/目/科/属汇总相对丰度，替代付费的 `usearch -sintax_summary` |
| `function_faprotax_plot.R` | ★新增。FAPROTAX 功能的全部可视化（7 张图） |
| `alpha_barplot_fix.R` | ★新增。原 `amplicon::alpha_barplot()` 调用了已废弃的 `summarise_()`，必然报错，此为等价替代 |
| `tax_maptree_treemap.R` | ★新增。树图备用实现（仅依赖 `treemap` 包） |
| `beta_div_gen.R` | 用 `vegan` 计算距离矩阵，替代付费的 `usearch -beta_div` |
| `alpha_rare_curve_gen.R` | 用 `vegan::rarecurve` 生成稀释曲线数据 |
| `alpha_rare_curve_custom.R` | 稀释曲线绘图 |

---

## 5. 每一步的结果说明

脚本共 12 个阶段，控制台会打印 `>>> [n/12]` 进度。下表为本次实际运行结果。

### 阶段 1：环境检查与预处理

- **做什么**：检查 vsearch/mafft/FastTree/Rscript/python/seqkit 是否存在；
  统一修正所有辅助脚本的 Windows 换行符（CRLF）；校验 `metadata.txt`。
- **日志**：`all required tools present` / `CRLF normalisation done`

### 阶段 2：元数据准备

- **做什么**：读取 `result/metadata.txt`，提取样本名与分组。
- **结果**：`9 samples: 1 2 3 4 5 6 7 8 9`，`groups: C H M`
- **中间文件**：`temp/samples.txt`（样本名列表）

### 阶段 3：双端拼接

- **命令**：`vsearch --fastq_mergepairs --fastq_minovlen 20 --fastq_maxdiffs 10`
- **结果**：**861,957 / 867,543 条（99.4%）拼接成功**
- **产物**：`temp/all.fq`（合并后的 FASTQ）、`temp/merged_<样本>.log`（每个样本的拼接统计）

### 阶段 4：质控、去噪与去嵌合体

| 子步骤 | 命令要点 | 结果 |
|---|---|---|
| 质控过滤 | `--fastq_minlen 200 --fastq_maxlen 600 --fastq_maxee_rate 0.01 --fastq_maxns 0` | **812,251** 条 |
| 去冗余 | `--derep_fulllength --minuniquesize 10` | **4,095** 条唯一序列 |
| UNOISE3 去噪 | `--cluster_unoise --minsize 10` | **1,599** 个 ASV |
| 去嵌合体 | `--uchime_ref`（RDP 参考） | **298** 个 ASV |

- **产物**：`result/raw/otus.fa`（298 条 ASV 代表序列，`>ASV_Uni_N`）
- **日志**：`temp/filter.log`、`temp/derep.log`、`temp/unoise.log`、`temp/uchime.log`

> **不切引物的原因**：本批数据引物已去除，read1 直接以 `TAGGGAATCTTCC…`
> 起始（341F 结合位点下游），拼接片段约 429 bp。若仍按原脚本切 29+18 bp，
> 会误切掉 47 bp 有效序列。

### 阶段 5：特征表构建与物种注释

| 项目 | 结果 |
|---|---|
| 特征表（`--usearch_global --id 0.97`） | 298 ASV × 9 样本 |
| SINTAX 注释（RDP，`--sintax_cutoff 0.1`） | 298/298 全部注释成功 |
| 去除非细菌/叶绿体/线粒体 | **全部保留 298 个 ASV**（叶绿体 0、线粒体 0、非特异 0） |

- **主要产物**
  - `result/otutab.txt` —— **最终特征表**（ASV × 样本，原始 counts）
  - `result/otus.fa` —— 过滤后的 ASV 序列
  - `result/otus.sintax` —— SINTAX 注释结果
  - `result/raw/otutab_nonBac.stat` —— 各样本过滤统计
  - `result/otutab.id`、`result/otus.tree` 等中间文件

### 阶段 6：等量抽平

- **命令**：`otutab_rare.R --depth 10000 --seed 1`
- **做什么**：按 10,000 reads/样本随机抽平，消除测序深度差异。
- **产物**
  - `result/otutab_rare.txt` —— 抽平后的特征表（后续 α/β 分析的基础）
  - `result/alpha/vegan.txt` —— 6 个 α 多样性指数（richness、chao1、ACE、shannon、simpson、invsimpson），以样本为行
  - 各样本 richness 范围：**54–150**

### 阶段 7：Alpha 多样性

| 产物 | 说明 |
|---|---|
| `result/alpha/alpha.txt` | usearch 计算的 22 个 α 多样性指数（含 chao1、shannon、simpson 等） |
| `result/alpha/alpha_rare.txt` | 稀释曲线数据（51 个深度点 × 9 样本） |
| `result/alpha/boxplot_{richness,chao1,ACE,shannon,simpson,invsimpson}.pdf` | 6 张箱线图，含 **Tukey HSD 显著性字母** |
| `result/alpha/barplot_richness.pdf` | richness 柱状图（均值 ± 标准差 + 显著性字母） |
| `result/alpha/alpha_rarefaction_curve.pdf` | **稀释曲线**（按分组均值 ± 标准误） |
| `result/alpha/otu_group_exist.txt` | 各组平均相对丰度 > 0.05% 的 ASV 列表 |
| `result/alpha/otu_group_exist.txt.C_H_M.vennDiagram.pdf` | **维恩图**（C/H/M 三组共有与特有 ASV） |
| `alpha_boxplot_TukeyHSD.txt` | 箱线图的多重比较统计表（生成于当前目录） |

**本次结果**：Shannon 指数 C 组最高（均值 ≈ 2.90），M 组居中（≈ 2.44），
H 组最低（≈ 1.80）；箱线图显著性字母显示 C 与 H 差异显著（a vs b）。
稀释曲线在 10,000 深度处仍缓慢上升，说明若需捕捉全部稀有 ASV 可加大测序量。

### 阶段 8：Beta 多样性

| 子步骤 | 命令 | 产物 |
|---|---|---|
| 构建系统树 | `mafft` + `FastTree -nt` | `result/otus.tree`（298 个叶节点） |
| 距离矩阵 | `beta_div_gen.R`（vegan） | `result/beta/bray_curtis.txt`、`jaccard.txt`、`euclidean.txt`、`unweighted_unifrac.txt` |
| PCoA | `beta_pcoa.R` | `bray_curtis.pcoa.pdf`（无标签）、`bray_curtis.pcoa.label.pdf`（带样本号） |
| CPCoA（约束排序） | `beta_cpcoa.R` | `bray_curtis.cpcoa.label.pdf` |
| 距离热图 | `sp_pheatmap.sh` | `bray_curtis.txt.pheatmap.pdf` |

**本次结果**：PCo1 解释 **54.16%**、PCo2 解释 **18.29%** 的变异，
C/M/H 三组在 PCoA 图上清晰分离。

> **说明**：本机 `picante` 只支持**非加权** UniFrac（加权版属 phyloseq 接口），
> 且 `FastTree` 输出的是无根树，脚本已自动按外群加根后再计算。

### 阶段 9：物种组成分析

**注释分级统计**（基于抽平表）：

| 级别 | 分类单元数 | 汇总表 |
|---|---|---|
| 门 Phylum | 5 | `result/tax/sum_p.txt` |
| 纲 Class | 8 | `result/tax/sum_c.txt` |
| 目 Order | 15 | `result/tax/sum_o.txt` |
| 科 Family | 29 | `result/tax/sum_f.txt` |
| 属 Genus | 49 | `result/tax/sum_g.txt` |

**产物**

| 文件 | 说明 |
|---|---|
| `result/taxonomy.txt` | 8 列物种注释表（Kingdom…Species） |
| `result/tax/sum_*.txt` | 各级别相对丰度表（`#OTUID` + 样本列，数值为百分比） |
| `result/tax/sum_{p,c,o,f,g}.stackplot.group.pdf` | 按分组均值的**堆叠柱状图**（Top 8 + Other） |
| `result/tax/sum_{p,c,o,f,g}.stackplot.sample.pdf` | 按单样本的堆叠柱状图 |
| `result/tax/sum_c.circlize.pdf` | 纲水平**弦图**（+ 图例） |
| `result/tax/tax_maptree.pdf` | Top 100 ASV **树图**（面积 = 平均相对丰度，颜色 = 门） |

**本次结果**：优势门为 **Proteobacteria（50.7%）、Firmicutes（45.1%）、
Bacteroidetes（3.8%）**。分组比较：C 组以 Firmicutes 为主（58.6%，
Proteobacteria 34.6%），M 组以 Proteobacteria 为主（71.1%，Firmicutes 23.1%），
H 组两者接近（Firmicutes 53.6%、Proteobacteria 46.2%）。

### 阶段 10：GreenGenes 有参比对

- **命令**：`vsearch --usearch_global --db db/gg/97_otus.fa --id 0.97`
- **结果**：生成 **99,322 个 GreenGenes 97% OTU** × 9 样本的 OTU 表
- **产物**：`result/gg/otutab.txt`、`result/gg/otutab.stat`

> 该步骤为有参比对结果，与 UNOISE3 ASV 结果互为参考。因参考库序列量大，
> 耗时较长（本次约 12 分钟），如需跳过可注释掉脚本第 12 节。

### 阶段 11：FAPROTAX 功能预测 ★新增

- **输入构造**：特征表 + SINTAX 物种注释 → `OTUID / 样本… / taxonomy` 格式
- **命令**：
  ```bash
  python temp/collapse_table.py \
      -i temp/fap_otu.txt -g db/script/FAPROTAX_1.2.7/FAPROTAX.txt \
      -o result/function/faprotax.txt -r result/function/faprotax_report.txt \
      -d taxonomy --omit_columns 0 --column_names_are_in first_data_line \
      --normalize_collapsed columns_before_collapsing
  ```
- **结果**：共注释出 **92 个功能组**，其中 **11 个丰度非零**：
  `chemoheterotrophy`、`aerobic_chemoheterotrophy`、`fermentation`、
  `nitrate_reduction`、`nitrate_respiration`、`nitrogen_respiration`、
  `human_gut`、`human_associated`、`mammal_gut`、
  `animal_parasites_or_symbionts`、`aromatic_compound_degradation`

**产物**

| 文件 | 说明 |
|---|---|
| `result/function/faprotax.txt` | 功能组 × 样本丰度表（92 行） |
| `result/function/faprotax_report.txt` | 每个功能组的匹配记录数报告 |
| `result/function/function_group_mean.txt` | 各功能组的分组均值 |
| `result/function/function_diff.txt` | 差异统计（均值、log2FC、p 值、FDR） |
| `result/function/function_heatmap_data.txt` | 热图所用数据（Top 30 功能） |
| `result/function/figure/function_stackplot.group.pdf` | 分组功能组成**堆叠柱状图** |
| `result/function/figure/function_stackplot.sample.pdf` | 单样本功能组成堆叠柱状图 |
| `result/function/figure/function_heatmap.pdf` | Top 30 功能**热图**（行 z-score，列按分组注释） |
| `result/function/figure/function_boxplot.pdf` | Top 12 功能**箱线图**（分面） |
| `result/function/figure/function_diff_bar.pdf` | 差异功能**分组柱状图**（Top 15） |
| `result/function/figure/function_volcano.pdf` | 差异功能**火山图**（log2FC vs -log10 FDR） |
| `result/function/figure/function_pca.pdf` | 功能谱 **PCA**（CLR 转换 + 68% 置信椭圆） |

**本次结果**：M 组相对 C 组， fermentation 与 gut 相关功能升高
（log2FC ≈ +2.16），而 aerobic_chemoheterotrophy 降低（log2FC ≈ -1.27）。
但每组仅 3 个重复，FDR 校正后无显著功能（最小 FDR = 0.41），
**结论应作为趋势参考而非确证**。

> **注意**：FAPROTAX 的功能组**并非互斥**。例如 `human_gut`、`mammal_gut`、
> `animal_parasites_or_symbionts` 描述的是同一批 taxa，因此丰度曲线完全相同。
> 脚本会自动打印 `NOTE: redundant (identical-profile) functional groups:`
> 提示这些冗余组，避免把它们当作独立证据解读。

### 阶段 12：差异分析 / STAMP / LEfSe

| 子步骤 | 命令/脚本 | 产物 |
|---|---|---|
| 差异丰度 | `compare.R`（edgeR，阈值 0.1） | `result/compare/M-C.txt` |
| 火山图 | `compare_volcano.R` | `result/compare/M-C.volcano.pdf` |
| 差异热图 | `compare_heatmap.sh` | `result/compare/M-C.heatmap.pdf` |
| 曼哈顿图 | `compare_manhattan.sh` | `M-C.manhattan.p.pdf`（门）、`M-C.manhattan.g.pdf`（属） |
| STAMP 输入 | `format2stamp.R` | `result/stamp/tax_{1..8}*.txt`（Kingdom…Species 分级表） |
| STAMP 图 | `compare_stamp.R` | `result/stamp/`（**本次无显著功能，按设计跳过**） |
| LEfSe 输入 | `format2lefse.R` | `result/lefse/LEfSe.txt`、`LEfSe2.txt` |

**本次结果**：M vs C 比较中，edgeR 判定 **1 个 ASV 富集、5 个 ASV 耗竭**、
63 个无差异。

> **STAMP 图为何没有生成？** 每组仅 3 个重复，属水平差异经 FDR 校正后
> 无一项达到 p < 0.05。脚本会明确打印
> `No feature passed the significance threshold (p < 0.05); STAMP figure skipped.`
> 后正常退出，这**属于预期行为，不是报错**。增加生物学重复后即可出图。

### 收尾：清理与校验

- 删除 `temp/` 下的临时 `.fq` / `.fa` 大文件
- 计算输入数据 MD5：`result/md5sum.txt`

---

## 6. 参数调整

所有可调参数集中在 `16s.sh` 顶部第 20–33 行：

```bash
WD="/home/cz/16s"        # 工作目录
DB="${WD}/db"            # 数据库目录
THREADS=8                # 线程数（vsearch / mafft）
CONDA_ENV="16s"          # conda 环境名
GROUP_COL="Group"        # 分组列名
COMPARE_PAIR="M-C"       # 差异比较组（实验组-对照组）
SEED=1                   # 随机种子
MIN_UNIQUE=10            # ASV 最小丰度阈值
RAREFY_DEPTH=10000       # 抽平深度
W=89; H=59               # 图片宽高（mm）
```

常见修改：

| 需求 | 修改 |
|---|---|
| 更换比较组 | `COMPARE_PAIR="H-C"` |
| 提高抽平深度 | `RAREFY_DEPTH=20000`（不可超过最小样本 reads 数 61,799） |
| 提高 ASV 灵敏度 | `MIN_UNIQUE=5` |
| 更换分组列 | `GROUP_COL="Treatment"`（需在 metadata 中存在） |
| 跳过 GreenGenes | 注释掉脚本 `# === 12. GreenGenes ===` 整节 |

---

## 7. 常见问题与注意事项

### 7.1 环境相关

**Q：为什么不能在 base 环境跑？**
A：base 的 R 为 4.5，而 `16s` 环境的 R 为 4.4。**跨版本混用已编译的 R 包
会直接崩溃**（例如 base 的 `cluster.so` 在 R 4.4 下报
`undefined symbol: R_chk_memcpy`）。脚本已固定使用 `16s` 环境，且不添加任何
外部 R 库路径。

**Q：为什么 FAPROTAX 报 `cannot import name 'NaN'`？**
A：`collapse_table.py` 使用 `from numpy import NaN`，该名称在 numpy ≥ 2 中已
移除。脚本会在 `temp/` 中生成一份打补丁的副本运行，**不改动原数据库文件**。

**Q：`amplicon` 装不上怎么办？**
A：该包已从 CRAN 下架，conda 也无打包。因它是纯 R 包，可从主 conda 的 R 库
直接复制到 `16s` 环境的库目录（见 3.2 第 3 步）。

### 7.2 数据相关

**Q：为什么脚本不切引物？**
A：数据中的引物已去除。脚本通过检测 read1 是否以 341F 结合位点下游起始来
判断。**若你的新数据保留了引物，需要自行加回 `--fastq_stripleft/right`。**

**Q：样本名是纯数字，会有问题吗？**
A：不会，但这是最容易踩的坑。R 的 `read.table()` 默认会把数值型列名改写为
`X1`、`X2`…，导致与 metadata 的样本名对不上，进而**静默产出空表**。所有
相关脚本均已加 `check.names=FALSE` 修正。**若你新增或修改脚本，请务必保留
该参数。**

**Q：ASV 标签为什么带 `centroid=` 前缀？**
A：`vsearch --cluster_unoise` 生成的代表序列首部为 `>centroid=Uni_N;size=…`。
脚本会把标签统一为 `ASV_Uni_N` 并**去掉 `;size=` 注释** —— 因为
`vsearch --otutabout` 输出的标签不带该后缀，若不一致，会导致非细菌过滤、
取序列、注释全部匹配失败。

### 7.3 结果解读

| 现象 | 原因 |
|---|---|
| 无显著差异功能/物种 | 每组仅 3 个重复，统计功效不足 |
| 多个功能组曲线完全相同 | FAPROTAX 功能组存在嵌套，非互斥 |
| CPCoA 提示 `Too few points to calculate an ellipse` | 每组仅 3 点，无法拟合椭圆，不影响排序结果 |
| 稀释曲线未饱和 | 测序深度不足以覆盖全部稀有 ASV |
| `alpha_boxplot_TukeyHSD.txt`、`Rplots.pdf` 出现在根目录 | 辅助脚本固定写入当前工作目录，属正常现象 |

### 7.4 重新运行

脚本可重复运行，结果文件会被直接覆盖。若改动了参数（如抽平深度、
比较组），建议先手动清理旧结果，避免残留文件造成误读：

```bash
# 保留唯一必需输入 metadata.txt
cp /home/cz/16s/result/metadata.txt /tmp/metadata.keep
rm -rf /home/cz/16s/result /home/cz/16s/temp /home/cz/16s/script
mkdir -p /home/cz/16s/result && cp /tmp/metadata.keep /home/cz/16s/result/metadata.txt

# 完整重跑（日志同时保存，便于排查）
bash /home/cz/16s/16s.sh 2>&1 | tee /home/cz/16s/run.log
```

**注意**：`result/metadata.txt` 是唯一必需输入，清理前务必先备份。

### 7.5 被修改脚本的备份

为打通流程而对 `db/script/` 下原有脚本所做的修改，均保留了 `.bak` 备份：

```
db/script/compare.R.bak
db/script/compare_heatmap.sh.bak
db/script/compare_manhattan.sh.bak
db/script/otu_mean.R.bak
db/script/otutab_filter_nonBac.R.bak
db/script/otutab_rare.R.bak
db/script/tax_maptree.R.bak
db/script/alpha_rare_curve_custom.R.bak
db/script/beta_div_gen.R.bak
```

对比差异：`diff db/script/compare.R.bak db/script/compare.R`

---

## 8. 数据库获取

因体积限制，以下数据库未随仓库提供。**若不运行 GreenGenes 比对（阶段 10），
仅需 RDP 库，而 RDP 库已随仓库提供，无需任何额外下载。**

### 8.1 GreenGenes 97% 参考序列（阶段 10 需要）

```bash
mkdir -p db/gg && cd db/gg

# 下载 gg_13_8_otus.tar.gz（约 90 MB）
wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz

# 解压并取出 97% 代表序列到 db/gg/97_otus.fa
tar xzf gg_13_8_otus.tar.gz
cp gg_13_8_otus/rep_set/97_otus.fasta 97_otus.fa
```

若下载源不可用，可改用备用镜像：
`http://210.75.224.110/db/gg/gg_13_8_otus.tar.gz`

**跳过该步骤的替代做法**：直接注释掉 `16s.sh` 中
`# === 12. GreenGenes reference mapping ===` 整节，其余分析不受影响。

### 8.2 RDP 16S 数据库

仓库已提供 `db/usearch/rdp_16s_v18.fa.gz`，脚本会自动解压。
如该文件损坏或需重新获取：

```bash
# USEARCH 官方提供的 RDP 数据库
wget -c https://www.drive5.com/sintax/rdp_16s_v18.fa.gz -O db/usearch/rdp_16s_v18.fa.gz
gunzip -k db/usearch/rdp_16s_v18.fa.gz    # 得到 rdp_16s_v18.fa
```

### 8.3 usearch 免费版二进制

仓库已提供 `db/usearch/usearch11.0.667_linux32`（2.9 MB）。
如需重新获取：

```bash
wget -c https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
gunzip usearch11.0.667_i86linux32.gz
mv usearch11.0.667_i86linux32 db/usearch/
chmod +x db/usearch/usearch11.0.667_linux32
```

> 该版本为 32 位免费授权版，仅支持 `-fastx_getseqs`、`-alpha_div` 等基础
> 子命令；`-unoise3`、`-sintax_summary` 需商业授权，本流程已用 vsearch
> 与自写 R 脚本替代。

### 8.4 FAPROTAX 数据库

仓库已提供 `db/script/FAPROTAX_1.2.7/`。官方地址：
<http://www.zoology.ubc.ca/louca/FAPROTAX/>

---

## 附：本次运行完整统计

| 阶段 | 数值 |
|---|---|
| 输入 read pairs | 867,543 |
| 拼接成功 | 861,957（99.4%） |
| 质控后 | 812,251 |
| 去冗余后 | 4,095 |
| 去噪 ASV | 1,599 |
| 去嵌合体后 ASV | **298** |
| 非细菌过滤后 | 298（去除 0） |
| 抽平深度 | 10,000 reads/样本 |
| 分类：门/纲/目/科/属 | 5 / 8 / 15 / 29 / 49 |
| FAPROTAX 功能组 | 92（11 个非零） |
| GreenGenes OTU | 99,322 |
| 生成图片 | 37 张 PDF |
| 运行耗时 | 约 25 分钟（含 GreenGenes 约 12 分钟） |
