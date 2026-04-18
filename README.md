# GibMut Primer Designer (v0.2.1)

[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[English](#english) | [简体中文](#简体中文)

---

## English

### Introduction
GibMut Primer Designer is a Python automation script for designing primers for **site-directed mutagenesis using Gibson Assembly**. It automatically generates optimal primer pairs that meet thermodynamic requirements for given gene sequences and mutation sites. It supports single-site mutations as well as collision detection and bridging primer design for close-proximity double mutations.

### Dependencies
Ensure you have the following Python libraries installed before running the script:
- `pandas`
- `biopython` (for sequence processing, melting temperature calculation, and sequence alignment)

Install the dependencies via pip:
```bash
pip install pandas biopython
```

### Design Logic
The core logic of this script ensures high amplification efficiency while generating homologous overlap regions required for Gibson Assembly:
1. **Candidate Binding Region Generation**: Uses sliding windows to extract 15-60bp binding fragments. Enforces a 3' GC clamp and uses the Nearest-Neighbor thermodynamic model (`Bio.SeqUtils.MeltingTemp.Tm_NN`) for Tm calculation with salt and primer concentration corrections.
2. **Best Primer Pair Selection**: Prioritizes primer combinations with a total length <= 59bp to reduce synthesis costs and mismatch rates. Strict parameters enforce dTm <= 3.0°C, binding length between 24-36bp, and Tm between 58-68°C. A multi-tier penalty system selects the optimal pair if no perfect match exists.
3. **Core Primer Design**: Strictly controls overlap length between 18-22bp, aiming for an overlap Tm ~5°C lower than the binding Tm.
4. **Double Mutation Collision Detection**: Uses local alignment to map single-mutation primers onto the full-length template. If two mutations are too close and their primers overlap, the script automatically merges them into a single block and designs a bridging primer.

### Usage Instructions
#### 1. Preparation
Place the script, your template sequence file (FASTA format), and the mutation list in the **same directory**.
The mutation list must be a `.csv` file without headers, one mutation per line. Example:
```csv
A123G
V124*
```

#### 2. Running the Script
Execute the following command in your terminal:
```bash
python "GibMutPrimerDesigner v0.2.1.py"
```

#### 3. Interactive Menu
- **Step 1:** Input the `.fa` or `.fasta` template filename (e.g., `template.fasta`).
- **Step 2:** Choose the operation mode:
  - **Mode 1 (Single Mutation):** Input the prepared CSV file. Generates `[filename]-single_mut_primers.csv` with detailed thermodynamic parameters.
  - **Mode 2 (Double Mutation):** Input the CSV generated from Mode 1. Evaluates primer collisions for all mutation pairs and generates `[filename]-double_mut_primers.csv`.

### Output Metrics
The resulting CSV includes:
- **Variant**: Mutation name (e.g., `A123G` or `A123G+V124*`).
- **Forward_Primer / Reverse_Primer**: Final 5' -> 3' sequences (including overlaps).
- **FP_Binding_Tm / RP_Binding_Tm**: Tm values of the template binding region.
- **Pair_dTm**: Difference in Tm between forward and reverse primers.
- **Overlap_Len / Overlap_Tm**: Length and Tm of the overlap region.
- **Warnings**: Alerts for primers >59bp, no GC clamp, abnormal Tm, or requirements for PAGE purification (>90bp).
- **Action**: The strategy taken (e.g., `Redesigned`, or `Use Single Primers Independently`).

### Important Notes
1. **Codon Optimization**: Built-in optimal codon tables for `human`, `ecoli`, and `yeast` are included. `ecoli` is currently hardcoded as the default in the `main()` function. Please modify the script manually if needed.
2. **Sequence Overhang Offset**: The script uses a hardcoded `overhang_len = 40` to offset the absolute coordinates of the mutation. If your FASTA starts directly from the start codon (CDS), change this variable to `0` in the code to prevent addressing errors.
3. **File Encoding**: Ensure your `.csv` and `.fa` files are pure text files without hidden characters.

---

## 简体中文

## 简介
GibMutPrimerDesigner 是一款用于**吉布森组装（Gibson Assembly）定点突变**的引物设计 Python 脚本。它能够自动化地为给定的基因序列和突变位点设计出符合热力学参数的最佳引物对，并支持单突变以及相近双突变的引物设计与引物碰撞检测。

## 依赖环境
在运行该脚本之前，请确保已安装以下 Python 库：
- `pandas`
- `biopython` (用于序列处理、熔解温度计算和序列比对)

可以通过以下命令安装依赖：
```bash
pip install pandas biopython
```

## 代码设计逻辑
该脚本的设计核心是在保证引物扩增效率的同时，生成符合 Gibson 组装要求的同源重叠区（Overlap）。其核心逻辑分为以下几个层级：

1. **候选结合区生成 (`get_all_candidates`)**：
   - 在突变位点两侧滑动窗口，提取长度在 15-60bp 的所有可能引物结合片段。
   - 检查 3' 端是否具有 GC 夹（GC clamp），以提升聚合酶结合稳定性。
   - 使用最近邻热力学模型（`Bio.SeqUtils.MeltingTemp.Tm_NN`）计算 Tm 值（加入盐浓度和引物浓度校正）。

2. **最佳引物对筛选 (`get_best_primer_pair`)**：
   - **最高优先级 (C0)**：优先过滤出使得单条引物预估总长度 $\le$ 59bp 的组合，以降低引物合成成本与错配率。
   - **参数约束**：筛选 dTm $\le$ 3.0°C、结合区长度在 24-36bp、Tm 在 58-68°C 的组合。
   - **多级打分 (Tier System)**：若无完美组合，则逐级放宽长度或 Tm 限制，最终通过罚分函数（Penalty Score）选出相对最优解，并在结果中输出相应警告（Warnings）。

3. **核心引物设计 (`design_core_primer`)**：
   - 严格控制 Overlap 区域长度在 18-22bp。
   - Overlap 的 Tm 值目标设定为比结合区 Tm 低约 5°C。
   - 优先选择能让正反向引物最大长度最小化（即 Overlap 尽量居中）的位点。

4. **双突变碰撞检测 (`design_double_mutations`)**：
   - 通过局部序列比对（Local Alignment）确定单突变引物在全长模板上的绝对覆盖范围（Footprint）。
   - 计算两两突变组合之间是否发生引物区域重叠（Collision）。
   - 如果不重叠（距离较远），建议分别独立使用单突变引物。
   - 如果发生重叠，则将两个突变合并为一个长 Block，重新设计能够同时跨越两个突变位点的桥接引物（Bridging Primer）。

## 使用说明

### 1. 准备工作
将脚本文件 `GibMutPrimerDesigner 0.2.1.py`、模板序列文件（FASTA 格式）和突变位点列表放在**同一个目录**下。

突变位点列表应为无表头的 `.csv` 文件，每行一个突变，格式例如：
```csv
A123G
V124*
```

### 2. 运行脚本
在终端或命令行中运行：
```bash
python "GibMutPrimerDesigner 0.2.1.py"
```

### 3. 交互式菜单
程序启动后会进行交互式引导：

- **第一步：提供 FASTA 模板文件**
  输入存放模板序列的 `.fa` 或 `.fasta` 文件名（需带后缀，例如 `template.fasta`）。

- **第二步：选择功能模式**
  - **模式 1：单突变引物设计**
    输入预先准备好的突变位点 `.csv` 文件名。程序会遍历所有位点，生成对应的正反向引物，并输出带有详细热力学参数的 `[你的文件名]-single_mut_primers.csv`。
  - **模式 2：双突变组合引物设计**
    输入**模式 1** 生成的单突变结果 CSV 文件名。程序会自动评估两两突变组合的引物碰撞情况，如果产生引物重叠则进行合并设计，最终生成 `[你的文件名]-double_mut_primers.csv`。

## 输出结果解析
生成的 CSV 结果包含以下核心列：
- **Variant**: 突变名称（如 `A123G` 或 `A123G+V124*`）。
- **Forward_Primer / Reverse_Primer**: 最终生成的 5' -> 3' 引物序列，已包含 Overlap 区与结合区。
- **FP_Binding_Tm / RP_Binding_Tm**: 模板结合区的 Tm 值。
- **Pair_dTm**: 正反向引物结合区 Tm 的差值（理想状态 $\le$ 3°C）。
- **Overlap_Len / Overlap_Tm**: 引物重叠区的长度与 Tm 值。
- **Warnings**: 警告信息。若引物总长超过 59bp、无 GC 夹或 Tm 异常，会在此提示。若长度超过 90bp，会提示需要 PAGE 纯化。
- **Action**: 该组合采取的策略（如重新设计 `Redesigned`，或双突变相距较远直接建议 `Use Single Primers Independently`）。

## 注意事项
1. **物种密码子优化**：目前脚本内置了 `human`、`ecoli`、`yeast` 的最优密码子表，但在 `main()` 函数中默认硬编码为 `species = "ecoli"`。若需更改物种密码子偏好性，请手动在脚本中修改。
2. **序列偏移 (Overhang Offset)**：脚本内部存在一个硬编码变量 `overhang_len = 40`，它会在计算突变绝对坐标时增加 40bp 的偏移量。如果你的 FASTA 是纯粹从起始密码子开始的 CDS 序列，这会导致寻址错误，请务必根据实际情况将其修改为 `0` 或序列前的实际标签/启动子长度。
3. **文件编码**：确保提供的 CSV 文件和 FASTA 文件为纯文本格式且不含不可见字符。
