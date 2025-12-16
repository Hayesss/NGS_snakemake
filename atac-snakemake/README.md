# ATAC-seq Analysis Pipeline (Snakemake)

基于原始 bash 脚本 (ngs00-atac-server) 改写的 Snakemake 自动化 ATAC-seq 分析流程。

## 目录结构

```
atac-snakemake/
├── Snakefile                    # 主 Snakemake 文件
├── config/
│   ├── config.yaml              # 配置文件
│   └── samples.tsv              # 样本信息表
├── workflow/
│   └── scripts/                 # Python 脚本
│       ├── atac-qc.py           # QC 统计和图表
│       ├── bam2bedshift.py      # Tn5 偏移校正
│       ├── countTable.py        # 读数计数
│       ├── count_table_merge.py # 合并计数矩阵
│       ├── generate_nucleosome_track.py  # 核小体定位
│       ├── normalizebedGraph.py # 深度归一化
│       ├── summit_filtration.py # Summit 过滤
│       └── summit_rename.py     # Summit 重命名
└── reference/                   # 参考文件 (需自行准备)
    ├── bowtie2_GRCh38_primary/  # Bowtie2 索引
    ├── GRCh38_chrom.sizes.tsv   # 染色体大小文件
    ├── hg38-blacklist.v2.bed    # ENCODE 黑名单
    └── ENCODE-GRCh38-V29-TSS.4k.bed  # TSS 注释文件
```

## 分析步骤

### 单样本处理 (rule all)

| 步骤 | 描述 | 输出目录 |
|------|------|----------|
| 01 | FastQC 原始数据质控 | 01-fastqc_raw/ |
| 02 | Cutadapt 接头去除 | 02-clean_fastq/ |
| 03 | FastQC 清洗后质控 | 03-fastqc_clean/ |
| 04 | Bowtie2 比对 | 04-align/ |
| 05 | 过滤 (去除 chrM, 低质量) | 05-filter/ |
| 06 | Picard 去重复 | 06-dedup/ |
| 07 | Tn5 偏移校正 (+4/-5 bp) | 07-shift/ |
| 08 | 深度归一化 + BigWig | 08-normalize/ |
| 09 | QC 统计指标 | 09-qc/ |
| 10 | MACS2 Peak Calling | 10-peaks/ |
| 11 | 核小体定位 | 11-nucleosome/ |
| 12 | 读数计数表 | 12-counts/ |

### 下游分析 (rule all_downstream)

| 步骤 | 描述 | 输出目录 |
|------|------|----------|
| DS01 | 合并重复样本 (BED) | 01-merge/ |
| DS02 | 分组 Peak Calling | 02-peaks/ |
| DS03 | Summit 整合 | 02-peaks/ |
| DS04 | 合并计数矩阵 | 04-counts/ |
| DS05 | 分组 BigWig 轨道 | 05-tracks/ |
| DS06 | 合并 BAM 文件 | 06-merge_bam/ |
| DS07 | HOMER Motif 分析 (可选) | 07-homer/ |
| DS08 | DeepTools 可视化 (可选) | 08-deeptools/ |
| DS09 | TF 足迹分析 (可选) | 09-footprint/ |

## 使用方法

### 1. 准备配置文件

编辑 `config/config.yaml`：
```yaml
# 设置原始数据路径
raw_data_dir: "/path/to/your/raw/data"

# 设置参考文件路径
reference:
  bowtie2_index: "reference/bowtie2_GRCh38_primary/genome"
  chrom_sizes: "reference/GRCh38_chrom.sizes.tsv"
  blacklist: "reference/hg38-blacklist.v2.bed"
  tss_bed: "reference/ENCODE-GRCh38-V29-TSS.4k.bed"
  effective_genome_size: 2913022398  # hg38
```

### 2. 准备样本信息表

编辑 `config/samples.tsv`：
```
sample_name	gsm_id	srr_id	group
HSC_4983_rep1	GSM1937376	SRR2920466	HSC
HSC_4983_rep2	GSM1937378	SRR2920468	HSC
MPP_4983_rep1	GSM1937377	SRR2920467	MPP
```

- `sample_name`: 样本唯一标识符
- `gsm_id`: GEO 样本编号 (原始数据子目录名)
- `srr_id`: SRA 运行编号 (fastq 文件前缀)
- `group`: 样本分组 (用于下游分析，如合并重复样本)

### 3. 运行 Pipeline

```bash
# 进入 pipeline 目录
cd /home/zhs/script/atac-snakemake

# 干运行 (查看将执行的任务)
snakemake -n

# 运行单样本分析 (使用 8 个核心)
snakemake --cores 8

# 运行下游分析
snakemake all_downstream --cores 8

# 生成 DAG 流程图
snakemake --dag | dot -Tpng > dag.png

# 生成 rulegraph
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

### 4. 集群运行 (可选)

使用 SLURM 集群：
```bash
snakemake --cluster "sbatch -p normal -n {threads} -t 4:00:00" \
          --jobs 100 --cores 200
```

使用 profile 配置：
```bash
snakemake --profile slurm --cores 200
```

## 依赖软件

### 必需工具

- **Snakemake** >= 7.0
- **Python** >= 3.8 (pandas, numpy, matplotlib, pysam)
- **FastQC** >= 0.11
- **MultiQC** >= 1.0
- **Cutadapt** >= 3.0
- **Bowtie2** >= 2.4
- **SAMtools** >= 1.10
- **Picard** >= 2.0
- **BEDTools** >= 2.29
- **MACS2** >= 2.2
- **UCSC tools** (genomeCoverageBed, bedGraphToBigWig)

### 可选工具 (下游分析)

- **HOMER** (Motif 分析)
- **DeepTools** (可视化)
- **RGT-HINT** (TF 足迹分析)

### 安装依赖

使用 Conda 安装：
```bash
# 创建环境
conda create -n atac-seq python=3.10
conda activate atac-seq

# 安装必需工具
conda install -c bioconda snakemake fastqc multiqc cutadapt bowtie2 samtools picard bedtools macs2
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-genomecoveragebed
conda install pandas numpy matplotlib pysam

# 安装可选工具
conda install -c bioconda homer deeptools
pip install RGT
```

## 参考文件准备

### Bowtie2 索引

```bash
# 下载参考基因组
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# 构建索引
bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa reference/bowtie2_GRCh38_primary/genome
```

### 染色体大小文件

```bash
# 从 UCSC 获取
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from hg38.chromInfo" > reference/GRCh38_chrom.sizes.tsv
```

### ENCODE 黑名单

```bash
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
mv hg38-blacklist.v2.bed reference/
```

### TSS 注释文件

需要准备 BED6 格式的 TSS 区域文件 (通常 TSS ± 2kb)。

## 输出文件说明

### QC 指标 (09-qc/*.QC.table.txt)

| 列名 | 说明 |
|------|------|
| TotalRawReads | 原始读数总数 |
| OverallAlignmentRate% | 总体比对率 |
| FinalMappedReads | 最终保留读数 |
| chrM% | 线粒体读数比例 |
| BlackListReads% | 黑名单区域读数比例 |
| Duplicate% | PCR 重复比例 |
| TSSEnrichScore | TSS 富集分数 |

### Peak 文件

- `*_peaks.narrowPeak`: MACS2 原始 peak 输出
- `*_peaks.narrowPeak.Q0.05.bed`: Q值过滤后的 peak
- `*_peaks.narrowPeak.Q0.05.rm_blacklist.bed`: 去除黑名单后的最终 peak

### 可视化文件

- `*.norm.bw`: 归一化后的 BigWig 轨道文件
- `*.nf.bw`: 核小体自由区域 (NFR) 轨道
- `*.neu.bw`: 核小体区域轨道

## 注意事项

1. **内存使用**: Picard MarkDuplicates 可能需要较大内存，建议为 TMP_DIR 提供足够空间
2. **磁盘空间**: SAM 文件较大，pipeline 中标记为 `temp()` 会自动删除
3. **线程设置**: 根据集群资源调整 `config.yaml` 中的 threads 设置
4. **调试**: 使用 `snakemake -n` 和 `snakemake --dag` 检查 workflow

## 引用

如果使用本 pipeline，请引用：
- 原始脚本: https://github.com/QuKunLab/ATAC-pipe
- Snakemake: https://snakemake.readthedocs.io/
- 相关工具的官方引用

## 许可证

MIT License

