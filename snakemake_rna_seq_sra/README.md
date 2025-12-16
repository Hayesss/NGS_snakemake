## 文件与目录
- `Snakefile`：原始单项目流程，读取 `config.yaml`（或 `config2.yaml`/`config3.yaml`），保持不变。
- `Snakefile_v2`：新增多项目流程，读取 `config_v2.yaml`，输出到 `v2_results/<project>/...`。
- `config.yaml|config2.yaml|config3.yaml`：单项目示例配置。
- `config_v2.yaml`：多项目集中配置，可一次运行多个 SRP。
- `environment.yaml`：流程所需 conda 环境。
- `SRP*_metadata.tsv`：示例 metadata（pysradb 导出 TSV，字段至少含 `run_accession`、`experiment_accession`、`library_layout`）。

## 环境准备
```bash
conda env create -f environment.yaml
conda activate rna_seq
```
若使用外部 snakemake 环境，确保 `STAR`、`samtools`、`fastp`、`deeptools` 可用，并使用 `--use-conda --conda-frontend conda` 启用规则级环境。

## Metadata 获取
```bash
pysradb metadata SRP310170 --saveto SRP310170_metadata.tsv
```
将文件名填入配置（单项目用 `metadata` 字段；多项目用 `projects.<ID>.metadata`）。数据目录需为 `data_dir/SRXxxxx/SRRxxxx.sra` 的结构，以匹配 `fasterq-dump` 输入。

## 使用方式
### 1) 保留的单项目流程（Snakefile）
编辑 `config.yaml`（或 `config2.yaml`/`config3.yaml`），然后：
```bash
snakemake --use-conda --conda-frontend conda -j 40 --resources limit_dump=2 limit_merge=2
```
需要切换项目时，可通过 `--configfile config2.yaml` 等方式指定。

### 2) 新增的多项目流程（Snakefile_v2）
- 在 `config_v2.yaml` 中设置公共参数（索引、线程、输出目录）与 `projects` 列表。
- 运行：
```bash
snakemake -s Snakefile_v2 --use-conda --conda-frontend conda -j 40 --resources limit_dump=2 limit_merge=2
```
- 只跑某个项目可临时注释掉 `config_v2.yaml` 中其他 `projects`，或修改 `out_dir` 隔离输出。
- 目标文件：`v2_results/<project>/03_merged_counts/<project>.tsv`、`02_read_align/*_ReadsPerGene.out.tab`、`04_bigwig/*.bw`。

## 常用说明
- 质控日志输出：`v2_results/<project>/log/fastp/`。
- STAR 日志：`v2_results/<project>/log/<sample>_Log.final.out`。
- 覆盖度图：`04_bigwig/*.bw` 由 `bamCoverage` 生成，可在 `config_v2.yaml.resources.bam_coverage` 调整参数。
- 若使用 gzip 压缩输入，请在 `align_and_count` 中添加 `--readFilesCommand zcat`（默认关闭）。

## 最小配置示例（单项目）
```yaml
star_index: "/home//data/public_data/human_reference/STAR_index"
gtf: "/home//data/public_data/human_reference/gencode.v41.annotation.gtf"
data_dir: "/home/data/public_data/pysradb_downloads/SRP310170"
star_threads: 20
fastp_threads: 8
metadata: "SRP310170_metadata.tsv"
```