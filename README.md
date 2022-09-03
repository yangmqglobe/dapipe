**dapipe** 是一个自动化的 ATAC-Seq 分析流程，在完成配置后，可以自动进行可复现、高并行的 ATAC-Seq 数据分析分析。

## 流程总览

![images/rulegraph.png](images/rulegraph.png)

## 依赖安装

使用 [mamba](https://mamba.readthedocs.io/) 自动安装所有依赖

```bash
mamba env create -f environment.yml
conda activate dapipe
```

## 运行流程

首先需要创建一个本次分析的配置文件

```bash
mkdir project
cd project
cp dapipe/config.yaml .
```

浏览检查配置文件内的所有配置项（所有配置项都是必须的），并按照注释进行更改。然后运行流程

```bash
snakemake --snakefile dapipe/Snakefile -j 10 --configfile config.yaml
```

其中的`-j`是允许流程使用的最大进程数，更多参数参考 [snakemake 官方文档](https://snakemake.readthedocs.io/en/stable/executing/cli.html)。
