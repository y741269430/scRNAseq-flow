# 10X Cell Ranger pipeline
## 目录
- 1.下载软件，基因组，gtf，部署环境变量
- 2.利用自己下载好的基因组，构建References基因组
- 3.运行 cellranger count

参考    
[Getting Started with Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ov)      
[单细胞分析流程之Cell Ranger](https://www.jianshu.com/p/3f01016b5302)     

---
## 1.下载软件，基因组，gtf，部署环境变量
[软件下载地址，里面有References基因组可供下载](https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links)     

解压cell ranger软件    
```bash
tar -xzvf cellranger-7.0.0.tar.gz
```
在环境变量中加入以下内容
```bash
vim ~/.bashrc
```
```bash
export PATH="$PATH:/home/jjyang/downloads/cellranger-7.0.0"
```

### 基因组，gtf下载
这里我选择下载的是[release_M27](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/)版本的基因组及gtf，即：  
GRCm39.genome.fa.gz        
gencode.vM27.annotation.gtf.gz     

解压并删除原gz
```bash
gzip -d GRCm39.genome.fa.gz 
gzip -d gencode.vM27.annotation.gtf.gz 
```

## 2.利用自己下载好的基因组，构建References基因组
第一步如果下载了References基因组，这一步就可以省略了。
```bash
nohup cellranger mkref --genome=mm39 --nthreads=60 \
  --fasta=/home/jjyang/downloads/genome/mm39_GRCm39/ucsc_fa/GRCm39.genome.fa \
  --genes=/home/jjyang/downloads/genome/mm39_GRCm39/gencode.vM27.annotation.gtf &
```
## 3.运行 cellranger count
```bash
nohup cellranger count --id=sample_text \
  --localcores=24 \
  --localmem=64 \
  --transcriptome=/home/jjyang/downloads/mm39/ \
  --fastqs=/home/jjyang/fastq/ &
```
