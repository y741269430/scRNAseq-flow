# scRNAseq-flow

## 目录 ####
- 0.加载以下包
- 1.前期准备工作
- 2.进行质控（这里绘制的是线粒体基因，红细胞基因，核糖体基因的小提琴图）
- 3.评估质量指标 
- 4.筛选 
- 5.合并的Seurat对象 
- 6.整合数据并降维聚类 
- 7.细胞注释 
- 8.拆分样本，进行双细胞去除 
- 9.合并样本 
- 10.画图 
- 11.差异表达分析 
- 12.GO 富集分析 

---
## 0.版本 ####
```r
R.Version()$version.string  # "R version 4.3.0 (2023-04-21 ucrt)"  
R.Version()$platform  # "x86_64-w64-mingw32"
BiocManager::version() # ‘3.17’
```
rtools 版本 4.3.5863

## 这是一个很无赖的安装R包方法，就是直接去官网找对应版本的页面  
- https://bioconductor.org/packages/3.17/BiocViews.html#___Software
- 假如要安装一个CRAN的R包，找一个特定版本（旧版本！！！），就直接在这个页面找  https://cran.r-project.org/src/contrib/Archive/ （/后面加包的名字就会进入旧版本的页面，例如）
https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.2.tar.gz 这个是直接下载该packages手动安装
- 这个是让它自动安装并且把依赖的packages也一并安装  
install.packages('https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.2.tar.gz', repos = NULL, type = 'source')  

#### 一些装包过程中报错的解决方案   
https://github.com/satijalab/seurat-object/issues/166  
https://github.com/bwlewis/irlba/issues/70

## 0.加载以下包 ####
```r
library(SingleCellExperiment) # 1.22.0
library(Seurat) # 4.4.0
library(SeuratObject) # 4.1.4
library(tidyverse) # 2.0.0
library(Matrix) # 1.6-1 # 1.6-1.1 # 1.6-2 这三个版本都可以
library(scales) # 1.2.1
library(cowplot) # 1.1.1
library(RCurl) # 1.98-1.12
library(clustree) # 0.5.0
library(SingleR) # 2.2.0
library(clusterProfiler) # 4.8.1
library(org.Mm.eg.db) # 3.17.0
library(Scillus) # 0.5.0
library(ggpubr) # 0.6.0
library(ggplot2) # 3.4.2
library(DoubletFinder) # 2.0.3
library(harmony) # 1.2.0
library(openxlsx) # 4.2.5.2
library(MySeuratWrappers) # 0.1.0
library(ggsci) # 3.0.0
library(data.table) # 1.14.8
library(scRNAtoolVis) # 0.0.7
```
---
## 1.前期准备工作 ####

```r
# 调整R内允许对象大小的限制（默认值是500*1024 ^ 2 = 500 Mb）
options(future.globals.maxSize = 500 * 1024 ^ 2)
```
#### 设置工作路径 ####
```r
# 读取文件的路径
readpath = 'F:/R work/mmbrain/'
# 保存文件的路径
path = 'F:/R work/mmbrain/results/'
```
#### 创建一个函数来处理数据并创建 Seurat 对象 ####
```r
process_and_create_seurat <- function(data_dir, project) {
      # 读取数据
      count_matrix <- Read10X(data.dir = data_dir, gene.column = 2)
      
      # # 获取基因Symbol并处理NA值
      # gene_symbols <- mapIds(org.Rn.eg.db, keys = count_matrix@Dimnames[[1]], column = "SYMBOL", keytype = "ENSEMBL")
      # gene_symbols <- ifelse(is.na(gene_symbols), paste0("NA-", seq_along(gene_symbols)), gene_symbols)
      # count_matrix@Dimnames[[1]] <- gene_symbols
      
      # 创建 Seurat 对象
      seurat_obj <- CreateSeuratObject(counts = count_matrix,   
                                       min.cells = 3,           # 仅保留在至少3个细胞中有表达的基因（过滤低表达基因）
                                       min.features = 200,      # 仅保留具有至少200个基因表达的细胞（过滤低质量细胞）
                                       project = project)       # 设置项目名称（用于结果的标识）
      return(seurat_obj)
    }
```
#### 设置项目名称（工作路径+项目名称=数据存放的路径） ####
```r
projects <- c("A_con", "B_tre")
```
#### 创建一个空的 Seurat 对象列表 ####
```r
seurat_objects <- list()

# 循环处理每个项目
for (project in projects) {
  # 构建数据目录路径
  data_dir <- file.path(readpath, project, 'filtered_feature_bc_matrix/')
  
  # 处理数据并创建 Seurat 对象
  seurat_obj <- process_and_create_seurat(data_dir, project)
  
  # 将 Seurat 对象添加到列表中
  seurat_objects[[project]] <- seurat_obj
}

# 检查矩阵行名是否为SYMBOL
head(rownames(seurat_objects[[1]]@assays$RNA))
```
---
## 2.进行质控（这里绘制的是线粒体基因，红细胞基因，核糖体基因的小提琴图） ####

#### 检索规则 ####
```r
# 查看线粒体基因
rownames(seurat_objects[[1]]@assays$RNA)[grep('^mt-', rownames(seurat_objects[[1]]@assays$RNA))]
# 查看红细胞基因
rownames(seurat_objects[[1]]@assays$RNA)[grep('^Hba-', rownames(seurat_objects[[1]]@assays$RNA))]
# 查看核糖体基因
rownames(seurat_objects[[1]]@assays$RNA)[grep('^Rp[sl]', rownames(seurat_objects[[1]]@assays$RNA))]

# 以下是红细胞的基因list，可以人为补充进去，作为筛选
Hb.genes_total <- c("Hbb-bt","Hbb-bs","Hba-a3","Hba-a2","Hba-a1",
                    "Hba-a2.1","Alas2","Tent5c","Fech","Bpgm")
```                
#### 作图 ####    
```r
Hb.genes <- list()
plot <- list()
for (i in 1:length(seurat_objects)){
  # 将每个细胞每个UMI的基因数目添加到metadata中
  seurat_objects[[i]]$log10GenesPerUMI <- log10(seurat_objects[[i]]$nFeature_RNA) / log10(seurat_objects[[i]]$nCount_RNA)
  
  # 计算线粒体基因比率
  seurat_objects[[i]]$mitoRatio <- PercentageFeatureSet(object = seurat_objects[[i]], pattern = "^mt-")
  seurat_objects[[i]]$mitoRatio <- seurat_objects[[i]]@meta.data$mitoRatio / 100
  
  # 计算红细胞基因比率
  Hb.genes[[i]] <- rownames(seurat_objects[[i]]@assays$RNA)[match(Hb.genes_total, 
                                                                  rownames(seurat_objects[[i]]@assays$RNA))]
  Hb.genes[[i]] <- Hb.genes[[i]][!is.na(Hb.genes[[i]])]
  
  seurat_objects[[i]]$percent.Hb <- PercentageFeatureSet(seurat_objects[[i]], features = Hb.genes[[i]])
  
  # 计算核糖体基因比率
  seurat_objects[[i]]$percent.ribo <- PercentageFeatureSet(object = seurat_objects[[i]], pattern = "^Rp[sl]")
  seurat_objects[[i]]$percent.ribo <- seurat_objects[[i]]@meta.data$percent.ribo
  
  # 为metadata添加细胞ID
  seurat_objects[[i]]@meta.data$cells <- rownames(seurat_objects[[i]]@meta.data)
  
  # 画图，设置cutoff
  a1 <- VlnPlot(seurat_objects[[i]], features = c("nFeature_RNA"), pt.size = 0) + NoLegend()
  a2 <- VlnPlot(seurat_objects[[i]], features = c("nCount_RNA"), pt.size = 0) + NoLegend()
  a3 <- VlnPlot(seurat_objects[[i]], features = c("mitoRatio"), pt.size = 0) + NoLegend()
  a4 <- VlnPlot(seurat_objects[[i]], features = c("percent.Hb"), pt.size = 0) + NoLegend()
  a5 <- VlnPlot(seurat_objects[[i]], features = c("percent.ribo"), pt.size = 0) + NoLegend()
  
  plot[[i]] <- plot_grid(a1, a2, a3, a4, a5, align = "v", nrow = 1)
}

p <- plot_grid(plotlist = plot, align = "v", ncol = 1)

ggplot2::ggsave(paste0(path, "QC_VlnPlot_five.pdf"), plot = p, 
                height = 8, width = 13, dpi = 300, limitsize = FALSE)
```
#### 复制并添加3列到metadata，用来画图的 #### 
```r
# UMI（或 nCount_RNA）表示每个基因的转录本数量
# nFeature表示检测到的基因数量（即细胞中有多少基因表达）

seurat_objects <- lapply(seurat_objects, function(x){
  x@meta.data$sample <- x@meta.data$orig.ident
  x@meta.data$nUMI <- x@meta.data$nCount_RNA
  x@meta.data$nGene <- x@meta.data$nFeature_RNA
  return(x)
})

# 取出未过滤之前的metadata，用来画图的
metadata0 <- lapply(seurat_objects, function(x){
  x <- x@meta.data
  return(x)
})
metadata <- Reduce(rbind, metadata0)

table(metadata$sample)

# # 可以把未过滤之前的数据保存一下，反正我是不存的。（可选）
# saveRDS(seurat_objects, paste0(path, "seurat_objects.rds"))
```
---
## 3.评估质量指标 ####
### 图1 
```r
# 细胞计数 (Cell counts per sample)
# 可视化每个样本的细胞计数
before <- metadata %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar(position = "dodge", show.legend = TRUE) +
  geom_text(stat = 'count', aes(label = after_stat(count)), position = position_dodge(0.9), vjust = -0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Cell counts per sample")

# 每个细胞的UMI计数 (UMI counts per cell)
# 可视化每个细胞的UMI计数
a1 <- metadata %>% 
  ggplot(aes(color = sample, x = nUMI, fill = sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10(breaks = c(500, 1000, 5000, 10000, 20000)) + 
  theme_classic() +
  xlab("nUMI") +
  ylab("Cell density") +
  geom_vline(xintercept = c(500, 1000, 5000, 10000, 20000), linetype = 'dotted') +
  ggtitle("UMI counts per cell")

# 复杂度 (Complexity)
# 通过可视化基因数与UMI的比率（log10基因数/UMI）来表示基因表达的整体复杂性
a2 <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = c(0.8, 0.85, 0.9, 0.95), linetype = 'dotted') +
  ggtitle("Complexity")

# 每个细胞检测到的基因数分布
a3 <- metadata %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill = sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) + 
  geom_vline(xintercept = c(500, 1000, 2500, 5000, 10000), linetype = 'dotted') +
  ggtitle("Genes per cell")

# 每个细胞检测到的基因数量的分布（箱线图）
a4 <- metadata %>% 
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Genes detected per cell across samples")

p <- plot_grid(a1, a2, a3, a4, align = "v", nrow = 2)

ggplot2::ggsave(paste0(path, "QC_four.pdf"), plot = p, 
                height = 8, width = 12, dpi = 300, limitsize = FALSE)
``` 
### 图2    
```r
# 线粒体基因计数占比 (Mitochondrial counts ratio)
# 可视化每个细胞检测到的线粒体基因表达分布
a1 <- metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.25) + 
  theme_classic() +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.3)) + 
  geom_vline(xintercept = c(0.001, 0.01, 0.1, 0.3), linetype = 'dotted') +
  ggtitle("Mitochondrial counts ratio")

# 红细胞基因计数占比 (HB counts ratio)
# 可视化每个细胞检测到的红细胞基因表达分布
a2 <- metadata %>% 
  ggplot(aes(color=sample, x=percent.Hb, fill=sample)) + 
  geom_density(alpha = 0.25) + 
  theme_classic() +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.5)) + 
  geom_vline(xintercept = c(0.001, 0.01, 0.1, 0.5), linetype = 'dotted') +
  ggtitle("Hb counts ratio")

# 核糖体基因计数占比 (RB counts ratio)
# 可视化每个细胞检测到的核糖体基因表达分布
a3 <- metadata %>% 
  ggplot(aes(color=sample, x=percent.ribo, fill=sample)) + 
  geom_density(alpha = 0.25) + 
  theme_classic() +
  scale_x_log10(breaks = c(0.1, 1, 10)) + 
  geom_vline(xintercept = c(0.1, 1, 10), linetype = 'dotted') +
  ggtitle("Ribosome counts ratio")

p <- plot_grid(a1, a2, a3, align = "v", nrow = 1)

ggplot2::ggsave(paste0(path, "QC_Geom_density_three.pdf"), plot = p, 
                height = 3, width = 15, dpi = 300, limitsize = FALSE)

# 检测到的UMI数对比基因数 (UMIs vs. genes detected)
# 可视化每个细胞中检测到的基因数（nFeature_RNA）与UMI数（nCount_RNA）之间的关系，
# 颜色代表线粒体基因计数占比（mitoRatio），并观察是否存在大量低基因数或低UMI数的细胞。

p <- metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

ggplot2::ggsave(paste0(path, "QC_UMIs_vs_genes_detected.pdf"), plot = p, 
                height = 8, width = 8, dpi = 300, limitsize = FALSE)
``` 
---
## 4.筛选 ####
```r 
seurat_filter <- lapply(seurat_objects, function(x){
  x <- subset(x,
              subset = (metadata$nCount_RNA >= 500) &
                (nFeature_RNA >= 500) & (nFeature_RNA <= 20000) &
                (log10GenesPerUMI > 0.85) & (mitoRatio < 0.03) )
})

# 取出过滤之后的metadata，用来画图的
metadata2 <- lapply(seurat_filter, function(x){ x <- x@meta.data })

metadata2 <- Reduce(rbind, metadata2)

after <- metadata2 %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) +
  geom_bar(position = "dodge", show.legend = TRUE) +
  geom_text(stat = 'count', aes(label = after_stat(count)), position = position_dodge(0.9), vjust = -0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Cell counts per sample")

p <- plot_grid(before, after)

ggplot2::ggsave(paste0(path, "QC_NCells.pdf"), plot = p, 
                height = 6, width = 8, dpi = 300, limitsize = FALSE)
``` 
---
## 5.合并的Seurat对象 ####
```r
seurat_merged <- merge(x = seurat_filter[[1]], 
                       y = c(seurat_filter[[2]]),
                       add.cell.id = projects)

# 基于表达丰度的基因筛选 (Gene-level filtering based on expression abundance)
# 提取计数
counts <- GetAssayData(object = seurat_merged, slot = "counts")

# 生成一个逻辑矩阵，其中每个基因在每个细胞中的表达计数是否大于0
nonzero <- counts > 0

# 计算每个基因在多少个细胞中被检测到（表达计数大于0的细胞数）
# 如果一个基因在10个或更多细胞中被检测到，则保留该基因
keep_genes <- Matrix::rowSums(nonzero) >= 10

# 过滤得到在至少10个细胞中有表达的基因
filtered_counts <- counts[keep_genes, ]

# 用过滤后的基因表达矩阵重新创建一个Seurat对象，并保留原有的元数据
seurat_merged2 <- CreateSeuratObject(filtered_counts, 
                                     meta.data = seurat_merged@meta.data)
```
---
## （选做）把线粒体基因，红细胞基因，核糖体基因全部从矩阵中删除，但是不去除细胞 （以下内容为选做，不做就直接保存seurat_merged2） 

#### 过滤以下基因： 

0. 仅保留在至少3个细胞中有表达的基因（过滤低表达基因）（读取数据的时候已经跑过了） 
0. 仅保留具有至少200个基因表达的细胞（读取数据的时候已经跑过了） 
1. 线粒体基因（以 MT- 开头的基因）
2. 核糖体基因（以 RPS 或 RPL 开头的基因）
3. 血红蛋白基因（明确的基因列表 HBA1, HBA2, HBB 等）
4. 其他不需要的基因，包括（可以自己在里面添加检索的规则）：
   - 长非编码RNA基因
   - 假基因或低表达基因（例如 FO, FP, AL, AC 等开头的基因）
   - miRNA（MIR 开头的基因）
   - 其他特定基因（C1orf, AP0, -AS1 和 Z98 含有的基因等）

```r
adult.genes <- rownames(seurat_merged2@assays$RNA)

# 通过检索相关字段获取基因，并展示出前100个
Gm.gene <- grep('^Gm',adult.genes, value = T); head(Gm.gene, 100)
mit.gene <- grep('^mt-',adult.genes, value = T); head(mit.gene, 100)
rib.gene <- grep('^Rp[sl]',adult.genes,value =T); head(rib.gene, 100)
Rik.gene <- grep('Rik',adult.genes, value = T); head(Rik.gene, 100)
p.gene <- grep('\\.',adult.genes, value = T); head(p.gene, 100)
num.gene <- grep("[0-9]{4,}", adult.genes, value = T); head(num.gene, 100)
A.gene <- grep('^A[A-Z][0-9]',adult.genes, value = T); head(A.gene, 100)
B.gene <- grep('^B[A-Z][0-9]',adult.genes, value = T); head(B.gene, 100)

Hb.genes_total <- c("Hbb-bt","Hbb-bs","Hba-a3","Hba-a2","Hba-a1",
                    "Hba-a2.1","Alas2","Tent5c","Fech","Bpgm")
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/a1.png" />  

```r
# 删除所有不需要的基因，如果你有想去除的基因list，则在后面补充  
clean.gene <- adult.genes[which(!(adult.genes %in% c(Gm.gene, mit.gene, rib.gene, Rik.gene, 
                                                     p.gene, num.gene, A.gene, B.gene, Hb.genes_total)))]
# 更新Seurat对象，仅保留筛选后的基因  
seurat_merged2 <- seurat_merged2[clean.gene, ]

# 保存  
saveRDS(seurat_merged2, paste0(path, "seurat_merged.rds"))   

# 保存之后清空，释放内存  
rm(list = ls())
gc()
```
---
## 6.整合数据并降维聚类 ####
#### 设置工作路径 ####
```r
readpath = 'F:/R work/mmbrain/'
path = 'F:/R work/mmbrain/results/'
```
#### 使用Harmony做整合分析 ####
```r
seurat_integrated <- readRDS(paste0(path, "seurat_merged.rds"))

seurat_integrated <- Seurat::NormalizeData(seurat_integrated, verbose = FALSE, normalization.method = 'LogNormalize')
seurat_integrated <- FindVariableFeatures(seurat_integrated, selection.method = "vst", nfeatures = 2000)
seurat_integrated <- ScaleData(seurat_integrated) 
seurat_integrated <- RunPCA(seurat_integrated, features = VariableFeatures(seurat_integrated))
seurat_integrated <- RunHarmony(seurat_integrated, 'orig.ident')
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20, reduction = 'harmony')
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:20, reduction = 'harmony')
seurat_integrated <- FindClusters(seurat_integrated, resolution = seq(from = 0.4,by = 0.2, length = 3))

Idents(seurat_integrated) <- seurat_integrated$RNA_snn_res.0.4
```    
#### 画图UMAP ####
```r
p1 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1) 

ggplot2::ggsave(paste0(path, "UMAP_1.pdf"), plot = p1, 
                height = 5, width = 7, dpi = 300, limitsize = FALSE)
```    
#### 画图UMAP 分开样本 ####
```r
p2 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = TRUE, split.by = 'sample',
              label.size = 3) + coord_equal(ratio = 1) 

ggplot2::ggsave(paste0(path, "UMAP_split1.pdf"), plot = p2, 
                height = 5, width = 9, dpi = 300, limitsize = FALSE)
```
---
## 7.细胞注释 ####  

接下来，你可以选择是： 
- 1.先做双细胞去除再做细胞注释； 
- 2.先做细胞注释，然后根据细胞注释去做双细胞去除。（这里我做第二种。）  
我一般用细胞marker去做注释，你也可以用singleR包注释  
https://cellxgene.cziscience.com/cellguide 这个网站可以查人类或小鼠的细胞marker  

```r
marker_ls <- list(Excitatory_neuron = c('Slc17a7','Slc17a6','Rorb', 'Sulf2', 'Cux2') ,
                  
                  Inhibitory_neuron = c('Gad1', 'Gad2', 'Slc32a1', 'Slc6a1', 'Nrxn3', 
                                        'Pvalb', 'Sst', 'Npy1r', 'Npy2r', 'Npy5r', 'Vip'),
                  
                  Astrocyte = c('Gja1','Gfap','Atp1b2','Agt','Atp13a4','Aldh1l1',
                                'Camsap1','Acsl3','Acsl6','Aqp4'),
                  
                  Endothelial_cell = c('Cldn5','Pecam1','Flt1','Fn1','Tek','Abcb1a',
                                       'Egfl7','Adgrl4','Eng','Adgrf5','Igfbp7'),
                  
                  Ependymal_cell = c('Dnah11','Tmem212','Acta2','Adamts20','Ak7',
                                     'Ak9','Armc3','Cfap65','Ccdc114','Ccdc146'),
                  
                  Microglial_cell = c('Cx3cr1','Tmem119','C1qa','C1qb','Aif1','Itgam',
                                      'C1qc','Ccr5','Csf1r','Ctss','Cyth4',
                                      'Adgre1','Tgfb1','Apbb1ip'),
                  
                  Macrophage = c('Apoe','Mrc1','Itgam','Ptprc','Dock2','Abca9','Aoah','Arhgap15',
                                 'Bank1','Cd40','Cd84','Cfh','Cmah'),
                  
                  Oligodendrocyte = c('Mog','Mbp','Olig2','Ptgds','Apod','Aspa',
                                      'Ccp110','Cldn11','Cntn2','Efnb3',
                                      'Opalin','Mobp'),
                  
                  OPCs = c('Pdgfra','Cspg4', 'Ptprz1','Cspg5','Olig1',
                           'Lhfpl3','Epn2','Serpine2','Luzp2','Bcan'))
```  
#### 输入细胞群的marker基因，判断其表达丰度
```r
vln_ls <- list()
for (i in 1:length(marker_ls)) {
  vln_ls[[i]] <- VlnPlot(seurat_integrated, features = marker_ls[[i]], stack=T, pt.size = 0) + 
    ggtitle(names(marker_ls)[i]) 
}

a1 <- plot_grid(plotlist = vln_ls, ncol = 1)

ggplot2::ggsave(paste0(path, "VlnPlot_marker_1.pdf"), plot = a1,
                height = 30, width = 12, dpi = 300, limitsize = FALSE)

Idents(seurat_integrated)
```  
#### 把细胞填入以下的Idents中（这里面是我瞎填的，仅供教学使用）
```r
seurat_integrated <- RenameIdents(seurat_integrated,
                                  "0" = "Neuron",
                                  "1" = "Astrocyte",
                                  "2" = "Ependymal_cell",
                                  "3" = "Oligodendrocyte",
                                  "4" = "Neuron",
                                  "5" = "Neuron",
                                  "6" = "Neuron",
                                  "7" = "OPCs",
                                  "8" = "Neuron",
                                  "9" = "Oligodendrocyte",
                                  "10" = "Neuron",
                                  "11" = "Neuron",
                                  "12" = "Microglial_cell",
                                  "13" = "Ependymal_cell",
                                  "14" = "Endothelial_cell",
                                  "15" = "OPCs",
                                  "16" = "Neuron",
                                  "17" = "Ependymal_cell",
                                  "18" = "Neuron",
                                  "19" = "Oligodendrocyte",
                                  "20" = "Ependymal_cell")
```                             
#### 修改metadata ####       
```r
# 将这个Idents赋值给celltype
seurat_integrated$celltype <- Idents(seurat_integrated)

# 设置因子的顺序
seurat_integrated$celltype <- factor(seurat_integrated$celltype,
                                     levels=c('Neuron', 'Astrocyte', 'Endothelial_cell',
                                              'Ependymal_cell', 'Microglial_cell', 'Oligodendrocyte', 'OPCs'))
# 把因子的顺序再赋值给Idents
Idents(seurat_integrated) <- seurat_integrated$celltype

# 设置样本顺序
table(seurat_integrated$sample)
seurat_integrated$sample <- factor(seurat_integrated$sample,
                                   levels=c("A_con", "B_tre"))
```
#### 第二次画图UMAP ####
```r
p1 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1) 

ggplot2::ggsave(paste0(path, "UMAP_2.pdf"), plot = p1, 
                height = 5, width = 7, dpi = 300, limitsize = FALSE)
```   
#### 画图UMAP 分开样本 ####
```r   
p2 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = TRUE, split.by = 'sample',
              label.size = 3) + coord_equal(ratio = 1) 

ggplot2::ggsave(paste0(path, "UMAP_split2.pdf"), plot = p2, 
                height = 5, width = 9, dpi = 300, limitsize = FALSE)

# 保存
saveRDS(seurat_integrated, paste0(path, "seurat_integrated.rds"))


# 保存之后清空，释放内存
rm(list = ls())
gc()
```
## 8.拆分样本，进行双细胞去除 ####

具体参考  
- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)  
- [What is the maximum number of cells that can be profiled?](https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled)
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/doublet_rate.png" width="600" />

#### 设置工作路径 ####

```r
readpath = 'F:/R work/mmbrain/'
path = 'F:/R work/mmbrain/results/'
```
#### 去除双细胞 ####
```r
seurat_integrated <- readRDS(paste0(path, "seurat_integrated.rds"))

# 根据样本进行拆分
seurat_split <- SplitObject(seurat_integrated, split.by = "sample")

# 这时候你可以把seurat_integrated清理掉，释放内存
rm(seurat_integrated); gc()
```    
#### 这里我对第一个样本进行双细胞去除（一共两个样本）。耗时大概20分钟一个样本吧，内存消耗不大（剩下另外一个样本，也要重复这个工作，记得把名称改了不然会覆盖掉）
```r
input_seurat <- seurat_split[[1]]
```
#### 继续往下跑，跑完第一个跑第二个 #### 
```r
seurat_1 <- paramSweep_v3(input_seurat, PCs = 1:10, sct = FALSE)
seurat_2 <- summarizeSweep(seurat_1, GT = F)
seurat_3 <- find.pK(seurat_2)

# 这里面mpK的值，我一般会保存在注释里面
mpK <- as.numeric(as.vector(seurat_3$pK[which.max(seurat_3$BCmetric)])); mpK # mpK = 0.005

annotations <- input_seurat$celltype
homotypic.prop <- modelHomotypic(annotations); homotypic.prop

# 这里面DoubleRate 计算出了0.169784的双细胞，也就是说，它要删去你16%的细胞OMG
DoubletRate = ncol(input_seurat)*8*1e-6; DoubletRate # 按每增加1000个细胞，双细胞比率增加千分之8来计算
DoubletRate = DoubletRate/5 # 但是因为去除双细胞去得太多了，我就人为地把这个值，再除了5，即去除8%的细胞

# 估计双细胞比例，根据细胞亚群数量与DoubletRate进行计算。
nExp_poi <- round(DoubletRate*length(input_seurat$celltype)); nExp_poi  # 721 #最好提供celltype，而不是seurat_clusters。

# 估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)); nExp_poi.adj # 474

gc() # 清缓存
seurat_singlet <- doubletFinder_v3(input_seurat, PCs = 1:10, pN = 0.25, 
                                   pK = mpK, nExp = nExp_poi, reuse.pANN = F, sct = F)
gc()
seurat_singlet <- doubletFinder_v3(seurat_singlet, PCs = 1:10, pN = 0.25, 
                                   pK = mpK, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
gc() # 再清缓存

# 查看以下列名
colnames(seurat_singlet@meta.data)

# 根据列名进行选择

seurat_singlet$DF_hi.lo <- seurat_singlet$DF.classifications_0.25_0.005_721
seurat_singlet$DF_hi.lo[which(seurat_singlet$DF_hi.lo == "Doublet" & seurat_singlet$DF.classifications_0.25_0.005_474 == "Singlet")] <- "Doublet-Low Confidience"
seurat_singlet$DF_hi.lo[which(seurat_singlet$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"

seurat_singlet$DF_hi.lo <- seurat_singlet$DF.classifications_0.25_0.005_474

table(seurat_singlet$DF_hi.lo)

b1 <- DimPlot(seurat_singlet, reduction = 'umap', group.by ="DF_hi.lo") + coord_equal(ratio = 1) 
```
#### 保存文件的时候切记不要覆盖旧文件 #### 
```r
ggplot2::ggsave(paste0(path, "UMAP_seurat_sample_1.pdf"), plot = b1,
                height = 5, width = 7, dpi = 300, limitsize = FALSE)

saveRDS(seurat_singlet, paste0(path, "seurat_sample_1.rds"))
```
---
## 9.合并样本 ####

#### 设置工作路径 ####
```r
readpath = 'F:/R work/mmbrain/'
path = 'F:/R work/mmbrain/results/'
```
#### 读取两个样本 #### 
```r
seurat_sample_1 <- readRDS(paste0(path, "seurat_sample_1.rds"))
seurat_sample_2 <- readRDS(paste0(path, "seurat_sample_2.rds"))

table(seurat_sample_1$DF_hi.lo)
Idents(seurat_sample_1) <- seurat_sample_1$DF_hi.lo
seurat_sample_1 <- subset(seurat_sample_1, idents = 'Singlet')

table(seurat_sample_2$DF_hi.lo)
Idents(seurat_sample_2) <- seurat_sample_2$DF_hi.lo
seurat_sample_2 <- subset(seurat_sample_2, idents = 'Singlet')
```
#### 合并两个样本 #### 
```r
seurat_integrated <- merge(seurat_sample_1, y = c(seurat_sample_2))
```
#### 合并完，又做一遍整合分析的流程 #### 
```r
seurat_integrated <- Seurat::NormalizeData(seurat_integrated, verbose = FALSE, normalization.method = 'LogNormalize')
seurat_integrated <- FindVariableFeatures(seurat_integrated, selection.method = "vst", nfeatures = 2000)
seurat_integrated <- ScaleData(seurat_integrated) 
seurat_integrated <- RunPCA(seurat_integrated, features = VariableFeatures(seurat_integrated))
seurat_integrated <- RunHarmony(seurat_integrated, 'orig.ident')
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20, reduction = 'harmony')
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:20, reduction = 'harmony')
seurat_integrated <- FindClusters(seurat_integrated, resolution = seq(from = 0.4,by = 0.2, length = 3))
```
#### 直接把celltype信息赋值给Idents #### 
```r
Idents(seurat_integrated) <- seurat_integrated$celltype
```
#### 画图UMAP ####
```r
p1 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1) 

ggplot2::ggsave(paste0(path, "UMAP_3_Singlet.pdf"), plot = p1, 
                height = 5, width = 7, dpi = 300, limitsize = FALSE)
```
#### 画图UMAP 分开样本 ####
```r
p2 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = TRUE, split.by = 'sample',
              label.size = 3) + coord_equal(ratio = 1) 

ggplot2::ggsave(paste0(path, "UMAP_split3_Singlet.pdf"), plot = p2, 
                height = 5, width = 9, dpi = 300, limitsize = FALSE)
```
#### 保存，以后读入这个文件进行下游分析即可 ####
```r
saveRDS(seurat_integrated, paste0(path, "seurat_integrated_2.rds"))
```
---
## 10.画图 ####

#### 桑基图柱形图饼图 统计细胞在样本之间的比例 ####
```r
# 细胞比例

Cellnum <- table(Idents(seurat_integrated), seurat_integrated$sample) 
Cellratio <- prop.table(Cellnum, margin = 2)*100#计算各组样本不同细胞群比例

# savedata <- as.data.frame.array(Cellnum)
# savedata <- as.data.frame.array(Cellratio)
#write.csv(savedata, paste0(path, 'Cellratio_celltype_count.csv'), row.names = T)

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))

c1 <- ggplot(Cellratio,
             aes(x = Var2, stratum = Var1, alluvium = Var1,
                 y = Freq,
                 fill = Var1, label = sprintf("%.1f%%", Freq))) +
  geom_flow(width = 0.7) +
  geom_stratum(alpha = 0.9, width = 0.7, colour = NA) +
  geom_text(stat = "stratum", size = 5, color="white") +
  labs(x = "Sample", y = "Ratio", fill = "celltype") +
  theme_classic() + 
  theme(panel.border = element_rect(fill=NA, color="white", size=0.5, linetype="solid"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")

# 细胞数量

Cellnum2 <- as.data.frame(Cellnum)
colourCount = length(unique(Cellnum2$Var1))

c2 <- ggplot(Cellnum2,
             aes(x = Var2, stratum = Var1, alluvium = Var1,
                 y = Freq,
                 fill = Var1, label = Freq)) +
  geom_bar(stat = "identity") + # 柱形图
  geom_text(position = position_stack(vjust = 0.8), color="white") +  # 柱形图
  
  # geom_flow(width = 0.7) +  # 桑基图
  # geom_stratum(alpha = 0.9, width = 0.7, colour = NA) + # 桑基图
  # geom_text(stat = "stratum", size = 5, color="white") + # 桑基图
  
  labs(x = "Sample", y = "Ratio", fill = "celltype") +
  theme_classic() + 
  theme(panel.border = element_rect(fill=NA, color="white", size=0.5, linetype="solid"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")

p <- plot_grid(c1, c2, nrow = 1)

ggplot2::ggsave(paste0(path, "细胞比例柱形图.pdf"), plot = p, 
                height = 9, width = 9, dpi = 300, limitsize = FALSE)


# 圆环图
data <- lapply(split(Cellnum2, Cellnum2$Var2), function(x){x <- x; return(x)})

data <- lapply(data, function(x){
  x$percent <- paste0(round(x$Freq, 4) / 100, '%')
  x$name <- paste0(x$Var1, ' (', x$percent, ')')
  return(x)
})

plot <- list()

for (i in 1:length(data)) {
  
  plot[[i]] <- ggdonutchart(data[[i]], 'Freq',
                            label = 'percent',
                            fill = 'Var1', 
                            color = "white", 
                            lab.font = c(3, "bold", "black"), 
                            font.family = "serif",
                            lab.pos = "in") +
    labs(x = NULL, y = NULL, fill = paste0(names(data)[i])) +
    theme_classic() + 
    theme(axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_blank(),
          text = element_text(size = 12,face = 'bold', family = "serif"),
          legend.position = 'left',
          plot.caption = element_text(hjust = 0.5, vjust = 7, size = 13), 
          plot.background = element_rect(fill = "transparent", 
                                         colour = NA)) +
    guides(fill = guide_legend(reverse = F))
}

cowplot::plot_grid(plotlist = plot, nrow = 1)

ggplot2::ggsave(paste0(path, '细胞比例饼图.pdf'),
                height = 7, width = 9, dpi = 300, limitsize = FALSE)

```

#### 通过FindAllMarkers去查找每个细胞群高表达的基因 ####
```r
# 设置min.pct = 0.25参数过滤掉那些在25%以下细胞中检测到的基因 
markers_label <- FindAllMarkers(seurat_integrated,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.70)

write.csv(markers_label, paste0(path, "markers_label.csv"), row.names = F)
```

####  热图 画前10个markers genes 的DoHeatmap ####
```r
all.genes <- rownames(seurat_integrated)
seurat_integrated <- ScaleData(seurat_integrated, features = all.genes)

top_10 <- as.data.frame(markers_label %>% group_by(cluster) %>%
                          top_n(n = 10, wt = avg_log2FC))

p <- DoHeatmap(seurat_integrated, features = top_10$gene, label = T, assay = "RNA", identity.legend = F) +
  scale_fill_gradientn(colors = c('white', 'grey', 'firebrick3'))

ggplot2::ggsave(paste0(path, "DoHeatmap_top10_label.pdf"), plot = p,
                height = 12, width = 20, dpi = 300, limitsize = FALSE)
```

#### Dotplot图 画marker基因dotplot ####
```r
DotPlot(object = seurat_integrated, features = marker, cols = c('white','firebrick')) + coord_flip()

# 或者用这个函数美化
scRNAtoolVis::jjDotPlot(object = seurat_integrated, 
                        id = 'celltype',
                        xtree = F,ytree = F,
                        gene = markers_cluster,
                        rescale = T,
                        aesGroName = 'sample',
                        #point.geom = F,
                        #tile.geom = T,
                        dot.col = c('white','firebrick'),
                        rescale.min = -2,
                        rescale.max = 2,
                        midpoint = 0) + coord_flip()

ggplot2::ggsave(paste0(path, "jjDotPlot.pdf"), 
                height = 6, width = 6, dpi = 300, limitsize = FALSE)
```

#### 统计 细胞群比例，umi比例 ####
```r
metadata <- seurat_integrated@meta.data

table(metadata$celltype)

t1 <- data.frame(table(metadata$celltype))
colnames(t1) <- c('celltype', 'No.nuclei')
metadata <- merge(metadata, t1, 'celltype')
metadata$a1 <- metadata$nCount_RNA/1000
metadata$b1 <- metadata$nFeature_RNA/1000

num = length(table(metadata$celltype))

a1 <- ggplot(t1, aes(x = celltype, y = No.nuclei, fill = celltype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = No.nuclei), vjust = -1, hjust = 0.5, color = "black") +
  xlab("Cell Type") +
  ylab("No.Nuclei") +
  coord_flip() + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none'); a1

a2 <- ggplot(metadata, aes(y = celltype, x = a1, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, color = 'black') +
  ylab("Cell Type") +
  xlab("No.UMIs per nuclei (x1000)") +
  scale_fill_discrete() +
  guides(fill = guide_legend(reverse = T)) + theme_classic() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none'); a2

a3 <- ggplot(metadata, aes(y = celltype, x = b1, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, color = 'black') +
  ylab("Cell Type") +
  xlab("No.genes per nuclei (x1000)") +
  scale_fill_discrete() +
  scale_x_continuous(breaks = c(0:6))+
  guides(fill = guide_legend(reverse = T)) + theme_classic() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none'); a3

markers_cluster <- c("Nrxn3", "Gja1","Flt1","Dnah11","Ctss", "Mog", "Pdgfra")

a4 <- Seurat::VlnPlot(seurat_integrated, split.by = 'celltype', 
                      markers_cluster, stack=T, pt.size = 0, flip = F)+
  ylab(NULL) +
  guides(fill = guide_legend(reverse = T)) + 
  scale_fill_manual(values = hue_pal()(num)); a4

ggarrange(ggarrange(a1, a2, a3, ncol = 3),
          a4, nrow = 2) 

ggplot2::ggsave(paste0(path, "long_merge.pdf"),
                height = 8, width = 8, dpi = 300, limitsize = FALSE)
```                    
 
#### UMAP图 画基因表达丰度图 ####
```r
input_gene <- c("Nrxn3", "Gja1","Flt1","Dnah11","Ctss", "Mog", "Pdgfra")

FeaturePlot(seurat_integrated, input_gene, split.by = 'sample') & theme(legend.position = "right")

# 或者用这个函数美化
scRNAtoolVis::featureCornerAxes(object = seurat_integrated, 
                                reduction = 'umap', pSize = 0.01,
                                features = input_gene, groupFacet = 'sample', axes  = 'one', ) + 
  coord_equal() +
  theme_bw() + 
  scale_colour_gradientn(colours = c('grey', 'blue'))

ggplot2::ggsave(paste0(path, "FeaturePlot.pdf"), 
                height = 25, width = 10, dpi = 300, limitsize = FALSE)
```              

#### 小提琴图 基因在不同细胞亚群中，样本之间的表达差异 ####
```r
Seurat::VlnPlot(seurat_integrated, split.by = 'sample', group.by = 'celltype',
                input_gene, stack=T, pt.size = 0, flip = T) +
  xlab(NULL) +
  guides(fill = guide_legend(reverse = T))

# 使用MySeuratWrappers::VlnPlot()会得到不一样的效果。。

```
#### 箱线图 基因在指定的细胞亚群中，样本之间的表达差异 添加p值 ####
```r
subsets_cell <- subset(seurat_integrated, ident = 'Endothelial_cell')

my_comparisons <- list(c("A_con", "B_tre"))

ps <- ggscplot(object = subsets_cell,
               features = input_gene,
               mapping = aes(x = sample, y = value)) +
  geom_boxplot(aes(fill = sample), outlier.colour = "grey90") +
  facet_feature(strip.col = NULL)+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylim(0, 8) +
  stat_summary(fun = mean, geom = "line", linetype = 2,
               aes(group = 0), color = "#E63863", size = 0.7) +
  stat_compare_means(method = "t.test", # t.test
                     comparisons = my_comparisons, 
                     line.size = 1, 
                     #step.increase = 0.5, 
                     label.y = c(5,7,6))
```
