# 02_DoubletFinder

## 一、DoubletFinder 总览
参考来源：[chris-mcginnis-ucsf/DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder?tab=readme-ov-file#pk-selection)     
DoubletFinder 的分析流程可分为以下四个步骤：

1.  **生成人工双细胞**：基于现有的单细胞RNA测序数据生成人工双细胞。
2.  **数据预处理**：对真实数据与人工双细胞的合并数据集进行预处理。
3.  **计算 pANN**：执行主成分分析（PCA），并利用主成分距离矩阵计算每个细胞的人工最近邻比例（pANN）。
4.  **筛选双细胞**：根据预期的双细胞数量，对 pANN 值进行排序并设定阈值进行筛选。

DoubletFinder 包含以下关键参数：

- **`seurat`** 经过完整预处理的 Seurat 对象（即已依次执行 `NormalizeData`, `FindVariableGenes`, `ScaleData`, `RunPCA` 及 `RunTSNE` 或 `RunUMAP` 等操作后生成的对象）。
- **`PCs`** 具有统计学意义的主成分数量，以范围形式指定（例如 `PCs = 1:10`）。
- **`pN`** 生成人工双细胞的数量，以合并后真实-人工数据总量的比例表示。默认值为 25%，基于 Cell Hashing 和 Demuxlet 数据集的 pN-pK 参数扫描 ROC 分析表明，DoubletFinder 的性能在很大程度上不受 pN 值选择的影响 [McGinnis, C. S., Murrow, L. M., & Gartner, Z. J. (2019). ](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471219300730%3Fshowall%3Dtrue)
- **`pK`** 用于计算 pANN 时所采用的 PC 邻域大小，以合并后真实-人工数据总量的比例表示。此参数未设默认值，因 pK 需根据每个单细胞RNA测序数据集进行调整。最佳 pK 值应按下文所述策略进行估算。
- **`nExp`** 用于最终判定双细胞/单细胞的 pANN 阈值。该值可基于 10X 或 Drop-Seq 仪器的细胞上样密度进行估算，并根据同型双细胞的预估比例进行相应调整。

### 输入数据要求
-   **错误做法**：如果同时使用两个不同处理的样本，例如野生型（WT）和突变型（mutant）细胞系的整合后的数据运行DoubletFinder，软件将会生成由WT和mutant细胞混合的人工双细胞。这类双细胞在您的实际实验数据中是不可能存在的，因此会严重歪曲分析结果。
-   **正确做法**：对于将**单个样本**拆分到多个10X通道进行测序所产生的数据，运行DoubletFinder是完全可以的。

### 数据预处理工作流程

为确保最佳效果，请确保输入数据已清除低质量的细胞簇。推荐采用以下工作流程：

1.  **初始质控过滤**
    根据RNA的nUMI（独特分子标识符计数）手动设定阈值过滤原始基因表达矩阵。

2.  **标准预处理**
    使用标准工作流程对过滤后的数据进行预处理。
    -   通常包括：`NormalizeData`, `FindVariableGenes`, `ScaleData`, `RunPCA` 等步骤。

3.  **识别并移除低质量细胞簇**
    在初步聚类后，识别并移除具有以下一个或多个特征的细胞簇：
    -   **(A) 低RNA nUMI**
    -   **(B) 高线粒体基因读数百分比**
    -   **(C) 缺乏有信息量的标记基因**

4.  **最终运行DoubletFinder**
    **移除低质量细胞簇后，请重新执行完整的预处理流程（即上述第2步），然后再运行DoubletFinder。**

为了最大限度地提高DoubletFinder预测的准确性，我们致力于寻找一种不依赖于已知真实数据的评估指标，该指标应能指示出在Cell Hashing和Demuxlet数据中能使AUC（曲线下面积）最大化的pK值。均值-方差归一化的双峰系数（BCmvn）成功实现了这一目标。该系数的特点是，在能最大化AUC的最佳pK值处，它会呈现一个单一且易于识别的峰值。

<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/BCmvn.png" width="500" />    

---
## 二、单样本去除双细胞分析流程
读取上一部分析中的rds [01_scRNAseq_QC](https://github.com/y741269430/scRNAseq-flow/blob/main/01_scRNAseq_QC.md)    

1. 加载R包
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
#remotes::install_github('https://github.com/ekernf01/DoubletFinder', force = T)
#remotes::install_github("lyc-1995/MySeuratWrappers")
```
2. 计算双细胞率的函数
```r
get_multiplet_rate <- function(cellnums) {
  # 每500个细胞对应0.4%的双细胞率
  base_rate <- 0.004  # 0.4%
  cells_per_increment <- 500
  # 计算倍数
  increment <- cellnums / cells_per_increment
  # 计算双细胞率
  rate <- base_rate * increment
  # 最大不超过8%
  if (rate > 0.08) {
    rate <- 0.08
  }
  return(rate)
}
```
3. 配置环境并读取rds
```r
# 调整R内允许对象大小的限制（默认值是500*1024 ^ 2 = 500 Mb）
options(future.globals.maxSize = 500 * 1024 ^ 2)

# 设置工作目录，批量创建文件夹储存结果
setwd(r"{D:\R work\GSE171169_RAW\}")

seurat_filter <- readRDS("1_QC_Files/seurat_filter.rds")

```
4. 预处理 Seurat object    
```r
# 取一个样本进行示范
input_data <- seurat_filter[[1]]

seu_data <-  input_data %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 0.2)

table(seu_data$seurat_clusters)
```
5. 计算最优 pK (no ground-truth)    
```r
sweep.res.list <- paramSweep_v3(seu_data, PCs = 1:10, sct = F)
sweep.res.list <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.res.list)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])); mpK # mpK = 0.005
```
6. 可视化
```r
p <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) + 
  labs(x = 'pK', 
       y = 'BCmvn',
       title = 'Mean-variance normalized \nbimodality coefficient (BCmvn)',
       subtitle = paste0('Optimal pK = ', mpK)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "blue", size = 12)) +
  geom_point(color = 'red', size = 1.5) + 
  geom_line(color = 'blue', linewidth = 0.5); p

ggplot2::ggsave("2_DoubletFinder/01_BCmvn_distributions.pdf", plot = p, height = 4, width = 6, dpi = 300)
ggplot2::ggsave("2_DoubletFinder/01_BCmvn_distributions.png", plot = p, height = 4, width = 6, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/2_DoubletFinder/01_BCmvn_distributions.png" width="500" />    

7. 计算理论双细胞比例与同源双细胞比例
```r
annotations <- seu_data@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) # 同源双细胞比例估算

cellnums <- ncol(seu_data)
rate <- get_multiplet_rate(cellnums)
cat("Cell numbers:", cellnums, "-> Cell multiplet rate:", percent(rate, 0.001))

# 估计理论双细胞比例，根据seurat_clusters分群数量与DoubletRate进行计算。
nExp_poi <- round(rate*ncol(seu_data)); nExp_poi

# 估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)); nExp_poi.adj

# 检查列数是否是10列，计算双细胞时，会在metadata后面添加第11列，12列，13列
ncol(seu_data@meta.data)
```
8. 使用`doubletFinder_v3`计算双细胞数量
```r
seu_data <- doubletFinder_v3(seu_data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, 
                             reuse.pANN = F, sct = F)

## 检查第11列是否是doubletFinder输出的表头
colnames(seu_data@meta.data)[[11]]

seu_data <- doubletFinder_v3(seu_data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, 
                             reuse.pANN = colnames(seu_data@meta.data)[[11]], sct = F)

colnames(seu_data@meta.data)

table(seu_data@meta.data[,12]) # nExp_poi
table(seu_data@meta.data[,13]) # nExp_poi.adj

# 基于两种检测结果进行置信度分级
seu_data$DF_hi.lo <- ifelse(
  # 高置信度双细胞：两种方法都认为是双细胞
  seu_data@meta.data[,12] == "Doublet" & 
    seu_data@meta.data[,13] == "Doublet",
  "Doublet-High Confidence",
  
  ifelse(
    # 低置信度双细胞：只有一种方法认为是双细胞
    seu_data@meta.data[,12] == "Doublet" | 
      seu_data@meta.data[,13] == "Doublet",
    "Doublet-Low Confidence",
    
    # 单细胞：两种方法都认为是单细胞
    "Singlet"
  )
)
```
9. 可视化评估双细胞的结果
```r
p1 <- DimPlot(seu_data, reduction = 'umap', group.by ="seurat_clusters") + coord_equal(ratio = 1) 
p2 <- DimPlot(seu_data, reduction = 'umap', group.by ="DF_hi.lo") + coord_equal(ratio = 1) 
p <- plot_grid(p1, p2, nrow = 1, align = "hv", axis = "tblr")
p

ggplot2::ggsave("2_DoubletFinder/02_UMAP_DoubletFinder.pdf", plot = p, height = 6, width = 12, dpi = 300)
ggplot2::ggsave("2_DoubletFinder/02_UMAP_DoubletFinder.png", plot = p, height = 6, width = 12, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/2_DoubletFinder/02_UMAP_DoubletFinder.png" width="700" />    

```r
p <- Seurat::VlnPlot(seu_data, group.by = "DF_hi.lo", 
                     features = c("nCount_RNA", "nFeature_RNA"), 
                     pt.size = 0, ncol = 2) 

ggplot2::ggsave("2_DoubletFinder/03_Violin_DoubletFinder.pdf", plot = p, height = 5, width = 10, dpi = 300)
ggplot2::ggsave("2_DoubletFinder/03_Violin_DoubletFinder.png", plot = p, height = 5, width = 10, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/2_DoubletFinder/03_Violin_DoubletFinder.png" width="500" />    

10. 双细胞数量的统计表
```r
doubletFinder_res <- as.data.frame.matrix(t(table(seu_data$DF_hi.lo)))
doubletFinder_res$Cell_Discard_Rate <- percent(c(1 - c(doubletFinder_res[,3] / ncol(seu_data) )), 0.01)
doubletFinder_res$Cell_Retention_Rate <- percent(c(doubletFinder_res[,3] / ncol(seu_data) ), 0.01)
doubletFinder_res

write.table(doubletFinder_res, '2_DoubletFinder/singlet.txt', row.names = F, quote = F, sep = '\t')
write.xlsx(doubletFinder_res, '2_DoubletFinder/singlet.xlsx', rowNames = F)
# 生成Markdown表格
txt <- knitr::kable(doubletFinder_res, format = "markdown", align = 'c')
write.table(txt, '2_DoubletFinder/markdown.txt', row.names = F, quote = F, col.names = F)
```
11. 保存去除双细胞后的结果。  
```r
Idents(seu_data) <- seu_data$DF_hi.lo
seu_data_singlet <- subset(seu_data, idents = 'Singlet')
saveRDS(seu_data_singlet, "2_DoubletFinder/seurat_singlet.rds")
```
