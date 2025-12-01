# 06_Cell_Annotation_by_SingleR

配置环境
```r
# 设置工作目录，批量创建文件夹储存结果
setwd(r"{D:\R work\GSE171169_RAW\}")

source('my_scRNAseq.R')

# 创建输出目录
if (!dir.exists("6_Cell_Annotation_by_SingleR")) {
  dir.create("6_Cell_Annotation_by_SingleR")
}
```

### 1. 加载数据集
```r
seurat_integrated <- readRDS('4_Cell_Clusters/seurat_integrated.rds')
refdata1 <- get(load('ref_Mouse_all.RData'))
refdata2 <- get(load('ref_Mouse_imm.RData'))

table(refdata1$label.main)
table(refdata2$label.main)
```

### 2. 运行singler分析
```r
singler = SingleR(test= GetAssayData(seurat_integrated, slot = 'data') , 
                  ref = list(Mouse_all = refdata1,
                             Mouse_imm = refdata2) ,
                  labels = list(refdata1$label.main, 
                                refdata2$label.main) , 
                  clusters = seurat_integrated$seurat_clusters)

table(singler$labels)
```
绘图
```r
plotScoreHeatmap(singler, show_colnames = TRUE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/6_Cell_Annotation_by_SingleR/01_plotScoreHeatmap.png" width="500" />   

### 3. 将注释结果整合到单细胞结果里面
```r
singler_data <- data.frame(
  seurat_clusters = rownames(singler),
  SingleR_celltypes = singler$labels,
  stringsAsFactors = FALSE
)

cluster_match <- match(as.character(seurat_integrated$seurat_clusters), 
                       as.character(singler_data$seurat_clusters))

seurat_integrated$SingleR_celltypes <- singler_data$SingleR_celltypes[cluster_match]

if(any(is.na(seurat_integrated$SingleR_celltypes))) {
  seurat_integrated$SingleR_celltypes[is.na(seurat_integrated$SingleR_celltypes)] <- "Unassigned"
}

Idents(seurat_integrated) <- seurat_integrated$SingleR_celltypes

saveRDS(seurat_integrated, '6_Cell_Annotation_by_SingleR/seurat_integrated_SR.rds')
```
可视化
```r
p1 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = TRUE, 
              label.size = 3, group.by = 'seurat_clusters') + 
  coord_equal(ratio = 1) +
  ggtitle("Seurat Clusters")

p2 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = TRUE, 
              label.size = 3, group.by = 'SingleR_celltypes') + 
  coord_equal(ratio = 1) +
  ggtitle("SingleR Cell Type Annotation")

p <- plot_grid(p1, p2 ,nrow = 1)

ggplot2::ggsave("6_Cell_Annotation_by_SingleR/02_Cell_Annotaions_UMAP.pdf", plot = p,
                height = 7, width = 12, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("6_Cell_Annotation_by_SingleR/02_Cell_Annotaions_UMAP.png", plot = p,
                height = 7, width = 12, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/6_Cell_Annotation_by_SingleR/02_Cell_Annotaions_UMAP.png" width="600" />   

---
目录树     
```r
fs::dir_tree("6_Cell_Annotation_by_SingleR", recurse = 2)
```
```bash
6_Cell_Annotation_by_SingleR
├── 01_plotScoreHeatmap.pdf
├── 01_plotScoreHeatmap.png
├── 02_Cell_Annotaions_UMAP.pdf
├── 02_Cell_Annotaions_UMAP.png
└── seurat_integrated_SR.rds
```

### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-26
- 修改日期：2025-11-26
