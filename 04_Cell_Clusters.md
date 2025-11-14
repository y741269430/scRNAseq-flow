# 04_Cell_Clusters    

1. 输入准备
```r
# 设置工作目录，批量创建文件夹储存结果
setwd(r"{D:\R work\GSE171169_RAW\}")
# 加载R包
source('my_scRNAseq.R')

# 创建输出目录
if (!dir.exists("4_Cell_Clusters")) {
  dir.create("4_Cell_Clusters")
}
```
2. 读取去除双细胞后的矩阵
```r
singlet_files <- list.files(path = "2_DoubletFinder", 
                            pattern = "seurat_singlet.rds$", 
                            recursive = TRUE, 
                            full.names = TRUE)

# 提取样本名作为列表名称
sample_names <- sapply(singlet_files, function(x) basename(dirname(x)))

# 批量读取
seurat_list <- lapply(singlet_files, readRDS)
names(seurat_list) <- sample_names

# 合并样本
seurat_merged <- merge(x = seurat_list[[1]], 
                       y = c(seurat_list[[2]],
                             seurat_list[[3]],
                             seurat_list[[4]]
                             ),
                       add.cell.id = sample_names)

seurat_merged$Sample <- seurat_merged$orig.ident

saveRDS(seurat_merged, '3_QC_stat/seurat_merged.rds')

rm(list = ls()); gc()
```
3. 降维聚类
```r
# 读取数据
seurat_merged <- readRDS('3_QC_stat/seurat_merged.rds')
seurat_merged@meta.data <- seurat_merged@meta.data[, c(1:8)]
seurat_merged$Sample <- seurat_merged$orig.ident

# 3.1 运行 ####
seurat_integrated <- seurat_merged %>%
  NormalizeData() %>%                                    # 默认使用"LogNormalize"方法，考虑测序深度差异并进行对数转换
  FindVariableFeatures() %>%                             # 识别在细胞间表达变化较大的基因，用于后续降维分析
  ScaleData() %>%                                        # 线性变换，使每个基因在所有细胞中的均值为0，方差为1，避免高表达基因主导分析
  RunPCA(dims = 1:10) %>%                                # 主成分分析，线性降维以保留关键生物学差异
  RunHarmony('orig.ident') %>% 
  RunUMAP(dims = 1:10, reduction = 'harmony') %>%        # 基于前20个主成分进行UMAP非线性降维，用于可视化细胞群体关系
  FindNeighbors(dims = 1:10, reduction = 'harmony') %>% 
  FindClusters(resolution = 0.2)

Idents(seurat_integrated) <- seurat_integrated$RNA_snn_res.0.2
saveRDS(seurat_integrated, '4_Cell_Clusters/seurat_integrated.rds')
```
4. UMAP 绘制
4.1 UMAP 多个样本整合
```r
# 4.1
# UMAP 多个样本整合
p1 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1); p1

ggplot2::ggsave("4_Cell_Clusters/01_UMAP_all_sample.pdf", plot = p1, 
                height = 5, width = 7, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("4_Cell_Clusters/01_UMAP_all_sample.png", plot = p1, 
                height = 5, width = 7, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/01_UMAP_all_sample.png" width="600" />    

4.2 UMAP 分开样本绘制
```r
# 4.2
# UMAP 分开样本绘制
p2 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = TRUE, split.by = 'Sample',
              label.size = 3) + coord_equal(ratio = 1); p2

ggplot2::ggsave("4_Cell_Clusters/02_UMAP_split_sample.pdf", plot = p2, 
                height = 5, width = 9, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("4_Cell_Clusters/02_UMAP_split_sample.png", plot = p2, 
                height = 5, width = 9, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/02_UMAP_split_sample.png" width="500" />    

5. 细胞亚群比例统计    
5.1 输出每个细胞亚群的细胞比例
```r
# 5.1
Cellnum <- table(Idents(seurat_integrated), seurat_integrated$Sample) 

Cellratio <- round(prop.table(Cellnum, margin = 2) * 100, 2)
Cellratio_save <- Cellratio %>% 
  as.data.frame.matrix() %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.2f%%", round(., 2)))) %>% 
  rownames_to_column("Clusters")

write.xlsx(Cellratio_save, '4_Cell_Clusters/Cellratio_clusters.xlsx', rownames = F)
write.table(Cellratio_save, '4_Cell_Clusters/Cellratio_clusters.txt', quote = F, sep = '\t', row.names = F)
#txt <- knitr::kable(Cellratio_save, format = "markdown", align = 'c')
#write.table(txt, "4_Cell_Clusters/Cellratio_clusters_markdown.txt", row.names = F, quote = F, col.names = F)
```

| Clusters | 05d_N1 | 05d_N2 | 14d_N1 | 14d_N2 |
|:--------:|:------:|:------:|:------:|:------:|
|    0     | 24.76% | 22.80% | 15.97% | 29.88% |
|    1     | 20.41% | 16.30% | 11.24% | 20.44% |
|    2     | 2.95%  | 3.28%  | 32.39% | 17.76% |
|    3     | 22.78% | 19.18% | 4.93%  | 3.74%  |
|    4     | 11.22% | 20.86% | 5.42%  | 2.98%  |
|    5     | 4.89%  | 4.15%  | 12.49% | 11.62% |
|    6     | 2.53%  | 2.84%  | 9.99%  | 4.10%  |
|    7     | 6.50%  | 4.87%  | 2.39%  | 2.14%  |
|    8     | 1.10%  | 1.06%  | 2.43%  | 4.10%  |
|    9     | 1.01%  | 1.06%  | 1.54%  | 2.14%  |
|    10    | 0.93%  | 1.25%  | 1.17%  | 1.11%  |
|    11    | 0.93%  | 2.34%  | 0.04%  | 0.00%  |

5.2 输出每个细胞亚群的细胞数量
```r
# 5.2
Cellnum <- table(Idents(seurat_integrated), seurat_integrated$Sample) 
Cellnum_save <- Cellnum %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("Clusters")

write.xlsx(Cellnum_save, '4_Cell_Clusters/Cellnum_clusters.xlsx', rownames = F)
write.table(Cellnum_save, '4_Cell_Clusters/Cellnum_clusters.txt', quote = F, sep = '\t', row.names = F)
#txt <- knitr::kable(Cellnum_save, format = "markdown", align = 'c')
#write.table(txt, "4_Cell_Clusters/Cellnum_clusters_markdown.txt", row.names = F, quote = F, col.names = F)
```

| Clusters | 05d_N1 | 05d_N2 | 14d_N1 | 14d_N2 |
|:--------:|:------:|:------:|:------:|:------:|
|    0     |  587   |  730   |  395   |  671   |
|    1     |  484   |  522   |  278   |  459   |
|    2     |   70   |  105   |  801   |  399   |
|    3     |  540   |  614   |  122   |   84   |
|    4     |  266   |  668   |  134   |   67   |
|    5     |  116   |  133   |  309   |  261   |
|    6     |   60   |   91   |  247   |   92   |
|    7     |  154   |  156   |   59   |   48   |
|    8     |   26   |   34   |   60   |   92   |
|    9     |   24   |   34   |   38   |   48   |
|    10    |   22   |   40   |   29   |   25   |
|    11    |   22   |   75   |   1    |   0    |

5.3 数据框转换
```r
# 5.3
Cellratio_df <- as.data.frame(Cellratio)
Cellnum_df <- as.data.frame(Cellnum)

colnames(Cellratio_df) <- c("CellType", "Sample", "Percentage")
colnames(Cellnum_df) <- c("CellType", "Sample", "Count")

colourCount <- length(unique(Cellratio_df$CellType))
color_palette <- scales::hue_pal()(colourCount)
```

5.4 桑基图+柱状图 展示细胞亚群的比例与数量
```r
# 5.4
c1 <- ggplot(Cellratio_df,
             aes(x = Sample, stratum = CellType, alluvium = CellType,
                 y = Percentage,
                 fill = CellType, 
                 label = sprintf("%.1f%%", Percentage))) +
  geom_flow(width = 0.7, alpha = 0.7) +
  geom_stratum(alpha = 0.9, width = 0.7) +
  geom_text(stat = "stratum", size = 3, color = "white", fontface = "bold") +
  labs(x = "", y = "Percentage (%)", fill = "Cell Type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", size = 0.5)) + 
  coord_equal(ratio = 0.05) +    # 保持图像xy轴的比例（如果输出的图比例不合适就把这句代码注释掉）
  scale_fill_manual(values = color_palette)

c2 <- ggplot(Cellnum_df,
             aes(x = Sample, y = Count, fill = CellType)) +
  geom_col(position = "stack") +
  geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3, fontface = "bold") +
  labs(x = "", y = "Cell Count", fill = "Cell Type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "black", size = 0.5)) +
  scale_fill_manual(values = color_palette) + 
  coord_equal(ratio = 0.0015) +   # 保持图像xy轴的比例（如果输出的图比例不合适就把这句代码注释掉）
  scale_y_continuous(labels = scales::comma)

p_combined <- plot_grid(c1, c2, nrow = 1, labels = "AUTO")
ggsave("4_Cell_Clusters/03_Cell_proportion_stratum.pdf", plot = p_combined, height = 6, width = 12, dpi = 300)
ggsave("4_Cell_Clusters/03_Cell_proportion_stratum.png", plot = p_combined, height = 6, width = 12, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/03_Cell_proportion_stratum.png" width="600" />    

5.5 细胞亚群比例圆环图
```r
# 5.5
create_donut_plots <- function(cell_data) {
  data_list <- split(cell_data, cell_data$Sample)
  create_single_donut <- function(data, sample_name) {
    data <- data[order(data$Percentage, decreasing = TRUE), ]
    ggdonutchart(data, "Percentage",
                 label = "Percentage", 
                 fill = "CellType", 
                 color = "white", 
                 lab.font = c(3, "bold", "black")) +
      labs(title = sample_name, fill = "Cell Type") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            legend.position = "bottom",
            legend.text = element_text(size = 8)) +
      guides(fill = guide_legend(ncol = 1))
  }

  plots <- lapply(names(data_list), function(sample_name) {
    create_single_donut(data_list[[sample_name]], sample_name)
    })
  
  combined_plot <- plot_grid(plotlist = plots, nrow = 1)

  ggsave("4_Cell_Clusters/04_Cell_proportion_donut.pdf", plot = combined_plot, height = 5, width = 8, dpi = 300)
  ggsave("4_Cell_Clusters/04_Cell_proportion_donut.png", plot = combined_plot, height = 5, width = 8, dpi = 300)
  
  return(combined_plot)
}

donut_plot <- create_donut_plots(Cellratio_df)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/04_Cell_proportion_donut.png" width="500" />    

5.6 柱形图 展示细胞亚群的比例与数量
```r
# 5.6
colourCount <- length(unique(Cellratio_df$Sample))

# 细胞亚群百分比柱形图
c1 <- ggplot(Cellratio_df, aes(x = Percentage, y = CellType, fill = Sample, 
                               label = sprintf("%.1f%%", Percentage))) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_text(position = position_dodge(0.9), size = 2.5, hjust = 0.1, angle = 0) + 
  labs(x = "Cell Ratio (%)", y = "Cell Type", fill = "Sample") + 
  theme_classic() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5)
  ) +
  scale_fill_manual(values = hue_pal()(colourCount))

# 细胞亚群数量柱形图
c2 <- ggplot(Cellnum_df, aes(x = Count, y = CellType, fill = Sample, label = Count)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_text(position = position_dodge(0.9), size = 2.5, hjust = 0.1, angle = 0) + 
  labs(x = "Cell Number", y = "Cell Type", fill = "Sample") +
  theme_classic() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5)
  ) +
  scale_fill_manual(values = hue_pal()(colourCount))

p <- plot_grid(c1, c2, nrow = 1, labels = "AUTO", label_size = 14)

ggsave("4_Cell_Clusters/05_Cell_proportion_barplot.pdf", plot = p, height = 6, width = 15, dpi = 300)
ggsave("4_Cell_Clusters/05_Cell_proportion_barplot.png", plot = p, height = 6, width = 15, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/05_Cell_proportion_barplot.png" width="600" />    

6. 通过FindAllMarkers去查找每个细胞亚群高表达的基因    
6.1 查找每个细胞亚群高表达的基因
```r
# 6.1
Cluster_markers <- FindAllMarkers(seurat_integrated,
                                  min.pct = 0.1,          # 设置min.pct = 0.1参数代表在细胞亚群中，基因在10%以上的细胞中有表达
                                  logfc.threshold = 0.25, # 设置logfc.threshold = 0.25过滤掉那些在不同细胞亚群之间平均表达的差异倍数低于0.25的基因
                                  only.pos = TRUE) %>% dplyr::select(cluster, everything())

write.xlsx(Cluster_markers, "4_Cell_Clusters/Cluster_markers.xlsx", rowNames = F)
write.table(Cluster_markers, "4_Cell_Clusters/Cluster_markers.txt", row.names = F, quote = F, sep = '\t')

head(Cluster_markers)
```

| cluster |   p_val   | avg_log2FC | pct.1 | pct.2 | p_val_adj |     gene      |
|:-------:|:---------:|:----------:|:-----:|:-----:|:---------:|:-------------:|
|    0    | 0.0000000 | 3.3690340  | 0.999 | 0.745 | 0.0000000 |     Lyz2      |
|    0    | 0.0000000 | 2.8526944  | 0.509 | 0.058 | 0.0000000 |     Chil3     |
|    0    | 0.0000000 | 2.7074128  | 0.879 | 0.194 | 0.0000000 |     Tgfbi     |
|    0    | 0.0000000 | 2.6115775  | 0.614 | 0.042 | 0.0000000 |     Thbs1     |

6.2 热图 可视化每个细胞亚群top10高差异倍数基因
```r
# 6.2
all.genes <- rownames(seurat_integrated)
#seurat_integrated <- ScaleData(seurat_integrated, features = all.genes)

top_10 <- as.data.frame(Cluster_markers %>% 
                          group_by(cluster) %>%
                          top_n(n = 10, wt = avg_log2FC))

write.xlsx(Cluster_markers, "4_Cell_Clusters/Cluster_markers_top_10.xlsx", rowNames = F)
write.table(Cluster_markers, "4_Cell_Clusters/Cluster_markers_top_10.txt", row.names = F, quote = F, sep = '\t')

p <- MySeuratWrappers::DoHeatmap(seurat_integrated, 
                                 features = top_10$gene, 
                                 label = T, 
                                 #assay = "RNA", 
                                 identity.legend = F, 
                                 size = 3.5) +
  scale_fill_gradientn(colors = c('white', 'grey', 'firebrick3'))

ggsave("4_Cell_Clusters/06_DoHeatmap_top10_label.pdf", plot = p, height = 15, width = 20, dpi = 300, limitsize = F)
ggsave("4_Cell_Clusters/06_DoHeatmap_top10_label.png", plot = p, height = 15, width = 20, dpi = 300, limitsize = F)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/06_DoHeatmap_top10_label.png" width="500" />    

6.3 Dotplot图 画marker基因散点图
```r
# 6.3
top2_genes <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 2) %>%
  pull(gene) %>%
  unique()

DotPlot(object = seurat_integrated, features = top2_genes, cols = c('white','firebrick')) + 
  coord_flip() # 翻转坐标轴

ggsave("4_Cell_Clusters/07_DotPlot_top2.pdf", height = 6, width = 6, dpi = 300, limitsize = FALSE)
ggsave("4_Cell_Clusters/07_DotPlot_top2.png", height = 6, width = 6, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/07_DotPlot_top2.png" width="400" />    

6.4 jjDotPlot图 画marker基因散点图
```r
# 6.4
scRNAtoolVis::jjDotPlot(object = seurat_integrated, 
                        gene = top2_genes,
                        id = 'seurat_clusters',
                        xtree = F,
                        ytree = F,
                        rescale = T,
                        #aesGroName = 'Sample',
                        #point.geom = F,
                        #tile.geom = T,
                        dot.col = c('white','firebrick'),
                        rescale.min = -2,
                        rescale.max = 2,
                        midpoint = 0) + 
  coord_flip() # 翻转坐标轴

ggsave("4_Cell_Clusters/08_jjDotPlot_top2.pdf", height = 6, width = 7, dpi = 300, limitsize = FALSE)
ggsave("4_Cell_Clusters/08_jjDotPlot_top2.png", height = 6, width = 7, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/08_jjDotPlot_top2.png" width="400" />    

---
```r
fs::dir_tree("4_Cell_Clusters", recurse = 2)
```
```bash
4_Cell_Clusters
├── 01_UMAP_all_sample.pdf
├── 01_UMAP_all_sample.png
├── 02_UMAP_split_sample.pdf
├── 02_UMAP_split_sample.png
├── 03_Cell_proportion_stratum.pdf
├── 03_Cell_proportion_stratum.png
├── 04_Cell_proportion_donut.pdf
├── 04_Cell_proportion_donut.png
├── 05_Cell_proportion_barplot.pdf
├── 05_Cell_proportion_barplot.png
├── 06_DoHeatmap_top10_label.pdf
├── 06_DoHeatmap_top10_label.png
├── 07_DotPlot_top2.pdf
├── 07_DotPlot_top2.png
├── 08_jjDotPlot_top2.pdf
├── 08_jjDotPlot_top2.png
├── Cellnum_clusters.txt
├── Cellnum_clusters.xlsx
├── Cellnum_clusters_markdown.txt
├── Cellratio_clusters.txt
├── Cellratio_clusters.xlsx
├── Cellratio_clusters_markdown.txt
├── cluster_markers.txt
├── cluster_markers.xlsx
└── seurat_integrated.rds
```
### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-10
- 修改日期：2025-11-10











