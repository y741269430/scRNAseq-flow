# 04_Cell_Clusters    

### 1. 输入准备
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
### 2. 读取去除双细胞后的矩阵
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
### 3. 降维聚类
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
  RunPCA(dims = 1:20) %>%                                # 主成分分析，线性降维以保留关键生物学差异
  RunHarmony('orig.ident') %>% 
  RunUMAP(dims = 1:20, reduction = 'harmony') %>%        # 基于前20个主成分进行UMAP非线性降维，用于可视化细胞群体关系
  RunTSNE(dims = 1:20, reduction = 'harmony') %>% 
  FindNeighbors(dims = 1:20, reduction = 'harmony') %>% 
  FindClusters(resolution = 0.8)

Idents(seurat_integrated) <- seurat_integrated$RNA_snn_res.0.8
table(Idents(seurat_integrated))

saveRDS(seurat_integrated, '4_Cell_Clusters/seurat_integrated_tmp.rds')
```
### 4. UMAP 绘制
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

### 5. 细胞亚群比例统计    
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
|    0     | 18.22% | 15.58% | 4.49%  | 3.12%  |
|    1     | 8.35%  | 7.87%  | 11.04% | 17.10% |
|    2     | 14.13% | 11.21% | 4.37%  | 12.51% |
|    3     | 14.80% | 10.49% | 5.90%  | 10.95% |
|    4     | 10.92% | 20.42% | 4.41%  | 2.49%  |
|    5     | 1.69%  | 1.72%  | 20.02% | 12.15% |
|    6     | 4.98%  | 4.18%  | 12.37% | 11.80% |
|    7     | 2.57%  | 2.90%  | 9.95%  | 4.05%  |
|    8     | 0.67%  | 0.94%  | 10.43% | 4.76%  |
|    9     | 6.54%  | 4.18%  | 2.43%  | 2.18%  |
|    10    | 2.24%  | 1.81%  | 4.16%  | 7.75%  |
|    11    | 1.10%  | 1.09%  | 2.43%  | 4.19%  |
|    12    | 2.53%  | 3.28%  | 0.24%  | 0.49%  |
|    13    | 1.01%  | 1.37%  | 3.23%  | 1.42%  |
|    14    | 1.86%  | 2.87%  | 0.57%  | 0.45%  |
|    15    | 2.70%  | 2.53%  | 0.28%  | 0.04%  |
|    16    | 0.97%  | 1.53%  | 1.33%  | 1.25%  |
|    17    | 0.97%  | 1.03%  | 1.33%  | 1.87%  |
|    18    | 1.05%  | 3.15%  | 0.04%  | 0.00%  |
|    19    | 2.02%  | 1.19%  | 0.24%  | 0.45%  |
|    20    | 0.67%  | 0.62%  | 0.73%  | 0.98%  |

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
|    0     |  432   |  499   |  111   |   70   |
|    1     |  198   |  252   |  273   |  384   |
|    2     |  335   |  359   |  108   |  281   |
|    3     |  351   |  336   |  146   |  246   |
|    4     |  259   |  654   |  109   |   56   |
|    5     |   40   |   55   |  495   |  273   |
|    6     |  118   |  134   |  306   |  265   |
|    7     |   61   |   93   |  246   |   91   |
|    8     |   16   |   30   |  258   |  107   |
|    9     |  155   |  134   |   60   |   49   |
|    10    |   53   |   58   |  103   |  174   |
|    11    |   26   |   35   |   60   |   94   |
|    12    |   60   |  105   |   6    |   11   |
|    13    |   24   |   44   |   80   |   32   |
|    14    |   44   |   92   |   14   |   10   |
|    15    |   64   |   81   |   7    |   1    |
|    16    |   23   |   49   |   33   |   28   |
|    17    |   23   |   33   |   33   |   42   |
|    18    |   25   |  101   |   1    |   0    |
|    19    |   48   |   38   |   6    |   10   |
|    20    |   16   |   20   |   18   |   22   |

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
create_donut_plots <- function(Cellratio_df) {
  data_list <- split(Cellratio_df, Cellratio_df$Sample)
  
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
            legend.position = "none", legend.text = element_text(size = 8)) +
      guides(fill = guide_legend(nrow = 1, byrow = T))
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

### 6. 通过FindAllMarkers去查找每个细胞亚群高表达的基因    
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
top1_genes <- Cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1) %>%
  pull(gene) %>%
  unique()

DotPlot(object = seurat_integrated, features = top1_genes, cols = c('white','firebrick')) + 
  coord_flip() # 翻转坐标轴

ggsave("4_Cell_Clusters/07_DotPlot_top1.pdf", height = 6, width = 8, dpi = 300, limitsize = FALSE)
ggsave("4_Cell_Clusters/07_DotPlot_top1.png", height = 6, width = 8, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/07_DotPlot_top1.png" width="400" />    

6.4 jjDotPlot图 画marker基因散点图
```r
# 6.4
scRNAtoolVis::jjDotPlot(object = seurat_integrated, 
                        gene = top1_genes,
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

ggsave("4_Cell_Clusters/08_jjDotPlot_top1.pdf", height = 6, width = 7, dpi = 300, limitsize = FALSE)
ggsave("4_Cell_Clusters/08_jjDotPlot_top1.png", height = 6, width = 7, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/08_jjDotPlot_top1.png" width="400" />    

### 7. Pearson相关性热图 使用每个细胞亚群的Top 5 marker基因
```r
# 6.5
Cluster_markers <- read.xlsx('4_Cell_Clusters/Cluster_markers.xlsx')

# 选择每个cluster的top marker基因
top_markers <- Cluster_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup()

# 获取这些基因在所有cluster中的表达矩阵
expression_matrix <- AverageExpression(seurat_integrated,
                                       assays = "RNA",
                                       features = unique(top_markers$gene),
                                       slot = "data")$RNA

# 转置矩阵
cluster_expression_profiles <- t(expression_matrix)

# 计算cluster表达谱之间的相关性
corr_matrix <- cor(cluster_expression_profiles, method = "pearson")

# 添加细胞亚群注释信息
tag <- top_markers[top_markers$gene %in% rownames(corr_matrix), ]

# 获取cluster信息并保持factor顺序
cluster_levels <- levels(Idents(seurat_integrated))  # 获取factor的levels顺序
cluster_info <- factor(tag$cluster, levels = cluster_levels)

# 创建注释信息
annotation_col <- data.frame(Seurat_Cluster = cluster_info)
rownames(annotation_col) <- rownames(corr_matrix)

# 按照factor levels顺序设置颜色
n_clusters <- length(cluster_levels)
cluster_colors <- colorRampPalette(brewer.pal(8, "Set1"))(n_clusters)
names(cluster_colors) <- cluster_levels  # 按照levels顺序命名

annotation_colors <- list(Seurat_Cluster = cluster_colors)

# 绘制热图（添加anno色块）
plot <- pheatmap(corr_matrix,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = F,
         display_numbers = F,  # 显示相关系数
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 6,
         annotation_col = annotation_col,      # 列注释
         annotation_colors = annotation_colors, # 注释颜色
         main = "Cluster Correlation Based on Marker Gene Expression Profiles",
         fontsize = 10,
         border_color = NA,
         fontsize_row = 8,
         fontsize_col = 8,
         silent = T)

plot_grid(plot$gtable)

ggsave("4_Cell_Clusters/09_Pearson_Cluster_markers.pdf", plot = plot, height = 9, width = 10, dpi = 300, limitsize = FALSE)
ggsave("4_Cell_Clusters/09_Pearson_Cluster_markers.png", plot = plot, height = 9, width = 10, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/4_Cell_Clusters/09_Pearson_Cluster_markers.png" width="400" />    

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
├── 07_DotPlot_top1.pdf
├── 07_DotPlot_top1.png
├── 08_jjDotPlot_top1.pdf
├── 08_jjDotPlot_top1.png
├── 09_Pearson_Cluster_markers.pdf
├── 09_Pearson_Cluster_markers.png
├── Cellnum_clusters.txt
├── Cellnum_clusters.xlsx
├── Cellratio_clusters.txt
├── Cellratio_clusters.xlsx
├── Cluster_markers.txt
├── Cluster_markers.xlsx
├── Cluster_markers_top_10.txt
├── Cluster_markers_top_10.xlsx
└── seurat_integrated_tmp.rds
```
### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-10
- 修改日期：2025-11-29



















