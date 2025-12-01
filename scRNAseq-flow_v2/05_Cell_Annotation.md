# 05_Cell_Annotation

## 一、数据准备     
1.1 配置环境
```r
# 设置工作目录，批量创建文件夹储存结果
#setwd(r"{D:\R work\GSE171169_RAW\}")
setwd('/Users/mac/R work/GSE171169_RAW')

# 加载R包
source('my_scRNAseq.R')

# 创建输出目录
if (!dir.exists("5_Cell_Annotation")) {
  dir.create("5_Cell_Annotation")
}
```
1.2 读取矩阵并重新降维
```r
# 读取数据
seurat_integrated <- readRDS('4_Cell_Clusters/seurat_integrated_tmp.rds')
Idents(seurat_integrated) <- seurat_integrated$RNA_snn_res.0.8
```

## 二、细胞注释      
一般有两种方法，一是基于先验知识，即现有的细胞markers进行手动注释；二是使用自动化注释。      

### 方法一：
2.1 从cellmarkers官网或文献中获取细胞markers      
```r
# 手动注释，基于小提琴图以及UMAP图注释细胞：
table(Idents(seurat_integrated))

marker_ls <- list('ab_T_cells' = c('Trbc2', 'Cd3e', 'Cd28'),
                  'gd_T_cells' = c('Tgfbi', 'Ms4a6d', 'Ccl9'),
                  'NK_cells' = c('Nkg7', 'Ncr1', 'Ctla2a'),
                  'B_cells' = c('Il17a', 'Tcrg-C1', 'Serpinb1a'),
                  'Macrophages_1' = c('Cd79a', 'Ly6d', 'Cd79b'),
                  'Macrophages_2' = c('Cd209a', 'H2-Aa', 'Cd74'),
                  'Microglia_1' = c('Fcrls', 'Ctsl', 'Apoe'),
                  'Microglia_2' = c('S100a9', 'S100a8', 'Cxcl2'),
                  'Neutrophils' = c('Fscn1', 'Ccl22', 'Il4i1'),
                  "cDCs" = c('Tmem119', 'P2ry12'),
                  'pDCs' = c('Ccr9', 'Siglech', 'Bst2'))

vln_ls <- list()
for (i in 1:length(marker_ls)) {
  vln_ls[[i]] <- VlnPlot(seurat_integrated, features = marker_ls[[i]], stack = T, pt.size = 0) + 
    ggtitle(names(marker_ls)[i]) 
}

a1 <- plot_grid(plotlist = vln_ls, ncol = 2)

ggplot2::ggsave("5_Cell_Annotation/01_VlnPlot_marker.pdf", plot = a1,
                height = 30, width = 12, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("5_Cell_Annotation/01_VlnPlot_marker.png", plot = a1,
                height = 30, width = 12, dpi = 300, limitsize = FALSE)

Idents(seurat_integrated) <- seurat_integrated$RNA_snn_res.0.8

p1 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1)

p2 <- DimPlot(seurat_integrated,
              reduction = "tsne", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/01_VlnPlot_marker.png" width="400" />

2.2 手动注释     
```r
# 2.2 手动进行细胞注释
seurat_integrated <- RenameIdents(seurat_integrated,
                                  "0" = "Microglia_1",
                                  "1" = "gd_T_cells",
                                  "2" = "gd_T_cells",
                                  "3" = "Macrophages_2",
                                  "4" = "cDCs",
                                  "5" = "ab_T_cells",
                                  "6" = "NK_cells",
                                  "7" = "Macrophages_1",
                                  "8" = "B_cells",
                                  "9" = "Microglia_2",
                                  "10" = "Macrophages_2",
                                  "11" = "Neutrophils",
                                  "12" = "Macrophages_2",
                                  "13" = "ab_T_cells",
                                  "14" = "gd_T_cells",
                                  "15" = "Microglia_1",
                                  "16" = "pDCs",
                                  "17" = "ab_T_cells",
                                  "18" = "Microglia_2",
                                  "19" = "Microglia_1",
                                  "20" = "Macrophages_2")

# 细胞身份排序
Idents(seurat_integrated) <- factor(Idents(seurat_integrated), 
                                    levels = c('ab_T_cells', 'gd_T_cells', 'NK_cells', 'B_cells',
                                               'Macrophages_1', 'Macrophages_2', 'Microglia_1', 'Microglia_2',
                                               'Neutrophils', 'cDCs', 'pDCs'))

seurat_integrated$celltype <- Idents(seurat_integrated)

# 设置样本顺序
table(seurat_integrated$Sample)
seurat_integrated$Sample <- factor(seurat_integrated$Sample,
                                   levels=c("05d_N1", "05d_N2", "14d_N1", "14d_N2"))
```

## 三、细胞注释结果可视化      
3.1 UMAP/TSNE       
```r
p3 <- DimPlot(seurat_integrated,
              reduction = "umap", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1)

p4 <- DimPlot(seurat_integrated,
              reduction = "tsne", 
              label = T, 
              label.size = 3) + coord_equal(ratio = 1)

pU <- plot_grid(p1, p3, nrow = 1)
pT <- plot_grid(p2, p4, nrow = 1)

ggplot2::ggsave("5_Cell_Annotation/02_Cell_Annotaions_UMAP.pdf", plot = pU,
                height = 7, width = 12, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("5_Cell_Annotation/02_Cell_Annotaions_UMAP.png", plot = pU,
                height = 7, width = 12, dpi = 300, limitsize = FALSE)

ggplot2::ggsave("5_Cell_Annotation/03_Cell_Annotaions_TSNE.pdf", plot = pT,
                height = 7, width = 12, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("5_Cell_Annotation/03_Cell_Annotaions_TSNE.png", plot = pT,
                height = 7, width = 12, dpi = 300, limitsize = FALSE)

saveRDS(seurat_integrated, '5_Cell_Annotation/seurat_integrated_anno.rds')
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/02_Cell_Annotaions_UMAP.png" width="600" />

<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/03_Cell_Annotaions_TSNE.png" width="600" />


3.2 绘制热图
```r
# 3.2 绘制热图

input_genes <- as.character(unlist(marker_ls))

p <- MySeuratWrappers::DoHeatmap(seurat_integrated, 
                                 features = input_genes, 
                                 label = T, 
                                 #assay = "RNA", 
                                 identity.legend = F, 
                                 size = 3.5) +
  scale_fill_gradientn(colors = c('white', 'grey', 'firebrick3'))

ggsave("5_Cell_Annotation/04_DoHeatmap_anno_markers.pdf", plot = p, height = 7, width = 9, dpi = 300, limitsize = F)
ggsave("5_Cell_Annotation/04_DoHeatmap_anno_markers.png", plot = p, height = 7, width = 9, dpi = 300, limitsize = F)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/04_DoHeatmap_anno_markers.png" width="500" />

3.3 绘制点图
```r
# 3.3 绘制点图

DotPlot(object = seurat_integrated, features = input_genes, cols = c('white','firebrick')) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1) )

ggsave("5_Cell_Annotation/05_DotPlot_anno_markers.pdf", height = 8, width = 6, dpi = 300, limitsize = FALSE)
ggsave("5_Cell_Annotation/05_DotPlot_anno_markers.png", height = 8, width = 6, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/05_DotPlot_anno_markers.png" width="300" />

## 四、批量绘制 UMAP      
4.1 使用markers基因 批量绘制 UMAP        
```r
# 创建输出目录
if (!dir.exists("5_Cell_Annotation/UMAP_markers")) {
  dir.create("5_Cell_Annotation/UMAP_markers")
}

plot <- list()
for (i in 1:length(marker_ls)) {
  reduction = 'umap'  # tsne
  
  plot[[i]] <- FeaturePlot(seurat_integrated, marker_ls[[i]], ncol = 1, reduction = reduction) & 
    theme(legend.position = "right")
  ggsave(paste0("5_Cell_Annotation/UMAP_markers/UMAP_markers_", i, ".pdf"), plot = plot[[i]],
         height = 8, width = 4, dpi = 300, limitsize = FALSE)
  ggsave(paste0("5_Cell_Annotation/UMAP_markers/UMAP_markers_", i, ".png"), plot = plot[[i]],
         height = 8, width = 4, dpi = 300, limitsize = FALSE)
}
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/UMAP_markers/UMAP_markers_1.png" width="300" />

4.2 使用markers基因 批量绘制 UMAP SplitBySample          
```r
plot <- list()
for (i in 1:length(marker_ls)) {
  reduction = 'umap'  # tsne
  
  plot[[i]] <- FeaturePlot(seurat_integrated, marker_ls[[i]], split.by = 'Sample', reduction = reduction) & 
    theme(legend.position = "right")
  ggsave(paste0("5_Cell_Annotation/UMAP_markers/UMAP_markers_SplieBySample_", i, ".pdf"), plot = plot[[i]],
         height = 8, width = 12, dpi = 300, limitsize = FALSE)
  ggsave(paste0("5_Cell_Annotation/UMAP_markers/UMAP_markers_SplieBySample_", i, ".png"), plot = plot[[i]],
         height = 8, width = 12, dpi = 300, limitsize = FALSE)
}
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/UMAP_markers/UMAP_markers_SplieBySample_1.png" width="600" />

## 五、细胞亚群比例统计      
5.1 输出每个细胞亚群的细胞比例       
```r
# 5.1 输出每个细胞亚群的细胞比例
# 5.1
Cellnum <- table(Idents(seurat_integrated), seurat_integrated$Sample) 

Cellratio <- round(prop.table(Cellnum, margin = 2) * 100, 2)
Cellratio_save <- Cellratio %>% 
  as.data.frame.matrix() %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.2f%%", round(., 2)))) %>% 
  rownames_to_column("Clusters")

write.xlsx(Cellratio_save, '5_Cell_Annotation/Cellratio_clusters_anno.xlsx', rownames = F)
write.table(Cellratio_save, '5_Cell_Annotation/Cellratio_clusters_anno.txt', quote = F, sep = '\t', row.names = F)
#txt <- knitr::kable(Cellratio_save, format = "markdown", align = 'c')
#write.table(txt, "5_Cell_Annotation/Cellratio_clusters_anno_markdown.txt", row.names = F, quote = F, col.names = F)
```

5.2 输出每个细胞亚群的细胞数量     
```r
# 5.2
Cellnum <- table(Idents(seurat_integrated), seurat_integrated$Sample) 
Cellnum_save <- Cellnum %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("Clusters")

write.xlsx(Cellnum_save, '5_Cell_Annotation/Cellnum_clusters_anno.xlsx', rownames = F)
write.table(Cellnum_save, '5_Cell_Annotation/Cellnum_clusters_anno.txt', quote = F, sep = '\t', row.names = F)
#txt <- knitr::kable(Cellnum_save, format = "markdown", align = 'c')
#write.table(txt, "5_Cell_Annotation/Cellnum_clusters_anno_markdown.txt", row.names = F, quote = F, col.names = F)
```

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
ggsave("5_Cell_Annotation/06_Cell_proportion_stratum.pdf", plot = p_combined, height = 6, width = 12, dpi = 300)
ggsave("5_Cell_Annotation/06_Cell_proportion_stratum.png", plot = p_combined, height = 6, width = 12, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/06_Cell_proportion_stratum.png" width="600" />

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
  
  ggsave("5_Cell_Annotation/07_Cell_proportion_donut.pdf", plot = combined_plot, height = 5, width = 8, dpi = 300)
  ggsave("5_Cell_Annotation/07_Cell_proportion_donut.png", plot = combined_plot, height = 5, width = 8, dpi = 300)
  
  return(combined_plot)
}

donut_plot <- create_donut_plots(Cellratio_df)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/07_Cell_proportion_donut.png" width="600" />

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

ggsave("5_Cell_Annotation/08_Cell_proportion_barplot.pdf", plot = p, height = 6, width = 15, dpi = 300)
ggsave("5_Cell_Annotation/08_Cell_proportion_barplot.png", plot = p, height = 6, width = 15, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/08_Cell_proportion_barplot.png" width="600" />

5.7 小提琴图 基因在不同细胞亚群中，样本之间的表达差异
```r
# 5.7 小提琴图 基因在不同细胞亚群中，样本之间的表达差异
v1 <- Seurat::VlnPlot(seurat_integrated, split.by = 'Sample', group.by = 'celltype',
                      input_genes, stack=T, pt.size = 0, flip = T) +
  xlab(NULL) +
  guides(fill = guide_legend(reverse = T))

ggsave("5_Cell_Annotation/09_VlnPlot_markers_SplieBySample.pdf", plot = v1, height = 9, width = 12, dpi = 300)
ggsave("5_Cell_Annotation/09_VlnPlot_markers_SplieBySample.png", plot = v1, height = 9, width = 12, dpi = 300)

# v2 <- MySeuratWrappers::VlnPlot(seurat_integrated, split.by = 'Sample', group.by = 'celltype',
#                                 input_genes, stack=T, pt.size = 0) +
#   xlab(NULL) +
#   guides(fill = guide_legend(reverse = T))
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/09_VlnPlot_markers_SplieBySample.png" width="600" />

## 六、Boxplot_for_QC         
6.1 细胞类型选择（cells or nuclei）    
```r
metadata <- seurat_integrated@meta.data

table(metadata$celltype)

t1 <- data.frame(table(metadata$celltype))
colnames(t1) <- c('celltype', 'cell_nums')
metadata <- merge(metadata, t1, 'celltype')
metadata$a1 <- metadata$nCount_RNA/1000
metadata$b1 <- metadata$nFeature_RNA/1000

metadata$Type = 'cells'

# 设置细胞类型变量
if (metadata$Type[1] == 'cells') {
  Type <- 'cells'
} else if (metadata$Type[1] == 'nuclei') {
  Type <- 'nuclei'
} else {
  stop("错误：metadata$Type必须是'cells'或'nuclei'中的一个")
}

num = length(table(metadata$celltype))
```

6.2 绘图    
```r
a1 <- ggplot(t1, aes(x = celltype, y = cell_nums, fill = celltype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = cell_nums), vjust = 0, hjust = -0.1, color = "black") +
  labs(x = "Cell Type", y = paste0("Number of ", Type)) +
  coord_flip() + 
  theme_classic() + 
  theme(axis.title.y = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_blank(),
        legend.position = 'none')

a2 <- ggplot(metadata, aes(y = celltype, x = a1, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, color = 'black', box.colour = 'white') +
  labs(x = paste0("No.UMIs per ", Type, " (x1000)"), y = "Cell Type") +
  scale_fill_discrete() +
  guides(fill = guide_legend(reverse = T)) + 
  theme_classic() +
  coord_cartesian(xlim = c(0, 15)) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')

a3 <- ggplot(metadata, aes(y = celltype, x = b1, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, color = 'black', box.colour = 'white') +
  labs(x = paste0("No.genes per ", Type, " (x1000)"), y = "Cell Type") +
  scale_fill_discrete() +
  scale_x_continuous(breaks = c(0:6))+
  guides(fill = guide_legend(reverse = T)) + 
  theme_classic() + 
  #coord_cartesian(xlim = c(0, 5)) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')

p <- plot_grid(a1, a2, a3, nrow = 1)

marker_ls_1 <- as.character(unlist(lapply(marker_ls, function(x){ x <- x[[1]] })))

a4 <- Seurat::VlnPlot(seurat_integrated, split.by = 'celltype',
                      marker_ls_1, stack = T, pt.size = 0, flip = F)+
  ylab(NULL) +
  guides(fill = guide_legend(reverse = T)) +
  scale_fill_manual(values = hue_pal()(num))

# a4 <- MySeuratWrappers::VlnPlot(seurat_integrated, marker_ls_1, stack = T, 
#                                 pt.size = 0, direction = "horizontal") +
#   labs(x = "", y = "") +
#   guides(fill = guide_legend(reverse = F)) + 
#   scale_fill_manual(values = hue_pal()(num)) +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1),
#         axis.text.x = element_text(angle = 0, hjust = 1))

pmerge <- ggarrange(p, a4, nrow = 2) 

ggsave("5_Cell_Annotation/10_Boxplot_for_QC.pdf", plot = pmerge, height = 8, width = 8, dpi = 300)
ggsave("5_Cell_Annotation/10_Boxplot_for_QC.png", plot = pmerge, height = 8, width = 8, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/5_Cell_Annotation/10_Boxplot_for_QC.png" width="600" />

---
目录树  
```r
fs::dir_tree("5_Cell_Annotation", recurse = 2)
```
```
5_Cell_Annotation
├── 01_VlnPlot_marker.pdf
├── 01_VlnPlot_marker.png
├── 02_Cell_Annotaions_UMAP.pdf
├── 02_Cell_Annotaions_UMAP.png
├── 03_Cell_Annotaions_TSNE.pdf
├── 03_Cell_Annotaions_TSNE.png
├── 04_DoHeatmap_anno_markers.pdf
├── 04_DoHeatmap_anno_markers.png
├── 05_DotPlot_anno_markers.pdf
├── 05_DotPlot_anno_markers.png
├── 06_Cell_proportion_stratum.pdf
├── 06_Cell_proportion_stratum.png
├── 07_Cell_proportion_donut.pdf
├── 07_Cell_proportion_donut.png
├── 08_Cell_proportion_barplot.pdf
├── 08_Cell_proportion_barplot.png
├── 09_VlnPlot_markers_SplieBySample.pdf
├── 09_VlnPlot_markers_SplieBySample.png
├── 10_Boxplot_for_QC.pdf
├── 10_Boxplot_for_QC.png
├── Cellnum_clusters_anno.txt
├── Cellnum_clusters_anno.xlsx
├── Cellratio_clusters_anno.txt
├── Cellratio_clusters_anno.xlsx
├── UMAP_markers
│   ├── UMAP_markers_1.pdf
│   ├── UMAP_markers_1.png
│   ├── ...
│   ├── UMAP_markers_SplieBySample_1.pdf
│   ├── UMAP_markers_SplieBySample_1.png
│   ├── ...
├── seurat_integrated_anno.rds
└── seurat_integrated_tmp.rds
```


### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-14
- 修改日期：2025-11-29

