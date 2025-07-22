# 查看metadata

加载函数
```r
library(ggplot2)
library(patchwork)
create_violin_plot <- function(df, feature_name, group_var = "orig.ident") {
  ggplot(df, aes_string(x = group_var, y = feature_name, fill = group_var)) +
    geom_violin() +
    theme_minimal() +
    labs(y = feature_name, x = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))  # 隐藏图例
}
```
加载数据
```r
metadata <- read.csv('seurat_metadata.csv')
```

#### 查看质控图1 #### 
```r
a1 <- create_violin_plot(metadata, "nFeature_RNA")
a2 <- create_violin_plot(metadata, "nCount_RNA")
a3 <- create_violin_plot(metadata, "mitoRatio")
a4 <- create_violin_plot(metadata, "percent.Hb")
a5 <- create_violin_plot(metadata, "percent.ribo")

# 合并所有图像
combined_plot <- (a1 + a2 + a3 + a4 + a5) + plot_layout(nrow = 2)
```

#### 查看质控图2 #### 
```r
# 每个细胞的UMI计数 (UMI counts per cell)
# 可视化每个细胞的UMI计数
a1 <- metadata %>% 
  ggplot(aes(color = sample, x = nCount_RNA, fill = sample)) + 
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
  ggtitle("Number of genes per cell")

# 每个细胞检测到的基因数量的分布（箱线图）
a4 <- metadata %>% 
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Genes detected per cell across samples")

plot_grid(a1, a2, a3, a4, align = "v", nrow = 2)
```
#### 查看细胞数量图 #### 
```r
metadata %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) +
  geom_bar(position = "dodge", show.legend = TRUE) +
  geom_text(stat = 'count', aes(label = after_stat(count)), position = position_dodge(0.9), vjust = -0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Cell counts per sample")
```

#### 查看UMAP ####

以sample做分面展示
```r
ggplot(metadata, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = factor(RNA_snn_res.0.5)), size = 1) + 
  labs(title = "UMAP Plot", x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
  theme_minimal() +  
  theme(
    plot.title = element_text(hjust = 0.5),  
    axis.text = element_text(size = 12),  
    legend.position = "right"  # 将图例放置在右侧
  ) + 
  facet_wrap(~ sample) + 
  scale_color_discrete(name = "Cluster")
```

不分面图
```r
ggplot(metadata, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = factor(RNA_snn_res.0.5)), size = 1) + 
  labs(title = "UMAP Plot", x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
  theme_minimal() +  
  theme(
    plot.title = element_text(hjust = 0.5),  
    axis.text = element_text(size = 12),  
    legend.position = "right"  # 将图例放置在右侧
  ) +
  scale_color_discrete(name = "Cluster")
```
