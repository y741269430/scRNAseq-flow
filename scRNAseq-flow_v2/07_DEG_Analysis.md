# 07_DEG_Analysis

### 1. 输入准备
```r
# 设置工作目录，批量创建文件夹储存结果
setwd(r"{D:\R work\GSE171169_RAW\}")
# 加载R包
source('my_scRNAseq.R')

# 创建输出目录
if (!dir.exists("7_DEG_Analysis")) {
  dir.create("7_DEG_Analysis")
}
```

### 2. 加载矩阵
```r
seurat_integrated <- readRDS('5_Cell_Annotation/seurat_integrated_anno.rds')
Idents(seurat_integrated)
```

### 3. 准备差异表达分析的比较组
```r
seurat_integrated$celltype.exp <- paste(seurat_integrated$Sample, 
                                        seurat_integrated$celltype, 
                                        sep = "-")

Idents(seurat_integrated) <- seurat_integrated$celltype.exp

# 定义需要进行差异表达分析的细胞类型组合
degmeta <- data.frame(table(seurat_integrated$celltype.exp))
degmeta$Sample <- str_split_fixed(degmeta$Var1, '-', n = 2)[,1]
degmeta$cluster <- str_split_fixed(degmeta$Var1, '-', n = 2)[,2]
degmeta_ls <- lapply(split(degmeta, degmeta$cluster), function(x){ x <- x; return(x)})

combinations <- lapply(degmeta_ls, function(cell_data){
  TRE <- cell_data[cell_data$Sample == "14d_N1", "Var1"]
  CTRL <- cell_data[cell_data$Sample == "05d_N1", "Var1"]
  
  if(length(TRE) > 0 & length(CTRL) > 0){
    data.frame(
      Treatment = TRE,
      Control = CTRL, 
      Celltype = unique(cell_data$cluster)
    )
  }
})

group <- do.call(rbind, combinations)
rownames(group) <- NULL

# 查看最终的分组表格
knitr::kable(group, format = "markdown", align = 'c')
```

|      Treatment       |       Control        |   Celltype    |
|:--------------------:|:--------------------:|:-------------:|
|  14d_N1-ab_T_cells   |  05d_N1-ab_T_cells   |  ab_T_cells   |
|    14d_N1-B_cells    |    05d_N1-B_cells    |    B_cells    |
|     14d_N1-cDCs      |     05d_N1-cDCs      |     cDCs      |
|  14d_N1-gd_T_cells   |  05d_N1-gd_T_cells   |  gd_T_cells   |
| 14d_N1-Macrophages_1 | 05d_N1-Macrophages_1 | Macrophages_1 |
| 14d_N1-Macrophages_2 | 05d_N1-Macrophages_2 | Macrophages_2 |
|  14d_N1-Microglia_1  |  05d_N1-Microglia_1  |  Microglia_1  |
|  14d_N1-Microglia_2  |  05d_N1-Microglia_2  |  Microglia_2  |
|  14d_N1-Neutrophils  |  05d_N1-Neutrophils  |  Neutrophils  |
|   14d_N1-NK_cells    |   05d_N1-NK_cells    |   NK_cells    |
|     14d_N1-pDCs      |     05d_N1-pDCs      |     pDCs      |


### 4. 差异表达分析    
```r
# 执行差异分析
p_cutoff = 0.05
log2fc_cutoff = 0.25

DEG_list <- lapply(1:nrow(group), function(i){
  cat("DEGs:", group$Celltype[i], "\n")
  
  markers <- FindMarkers(
    object = seurat_integrated,     # 输入要计算的矩阵
    ident.1 = group$Treatment[i],
    ident.2 = group$Control[i],
    min.pct = 0.1,
    logfc.threshold = log2fc_cutoff,
    test.use = "wilcox"
  )

  markers$gene <- rownames(markers)
  markers$cluster <- group$Celltype[i]
  markers$comparison <- paste(str_split_fixed(group$Treatment[1], n = 2, '-')[,1], 
                              "vs", str_split_fixed(group$Control[1], n = 2, '-')[,1])
  
  return(markers)
})

names(DEG_list) <- group$Celltype

merged_DEG_list <- lapply(DEG_list, function(x){
  id_df <- bitr(x$gene, 
                fromType = 'SYMBOL', 
                toType = c('ENTREZID', 'GENENAME'), 
                OrgDb = 'org.Mm.eg.db')
  merged <- merge(id_df, x, by.x = "SYMBOL", by.y = "gene", all.y = TRUE)
  merged <- merged[order(merged$avg_log2FC, decreasing = TRUE), ]
  merged <- merged[!duplicated(merged$ENTREZID) | is.na(merged$ENTREZID), ]
  
  colnames(merged)[colnames(merged) == "SYMBOL"] <- "gene"
  
  return(merged)
})

merged_DEG_list_p <- lapply(merged_DEG_list, function(x){
  x <- subset(x, p_val < p_cutoff)
  return(x)
})

write.xlsx(merged_DEG_list, file = "7_DEG_Analysis/DEG_list_log_025.xlsx")
write.xlsx(merged_DEG_list_p, file = "7_DEG_Analysis/DEG_list_log_025_p_005.xlsx")
save(merged_DEG_list, merged_DEG_list_p, file = "7_DEG_Analysis/DEG_list.RData")
```

| 变量名 (Variable) | 英文解释 | 中文解释 |
| :---: | :---: | :---: |
| **avg_log2FC** | Log fold-change of the average expression between the two groups. Positive values indicate higher expression in the first group. | 两组间平均表达差异倍数的对数值。正值表示该基因在第一组中的表达更高。 |
| **pct.1** | The percentage of cells where the gene is detected in the first group. | 第一组中检测到该基因表达的细胞百分比。 |
| **pct.2** | The percentage of cells where the gene is detected in the second group. | 第二组中检测到该基因表达的细胞百分比。 |
| **p_val** | Raw p-value from the hypothesis test. | 假设检验后得到的原始P值。 |
| **p_val_adj** | Adjusted p-value, based on Bonferroni correction using all genes in the dataset. | 基于Bonferroni多重检验校正后得到的校正P值。 |

### 5. 可视化 
```r
# 5.1 绘制差异基因火山图
data1 <- Reduce(rbind, merged_DEG_list)

a1 <- jjVolcano(diffData = data1, 
                pSize = 1.5, 
                log2FC.cutoff = log2fc_cutoff,
                tile.col = brewer.pal(11, "RdBu"),
                topGeneN = 5, 
                aesCol = c("#bbe0f2", "#f894af"),
                polar = F) +                         # 是否为环形火山图
  geom_hline(yintercept = 0.5, linetype = 'dotted') +
  geom_hline(yintercept = -0.5, linetype = 'dotted') + 
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank()) + 
  ggtitle(unique(data1$comparison))

ggplot2::ggsave("7_DEG_Analysis/01_DEG_volcano.pdf", plot = a1,
                height = 8, width = 14, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("7_DEG_Analysis/01_DEG_volcano.png", plot = a1,
                height = 8, width = 14, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/7_DEG_Analysis/01_DEG_volcano.png" width="600" />   

```r
# 5.2 统计每个细胞类型的上下调基因数量（以p_val < 0.05和|log2FC| >= 0.25为标准）

p_cutoff = 0.05
log2fc_cutoff = 0.25

deg_stats <- lapply(names(merged_DEG_list_p), function(celltype) {
  df <- merged_DEG_list_p[[celltype]]

  up_genes <- sum(df$p_val < p_cutoff & df$avg_log2FC >= log2fc_cutoff, na.rm = TRUE)
  down_genes <- sum(df$p_val < p_cutoff & df$avg_log2FC <= -log2fc_cutoff, na.rm = TRUE)
  
  data.frame(celltype = celltype, up = up_genes, down = down_genes)
})

deg_stats_df <- do.call(rbind, deg_stats)

knitr::kable(deg_stats_df, format = "markdown", align = 'c')
```

|   celltype    | up  | down |
|:-------------:|:---:|:----:|
|  ab_T_cells   | 213 | 1131 |
|    B_cells    | 105 | 492  |
|     cDCs      | 307 | 242  |
|  gd_T_cells   | 424 | 376  |
| Macrophages_1 | 302 | 148  |
| Macrophages_2 | 341 | 386  |
|  Microglia_1  | 186 | 566  |
|  Microglia_2  | 417 | 309  |
|  Neutrophils  | 132 | 346  |
|   NK_cells    | 398 | 175  |
|     pDCs      | 177 | 281  |

```r
deg_stats_long <- deg_stats_df %>% 
  pivot_longer(cols = c(up, down), names_to = "regulation", values_to = "count")

# 绘制分组柱形图
a2 <- ggplot(deg_stats_long, aes(x = celltype, y = count, fill = regulation)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = count), 
            position = position_dodge(0.7), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("up" = "#e74c3c", "down" = "#3498db")) +
  labs(title = "Statistics of DEG in each cell types",
       subtitle = paste("Criteria: p_val <",p_cutoff, "& |log2FC| ≥", log2fc_cutoff),
       x = "Cell types", 
       y = "DEG numbers",
       fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot2::ggsave("7_DEG_Analysis/02_DEG_counts.pdf", plot = a2,
                height = 6, width = 8, dpi = 300, limitsize = FALSE)
ggplot2::ggsave("7_DEG_Analysis/02_DEG_counts.png", plot = a2,
                height = 6, width = 8, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/7_DEG_Analysis/02_DEG_counts.png" width="600" />   

```r
fs::dir_tree("7_DEG_Analysis", recurse = 2)
```
```bash
7_DEG_Analysis
├── 01_DEG_volcano.pdf
├── 01_DEG_volcano.png
├── 02_DEG_counts.pdf
├── 02_DEG_counts.png
├── DEG_list.RData
├── DEG_list_log_025.xlsx
└── DEG_list_log_025_p_005.xlsx
```

### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-29
- 修改日期：2025-11-29
