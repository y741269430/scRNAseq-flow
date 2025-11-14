# 01_scRNAseq_QC   
读取——质控——过滤
1. 配置环境
```r
# 调整R内允许对象大小的限制（默认值是500*1024 ^ 2 = 500 Mb）
options(future.globals.maxSize = 500 * 1024 ^ 2)

# 设置工作目录，批量创建文件夹储存结果
setwd(r"{D:\R work\GSE171169_RAW\}")

source('my_scRNAseq.R')

# 创建输出目录
if (!dir.exists("1_QC_Files")) {
  dir.create("1_QC_Files")
}

# 查看RawData下的文件夹名称，后续每个样本都是命名为该名称
projects <- list.files("RawData"); projects
```

2. 创建函数批量读取单细胞矩阵    
```r
process_multiple_projects <- function(projects, base_dir = 'RawData') {

  process_single_project <- function(project) {
    data_dir <- file.path(base_dir, project)
    
    if (!dir.exists(data_dir)) {
      warning(paste("目录不存在:", data_dir))
      return(NULL)
    }

    count_matrix <- Read10X(data.dir = data_dir, gene.column = 2)
    # #如果count_matrix的行名的ENSEMBL ID，则可以通过以下命令转换为SYMBOL，即Gene name，并处理NA值
    # gene_symbols <- mapIds(org.Mm.eg.db, keys = count_matrix@Dimnames[[1]], column = "SYMBOL", keytype = "ENSEMBL")
    # gene_symbols <- ifelse(is.na(gene_symbols), paste0("NA-", seq_along(gene_symbols)), gene_symbols)
    # count_matrix@Dimnames[[1]] <- gene_symbols
    
    seurat_obj <- CreateSeuratObject(
      counts = count_matrix,   
      min.cells = 3, 
      min.features = 200,
      project = project
    ) 
    return(seurat_obj)
  }

  seurat_list <- lapply(projects, process_single_project)

  names(seurat_list) <- projects
  seurat_list <- Filter(Negate(is.null), seurat_list)
  
  return(seurat_list)
}

seurat_objects <- process_multiple_projects(projects)

saveRDS(seurat_objects, "1_QC_Files/seurat_objects.rds")
```
3. 生成质控指标    
```r
# 检查矩阵行名，是ENSEMBL还是gene name（SYMBOL）
head(rownames(seurat_objects[[1]]@assays$RNA))

# 定义过滤的基因列表
mt_genes <- grep('^mt-', rownames(seurat_objects[[1]]), value = TRUE)
Hb_genes_total <- c("Hbq1a","Hbb-y","Hbb-bh1","Hbb-bs","Hba-x","Hba-a1","Hba-a2","Hba-a3", 
                    "Hbq1b","Hbb-bt","Hbb-bh2","Hbb-bh3","Hba-ps4","Hbb-bh0","Hba-ps3")

# 批量生成对应指标
seurat_objects <- lapply(seurat_objects, function(obj) {
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)         # 计算log10GenesPerUMI
  obj$mitoRatio <- PercentageFeatureSet(obj, features = mt_genes) / 100           # 计算线粒体基因比率
  hb_genes <- intersect(Hb_genes_total, rownames(obj))                            # 计算红细胞基因比率
  obj$percent.Hb <- if (length(hb_genes) > 0) { 
    PercentageFeatureSet(obj, features = hb_genes) 
    } else { 0 }
  obj$percent.ribo <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")              # 计算核糖体基因比率
  obj$cells <- rownames(obj@meta.data)
  return(obj)
})

# 合并metadata
combined_meta <- do.call(rbind, lapply(seq_along(seurat_objects), function(i) {
  data.frame(Sample = ifelse(!is.null(names(seurat_objects)), names(seurat_objects)[i], paste0("Sample", i)),
             seurat_objects[[i]]@meta.data)
}))
```
4. 质控1    
```r
# 创建质控图
plots <- list(
  ggplot(combined_meta, aes(Sample, nFeature_RNA, fill = Sample)) + geom_violin() + labs(title = "Gene Count"),
  ggplot(combined_meta, aes(Sample, nCount_RNA, fill = Sample)) + geom_violin() + labs(title = "UMI Count"),
  ggplot(combined_meta, aes(Sample, mitoRatio, fill = Sample)) + geom_violin() + labs(title = "MT Ratio"),
  ggplot(combined_meta, aes(Sample, percent.Hb, fill = Sample)) + geom_violin() + labs(title = "Hb Percent"),
  ggplot(combined_meta, aes(Sample, percent.ribo, fill = Sample)) + geom_violin() + labs(title = "Ribo Percent")
)

# 应用统一主题
plots <- lapply(plots, function(p) {
  p + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
})

p <- plot_grid(plotlist = plots, nrow = 2)
ggsave("1_QC_Files/01_QC_VlnPlot.pdf", p, height = 5, width = 7, dpi = 300)
ggsave("1_QC_Files/01_QC_VlnPlot.png", p, height = 5, width = 7, dpi = 300)

# 细胞计数 (Cell Numbers before Filter)
cell_nums_before_filter <- data.frame(table(combined_meta$Sample))

# 细胞计数可视化
before <- combined_meta %>%
  ggplot(aes(x = Sample, fill = Sample)) +
  geom_bar(position = "dodge", show.legend = TRUE) +
  geom_text(stat = 'count', aes(label = after_stat(count)), position = position_dodge(1), vjust = -0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Cell Numbers before Filter"); before
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/1_QC_Files/01_QC_VlnPlot.png" width="500" />    

5. 质控2    
```r
# 每个细胞的UMI计数 (UMI counts per cell)
a1 <- combined_meta %>% 
  ggplot(aes(color = Sample, x = nCount_RNA, fill = Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10(breaks = c(500, 1000, 5000, 10000, 20000)) + 
  theme_classic() +
  xlab("nUMI") +
  ylab("Cell density") +
  geom_vline(xintercept = c(500, 1000, 5000, 10000, 20000), linetype = 'dotted') +
  ggtitle("UMI counts per cell")

# 复杂度 (Complexity)
# 通过可视化基因数与UMI的比率（log10基因数/UMI）来表示基因表达的整体复杂性
a2 <- combined_meta %>%
  ggplot(aes(x=log10GenesPerUMI, color = Sample, fill = Sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = c(0.8, 0.85, 0.9, 0.95), linetype = 'dotted') +
  ggtitle("Complexity")

# 每个细胞检测到的基因数分布
a3 <- combined_meta %>% 
  ggplot(aes(color=Sample, x=nFeature_RNA, fill = Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) + 
  geom_vline(xintercept = c(500, 1000, 2500, 5000, 10000), linetype = 'dotted') +
  ggtitle("Number of genes per cell")

# 每个细胞检测到的基因数量的分布（箱线图）
a4 <- combined_meta %>% 
  ggplot(aes(x=Sample, y=log10(nFeature_RNA), fill=Sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Genes detected per cell across Samples")

p <- plot_grid(a1, a2, a3, a4, align = "v", nrow = 2)

ggsave("1_QC_Files/02_QC_Four.pdf", p, height = 8, width = 12, dpi = 300)
ggsave("1_QC_Files/02_QC_Four.png", p, height = 8, width = 12, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/1_QC_Files/02_QC_Four.png" width="600" />    

6. 质控3    
```r
# 线粒体基因计数占比 (Mitochondrial counts ratio)
# 可视化每个细胞检测到的线粒体基因表达分布
a1 <- combined_meta %>% 
  ggplot(aes(color=Sample, x=mitoRatio, fill=Sample)) + 
  geom_density(alpha = 0.25) + 
  theme_classic() +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.3)) + 
  geom_vline(xintercept = c(0.001, 0.01, 0.1, 0.3), linetype = 'dotted') +
  ggtitle("Mitochondrial counts ratio")

# 红细胞基因计数占比 (HB counts ratio)
# 可视化每个细胞检测到的红细胞基因表达分布
a2 <- combined_meta %>% 
  ggplot(aes(color=Sample, x=percent.Hb, fill=Sample)) + 
  geom_density(alpha = 0.25) + 
  theme_classic() +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.5)) + 
  geom_vline(xintercept = c(0.001, 0.01, 0.1, 0.5), linetype = 'dotted') +
  ggtitle("Hb counts ratio")

# 核糖体基因计数占比 (RB counts ratio)
# 可视化每个细胞检测到的核糖体基因表达分布
a3 <- combined_meta %>% 
  ggplot(aes(color=Sample, x=percent.ribo, fill=Sample)) + 
  geom_density(alpha = 0.25) + 
  theme_classic() +
  scale_x_log10(breaks = c(0.1, 1, 10)) + 
  geom_vline(xintercept = c(0.1, 1, 10), linetype = 'dotted') +
  ggtitle("Ribosome counts ratio")

p <- plot_grid(a1, a2, a3, align = "v", nrow = 2)

ggplot2::ggsave("1_QC_Files/03_QC_Density.pdf", plot = p, height = 6, width = 8, dpi = 300)
ggplot2::ggsave("1_QC_Files/03_QC_Density.png", plot = p, height = 6, width = 8, dpi = 300)
```    
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/1_QC_Files/03_QC_Density.png" width="600" />    

```r
# 检测到的UMI数对比基因数 (UMIs vs. genes detected)
# 可视化每个细胞中检测到的基因数（nFeature_RNA）与UMI数（nCount_RNA）之间的关系，
# 颜色代表线粒体基因计数占比（mitoRatio），并观察是否存在大量低基因数或低UMI数的细胞。

p <- combined_meta %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Sample)

ggplot2::ggsave("1_QC_Files/04_QC_UMIs_vs_genes.pdf", plot = p, height = 8, width = 8, dpi = 300)
ggplot2::ggsave("1_QC_Files/04_QC_UMIs_vs_genes.png", plot = p, height = 8, width = 8, dpi = 300)
```    
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/1_QC_Files/04_QC_UMIs_vs_genes.png" width="600" />     
  
7. 细胞过滤     
```r    
# 细胞过滤，具体情况具体分析，我这里的指标是根据文献的指标而定的。
seurat_filter <- lapply(seurat_objects, function(x){
  x <- subset(x, subset = 
                #(metadata$nCount_RNA >= 500) &
                (nFeature_RNA >= 200) & 
                (nFeature_RNA <= 5000) &
                #(log10GenesPerUMI > 0.85) & 
                (mitoRatio < 0.05) )
})

# 统计
filter_meta <- do.call(rbind, lapply(seq_along(seurat_filter), function(i) {
  data.frame(Sample = ifelse(!is.null(names(seurat_filter)), names(seurat_filter)[i], paste0("Sample", i)),
             seurat_filter[[i]]@meta.data)
}))

# 细胞计数 (Cell Numbers before Filter)
cell_nums_after_filter <- data.frame(table(filter_meta$Sample))
```
8. 可视化过滤前后的细胞计数    
```r    
after <- filter_meta %>%
  ggplot(aes(x = Sample, fill = Sample)) +
  geom_bar(position = "dodge", show.legend = TRUE) +
  geom_text(stat = 'count', aes(label = after_stat(count)), position = position_dodge(1), vjust = -0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Cell Numbers after Filter")

p <- plot_grid(before, after); p

ggplot2::ggsave("1_QC_Files/05_QC_NCells.pdf", plot = p, height = 4, width = 8, dpi = 300)
ggplot2::ggsave("1_QC_Files/05_QC_NCells.png", plot = p, height = 4, width = 8, dpi = 300)
```    
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/1_QC_Files/05_QC_NCells.png" width="500" />      

9. 统计过滤前后的细胞数量，并保存结果    
```r    
cell_nums <- cbind(cell_nums_before_filter, cell_nums_after_filter[,2])
cell_nums$Vaild <- percent(cell_nums[,3]/cell_nums[,2])
cell_nums$Filter <- percent(1 - cell_nums[,3]/cell_nums[,2])

colnames(cell_nums)[1:3] <- c('Sample','Cell_nums_before_filter','Cell_nums_after_filter')

write.table(cell_nums, '1_QC_Files/cell_nums.txt', row.names = F, quote = F, sep = '\t')
write.xlsx(cell_nums, '1_QC_Files/cell_nums.xlsx', rowNames = F)
# 生成Markdown表格
txt <- knitr::kable(cell_nums, format = "markdown", align = 'c')
write.table(txt, '1_QC_Files/markdown.txt', row.names = F, quote = F, col.names = F)

saveRDS(seurat_filter, "1_QC_Files/seurat_filter.rds")
```    
|Sample | Cell_nums_before_filter| Cell_nums_after_filter|Vaild  |Filter |
|:------|:-----------------------:|:----------------------:|:------|:------|
|05d_N1 |                    2542|                   2418|95.12% |4.88%  |
|05d_N2 |                    3484|                   3289|94.40% |5.60%  |
|14d_N1 |                    2583|                   2524|97.72% |2.28%  |
|14d_N2 |                    2313|                   2288|98.92% |1.08%  |

---
目录树      
```r
fs::dir_tree("1_QC_Files", recurse = 2)
```
```
1_QC_Files
├── 01_QC_VlnPlot.pdf
├── 01_QC_VlnPlot.png
├── 02_QC_Four.pdf
├── 02_QC_Four.png
├── 03_QC_Density.pdf
├── 03_QC_Density.png
├── 04_QC_UMIs_vs_genes.pdf
├── 04_QC_UMIs_vs_genes.png
├── 05_QC_NCells.pdf
├── 05_QC_NCells.png
├── cell_nums.txt
├── cell_nums.xlsx
├── markdown.txt
├── seurat_filter.rds
└── seurat_objects.rds
```

### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-08
- 修改日期：2025-11-08


