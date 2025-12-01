# 08_DEG_Enrichment

### 1. 输入准备
```r
# 设置工作目录，批量创建文件夹储存结果
setwd(r"{D:\R work\GSE171169_RAW\}")
# 加载R包
source('my_scRNAseq.R')

# 创建输出目录
if (!dir.exists("8_DEG_Enrichment")) {
  dir.create("8_DEG_Enrichment")
}
if (!dir.exists("8_DEG_Enrichment/GO_Enrichment")) {
  dir.create("8_DEG_Enrichment/GO_Enrichment")
}
if (!dir.exists("8_DEG_Enrichment/KEGG_Enrichment")) {
  dir.create("8_DEG_Enrichment/KEGG_Enrichment")
}
```
### 2. 加载数据
```r
load("7_DEG_Analysis/DEG_list.RData")
data1 <- Reduce(rbind, merged_DEG_list)
```

### 3.1 执行GO BP 的富集分析（CC和MF也是同理）
```r
# 执行GO BP 的富集分析（CC和MF也是同理）
enrich_go_bp_list <- lapply(names(merged_DEG_list_p), function(celltype) {
  cat("正在进行GO富集分析(Biological Process):", celltype, "\n")
  df <- merged_DEG_list_p[[celltype]]
  result <- clusterProfiler::enrichGO(
    gene = df$ENTREZID,
    ont = "BP",          # CC MF
    OrgDb = 'org.Mm.eg.db',
    keyType = 'ENTREZID',
    readable = TRUE,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
})

names(enrich_go_bp_list) <- names(merged_DEG_list_p)
enrich_go_bp_list <- enrich_go_bp_list[!sapply(enrich_go_bp_list, is.null)]
```

### 3.2 获取每个通路中上调和下调的基因
```r
GO_BP_updown <- list()

for(celltype in names(enrich_go_bp_list)) {
  cat("正在处理:", celltype, "\n")

  GO_BP_updown[[celltype]] <- add_regulation_direction(
    enrich_go_bp_list[[celltype]], 
    merged_DEG_list_p[[celltype]]
  )
}

# 批量保存
for (i in 1:length(GO_BP_updown)) {
  write.xlsx(GO_BP_updown[[i]]@result, paste0('8_DEG_Enrichment/GO_Enrichment/GO_BP-', names(GO_BP_updown)[[i]],'.xlsx'))
}
```

### 3.3 将所有celltype的@result提取到一个list中
```r
GO_BP_results_list <- lapply(GO_BP_updown, function(x) { return(x@result) })
write.xlsx(GO_BP_results_list, file = '8_DEG_Enrichment/GO_BP_results_list.xlsx')

save(GO_BP_updown, GO_BP_results_list, file = '8_DEG_Enrichment/GO_BP.RData')
```

### 4.1 执行KEGG 的富集分析
```r
# 执行KEGG 的富集分析
enrich_kegg_list <- lapply(names(merged_DEG_list_p), function(celltype) {
  cat("正在进行GO富集分析(Biological Process):", celltype, "\n")
  df <- merged_DEG_list_p[[celltype]]
  result <- my_enrichKEGG(
    gene = df$ENTREZID,
    organism  = 'Mmu',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
})

names(enrich_kegg_list) <- names(merged_DEG_list_p)

enrich_kegg_list <- enrich_kegg_list[!sapply(enrich_kegg_list, is.null)]

enrich_kegg_list <- lapply(enrich_kegg_list, function(x) {
  setReadable(x, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
})
```


### 4.2 获取每个通路中上调和下调的基因
```r
KEGG_updown <- list()

for(celltype in names(enrich_kegg_list)) {
  cat("正在处理:", celltype, "\n")
  
  KEGG_updown[[celltype]] <- add_regulation_direction(
    enrich_kegg_list[[celltype]], 
    merged_DEG_list_p[[celltype]]
  )
}

# 批量保存
for (i in 1:length(KEGG_updown)) {
  write.xlsx(KEGG_updown[[i]]@result, paste0('8_DEG_Enrichment/KEGG_Enrichment/KEGG-', names(KEGG_updown)[[i]],'.xlsx'))
}
```

### 4.3 将所有celltype的@result提取到一个list中
```r
KEGG_results_list <- lapply(KEGG_updown, function(x) { return(x@result) })
write.xlsx(KEGG_results_list, file = '8_DEG_Enrichment/KEGG_results_list.xlsx')

save(KEGG_updown, KEGG_results_list, file = '8_DEG_Enrichment/KEGG.RData')
```

### 5 GO可视化    
```r
# GO BP
load('8_DEG_Enrichment/GO_BP.RData')

# 构建可视化所需的表格
enrichment_df <- lapply(GO_BP_results_list, function(df) {
  df_sorted <- df[order(df$pvalue, decreasing = FALSE), ]
  df_processed <- df_sorted %>%
    mutate(
      GeneRatio_value = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      EnrichmentFactor = (as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)) / 
                            as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 1))),
      total_count = up_count + down_count
    )
  return(df_processed)
})
```
```r
# 批量绘制dotplot
plot <- list()
for (i in 1:length(enrichment_df)) {
  plot[[i]] <- mysc_dotplot(enrichment_df[[i]], title = paste0("Enrichment Analysis Dotplot in \n", names(enrichment_df)[i]))
  
  ggsave(paste0('8_DEG_Enrichment/GO_Enrichment/01_GOBP_top10_dotplot_',names(enrichment_df)[i],'.pdf'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
  
  ggsave(paste0('8_DEG_Enrichment/GO_Enrichment/01_GOBP_top10_dotplot_',names(enrichment_df)[i],'.png'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
}
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/8_DEG_Enrichment/GO_Enrichment/01_GOBP_top10_dotplot_ab_T_cells.png" width="600" />   
```r
# 批量绘制barplot
plot <- list()
for (i in 1:length(enrichment_df)) {
  plot[[i]] <- mysc_barplot(enrichment_df[[i]], title = paste0("Enrichment Analysis Barplot in \n", names(enrichment_df)[i]))
  
  ggsave(paste0('8_DEG_Enrichment/GO_Enrichment/02_GOBP_top10_barplot_',names(enrichment_df)[i],'.pdf'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
  
  ggsave(paste0('8_DEG_Enrichment/GO_Enrichment/02_GOBP_top10_barplot_',names(enrichment_df)[i],'.png'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
}
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/8_DEG_Enrichment/GO_Enrichment/02_GOBP_top10_barplot_ab_T_cells.png" width="600" />   

### 6 KEGG可视化    
```r
# KEGG
load('8_DEG_Enrichment/KEGG.RData')

# 构建可视化所需的表格
enrichment_df <- lapply(KEGG_results_list, function(df) {
  df_sorted <- df[order(df$pvalue, decreasing = FALSE), ]
  df_processed <- df_sorted %>%
    mutate(
      GeneRatio_value = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      EnrichmentFactor = (as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)) / 
                            as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 1))),
      total_count = up_count + down_count
    )
  return(df_processed)
})

for (i in 1:length(enrichment_df)) {
  enrichment_df[[i]]$Description <- str_split_fixed(enrichment_df[[i]]$Description, ' - Mus', n = 2)[,1]
}

# 批量绘制dotplot
plot <- list()
for (i in 1:length(enrichment_df)) {
  plot[[i]] <- mysc_dotplot(enrichment_df[[i]], title = paste0("Enrichment Analysis Dotplot in \n", names(enrichment_df)[i]))
  
  ggsave(paste0('8_DEG_Enrichment/KEGG_Enrichment/01_KEGG_top10_dotplot_',names(enrichment_df)[i],'.pdf'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
  
  ggsave(paste0('8_DEG_Enrichment/KEGG_Enrichment/01_KEGG_top10_dotplot_',names(enrichment_df)[i],'.png'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
}
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/8_DEG_Enrichment/KEGG_Enrichment/01_KEGG_top10_dotplot_ab_T_cells.png" width="600" />   
```r
# 批量绘制barplot
plot <- list()
for (i in 1:length(enrichment_df)) {
  plot[[i]] <- mysc_barplot(enrichment_df[[i]], title = paste0("Enrichment Analysis Barplot in \n", names(enrichment_df)[i]))
  
  ggsave(paste0('8_DEG_Enrichment/KEGG_Enrichment/02_KEGG_top10_barplot_',names(enrichment_df)[i],'.pdf'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
  
  ggsave(paste0('8_DEG_Enrichment/KEGG_Enrichment/02_KEGG_top10_barplot_',names(enrichment_df)[i],'.png'), 
         plot = plot[[i]], height = 5, width = 7, dpi = 300, limitsize = FALSE)
}
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/8_DEG_Enrichment/KEGG_Enrichment/02_KEGG_top10_barplot_ab_T_cells.png" width="600" />   

---
目录树  
```r
fs::dir_tree("8_DEG_Enrichment", recurse = 2)
```
```bash
8_DEG_Enrichment
├── GO_BP.RData
├── GO_BP_results_list.xlsx
├── GO_Enrichment
│   ├── 01_GOBP_top10_dotplot_ab_T_cells.pdf
│   ├── 01_GOBP_top10_dotplot_ab_T_cells.png
│   ├── 01_GOBP_top10_dotplot_B_cells.pdf
│   ├── 01_GOBP_top10_dotplot_B_cells.png
│   ├── ...
│   ├── 02_GOBP_top10_barplot_ab_T_cells.pdf
│   ├── 02_GOBP_top10_barplot_ab_T_cells.png
│   ├── 02_GOBP_top10_barplot_B_cells.pdf
│   ├── 02_GOBP_top10_barplot_B_cells.png
│   ├── ...
│   ├── GO_BP-ab_T_cells.xlsx
│   ├── GO_BP-B_cells.xlsx
│   ├── ...
├── KEGG.RData
├── KEGG_Enrichment
│   ├── 01_KEGG_top10_dotplot_ab_T_cells.pdf
│   ├── 01_KEGG_top10_dotplot_ab_T_cells.png
│   ├── 01_KEGG_top10_dotplot_B_cells.pdf
│   ├── 01_KEGG_top10_dotplot_B_cells.png
│   ├── ...
│   ├── 02_KEGG_top10_barplot_ab_T_cells.pdf
│   ├── 02_KEGG_top10_barplot_ab_T_cells.png
│   ├── 02_KEGG_top10_barplot_B_cells.pdf
│   ├── 02_KEGG_top10_barplot_B_cells.png
│   ├── ...
│   ├── KEGG-ab_T_cells.xlsx
│   ├── KEGG-B_cells.xlsx
│   ├── ...
└── KEGG_results_list.xlsx
```

### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-30
- 修改日期：2025-12-01

