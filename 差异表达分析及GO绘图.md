# 差异表达分析及GO绘图

## 目录 ####
- 1.差异表达分析
- 1.1 制作一个contrast的list，用来做循环
- 1.2 提取差异基因
- 1.3 差异基因火山图
- 1.4 计算差异基因数量柱形图
- 1.5 探讨FindMarkers是如何进行log2FC计算的
- 2.GO 富集分析 

## 1.差异表达分析 ####    
参考 以下这篇文章的差异表达基因选取是 +/- 0.25    
[Nat Neurosci (2024)](https://www.nature.com/articles/s41593-024-01791-4#Sec2)  Analysis of normal aging brains 部分    

```r
# 加载文件
path = 'F:/R work/mmbrain/results/'
seurat_integrated <- readRDS(paste0(path, "seurat_integrated_2.rds"))

# 将细胞类型与样本名称进行拼接，即：celltype_sample
seurat_integrated$celltype.exp <- paste(Idents(seurat_integrated), seurat_integrated$sample, sep = "_")

Idents(seurat_integrated) <- seurat_integrated$celltype.exp

table(seurat_integrated@meta.data$celltype.exp)
```
#### 1.1 制作一个contrast的list，用来做循环
```r
# 定义需要进行差异表达分析的细胞类型组合
degmeta <- data.frame(table(seurat_integrated@meta.data$celltype.exp))

# 为什么这里有这句命令，假设我们有2个样本，那么这里的`each`就填2，意思就是每个一celltype，用一个数字来表示它。  
degmeta$cluster <- rep(1:length(unique(seurat_integrated$celltype)), each = 2)

```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/degmeta.png" width="500" />

```r
# 分割每个样本进行DEG
degmeta_ls <- lapply(split(degmeta, degmeta$cluster), function(x){ x <- x; return(x)})

# 定义需要进行差异表达分析的细胞类型组合
combinations <- list()

for(i in 1:length(degmeta_ls)){
  combinations[[i]] <- list(c(degmeta_ls[[i]][2,1], degmeta_ls[[i]][1,1]))
}

names(combinations) <- names(degmeta_ls)

head(combinations[[1]])
```
#### 运行以下循环
```r
# 创建一个空的列表来保存结果
ALL <- list()

# 循环遍历每个细胞类型组合进行差异表达分析
for(j in 1:length(unique(seurat_integrated$celltype)) ){
  for (i in 1:length(unique(combinations[[1]])) ) {
    
    # j 是 length(unique(seurat_integrated$celltype)) 个细胞亚群
    # i 是 length(unique(combinations[[1]])) 个对比
    
    ident_1 <- combinations[[j]][[i]][1]
    ident_2 <- combinations[[j]][[i]][2]
    
    # avg_logFC: log fold-chage of the average expression between the two groups. 
    #            Positive values indicate that the gene is more highly expressed in the first group (前比后，高)
    # pct.1: The percentage of cells where the gene is detected in the first group
    # pct.2: The percentage of cells where the gene is detected in the second group
    # p_val_adj: Adjusted p-value, based on bonferroni correction using all genes in the dataset
    
    # p_val：假设检验后得到的原始P值
    # avg_logFC：两组之间平均表达差异倍数的对数值。正值表示该基因在第一组中的表达更高。
    # pct.1：第一组中检测到表达该基因的细胞所占的百分比
    # pct.2：第二组中检测到表达该基因的细胞所占的百分比
    # p_val_adj：bonferroni多重检验校正后得到的校正后的P值。
    
    markers <- FindMarkers(
      object = seurat_integrated, # 记得改
      ident.1 = ident_1,
      ident.2 = ident_2,
      min.cells.group = 1,
      logfc.threshold = 0.25
    )
    
    # 根据组合名称添加到DEG列表中
    ALL[[paste(ident_1, ident_2, sep = "-vs-")]] <- markers
  }
}

# 这时候你可以改一下list的名字，不然名称太长，xlsx保存不了
names(ALL)
names(ALL) <- c("Astro","Endo","Epen","Micro","Neuron","Oligo","OPCs")

ALL <- lapply(ALL, function(x){
  x$SYMBOL <- rownames(x)
  x <- x[order(x$avg_log2FC, decreasing = T), ]
  return(x)
})
```
#### 1.2 提取差异基因
```r
DEG <- lapply(ALL, function(x){ 
  x <- subset(x, p_val < 0.05 & abs(avg_log2FC) >= 0.5) 
  x <- x[order(x$avg_log2FC, decreasing = T), ]
  return(x)
})

save(ALL, DEG, file = paste0(path, "DEG_list.RData"))

write.xlsx(DEG, file = paste0(path, "DEG_list_log_05.xlsx"))
```
#### 1.3 差异基因火山图 ####
```r
path = 'F:/R work/mmbrain/results/'
load(paste0(path, "DEG_list.RData"))

ALL2 <- lapply(ALL, function(x){
  x$cluster <- NA
  x$gene <- x$SYMBOL
  return(x)
})

for (i in 1:length(ALL2)) { ALL2[[i]]$cluster <- names(ALL2)[i] }

gap <- length(combinations[[1]])
data1 <- Reduce(rbind, ALL2[c(seq(1, length(ALL2), gap))])

a1 <- jjVolcano(diffData = data1, pSize = 1.5, log2FC.cutoff = 0.5,
                tile.col = corrplot::COL2('RdBu', 15)[c(4:12)], 
                topGeneN = 10,aesCol = c("#bbe0f2", "#f894af"),
                polar = F) +
  geom_hline(yintercept = 0.5, linetype = 'dotted') +
  geom_hline(yintercept = -0.5, linetype = 'dotted') + NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank())

ggplot2::ggsave(paste0(path, "火山图.pdf"), plot = a1,
                height = 12, width = 13, dpi = 300, limitsize = FALSE)
```

#### 1.4 计算差异基因数量柱形图 ####
```r

# 统计差异基因列表的每个dataframe行数
deg_num <- list()
for (i in 1:length(DEG)) { deg_num[[i]] <- length(DEG[[i]]$SYMBOL) }

names(deg_num) <- names(ALL)

# 调整柱形图展示的顺序
cluster_names <- factor(c('Neuron', 'Astrocyte', 'Endothelial_cell',
                           'Ependymal_cell', 'Microglial_cell', 'Oligodendrocyte', 'OPCs'), 
                        levels = c('Neuron', 'Astrocyte', 'Endothelial_cell',
                                    'Ependymal_cell', 'Microglial_cell', 'Oligodendrocyte', 'OPCs'))

# 需要展示的对比
cpname <- c('Tre1 vs CTRL', 'Tre2 vs CTRL')

deg_num2 <- as.data.frame(unlist(deg_num))
# 重复cpname多少次（细胞亚群的数量）
deg_num2$name <- rep(cpname, length(cluster_names))
# 重复cluster_names多少次（cpname的数量）
deg_num2$name <- factor(deg_num2$name, levels = cpname)

deg_num2$celltype <- rep(cluster_names, each = length(cpname))

colnames(deg_num2)[1] <- 'Value'
rownames(deg_num2) <- 1:nrow(deg_num2)

head(deg_num2)

ggplot(deg_num2, aes(x=celltype, y=Value, fill=celltype, label = Value)) +
  geom_bar(stat="identity") + # 使用实际值作为高度
  geom_text(position = position_stack(vjust = 0.8), color="black") + 
  facet_wrap(~name) + # 每个name值一个分面
  labs(title="DEG Numbers by Cell Type and Condition Comparison",
       x=NULL,
       y=NULL) + 
  theme_bw() + 
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c(hue_pal()( length(cluster_names) ))) # 颜色为细胞亚群数量

ggplot2::ggsave(paste0(path, "差异表达基因数量统计图.pdf"),
                height = 6, width = 8, dpi = 300, limitsize = FALSE)
```

#### 1.5 探讨FindMarkers是如何进行log2FC计算的 ####    
参考[关于 FindMarkers与AverageExpression 两个函数的差异](http://www.bio-info-trainee.com/8707.html)    
该文章其中一条函数挺有用的，可以用来观察某个基因的表达量，在某个细胞中的频数    
```r
hist(as.matrix(seurat_integrated@assays$RNA@data)['Ccr5', Idents(seurat_integrated)=='BAMs'])
```

---
## 12.GO 富集分析 ####
### 差异表达基因分开上下调 ####    
```r
rna <- c(lapply(DEG, function(x){x <- x[x$avg_log2FC > 0,]}),
          lapply(DEG, function(x){x <- x[x$avg_log2FC < 0,]}))

names_up <- character()
names_down <- character()

for (i in 1:length(DEG)) {names_up <- c(names_up, paste0('UP ', names(DEG)[i]))}
for (i in 1:length(DEG)) {names_down <- c(names_down, paste0('DOWN ', names(DEG)[i]))}

names(rna) <- c(names_up, names_down)

rna <- lapply(rna, function(x){x <- x$SYMBOL})
```

### 差异表达基因不分开上下调 ####    
```r
rna <- lapply(DEG, function(x){x <- x$SYMBOL})
```

### GO #### 
```r
BP <- clusterProfiler::compareCluster(rna, fun = "enrichGO", ont = "BP", 
                                      OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', readable = T)

CC <- clusterProfiler::compareCluster(rna, fun = "enrichGO", ont = "CC", 
                                      OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', readable = T)

MF <- clusterProfiler::compareCluster(rna, fun = "enrichGO", ont = "MF", 
                                      OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', readable = T)

save(BP, CC, MF, file = paste0(path, "DEG_GO.RData"))
```

#### 假如有两个组进行比较 ####
比如说 tre2 vs con、 tre1 vs con， 以下函数可以同时展示两个组中，-log10pvalue前20的通路    

```r
MultiPathway <- function(x, num = 20) {
  # 分组并提取每组前20个结果
  GO_2 <- lapply(split(x@compareClusterResult, x@compareClusterResult$Cluster), function(x) {
    head(x, num)
  })
  
  # 合并所有结果
  GO_3 <- Reduce(rbind, GO_2)
  
  # 提取基因比例中的总基因数
  GO_3$rati <- as.numeric(str_split_fixed(GO_3$GeneRatio, '/', n = 2)[, 2])
  
  # 计算比例
  GO_3$ratio <- (GO_3$Count / GO_3$rati) * 100
  
  # 重命名列
  colnames(GO_3)[1] <- 'group'
  
  # 计算 -log10(pvalue)
  enrich2 <- GO_3 %>%
    mutate(log10pvalue = -log10(pvalue))
  
  # 处理描述字段，使其适合显示在图中
  enrich2 <- enrich2 %>% mutate(Description = str_wrap(Description, width = 25))
  
  # 按 -log10(pvalue) 排序
  enrich2 <- enrich2 %>% arrange(desc(log10pvalue))
  
  # 根据group数量取颜色
  color <- c(hue_pal()(length(unique(GO_3$group))))
  
  # 创建柱形图
  p <- ggplot(enrich2, aes(x = reorder(Description, log10pvalue), y = log10pvalue, fill = group)) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
    coord_flip() +
    labs(x = "Pathway Description", y = '-log10(pvalue)') +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) +
    geom_text(aes(label = paste0(round(ratio, 2), "%")), vjust = 0.5, size = 3) +
    ggplot2::scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) + 
    scale_fill_manual(values = color)
  
  return(p)
}

MultiPathway(BP)
```
