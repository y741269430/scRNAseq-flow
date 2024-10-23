# 差异表达分析及GO绘图

## 目录 ####
- 1.差异表达分析 
- 2.GO 富集分析 

## 1.差异表达分析 ####
```r
# 加载文件
path = 'F:/R work/mmbrain/results/'
seurat_integrated <- readRDS(paste0(path, "seurat_integrated_2.rds"))

# 将细胞类型与样本名称进行拼接，即：celltype_sample
seurat_integrated$celltype.exp <- paste(Idents(seurat_integrated), seurat_integrated$sample, sep = "_")

Idents(seurat_integrated) <- seurat_integrated$celltype.exp

table(seurat_integrated@meta.data$celltype.exp)
```
#### 制作一个contrast的list，用来做循环
```r
# 定义需要进行差异表达分析的细胞类型组合
degmeta <- data.frame(table(seurat_integrated@meta.data$celltype.exp))

# 为什么这里有这句命令，假设我们有2个样本，那么这里的`each`就填2，意思就是每个一celltype，用一个数字来表示它。  
degmeta$cluster <- rep(1:length(unique(seurat_integrated$celltype)), each = 2)

```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/degmeta.png" width="600" />

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
#### 提取差异基因
```r
DEG <- lapply(ALL, function(x){ 
  x <- subset(x, p_val < 0.05 & abs(avg_log2FC) >= 0.5) 
  x <- x[order(x$avg_log2FC, decreasing = T), ]
  return(x)
})

save(ALL, DEG, file = paste0(path, "DEG_list.RData"))

write.xlsx(DEG, file = paste0(path, "DEG_list_log_05.xlsx"))
```
#### 差异基因火山图 ####
```r
path = 'F:/R work/mmbrain/results/'
load(paste0(path, "DEG_list.RData"))

ALL <- lapply(ALL, function(x){
  x$cluster <- NA
  x$gene <- x$SYMBOL
  return(x)
})

for (i in 1:length(ALL)) { ALL[[i]]$cluster <- names(ALL)[i] }

data1 <- Reduce(rbind, ALL[c(seq(1, length(ALL), gap))])

a1 <- jjVolcano(diffData = data1, pSize = 1.5, log2FC.cutoff = 0.5,
                tile.col = corrplot::COL2('RdBu', 15)[c(4:12)], 
                topGeneN = 10,aesCol = c("#bbe0f2", "#f894af"),
                polar = F) +
  geom_hline(yintercept = 0.5, linetype = 'dotted') +
  geom_hline(yintercept = -0.5, linetype = 'dotted') + NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank())

ggplot2::ggsave(paste0(path, "Volcano.pdf"), plot = a1,
                height = 12, width = 13, dpi = 300, limitsize = FALSE)
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

### 差异表达基因不分开上下调 #### 
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
