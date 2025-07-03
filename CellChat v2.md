# CellChat    

## 目录
- 1.读取RDS，选择样本进行cellchat
- 2.Extract the inferred cellular communication network as a data frame
- 3.Compare the total number of interactions and interaction strength
- 4.Compare the number of interactions and interaction strength among different cell populations
- 5.Identify altered signaling with distinct interaction strength      
- 6.Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
- 7.导出通路里面的基因
- 8.查看基因所在通路的和弦图
- 9.Visually compare cell-cell communication using Chord diagram


参考    
[单细胞分析之细胞交互-3：CellChat](https://www.jianshu.com/p/b3d26ac51c5a)     
[CellChat Github](https://github.com/jinworks/CellChat)      
[Inference and analysis of cell-cell communication using CellChat](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html)      
[Comparison analysis of multiple datasets using CellChat](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html)     
---

## 1.读取RDS，选择样本进行cellchat   
```r
readpath <- 'results/'
seurat_integrated <- readRDS(paste0(readpath, "seurat_scRNAseq.rds"))

seurat_integrated2 <- subset(seurat_integrated, subset = c((sample == 'C') | (sample == 'S') | (sample == 'R')))

path = 'results/cellchat/'

table(seurat_integrated$sample)

cell_list <- list(C = subset(seurat_integrated, subset = (sample == 'C')),
                  S = subset(seurat_integrated, subset = (sample == 'S')),
                  R = subset(seurat_integrated, subset = (sample == 'R')) )
                  
# Part I: Data input & processing and initialization of CellChat object ####
cell_list <- lapply(cell_list, function(x){
  x <- createCellChat(x@assays$RNA@data, 
                      meta = x@meta.data, 
                      group.by = "celltype")
  x@DB <- CellChatDB.mouse
  return(x)
})

future::nbrOfWorkers()
future::nbrOfFreeWorkers()
future::plan("multicore", workers = 64) # do parallel
future::nbrOfWorkers()
```
使用以下脚本批量运行
```r
cell_list <- lapply(cell_list, function(x){
  # 1.3 Preprocessing the expression data for cell-cell communication analysis
  x <- CellChat::subsetData(x) # This step is necessary even if using the whole database
  x <- CellChat::identifyOverExpressedGenes(x)
  x <- CellChat::identifyOverExpressedInteractions(x)
  x <- CellChat::projectData(x, PPI.mouse)
  
  # Part II: Inference of cell-cell communication network ####
  # 2.1 Compute the communication probability and infer cellular communication network (慢) ####
  x <- CellChat::computeCommunProb(x, raw.use = TRUE, population.size = TRUE) 
  x <- CellChat::filterCommunication(x, min.cells = 10)
  
  # 2.2 Extract the inferred cellular communication network as a data frame ####
  # all the inferred cell-cell communications at the level of ligands/receptors
  # df.net <- CellChat::subsetCommunication(cellchat)
  # write.csv(df.net, paste0(path, "all_cellchat_cell_chat_D.csv"))
  
  # access the the inferred communications at the level of signaling pathways
  # df.net1 <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  # write.csv(df.net1, paste0(path, "all_cellchat_significant_cell_chat_D.csv"))
  
  # 2.3 Infer the cell-cell communication at a signaling pathway level ####
  x <- CellChat::computeCommunProbPathway(x)
  
  # 2.4 Calculate the aggregated cell-cell communication network ####
  x <- CellChat::aggregateNet(x)
  x <- CellChat::netAnalysis_computeCentrality(x, slot.name = "netP") 
})

cellchat <- mergeCellChat(cell_list, add.names = names(cell_list), cell.prefix = T)

save(cell_list, cellchat, file = paste0(path, "cellchat_all.RData"))
```
加载数据准备画图  
```r
load(paste0(path, "cellchat_all.RData"))
col_sample <- hue_pal()(3)

levels(cellchat@idents$joint)
names(cell_list)
```
## 2. Extract the inferred cellular communication network as a data frame
```r
df.net <- CellChat::subsetCommunication(cellchat)
```
## 3. Compare the total number of interactions and interaction strength
```r
a1 <- compareInteractions(cellchat, show.legend = F, group = c(1:3), measure = 'count', color.use = col_sample[1:3])
a2 <- compareInteractions(cellchat, show.legend = F, group = c(1:3), measure = 'weight', color.use = col_sample[1:3])

a1 <- data.frame(a1$data)
a2 <- data.frame(a2$data)

data <- cbind(a1, 'weight' = a2$count)

ggplot(data, aes(x = dataset)) +
  geom_bar(aes(y = count, fill = dataset), stat = "identity", alpha = 0.8, width = 0.4, position = position_nudge(x = -0.2)) +
  geom_text(aes(y = count, label = count), vjust = -0.5, position = position_nudge(x = -0.2), size = 3) +
  geom_bar(aes(y = weight * (max(count) / max(weight)), color = dataset), stat = "identity", fill = NA, size = 1, width = 0.4, position = position_nudge(x = 0.2), alpha = 0.8) +
  geom_text(aes(y = weight * (max(count) / max(weight)), label = weight), vjust = -0.5, position = position_nudge(x = 0.2), size = 3) +
  geom_hline(yintercept = 0, color = "black", size = 1) +  # 在 y = 0 处添加 x 轴线
  scale_y_continuous(
    name = "Number of inferred interactions",
    sec.axis = sec_axis(~ . / (max(data$count) / max(data$weight)), name = "Inteaction strength")) +
  scale_fill_manual(values = col_sample) +
  scale_color_manual(values = col_sample) +
  labs(x = NULL, title = "Cell chat") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 10, face = "bold", vjust = 1.3),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none")

ggplot2::ggsave(paste0(path, "compareInteractions.pdf"),
                height = 6, width = 8, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/compareInteraction.jpg" width="400" />     

## 4. Compare the number of interactions and interaction strength among different cell populations    
4.1 The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where `red` (or `blue`) colored edges represent increased(or decreased) signaling in the second dataset compared to the first one. (比较不同细胞群之间的相互作用数量和强度-线图 Export as A4 (`weight` 可改成 `count`))     

```r
par(mfrow = c(1,3))
netVisual_diffInteraction(cellchat, comparison = c(1,2), weight.scale = T, measure = 'weight', title.name = 'S vs C')
netVisual_diffInteraction(cellchat, comparison = c(1,3), weight.scale = T, measure = 'weight', title.name = 'R vs C')
netVisual_diffInteraction(cellchat, comparison = c(3,2), weight.scale = T, measure = 'weight', title.name = 'S vs R')
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/netVisual_diffInteraction.jpg" width="600" />       

    
4.2 CellChat can also show differential number of interactions or interaction strength in greater details using a heatmap. The top colored bar plot represents the sum of each column of the absolute values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of each row of the absolute values (outgoing signaling). Therefore, the bar height indicates the degree of change in terms of the number of interactions or interaction strength between the two conditions. In the colorbar, `red` (or `blue`) represents increased (or decreased) signaling in the second dataset compared to the first one. (比较不同细胞群之间的相互作用数量和强度-热图 Export as A4 (`weight` 可改成 `count`))      

```r
b1 <- netVisual_heatmap(cellchat, measure = 'weight', comparison = c(1,2), title.name = 'S vs C')
b2 <- netVisual_heatmap(cellchat, measure = 'weight', comparison = c(1,3), title.name = 'R vs C')
b3 <- netVisual_heatmap(cellchat, measure = 'weight', comparison = c(3,2), title.name = 'S vs R')

b1+b2+b3 # pdf(4:11)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/netVisual_heatmap.jpg" width="600" />     

## 5.Identify altered signaling with distinct interaction strength    
5.1 Compare the overall information flow of each signaling pathway or ligand-receptor pair (比较每个信号通路的整体信息流)     
    
```r
h1 <- rankNet(cellchat, comparison = c(1,2), mode = 'comparison', stacked = T, do.stat = T, color.use = col_sample[c(1,2)], title = 'S vs C')
h2 <- rankNet(cellchat, comparison = c(1,3), mode = 'comparison', stacked = T, do.stat = T, color.use = col_sample[c(1,3)], title = 'R vs C') 
h3 <- rankNet(cellchat, comparison = c(3,2), mode = 'comparison', stacked = T, do.stat = T, color.use = col_sample[c(3,2)], title = 'S vs R')

h1+h2+h3 # pdf 14:20
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/rankNet.jpg" width="600" />    

将这些差异的通路取交集看看    
```r
a1 <- h1$data[h1$data$pvalues < 0.05, ]
a2 <- h2$data[h2$data$pvalues < 0.05, ]
a3 <- h3$data[h3$data$pvalues < 0.05, ]

data_ls <- list('S vs C' = a1$name, 
                'R vs C' = a2$name, 
                'S vs R' = a3$name)

ggvenn(data_ls, 
       fill_color = col_sample,
       stroke_color = 'white',
       text_size = 6)

ggsave(paste0(path, 'Venn.pdf'), height = 7, width = 10, dpi = 300, limitsize = FALSE)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/ggvenn.jpg" width="250" />     

5.2 Compare outgoing (or incoming) signaling patterns associated with each cell population
```r
# 定义函数
netAnalysis_signalingRole_heatmap_modified <- function(object, signaling = NULL,
                                                       pattern = c("outgoing", "incoming", "all"), slot.name = "netP",
                                                       color.use = NULL, color.heatmap = "BuGn", 
                                                       title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
                                                       cluster.rows = FALSE, cluster.cols = FALSE,
                                                       show_column_dend = FALSE, show_row_dend = FALSE,
                                                       col_order = NULL,show_annot = FALSE, name = "Heatmap1"){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"), 
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                             name = color.heatmap))))(100)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  mat[is.na(mat)] <- 0
  col_order0 <- seq_len(ncol(mat))
  if(is.null(col_order))
    col_order <- col_order0
  if(show_annot){
    ht1 = Heatmap(mat[,col_order], col = color.heatmap.use, 
                  na_col = "white", name = name(), bottom_annotation = col_annotation, 
                  top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                  cluster_columns = cluster.cols, row_names_side = "left", 
                  show_column_dend = show_column_dend, show_row_dend = show_row_dend,
                  row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                  column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                                                                             "cm"), height = unit(height, "cm"), column_title = title, 
                  column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                              fontface = "plain"), title_position = "leftcenter-rot", 
                                              border = NA, at = legend.break,
                                              legend_height = unit(20, "mm"),
                                              labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm")))
    
  }else{
    
    ht1 = Heatmap(mat[,col_order], col = color.heatmap.use, na_col = "white" ,
                  cluster_rows = cluster.rows, 
                  name = name,
                  cluster_columns = cluster.cols, row_names_side = "left", 
                  show_column_dend = show_column_dend, show_row_dend = show_row_dend,
                  row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                  column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"), column_title = title, 
                  column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),
                                              title_position = "leftcenter-rot", 
                                              border = NA, at = legend.break, legend_height = unit(20, "mm"),
                                              labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))
  }
  return(ht1)
  
}
```
```r
i = 1
col_sample <- hue_pal()(3)
# 读取上述rankNet结果
h1 <- rankNet(cellchat, comparison = c(1,2), mode = 'comparison', stacked = T, do.stat = T, color.use = col_sample[c(1,2)], title = 'S vs C')

# 设置阈值
top_paths <- h1$data[h1$data$pvalues < 0.05, ]

# 设置前后20行
top_paths2 <- c(tail(top_paths,20)$name,head(top_paths,20)$name)

# 绘图
ht1 = netAnalysis_signalingRole_heatmap_modified(cell_list[[i]],
                                                 pattern = "all",
                                                 name = "C",
                                                 signaling = top_paths2,
                                                 title = 'C(test)',
                                                 width = 8.5, height = 16,
                                                 color.heatmap = "GnBu",
                                                 cluster.rows = T,
                                                 show_column_dend = FALSE, 
                                                 show_row_dend = F,
                                                 cluster.cols = F, 
                                                 show_annot = FALSE,
                                                 font.size = 12)

ht2 = netAnalysis_signalingRole_heatmap_modified(cell_list[[i+1]], 
                                                 pattern = "all", 
                                                 name = "S",
                                                 signaling = top_paths2,
                                                 title = 'S(test)',
                                                 width = 8.5, height = 16,
                                                 color.heatmap = "GnBu",
                                                 cluster.rows = F,
                                                 show_column_dend = FALSE,
                                                 show_row_dend = FALSE,
                                                 font.size = 12)
ht1+ht2

pdf(paste0(path, "netAnalysis_signalingRole_heatmap.pdf"), width = 10, height = 10)
ht1 + ht2
dev.off()
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/netAnalysis_signalingRole_heatmap_modified.jpg" width="400" />   

## 6.Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
```r
netVisual_bubble(cellchat, 
                 sources.use = c(1),          #发射端细胞
                 targets.use = c(1:8),        #接收端细胞
                 angle.x = 45, 
                 comparison = c(1:3),         #组别
                 thresh = 0.01,               #阈值默认0.05
                 min.quantile = 0.7, 
                 remove.isolate = F,
                 #signaling = c("CCL","CXCL") #选择展示的通路
                 color.text = col_sample)     #配色

# 挑选配受体展示
pairLR.use <- as.data.frame(c('PSAP_GPR37L1', 'NRXN3_NLGN3', 'NRXN3_NLGN1', 'NRXN1_NLGN1', 'NEGR1_NEGR1'))
colnames(pairLR.use) <- 'interaction_name'

netVisual_bubble(cellchat,
                 sources.use = c(1),
                 targets.use = c(1:8),
                 comparison = c(1:3),
                 angle.x = 45,
                 min.quantile = 0.7,
                 remove.isolate = F,
                 color.text = col_sample,
                 pairLR.use = pairLR.use)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/netVisual_bubble.jpg" width="400" />

6.1 Identify dysfunctional signaling by comparing the communication probabities     
Moreover, CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. This can be done by specifying `max.dataset` and `min.dataset` in the function `netVisual_bubble`. The increased signaling means these signaling have higher communication probability (strength) in the second dataset compared to the first dataset. The ligand-receptor pairs shown in the bubble plot can be accessed via `gg1$data`.     
```r
gg1 <- netVisual_bubble(cellchat, 
                        sources.use = 1, 
                        targets.use = c(1:3),  
                        comparison = c(1,2), 
                        max.dataset = 2, 
                        title.name = "Increased signaling in S", 
                        angle.x = 45, 
                        remove.isolate = T)

gg2 <- netVisual_bubble(cellchat, 
                        sources.use = 1, 
                        targets.use = c(1:3),  
                        comparison = c(1,2), 
                        max.dataset = 1, 
                        title.name = "Decreased signaling in S", 
                        angle.x = 45, 
                        remove.isolate = T)

gg1 + gg2
```

## 7.导出通路里面的基因
```r
pathways.show = c('PSAP', 'NEGR', 'L1CAM')

gene <- list()

for (i in pathways.show) { gene[[i]] <- extractEnrichedLR(cellchat, signaling = i, geneLR.return = TRUE, enriched.only = T)$geneLR }

genes <- Reduce(c, gene)

Seurat::VlnPlot(seurat_integrated2, split.by = 'sample', group.by = 'celltype',
                genes[c(1:25)], pt.size = 0, flip = T, stack = T)+
  ylab(NULL) +
  xlab(NULL) +
  guides(fill = guide_legend(reverse = T))

# 或
plotGeneExpression(cellchat, signaling = "PSAP", split.by = "datasets", colors.ggplot = T, type = "violin")

```
## 8.查看基因所在通路的和弦图
`netVisual_chord_gene`和`netVisual_chord_cell`是不太一样的注意区别    
```r
pathways.show = c('PSAP', 'NEGR', 'L1CAM')

netVisual_chord_gene(cell_list[[1]], sources.use = c(1:8), targets.use = c(1:8), 
                     lab.cex = 0.5, legend.pos.y = 30, 
                     signaling = pathways.show[1], 
                     title.name = paste0(pathways.show[1], ' - C'))

netVisual_chord_gene(cell_list[[2]], sources.use = c(1:8), targets.use = c(1:8), 
                     lab.cex = 0.5, legend.pos.y = 30, 
                     signaling = pathways.show[1], 
                     title.name = paste0(pathways.show[1], ' - S'))

netVisual_chord_gene(cell_list[[3]], sources.use = c(1:8), targets.use = c(1:8), 
                     lab.cex = 0.5, legend.pos.y = 30, 
                     signaling = pathways.show[1], 
                     title.name = paste0(pathways.show[1], ' - R'))
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/netVisual_chord_gene.jpg" width="250" />

## 9.Visually compare cell-cell communication using Chord diagram
```r
pathways.show <- 'PTN'

par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(cell_list)) {
  netVisual_aggregate(cell_list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cell_list)[i]))
}
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/netVisual_aggregate_select.jpg" width="600" />

```r
par(mfrow = c(1,3), xpd=TRUE)
ht <- list()
for (i in 1:length(cell_list)) {
  ht[[i]] <- netVisual_heatmap(cell_list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(cell_list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm"))
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/netVisual_heatmap_select.jpg" width="600" />





                 
