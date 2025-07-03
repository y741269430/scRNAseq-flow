# Cellchat    

## 目录

参考    
[单细胞分析之细胞交互-3：CellChat](https://www.jianshu.com/p/b3d26ac51c5a)     

---

## 1.读取RDS，选择样本进行cellchat    
```r
readpath <- 'results/'
seurat_integrated <- readRDS(paste0(readpath, "seurat_scRNAseq.rds"))

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
## 2.画图  
```r
load(paste0(path, "cellchat_all.RData"))
col_sample <- hue_pal()(3)

levels(cellchat@idents$joint)
names(cell_list)

# 导出互作信息
df.net <- CellChat::subsetCommunication(cellchat)
```

比较不同细胞群之间的相互作用数量和强度-线图 Export as A4 (weight 可改成 count)
```r
par(mfrow = c(1,3))
netVisual_diffInteraction(cellchat, comparison = c(1,2), weight.scale = T, measure = 'weight', title.name = 'S vs C')
netVisual_diffInteraction(cellchat, comparison = c(1,3), weight.scale = T, measure = 'weight', title.name = 'R vs C')
netVisual_diffInteraction(cellchat, comparison = c(3,2), weight.scale = T, measure = 'weight', title.name = 'S vs R')
```
比较不同细胞群之间的相互作用数量和强度-热图 Export as A4 (weight 可改成 count)
```r
b1 <- netVisual_heatmap(cellchat, measure = 'weight', comparison = c(1,2), title.name = 'S vs C')
b2 <- netVisual_heatmap(cellchat, measure = 'weight', comparison = c(1,3), title.name = 'R vs C')
b3 <- netVisual_heatmap(cellchat, measure = 'weight', comparison = c(3,2), title.name = 'S vs R')

b1+b2+b3 # pdf(4:11)
```
比较每个信号通路的整体信息流
```r
h1 <- rankNet(cellchat, comparison = c(1,2), mode = 'comparison', stacked = T, do.stat = T, color.use = col_sample[c(1,2)], title = 'S vs C')
h2 <- rankNet(cellchat, comparison = c(1,3), mode = 'comparison', stacked = T, do.stat = T, color.use = col_sample[c(1,3)], title = 'R vs C') 
h3 <- rankNet(cellchat, comparison = c(3,2), mode = 'comparison', stacked = T, do.stat = T, color.use = col_sample[c(3,2)], title = 'S vs R')

h1+h2+h3 # pdf 14:20
```
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


