# 基因表达UMAP图
使用UMAP对基因进行展示，尤其是当某些基因表达在非目标细胞中，我们把它们的表达量设置为0   


### 假如有3个基因需要展示  'Slc17a6', 'Slc32a1', 'Gad2'， 先选一个基因看看表达情况  
```r
input_gene = 'Slc17a6'
p1 <- scRNAtoolVis::FeatureCornerAxes(object = seurat_integrated, 
                                      reduction = 'umap',
                                      pSize = 0.5,
                                      features = input_gene,
                                      groupFacet = 'sample') + 
  coord_equal() + theme_bw() + 
  scale_colour_gradientn(colours = c('grey', 'firebrick3')); p1
```

### 假如它表达在非神经元细胞上，那我们把它表达量设置为0    
```r
# 标记神经元细胞（假设已知神经元的细胞群体标记）
object_cells <- WhichCells(object = seurat_integrated, ident = c('Glut'))

# 1. 提取 Slc17a6 的表达数据
expr <- FetchData(seurat_integrated, vars = input_gene)

# 2. 找到非神经元细胞，并将其表达值设置为0
expr[!Cells(seurat_integrated) %in% object_cells, input_gene] <- 0

# 3. 将修改后的表达矩阵重新赋值回 Seurat 对象
seurat_integrated[['RNA']]@data[input_gene, ] <- expr[,1]

# 4. 绘制图像
p2 <- scRNAtoolVis::FeatureCornerAxes(object = seurat_integrated, 
                                      reduction = 'umap',
                                      pSize = 0.5,
                                      features = input_gene,
                                      groupFacet = 'sample') + 
  coord_equal() + theme_bw() + 
  scale_colour_gradientn(colours = c('grey', 'firebrick3')); p2
```

### 比较看看前后有什么区别    
```r
p = plot_grid(p1, p2, ncol = 1)
ggplot2::ggsave(paste0(path, 'FeatureCornerAxes_', input_gene, '.pdf'), plot = p, 
                height = 8, width = 8, dpi = 300, limitsize = FALSE)
```
