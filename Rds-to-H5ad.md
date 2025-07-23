# Rds-to-H5ad

## 转换rds 到 h5ad    
```r
readpath = 'F:/R work/单细胞结果/'
seurat <- readRDS(paste0(readpath, "seurat.rds"))

# 这一步需要把细胞群转成字符串，因为h5ad文件会把因子转换成数值型，导致细胞群变成数字。
seurat$celltype <- as.character(seurat$celltype)
SaveH5Seurat(seurat, filename = paste0(readpath, "seurat.h5Seurat"))
Convert(paste0(readpath, "seurat.h5Seurat"), dest = "h5ad")
```
## 打开python 验证     
```python
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rc_context

adata = sc.read_h5ad(r'F:\R work\单细胞结果\seurat.h5ad')
```
查看基本信息
```python
adata
adata.obs
adata.obsm
```
指定你想要的因子顺序
```python
import pandas as pd

celltype_order = ['D1_MSN', 'D2_MSN','IN','Astro', 'Oligo','Micro', 'BAMs','Endo', 'OPCs']

# 转换为 category 类型并设置顺序
adata.obs['celltype'] = pd.Categorical(
    adata.obs['celltype'],
    categories=celltype_order,
    ordered=True  # 可选：是否为有序类别
)

print(adata.obs['celltype'].cat.categories)
print(adata.obs['celltype'].head())
```
画UMAP
```python
sc.pl.umap(adata, 
           color='celltype', 
           frameon=False, 
           title='UMAP projection of the dataset', 
           legend_loc='right margin',  
           legend_fontsize=12,         
           legend_fontoutline=2,
           palette='Paired',
           show=False)                 

#plt.savefig(r'F:\R work\单细胞结果\1_results_NAc\umap_plot2.pdf', format='pdf', bbox_inches='tight')
```
plt.show()



