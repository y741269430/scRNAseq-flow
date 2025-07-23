# Rds-to-H5ad

## 1. 转换rds 到 h5ad    
```r
readpath = 'F:/R work/单细胞结果/'
seurat <- readRDS(paste0(readpath, "seurat.rds"))

# 这一步需要把细胞群转成字符串，因为h5ad文件会把因子转换成数值型，导致细胞群变成数字。
seurat$celltype <- as.character(seurat$celltype)
SaveH5Seurat(seurat, filename = paste0(readpath, "seurat.h5Seurat"))
Convert(paste0(readpath, "seurat.h5Seurat"), dest = "h5ad")
```
## 2. 打开python 验证     
```python
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rc_context

os.getcwd()

outdir = 'F:/R work/单细胞结果/1_results_NAc/'
adata = sc.read_h5ad(outdir + 'seurat_NAc_early.h5ad')
```
查看基本信息
```python
adata
adata.obs
adata.obsm
```
指定你想要的因子顺序
```python
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

#plt.savefig(outdir + '00-UMAP_clusters.pdf', format='pdf', bbox_inches='tight')

plt.show()
```
## 3.使用python绘图   
创建细胞-marker字典
```python
marker_genes_dict = {
    'D1_MSN': ['Drd1'],
    'D2_MSN': ['Drd2'],
    'IN': ['Resp18'],
    'Astro': ['Gja1'],
    'Oligo': ['Tmem119'],
    'Micro': ['Ms4a7'],
    'BAMs': ['Mog'],
    'Endo': ['Pdgfra'],
    'OPCs': ['Cldn5'],
}
```
1.画散点图
```python
sc.pl.dotplot(adata, marker_genes_dict, 'celltype', dendrogram=True)
pl.savefig(outdir + "01-Dotplot_markers.pdf", format='pdf')
```
2.画小提琴图
```python
ax = sc.pl.stacked_violin(adata, marker_genes_dict, groupby='celltype', swap_axes=False, dendrogram=True)
pl.savefig(outdir + "02-stacked-violin.pdf", format='pdf') 
```
3.画矩阵热图
```python
sc.pl.matrixplot(adata, marker_genes_dict, 'celltype', dendrogram=True, cmap='Blues', standard_scale='var', colorbar_title='column scaled\nexpression')
pl.savefig(outdir + "03-matrixplot.pdf", format='pdf')   
```
4.画矩阵热图（标准化）
```python
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.pl.matrixplot(adata, marker_genes_dict, 'celltype', dendrogram=True, colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
#pl.savefig(outdir + "04-matrixplot-scaled.pdf", format='pdf')   
```
5.组合图
```python
import matplotlib.pyplot as pl
fig, (ax1, ax2, ax3, ax4) = pl.subplots(1, 4, figsize=(25,4), gridspec_kw={'wspace':0.9})
ax1_dict = sc.pl.dotplot(adata, marker_genes_dict, groupby='celltype', ax=ax1, show=False)
ax2_dict = sc.pl.stacked_violin(adata, marker_genes_dict, groupby='celltype', ax=ax2, show=False)
ax3_dict = sc.pl.matrixplot(adata, marker_genes_dict, groupby='celltype', ax=ax3, show=False, cmap='viridis')
ax4_dict = sc.pl.matrixplot(adata, marker_genes_dict, groupby='celltype', ax=ax4, colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

#pl.savefig(outdir + "05-plot_combined.pdf", format='pdf') 
plt.show()
```
