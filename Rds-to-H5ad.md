# Rds-to-H5ad

## 转换rds 到 h5ad    
```r
readpath = 'F:/R work/单细胞结果/'
seurat <- readRDS(paste0(readpath, "seurat.rds"))

SaveH5Seurat(seurat, filename = paste0(readpath, "seurat.h5Seurat"))
Convert(paste0(readpath, "seurat.h5Seurat"), dest = "h5ad")
```
## 打开python 验证     
```python
import scanpy as sc
import pandas as pd
import os
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as pl

adata = sc.read_h5ad(r'F:\R work\单细胞结果\seurat.h5ad')
```

查看基本信息
```python
adata
adata.obs
adata.obsm
```
