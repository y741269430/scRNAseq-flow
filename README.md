# scRNAseq-flow

    library(SingleCellExperiment)
    library(Seurat)
    library(tidyverse)
    library(Matrix)
    library(scales)
    library(cowplot)
    library(RCurl)
    library(clustree)
    library(SingleR)
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(Scillus)
    library(ggpubr)
    library(DoubletFinder)
    library(openxlsx)
    library(MySeuratWrappers)
    library(ggsci)
    library(data.table)

---
#### 1.前期准备工作 ####

#### 调整R内允许对象大小的限制（默认值是500*1024 ^ 2 = 500 Mb）
    options(future.globals.maxSize = 500 * 1024 ^ 2)

#### 设置工作路径 ####
    readpath = 'F:/R work/mmbrain/'
    path = 'F:/R work/mmbrain/results/'

#### 创建一个函数来处理数据并创建 Seurat 对象 ####
    process_and_create_seurat <- function(data_dir, project) {
      # 读取数据
      count_matrix <- Read10X(data.dir = data_dir, gene.column = 2)
      
      # # 获取基因Symbol并处理NA值
      # gene_symbols <- mapIds(org.Rn.eg.db, keys = count_matrix@Dimnames[[1]], column = "SYMBOL", keytype = "ENSEMBL")
      # gene_symbols <- ifelse(is.na(gene_symbols), paste0("NA-", seq_along(gene_symbols)), gene_symbols)
      # count_matrix@Dimnames[[1]] <- gene_symbols
      
      # 创建 Seurat 对象
      seurat_obj <- CreateSeuratObject(counts = count_matrix,   
                                       min.cells = 3,           # 仅保留在至少3个细胞中有表达的基因（过滤低表达基因）
                                       min.features = 200,      # 仅保留具有至少200个基因表达的细胞（过滤低质量细胞）
                                       project = project)       # 设置项目名称（用于结果的标识）
      
      
      return(seurat_obj)
    }

#### 设置项目名称（工作路径+项目名称=数据存放的路径） ####
    projects <- c("A_con", "B_tre")

#### 创建一个空的 Seurat 对象列表 ####
    seurat_objects <- list()
    
    # 循环处理每个项目
    for (project in projects) {
      # 构建数据目录路径
      data_dir <- file.path(path0, project, 'filtered_feature_bc_matrix/')
      
      # 处理数据并创建 Seurat 对象
      seurat_obj <- process_and_create_seurat(data_dir, project)
      
      # 将 Seurat 对象添加到列表中
      seurat_objects[[project]] <- seurat_obj
    }
    
    # 检查矩阵行名是否为SYMBOL
    head(rownames(seurat_objects[[1]]@assays$RNA))
