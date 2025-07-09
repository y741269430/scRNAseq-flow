# pyscenic预测转录因子调控网络
## 目录 ####
- 1.从seurat中选择基因
- 2.把seurat转成loom
- 3.输入loom构建调控网络
- 4.提取 out_SCENIC.loom 信息
- 5.读取rds
- 6.选择转录因子进行展示
- 7.查看不同单细胞亚群的转录因子活性平均值
- 8.挑选各个单细胞亚群特异性的转录因子
- 9.绘制每个cluster top5的转录因子

---
参考
[Python版SCENIC转录因子分析（四）一文就够了](https://cloud.tencent.com/developer/article/2228252)    

## 0.数据库下载 ####
从此下载[Welcome to the cisTarget resources website!](https://resources.aertslab.org/cistarget/)    
该网站主要下载3个部分`databases`,`motif2tf`,`tf_lists`    
1.`tf_lists`里面，下载`allTFs_mm.txt`即可（直接wget）    
```bash
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt
```

2.`motif2tf`里面，小鼠只有`motifs-v9-nr.mgi-m0.001-o0.0.tbl` 和 `motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl`    
- v9: Annotations based on the 2017 cisTarget motif collection. Use these files if you are using the mc9nr databases.    
- v10: Annotations based on the 2022 SCENIC+ motif collection. Use these files if you are using the mc_v10_clust databases.
- 二选一下载(下载速度有点慢，直接IDM下了再传服务器吧)    
```bash
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
```

3.`databases`这里就比较复杂，v10版本分为[gene_based](https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/)和[region_based](https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/), 官网是这样解释的：    
- `refseq_r80`: Gene-based databases based on RefSeq 80 genes. To be used with (py)SCENIC and motif enrichment in gene sets.    
- `screen`: Region-based databases based on ENCODE SCREEN regions. To be used with pycisTarget/SCENIC+ and motif enrichment in region sets.    
- 两者差距太大了一个200多兆一个17G，这里我只选了`gene_based`
 
4.`gene_based`这里又细分了`scores`or`rankings`以及`TSS+/-10kb`or`500bpUp100Dw`，4种组合4种结果任君选择。。。    
Select motif database:    
- `scores`: Matrix containing motifs as rows and genes as columns and cluster-buster CRM scores as values. To be used with DEM.    
- `rankings`: Matrix containing motifs as rows and genes as columns and ranking position for each gene and motif (based on CRM scores) as values. To be used with cisTarget (R).    
The search space around the TSS of the gene in which the motif is scored is indicated in the database name:    
- `TSS+/-10kb`: 10kb around the TSS (total: 20kb).    
- `500bpUp100Dw`: 500bp upstream of TSS, and 100bp downstream.
- 这个官网教程[tutorial](https://pyscenic.readthedocs.io/en/latest/tutorial.html)推荐用`rankings`，至于距离的话`TSS+/-10kb`能更好地覆盖大多数 TF 结合位点，如研究启动子近端调控的话就选`500bpUp100Dw`    
我自己就下载`mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`(下载速度有点慢，直接IDM下了再传服务器吧)    
```bash
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather
```

## 1.从seurat中选择基因 ####

```r
readpath = './results/'
seurat_integrated <- readRDS(paste0(readpath, "seurat.rds"))

seurat_integrated <- subset(seurat_integrated, subset = 
                              (sample == 'C')| 
                              (sample == 'S')| 
                              (sample == 'R')) 

seurat.data <- subset(seurat_integrated, subset = (celltype == 'Endo'))

saveRDS(seurat.data, paste0(readpath, "seurat_Endo.rds"))

# 输入基因(只是示范)
load('./results/DEG.RData')
gene <- unique(c(ALL[[1]]$SYMBOL, ALL[[2]]$SYMBOL, ALL[[3]]$SYMBOL))

# 取出矩阵
matrix_cell <- seurat.data@assays$RNA@counts
# 选择基因
matrix_cell <- matrix_cell[rownames(matrix_cell) %in% gene, ]
matrix_cell <- as.data.frame(matrix_cell)

dir.create('./results/scenic/')

# 保存
write.csv(matrix_cell, './results/scenic/subset_scienic.csv')
```
## 2.把seurat转成loom ####
```
vim change.py
```
```python
import os
import sys
import loompy as lp
import numpy as np
import scanpy as sc

# 检查是否正确传入文件路径参数
if len(sys.argv) != 3:
    print("Usage: python <script_name> <input_csv_path> <output_loom_path>")
    sys.exit(1)

input_csv_path = sys.argv[1]  # 读取的CSV文件路径
output_loom_path = sys.argv[2]  # 生成的LOOM文件路径

# 读取 CSV 文件
x = sc.read_csv(input_csv_path)

# 转置数据
x = x.T  # 转置操作，细胞变基因，基因变细胞

# 创建 Loom 文件
row_attrs = {"Gene": np.array(x.var_names)}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create(output_loom_path, x.X.transpose(), row_attrs, col_attrs)

# 读取 Loom 文件
adata = sc.read_loom(output_loom_path)

# 打印基本信息到文件
with open("loom_info.txt", "w") as f:
    f.write(f"Shape: {adata.shape}\n")  # 表达矩阵的形状 (细胞数, 基因数)
    f.write(f"Var (genes): {adata.var_names[:5]}\n")  # 前5个基因名称
    f.write(f"Obs (cells): {adata.obs_names[:5]}\n")  # 前5个细胞ID
    f.write(f"Data file {output_loom_path} created successfully.\n")

print("Information has been written to loom_info.txt.")
```
```bash
python change.py ./results/scenic/subset_scienic.csv ./results/scenic/scrna.loom
```
## 3.输入loom构建调控网络 ####
```
vim scenic.py
```
```python
import os
import subprocess
import sys
import pyscenic

# 获取运行时输入的参数
if len(sys.argv) != 3:
    print("Usage: python <script_name> <input_loom_path> <output_dir>")
    sys.exit(1)

input_loom = sys.argv[1]  # 输入的 Loom 文件路径
output_dir = sys.argv[2]  # 输出文件的生成路径

# 确保输出路径存在
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 设置其他文件路径
tfs = '/home/jjyang/downloads/SCENIC_db/allTFs_mm.txt'
feather = '/home/jjyang/downloads/SCENIC_db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
tbl = '/home/jjyang/downloads/SCENIC_db/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'

# 创建一个日志文件记录每个步骤
log_file = os.path.join(output_dir, "scenic_script_log.txt")

def log_step(step_name):
    with open(log_file, 'a') as f:
        f.write(f"Starting {step_name}...\n")
        f.flush()

def run_step(command, step_name):
    log_step(step_name)
    try:
        # 使用 subprocess 运行命令
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(log_file, 'a') as f:
            f.write(f"{step_name} completed successfully.\n")
        return result
    except subprocess.CalledProcessError as e:
        with open(log_file, 'a') as f:
            f.write(f"{step_name} failed. Error: {e.stderr.decode()}\n")
        print(f"{step_name} failed. Check {log_file} for details.")
        sys.exit(1)

# 1. GRN step
grn_command = f"nohup pyscenic grn --num_workers 20 --output {os.path.join(output_dir, 'adj.sample.tsv')} --method grnboost2 {input_loom} {tfs} &"
run_step(grn_command, "GRN step")

# 2. Contextualization step
ctx_command = f"nohup pyscenic ctx {os.path.join(output_dir, 'adj.sample.tsv')} {feather} --annotations_fname {tbl} --expression_mtx_fname {input_loom} --mode 'dask_multiprocessing' --output {os.path.join(output_dir, 'reg.csv')} --num_workers 30 --mask_dropouts &"
run_step(ctx_command, "Contextualization step")

# 3. AUCell step
aucell_command = f"nohup pyscenic aucell {input_loom} {os.path.join(output_dir, 'reg.csv')} --output {os.path.join(output_dir, 'out_SCENIC.loom')} --num_workers 15"
run_step(aucell_command, "AUCELL step")

print("Script finished successfully.")
```
```bash
python scenic.py ./results/scenic/scrna.loom ./results/scenic/output
```
## 4.提取 out_SCENIC.loom 信息 ####
```r
loom <- open_loom('./results/scenic/output/out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)
```
## 5.读取rds ####
```r
seurat.data <- readRDS(paste0('F:/output/', "seurat_Endo.rds"))
seurat.data$seurat_clusters <- seurat.data$RNA_snn_res.0.4

sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]

dim(sub_regulonAUC)
seurat.data

#确认是否一致
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$seurat_clusters)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4]
```
## 6.选择转录因子进行展示 ####
```r
regulonsToPlot = c('Atf1(+)','Bach1(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))

p1 = DotPlot(seurat.data, features = unique(regulonsToPlot), group.by = 'seurat_clusters') + RotatedAxis()
p2 = RidgePlot(seurat.data, features = regulonsToPlot, ncol = 2, group.by = 'seurat_clusters') 
p3 = VlnPlot(seurat.data, features = regulonsToPlot, pt.size = 0, group.by = 'seurat_clusters')
p4 = FeaturePlot(seurat.data,features = regulonsToPlot)

wrap_plots(p1,p2,p3,p4)
```
## 7.查看不同单细胞亚群的转录因子活性平均值 ####
```r
# Split the cells by cluster:
selectedResolution <- "seurat_clusters" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,1])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)

regulonActivity_byGroup_Scaled = na.omit(regulonActivity_byGroup_Scaled)

Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
```
## 8.挑选各个单细胞亚群特异性的转录因子 ####
```r
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), 1]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
```
## 9.绘制每个cluster top5的转录因子 ####
```r
rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)
```









