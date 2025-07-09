# pyscenic预测转录因子调控网络
## 目录 ####
- 1.从seurat中选择基因
- 2.把seurat转成loom
- 3.输入loom构建调控网络

## 1.从seurat中选择基因 ####
```r
# 输入基因
gene = c()
# 取出细胞亚群
subsets_cell <- subset(seurat_integrated, idents = c('celltype'))
# 取出矩阵
matrix_cell <- subsets_cell@assays$RNA@counts
# 选择基因
matrix_cell <- matrix_cell[rownames(matrix_cell) %in% gene, ]
matrix_cell <- as.data.frame(matrix_cell)
# 保存
write.csv(matrix_cell, 'test.csv')
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
python change.py test.csv scrna.loom
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
tfs = '/home/jjyang/downloads/SCENIC_db/v1/mm_mgi_tfs.txt'
feather = '/home/jjyang/downloads/SCENIC_db/v1/mm9-tss-centered-10kb-7species.mc9nr.feather'
tbl = '/home/jjyang/downloads/SCENIC_db/v1/motifs-v9-nr.mgi-m0.001-o0.0.tbl'

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
aucell_command = f"nohup pyscenic aucell {input_loom} {os.path.join(output_dir, 'reg.csv')} --output {os.path.join(output_dir, 'out_SCENIC.loom')} --num_workers 3"
run_step(aucell_command, "AUCELL step")

print("Script finished successfully.")
```
```bash
python scenic.py scrna.loom ./output
```


