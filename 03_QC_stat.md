# 03_QC_stat
统计以下质控指标       
|     指标英文名称      |                中文含义                 |
|:---------------------:|:---------------------------------------:|
|        Sample         |                 样本名                  |
|  Mean_Reads_per_Cell  |      平均读数/细胞 - 衡量测序深度       |
| Median_Genes_per_Cell |   中位数基因数/细胞 - 衡量检测灵敏度    |
|      Total_Cells      |    总细胞数 - 样本中检测到的细胞总数    |
|      Total_Reads      |     总读数 - 样本所有细胞的UMI总数      |
|  Mean_Genes_per_Cell  |  平均基因数/细胞 - 检测到的基因平均数   |
|       Min_Reads       |     最小读数 - 单个细胞的最小UMI数      |
|       Max_Reads       |     最大读数 - 单个细胞的最大UMI数      |
|       Min_Genes       | 最小基因数 - 单个细胞检测到的最少基因数 |
|       Max_Genes       | 最大基因数 - 单个细胞检测到的最多基因数 |


配置环境
```r
# 设置工作目录，批量创建文件夹储存结果
setwd(r"{D:\R work\GSE171169_RAW\}")

source('my_scRNAseq.R')

if (!dir.exists("3_QC_stat")) {
  dir.create("3_QC_stat")
}
```

1. 细胞质量指标统计_01_原始矩阵
```r
# 01_原始矩阵 ####
seurat_objects <- readRDS("1_QC_Files/seurat_objects.rds")

combined_meta <- do.call(rbind, lapply(seq_along(seurat_objects), function(i) {
  data.frame(Sample = ifelse(!is.null(names(seurat_objects)), 
                             names(seurat_objects)[i], paste0("Sample", i)),
             seurat_objects[[i]]@meta.data)
}))

qc_results1 <- generate_qc_report(combined_meta, "3_QC_stat/细胞质量指标统计_01_原始矩阵")

print(qc_results1)
knitr::kable(qc_results1, format = "markdown", align = 'c')
```
| Sample | Mean_Reads_per_Cell | Median_Genes_per_Cell | Total_Cells | Total_Reads | Mean_Genes_per_Cell | Min_Reads | Max_Reads | Min_Genes | Max_Genes |
|:------:|:-------------------:|:---------------------:|:-----------:|:-----------:|:-------------------:|:---------:|:---------:|:---------:|:---------:|
| 05d_N1 |       4220.10       |        1346.5         |    2542     |  10727492   |       1421.71       |    500    |   23352   |    270    |   4503    |
| 05d_N2 |       3904.03       |        1203.0         |    3484     |  13601657   |       1340.54       |    504    |   22986   |    205    |   4667    |
| 14d_N1 |       3594.61       |        1138.0         |    2583     |   9284876   |       1295.52       |   1179    |   24107   |    305    |   4600    |
| 14d_N2 |       4650.55       |        1501.0         |    2313     |  10756714   |       1542.66       |   1317    |   24921   |    341    |   5151    |    


2. 细胞质量指标统计_02_初次过滤
```r
# 02_初次过滤 ####
seurat_filter <- readRDS("1_QC_Files/seurat_filter.rds")

combined_meta_filter <- do.call(rbind, lapply(seq_along(seurat_filter), function(i) {
  data.frame(Sample = ifelse(!is.null(names(seurat_filter)), 
                             names(seurat_filter)[i], paste0("Sample", i)),
             seurat_filter[[i]]@meta.data)
}))

qc_results2 <- generate_qc_report(combined_meta_filter, "3_QC_stat/细胞质量指标统计_02_初次过滤")

print(qc_results2)
knitr::kable(qc_results2, format = "markdown", align = 'c')
```
| Sample | Mean_Reads_per_Cell | Median_Genes_per_Cell | Total_Cells | Total_Reads | Mean_Genes_per_Cell | Min_Reads | Max_Reads | Min_Genes | Max_Genes |
|:------:|:-------------------:|:---------------------:|:-----------:|:-----------:|:-------------------:|:---------:|:---------:|:---------:|:---------:|
| 05d_N1 |       4355.45       |         1404          |    2418     |  10531470   |       1459.72       |    504    |   23352   |    281    |   4503    |
| 05d_N2 |       4037.87       |         1252          |    3289     |  13280539   |       1379.08       |    513    |   22986   |    229    |   4667    |
| 14d_N1 |       3626.28       |         1150          |    2524     |   9152732   |       1304.48       |   1179    |   24107   |    305    |   4600    |
| 14d_N2 |       4667.62       |         1506          |    2288     |  10679524   |       1547.62       |   1317    |   24921   |    341    |   4558    |


3. 细胞质量指标统计_03_双细胞去除
```r
# 03_双细胞去除 ####
singlet_files <- list.files(path = "2_DoubletFinder", 
                            pattern = "seurat_singlet.rds$", 
                            recursive = TRUE, 
                            full.names = TRUE)

# 提取样本名作为列表名称
sample_names <- sapply(singlet_files, function(x) basename(dirname(x)))

# 批量读取
seurat_list <- lapply(singlet_files, readRDS)
names(seurat_list) <- sample_names

# 合并样本
seurat_merged <- merge(x = seurat_list[[1]], 
                       y = c(seurat_list[[2]],
                             seurat_list[[3]],
                             seurat_list[[4]]),
                       add.cell.id = sample_names)
rm(seurat_list); gc()

seurat_merged$Sample <- seurat_merged$orig.ident

qc_results3 <- generate_qc_report(seurat_merged@meta.data, "3_QC_stat/细胞质量指标统计_03_双细胞去除")

print(qc_results3)
knitr::kable(qc_results3, format = "markdown", align = 'c')
```
| Sample | Mean_Reads_per_Cell | Median_Genes_per_Cell | Total_Cells | Total_Reads | Mean_Genes_per_Cell | Min_Reads | Max_Reads | Min_Genes | Max_Genes |
|:------:|:-------------------:|:---------------------:|:-----------:|:-----------:|:-------------------:|:---------:|:---------:|:---------:|:---------:|
| 05d_N1 |       4257.85       |         1380          |    2371     |  10095355   |       1437.52       |    504    |   23352   |    281    |   4503    |
| 05d_N2 |       3785.71       |         1227          |    3202     |  12121851   |       1327.14       |    513    |   22986   |    229    |   4321    |
| 14d_N1 |       3552.80       |         1135          |    2473     |   8786067   |       1285.50       |   1179    |   22378   |    305    |   4310    |
| 14d_N2 |       4564.51       |         1496          |    2246     |  10251896   |       1525.12       |   1317    |   18218   |    341    |   3835    |

---
4. 细胞质量指标统计_整合指标并绘图    
```r
# 1. 数据整合
qc_results1$Stage <- "01_Raw_Data"
qc_results2$Stage <- "02_Filtered_Data"
qc_results3$Stage <- "03_Final_Data"

combined_qc <- bind_rows(qc_results1, qc_results2, qc_results3) %>%
  mutate(Stage = factor(Stage, levels = c("01_Raw_Data", "02_Filtered_Data", "03_Final_Data")))

write.xlsx(combined_qc, '3_QC_stat/combined_qc.xlsx')
write.table(combined_qc, '3_QC_stat/combined_qc.txt', quote = F, sep = '\t', row.names = F)
```
可视化图表1
```r
## 图1：平均读数/细胞的变化趋势
p1 <- ggplot(combined_qc, aes(x = Stage, y = Mean_Reads_per_Cell, 
                              group = Sample, color = Sample)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Mean Reads per Cell Trend", 
       x = "", y = "Mean Reads per Cell") +
  theme_bw() +
  theme(legend.position = "top")

## 图2：中位数基因数/细胞的变化趋势
p2 <- ggplot(combined_qc, aes(x = Stage, y = Median_Genes_per_Cell, 
                              group = Sample, color = Sample)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Median Genes per Cell Trend", 
       x = "", y = "Median Genes per Cell") +
  theme_bw() +
  theme(legend.position = "top")

## 图3：总细胞数变化
p3 <- ggplot(combined_qc, aes(x = Stage, y = Total_Cells, fill = Sample)) +
  geom_col(position = position_dodge(0.9)) +
  geom_text(aes(label = scales::comma(Total_Cells)), 
            position = position_dodge(0.9), 
            vjust = -0.1, size = 3, fontface = "bold") +
  labs(title = "Total Cells Count", x = "Processing Stage", y = "Total Cells") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma)

trend_plots <- (p1 | p2) / p3
print(trend_plots)

ggsave("3_QC_stat/QC_Combined.pdf", trend_plots, width = 12, height = 8, dpi = 300)
ggsave("3_QC_stat/QC_Combined.png", trend_plots, width = 12, height = 8, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/3_QC_stat/QC_Combined.png" width="600" />

可视化图表2
```r
## 图4：样本间比较（03_Final_Data）
final_data <- combined_qc %>% filter(Stage == "03_Final_Data")

p4 <- ggplot(final_data, aes(x = Sample, y = Mean_Reads_per_Cell, fill = Sample)) +
  geom_col() +
  geom_text(aes(label = round(Mean_Reads_per_Cell, 0)), vjust = -0.1) +
  labs(title = "Final Data: Mean Reads per Cell", x = "", y = "Mean Reads per Cell") +
  theme_bw()

p5 <- ggplot(final_data, aes(x = Sample, y = Median_Genes_per_Cell, fill = Sample)) +
  geom_col() +
  geom_text(aes(label = round(Median_Genes_per_Cell, 0)), vjust = -0.1) +
  labs(title = "Final Data: Median Genes per Cell", x = "", y = "Median Genes per Cell") +
  theme_bw()

final_comparison <- p4 / p5
print(final_comparison)

ggsave("3_QC_stat/QC_Final_Comparison.pdf", final_comparison, width = 10, height = 6, dpi = 300)
ggsave("3_QC_stat/QC_Final_Comparison.png", final_comparison, width = 10, height = 6, dpi = 300)
```
<img src="https://github.com/y741269430/scRNAseq-flow/blob/main/img/3_QC_stat/QC_Final_Comparison.png" width="600" />

---
目录树     
```r
fs::dir_tree("3_QC_stat", recurse = 2)
```
```bash
3_QC_stat
├── QC_Combined.pdf
├── QC_Combined.png
├── QC_Final_Comparison.pdf
├── QC_Final_Comparison.png
├── combined_qc.txt
├── combined_qc.xlsx
├── 细胞质量指标统计_01_原始矩阵.txt
├── 细胞质量指标统计_01_原始矩阵.xlsx
├── 细胞质量指标统计_02_初次过滤.txt
├── 细胞质量指标统计_02_初次过滤.xlsx
├── 细胞质量指标统计_03_双细胞去除.txt
└── 细胞质量指标统计_03_双细胞去除.xlsx
```

### 联系方式    
- 作者：JJYang
- 邮箱：y741269430@163.com
- 创建日期：2025-11-08
- 修改日期：2025-11-08
