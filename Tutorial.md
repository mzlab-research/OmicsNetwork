---
title: "OmicsNetWork"
author: "denghaoke"
date: "2025-09-18"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# R包使用教程：网络分析与可视化
```{r cars}
library(OmicsNetwork)
```

## I. 网络流程分析功能

### 1.1 特定组的稳定子网络分析 - stable_subnetwork函数

```{r pressure1, echo=TRUE,eval=FALSE}
# 加载示例数据
data("metabolite_data")

# 运行分析（为加速计算，nBoots设置为5）
subnet_result <- stable_subnetwork(
  metabolite_data$count_table,
  metabolite_data$samplelist,
  group_name = 'T',
  metabolite_data$annotation_table,
  nBoots = 5,
  stability_threshold = 0.4,
  bootnet_R_threshold = 0.6,
  species = 'hsa'
)

# 保存结果数据和图像
pipline_save(
  Network = subnet_result,
  outdir = "./",
  R_threshold = 0.6
)
```

### 1.2 实验组与对照组的稳定差异网络分析 - differential_network函数
```{r pressure2, echo=TRUE,eval=FALSE}
# 加载示例数据
data("metabolite_data")

# 运行分析（为加速计算，nBoots设置为5）
diffnet_result <- differential_network(
  count_table = metabolite_data$count_table,
  samplelist = metabolite_data$samplelist,
  annotation_table = metabolite_data$annotation_table,
  diff_table = metabolite_data$diff_table,
  compare_group = 'T:N',
  nBoots = 5,
  bootnet_R_threshold = 0.3,
  stability_threshold = 0.4,
  species = 'hsa'
)

# 保存结果数据和图像
pipline_save(
  Network = diffnet_result,
  outdir = "./",
  R_threshold = 0.6
)
```

### 1.3 实验组与对照组的稳定多重网络分析 - multiplex_network函数
```{r pressure3, echo=TRUE,eval=FALSE}
# 加载示例数据
data("protein_data")

# 运行分析（为加速计算，nBoots设置为5）
multinet_result <- multiplex_network(
  count_table = protein_data$count_table,
  samplelist = protein_data$samplelist,
  annotation_table = protein_data$annotation_table,
  diff_table = protein_data$diff_table,
  Interaction_table = protein_data$Interaction_table,
  compare_group = 'T:N',
  nBoots = 5,
  stability_threshold = 0.4,
  species = 'hsa'
)

# 保存结果数据和图像
pipline_save(
  Network = multinet_result,
  outdir = "./",
  R_threshold = 0.6
)
```

## II. 独立分析模块功能

### 2.1 稳定性分析 - run_corStability函数
```{r pressure4, echo=TRUE}
# 创建示例网络数据
data("stable_subnetwork_result")
data("metabolite_data")
filter_table <- stable_subnetwork_result@PreCor@filter_table
p_filter_table <- stable_subnetwork_result@PreCor@precor$edges |> dplyr::select(from, to, padj)

# 执行稳定性分析（为加速计算，nBoots设置为5）
StableNetwork <- run_corStability(
  count_table = filter_table, 
  annotation_table = metabolite_data$annotation_table,
  p_filter_table = p_filter_table,
  group_name = 'T',
  nBoots = 5,  
  stability_threshold = 0.4,
  bootnet_R_threshold = 0.6
)

# 可视化
network_show(Network = StableNetwork, plot_type = 'overall_network', show_node_legend = TRUE)
```

### 2.2 子网络聚类分析 - run_cluster函数
```{r pressure5, echo=TRUE}
# 创建示例网络数据
data("stable_subnetwork_result")
nodes <- stable_subnetwork_result@StableNetwork@bootnet_result_filter@bootnet_node
edges <- stable_subnetwork_result@StableNetwork@bootnet_result_filter@bootnet_edge

# 执行子网络聚类分析
subNetwork_results <- run_cluster(nodes = nodes, edges = edges, group_name = 'T')

# 可视化
# 总体子网络划分图
network_show(Network = subNetwork_results, plot_type = 'overall_cluster_network', show_node_legend = TRUE)

# 子网络显示：以子网络1为例
network_show(
  Network = subNetwork_results, 
  plot_type = 'sub_network', 
  subnetwork_name = c('subnet_1'),
  node_colortype = 'Class', 
  show_node_name = TRUE, 
  show_edge_legend = TRUE,
  show_node_legend = TRUE, 
  centrality_scatterplot = FALSE
)
```

### 2.3 分子互作网络分析 

#### 2.3.1 支持基于STRING数据库的小鼠（mmu）和人(hsa)的蛋白互作网络构建 - run_interaction_table和run_interaction_network函数
```{r pressure6, echo=TRUE}
# 整体互作表生成
interaction_table <- run_interaction_table(
  node_list = c("TP53", "BRCA1", "EGFR"),
  species = 'hsa',
  score_threshold = 600
)

Diff_anno <- data.frame(
  feature_ID = c("TP53", "BRCA1", "EGFR"),
  FC = c(1.3, 0.3, 1.4),
  State = c("Up", "Down", "Up"),
  Class = c("p", "p", "p")
)

# 差异特征互作网络构建
interaction_network <- run_interaction_network(
  Interaction_table = interaction_table,
  Diff_anno = Diff_anno,
  compare_group = "treatment:control"
)

# 可视化
network_show(
  Network = interaction_network,
  plot_type = 'interaction_network',
  show_node_name = TRUE
)
```


#### 2.3.2 也可以直接导入互作表格，支持其他类型的相互作用关系 - run_interaction_network函数
```{r pressure7, echo=TRUE}
data("protein_data")
data("multiplex_network_result")
interaction_data <- protein_data$Interaction_table
diff_anno <- multiplex_network_result@Diff_anno

interaction_network <- run_interaction_network(
  Interaction_table = interaction_data,
  Diff_anno = diff_anno,
  compare_group = 'T:N'
)

# 可视化
network_show(
  Network = interaction_network,
  plot_type = 'interaction_network',
  show_node_name = TRUE
)
```

### 2.4 节点布局一致的条件特异性相关性网络 - run_conditional_network函数
```{r pressure8, echo=TRUE}
# 创建示例网络数据
data("differential_network_result")
count_table = differential_network_result@PreData@count_table
samplelist = differential_network_result@PreData@samplelist
Diff_anno = differential_network_result@Diff_anno
node_list = differential_network_result@node_list

# 运行条件网络分析
Conditional_network <- run_conditional_network(
  count_table = count_table, 
  samplelist = samplelist,
  compare_group = 'T:N',
  Diff_anno = Diff_anno, 
  node_list = node_list,
  nBoots = 5, 
  stability_threshold = 0.4,
  bootnet_R_threshold = 0.3
)

# 可视化
# 实验组
network_show(
  Network = Conditional_network,
  plot_type = 'case_overall_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean'
)

# 对照组
network_show(
  Network = Conditional_network,
  plot_type = 'control_overall_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean'
)
```

### 2.5 节点布局一致的条件特异性复合网络 - run_conditional_multiplexnetwork函数
```{r pressure9, echo=TRUE}
# 创建示例网络数据
data("multiplex_network_result")
Interaction_network = multiplex_network_result@Interaction_network
Conditional_network = multiplex_network_result@Conditional_network

# 运行条件多重网络分析
Conditional_multiplexnetwork <- run_conditional_multiplexnetwork(
  Interaction_network = Interaction_network,
  Conditional_network = Conditional_network,
  compare_group = 'T:N'
)

# 可视化
# 实验组
network_show(
  Network = Conditional_multiplexnetwork,
  plot_type = 'case_multi_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean'
)

# 对照组
network_show(
  Network = Conditional_multiplexnetwork,
  plot_type = 'control_multi_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean'
)
```

### 2.6 差异网络分析 - run_diff_network函数

#### 2.6.1 使用条件网络对象执行差异分析
```{r pressure10, echo=TRUE}
# 创建示例网络数据
data("differential_network_result")
Conditional_network = differential_network_result@Conditional_network

# 运行差异网络分析
Differential_network <- run_diff_network(
  Conditional_network = Conditional_network,
  edge_FC_threshold = 1.2,
  edge_p_threshold = 0.05,
  compare_group = 'T:N'
)

# 可视化
network_show(
  Network = Differential_network,
  plot_type = 'diff_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Log2FC'
)
```

#### 2.6.2 使用多重网络对象执行差异分析
```{r pressure11, echo=TRUE}
# 创建示例网络数据
data("multiplex_network_result")
Conditional_network = multiplex_network_result@Conditional_network
Conditional_multiplexnetwork = multiplex_network_result@Conditional_multiplexnetwork

Differential_multiplexnetwork <- run_diff_network(
  Conditional_network = Conditional_network,
  Conditional_multiplexnetwork = Conditional_multiplexnetwork,
  edge_FC_threshold = 1.2,
  edge_p_threshold = 0.05,
  compare_group = 'T:N'
)

# 可视化
network_show(
  Network = Differential_multiplexnetwork,
  plot_type = 'diff_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Log2FC'
)
```

### 2.7 网络节点富集分析 

#### 2.7.1 总网络富集分析 - run_enrichment函数
```{r pressure12, echo=TRUE}
# 创建示例网络数据
data("stable_subnetwork_result")
subnetworks_nodes <- stable_subnetwork_result@StableNetwork@bootnet_result_filter@bootnet_node

# 运行富集分析
Enrichment <- run_enrichment(
  annotation_table_select = subnetworks_nodes,
  species = 'hsa',
  group_name = 'T'
)

# 可视化
network_show(Network = Enrichment, plot_type = 'enrichment')
```

#### 2.7.2 子网络富集分析
```{r pressure13, echo=TRUE}
# 创建示例网络数据
data("stable_subnetwork_result")
subnetworks_nodes <- lapply(
  stable_subnetwork_result@SubNetwork@subnetworks, 
  function(subnet) {
    subnet$nodes
  }
)

# 运行富集分析
Subnet_Enrichment <- run_enrichment(
  annotation_table_select = subnetworks_nodes,
  species = 'hsa',
  group_name = 'T'
)

# 可视化（示例：subnet_1）
network_show(
  Network = Subnet_Enrichment,
  plot_type = 'subnetwork_enrichment',
  subnetwork_name = c("subnet_1")
)
```

#### 2.7.3 差异总网络富集分析（复合总网络适用） - run_diff_enrichment函数
```{r pressure14, echo=TRUE}
data("differential_network_result")
Differential_network = differential_network_result@Differential_network

Enrichment <- run_diff_enrichment(
  Differential_network = Differential_network,
  species = 'hsa',
  enrichment_p_threshold = 0.05,
  compare_group = 'T:N'
)

# 可视化
network_show(Network = Enrichment, plot_type = 'enrichment')
```

#### 2.7.4 差异子网络富集分析（复合子网络适用） - run_diffsubnet_enrichment函数
```{r pressure15, echo=TRUE}
Differential_subnetwork = differential_network_result@Differential_subnetwork

DiffSubnet_Enrichment <- run_diffsubnet_enrichment(
  Differential_subnetwork = Differential_subnetwork,
  species = 'hsa',
  enrichment_p_threshold = 0.05,
  compare_group = 'T:N'
)

# 可视化（示例：subnet_1）
network_show(
  Network = DiffSubnet_Enrichment,
  plot_type = 'differential_subnetwork_enrichment',
  subnetwork_name = c("subnet_1")
)
```

# III. 绘图函数

## 3.1 整体流程结果数据和图像保存函数

### pipline_save函数，使用方法参考第I部分

## 3.2 可视化函数

### network_show函数，使用方法参考第II部分

