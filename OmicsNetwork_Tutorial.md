# OmicsNetWork: Multi-Omics Data Network Analysis R Package

OmicsNetWork is a comprehensive toolkit specifically designed for multi-omics data network analysis, aimed at helping researchers identify key functional modules and regulatory networks from complex biological data. This package provides a one-stop solution from **data filtering**, **molecular interaction network construction**, **condition-specific correlation network construction**, **correlation stability assessment**, **multiple network integration**, **correlation difference testing**, to **functional enrichment analysis**, with particular emphasis on network structure stability and biological significance mining.

Its core advantage lies in integrating various analytical methods such as **condition-specific correlation network analysis**, **differential network comparison**, and **multiple network integration**, and provides rich visualization functions, enabling researchers to intuitively explore and understand the hidden biological patterns and mechanisms in omics data.

## Table of Contents

[TOC]

## Usage Conditions
* R version >= 3.4.3

```R
# Install the package
devtools::install_local("OmicsNetwork_0.0.0.9000.tar.gz")
# Load the package
library(OmicsNetwork)
```
## I. Network Pipeline Analysis Functions
### 1.1 Stable Subnetwork Analysis for Specific Group - stable_subnetwork Function
Stable Subnetwork Analysis is a function that executes a complete stable subnetwork analysis pipeline, focusing on identifying and parsing highly connected stable network modules under a single experimental condition. This function implements end-to-end analysis from raw data to functional annotation.

#### 1.1.1 Analysis Pipeline
**Main Pipeline:**
  
  This function performs a complete stable subnetwork analysis pipeline for a specific experimental group (e.g., treatment group 'T'), including the following analyses:
  
* **Data Preparation:** Check column names, filter molecular quantitative values and sample table for the group (T) based on the group name.
* **Correlation Calculation:** Calculate inter-molecular associations using the specified correlation method (e.g., spearman). If there are more than 1000 molecules, select the top 1000 molecules by connectivity for subsequent analysis.
* **Network Stability Assessment:** Evaluate network connection stability through Bootstrap resampling, identifying network connections that stably exist across different resamplings.
* **Molecule Selection:** Filter important molecules based on correlation stability and significance indicators.
* **Subnetwork Identification:** Extract subnetwork modules from the total network, partitioning modules based on connection tightness.
* **Enrichment Analysis:** Perform functional enrichment analysis on the overall network and each subnetwork.
  
**Main Parameter Description:**

- count_table: Molecular expression table, must contain a feature_ID column.
- samplelist: Sample information table, containing sample and group columns.
- group_name: Analysis group name (e.g., 'T').
- annotation_table: Node annotation information table (feature_ID, Class, KEGG.ID).
- filter_num: Maximum molecule limit for the network, default 1000 (if exceeded, filter by connectivity).
- nBoots: Number of Bootstrap resampling iterations (default 500, set to 5 in the example for faster demonstration).
- stability_threshold: Stability threshold, used to filter stable edges based on confidence interval range (CIrange).
- bootnet_R_threshold: Correlation coefficient threshold for filtering weak correlation edges.
- bootnet_p_threshold: Correlation p-value threshold, default 0.05.
- species: Species abbreviation (e.g., "hsa" for human, "mmu" for mouse).

```R
# Load small example data
data("metabolite_data")
# Or read data
metabolite_data=readRDS("metabolite_data.rds")

# Run analysis (nBoots set to 5 for faster computation)
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

# Save result data and images. If the data volume is too large, set the StabilityTest parameter to FALSE to speed up saving by not outputting stability test results.
pipline_save(
  Network = subnet_result,
  outdir = "./",
  R_threshold = 0.6
)
```

#### 1.1.2 Partial Result Display

```R
# Network clustering and function
cluster_anno_plot1<-network_show(Network=subnet_result,plot_type="overall_cluster_network", 
                                 add_enrichement=TRUE,show_node_legend =TRUE, palette_style = "morandi")

cluster_anno_plot1@plot
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Network_Clustering_T.png)

```R
# View subnetwork 1 (node size represents betweenness)
subnet_betweennessplot1<-network_show(Network=subnet_result,plot_type="sub_network",
                                      subnetwork_name=c("subnet_1"),node_colortype="Class",show_node_name = TRUE,
                                      show_edge_legend=TRUE,show_node_legend=TRUE,add_Centrality="betweenness",
                                      centrality_scatterplot=FALSE,node_name_size = 3, palette_style = "morandi")

subnet_betweennessplot1@plot

```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Betweenness_subnet_1_T.png)

```R
# Topological analysis of subnetwork 1
subnet_centrality<-network_show(Network=subnet_result,plot_type="sub_network"
                                ,subnetwork_name=c("subnet_1"), node_colortype="Class",show_node_name = TRUE,
                                show_edge_legend=TRUE,node_name_size=3,show_node_legend=TRUE,
                                add_Centrality=c("betweenness","degree","eigenvector"), palette_style = "morandi")

subnet_centrality@plot$subnet_1

```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Centrality_subnet_1_T.png)

### 1.2 Stable Differential Network Analysis for Experimental vs Control Groups - differential_network Function

**Function Overview:**
  differential_network is a function that performs comprehensive differential network analysis, implementing a complete analysis pipeline from data preparation to functional enrichment. This function not only identifies expression differences at the node level but, more importantly, reveals condition-specific changes in network topology, providing in-depth insights into the functional restructuring of biological processes.

#### 1.2.1 Analysis Pipeline
**Main Pipeline:**
  
This function compares two groups of samples (e.g., treatment T vs control N) and performs the following analyses:
  
* **Data Preparation:** Check column names, filter molecular quantitative values and sample table for the two groups based on the comparison scheme (T:N).
* **Differential Analysis:** Identify molecules with significant expression differences between the two groups, filtering differentially expressed molecules based on fold change (FC) and statistical significance (p-value/q-value) thresholds.
* **Conditional Network Construction:** Build correlation networks for the differential molecules separately for the experimental and control groups, calculating inter-molecular associations using the specified correlation method (e.g., spearman).
* **Network Stability Assessment:** Evaluate the connection stability of the two group networks through Bootstrap resampling, identifying network connections that stably exist across different resamplings.
* **Differential Network Analysis:** Systematically compare the two group networks, identifying different types of differential connections: significantly enhanced connections, significantly weakened connections, and group-specific connections (connections present only in one group).
* **Differential Subnetwork Identification:** Extract subnetwork modules from the differential network, partitioning modules based on connection tightness.
* **Enrichment Analysis:** Perform functional enrichment analysis separately for the overall differential network and for the molecular lists involved in each correlation difference type (enhanced/weakened/specific, etc.) within the differential subnetworks.

**Main Parameter Description:**
  
- count_table: Molecular expression table, must contain a feature_ID column.
- samplelist: Sample information table, containing sample and group columns.
- compare_group: Comparison group format (e.g., "T:N").
- annotation_table: Node annotation information table (feature_ID, Class, KEGG.ID).
- diff_table: Pre-calculated differential analysis result table; if not provided, it will be calculated automatically.
- nBoots: Number of Bootstrap resampling iterations (default 25, set to 5 in the example for faster demonstration).
- bootnet_R_threshold: Correlation coefficient threshold for filtering weak correlation edges.
- stability_threshold: Stability threshold for filtering stable edges.
- edge_FC_threshold: Fold change threshold for differential connections, default 1.2.
- edge_p_threshold: Significance threshold for differential connections, default 0.05.
- species: Species abbreviation (e.g., "hsa" for human, "mmu" for mouse).

```R
# Load small example data
data("metabolite_data")
# Or read data
metabolite_data=readRDS("metabolite_data.rds")

# Run analysis (nBoots set to 5 for faster computation)
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

# Save result data and images. If the data volume is too large, set the StabilityTest parameter to FALSE to speed up saving by not outputting stability test results.
pipline_save(
  Network = diffnet_result,
  outdir = "./",
  R_threshold = 0.3
)
```

#### 1.2.2 Partial Result Display

```R
# read data
diffnet_result=readRDS("differential_network_result.rds")

# Conditional differential network
diff_net_top_plot<-network_show(Network=diffnet_result,
  plot_type="differential_network",
  node_colortype="Log2FC",focus=c("all"),
  show_edge_legend = TRUE,show_node_legend = TRUE,
  show_node_name = FALSE,node_size=5,node_name_size=3, palette_style = "morandi"
 )

diff_net_top_plot@plot
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Differential_network_1_T-vs-N.png)

```R
# read data
diffnet_result=readRDS("differential_network_result.rds")

# Conditional differential subnetwork
diff_subnet_top_plot<-network_show(Network=diffnet_result,
                                   plot_type="differential_subnetwork",
                                   node_colortype="Log2FC",focus=c("all"),subnetwork_name = c("subnet_1"),
                                   show_edge_legend = TRUE,show_node_legend = TRUE,
                                   show_node_name = TRUE,node_size=5,node_name_size=3, palette_style = "morandi"
)

diff_subnet_top_plot@plot
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Differential_subnetwork_1_T-vs-N.png)

### 1.3 Stable Multiplex Network Analysis for Experimental vs Control Groups - multiplex_network Function
**Function Overview:**
  Multiplex Network Analysis is a function that performs comprehensive multiplex network analysis by integrating molecular interaction information and expression correlation data to construct multi-level biological networks and perform systematic comparative analysis.

#### 1.3.1 Analysis Pipeline
**Main Pipeline:**
  
This function integrates molecular interaction information with expression correlation to perform multiplex network analysis:
  
* **Data Preparation:** Check column names, filter molecular quantitative values and sample table for the two groups based on the comparison scheme (T:N).
* **Differential Analysis:** Identify molecules with significant expression differences between the two groups, filtering differentially expressed molecules based on fold change (FC) and statistical significance (p-value/q-value) thresholds.
* **Molecular Interaction Network Construction:** Build a differential protein-protein interaction network based on the STRING database; also supports direct import of other types of molecular interaction tables.
* **Conditional Correlation Networks:** Build correlation networks for the differential molecules separately for the experimental and control groups, calculating inter-molecular associations using the specified correlation method (e.g., spearman).
* **Stability Assessment:** Evaluate the connection stability of the two group networks through Bootstrap resampling, identifying network connections that stably exist across different resamplings.
* **Network Integration:** Integrate the interaction network with the correlation networks to form a multiplex network.
* **Differential Multiplex Network Analysis:** Identify differences in the multiplex network between the experimental and control groups.
* **Functional Analysis:** Perform functional enrichment analysis separately for the overall differential multiplex network and for the molecular lists involved in each correlation difference type (enhanced/weakened/specific, etc.) within the differential subnetworks.

**Main Parameter Description:**
  
- count_table: Molecular expression table, must contain a feature_ID column.
- quantitative_table: Quantitative data table for differential analysis; if not provided, count_table is used.
- samplelist: Sample information table, containing sample and group columns.
- compare_group: Comparison group format (e.g., "T:N").
- annotation_table: Node annotation information table (feature_ID, Class, KEGG.ID).
- diff_table: Pre-calculated differential analysis result table.
- Interaction_table: Molecular interaction table; if not provided, it will be retrieved from the STRING database.
- score_threshold: Interaction score threshold, default 600.
- nBoots: Number of Bootstrap resampling iterations (default 50, set to 5 in the example for faster demonstration).
- bootnet_R_threshold: Correlation coefficient threshold for filtering weak correlation edges.
- stability_threshold: Stability threshold for filtering stable edges.
- edge_FC_threshold: Fold change threshold for differential connections, default 1.2.
- edge_p_threshold: Significance threshold for differential connections, default 0.05.
- species: Species abbreviation (e.g., "hsa" for human, "mmu" for mouse).

```R
# Load small example data
data("protein_data")
# Or read data
protein_data=readRDS("protein_data.rds")

# Run analysis (nBoots set to 5 for faster computation)
multinet_result <- multiplex_network(
  count_table = protein_data$count_table,
  samplelist = protein_data$samplelist,
  annotation_table = protein_data$annotation_table,
  diff_table = protein_data$diff_table,
  Interaction_table = protein_data$Interaction_table,
  compare_group = 'T:N',
  nBoots = 5,
  bootnet_R_threshold = 0.3,
  stability_threshold = 0.4,
  species = 'hsa'
)

# Save result data and images. If the data volume is too large, set the StabilityTest parameter to FALSE to speed up saving by not outputting stability test results.
pipline_save(
  Network = multinet_result,
  outdir = "./",
  R_threshold = 0.3
)
```

#### 1.3.2 Partial Result Display

```R
# read data
multinet_result=readRDS("multiplex_network_result.rds")

# Conditional differential multiplex subnetwork
diff_multisubnet_top_plot<-network_show(Network=multinet_result,
                                        plot_type="differential_network",node_colortype="Log2FC",focus=c("all"),
                                        subnetwork_name = c("subnet_1"),show_edge_legend = TRUE,show_node_legend = TRUE,
                                        show_node_name = TRUE,node_size=5,node_name_size=3, palette_style = "morandi")

diff_multisubnet_top_plot@plot
```
![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Differential_multiplex_network_T-vs-N.png)

## II. Independent Analysis Module Functions
### 2.1 Stability Analysis - run_corStability Function
**Function Overview:**
  Stability analysis is used to assess the robustness of correlation networks across different data subsets. Using the Bootstrap resampling method, it calculates stability indicators for each edge in the network over multiple resamplings, filtering out statistically stable network connections.

**Main Parameter Description:**
  
- count_table: Feature expression table, must contain a feature_ID column.
- annotation_table: Node annotation information table (feature_ID, Class, KEGG.ID).
- group_name: Analysis group name.
- nBoots: Number of Bootstrap resampling iterations (default 50, set to 5 in the example for faster demonstration).
- stability_threshold: Stability threshold, used to filter stable edges based on confidence interval range (CIrange).
- bootnet_R_threshold: Correlation coefficient threshold for filtering weak correlation edges.
- cor_method: Correlation calculation method ("spearman" or "cor").

**Usage Example:**
  
  ```R
# Create example network data
data("stable_subnetwork_result")
data("metabolite_data")

# Or read data
metabolite_data=readRDS("metabolite_data.rds")
stable_subnetwork_result=readRDS("stable_subnetwork_result.rds")

filter_table <- stable_subnetwork_result@PreCor@filter_table
annotation_table<-metabolite_data$annotation_table

# Perform stability analysis (nBoots set to 5 for faster computation)
StableNetwork <- run_corStability(
  count_table = filter_table, 
  annotation_table = annotation_table,
  p_filter_table = p_filter_table,
  group_name = 'T',
  nBoots = 5,
  stability_threshold = 0.4,
  bootnet_R_threshold = 0.6
)

# Visualization: Stable correlation network
network_show(Network = StableNetwork, plot_type = 'overall_network', show_node_legend = TRUE, palette_style = "morandi")
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_corStability_T.png)

### 2.2 Subnetwork Clustering Analysis - run_cluster Function
**Function Overview:**
  Subnetwork clustering analysis identifies functional modules within a network through community detection algorithms, decomposing complex biological networks into biologically meaningful subnetworks. It supports fast greedy and Louvain algorithms and automatically optimizes clustering resolution.

**Main Parameter Description:**
  
- nodes: Node data frame, must contain "node" and "Class" columns.
- edges: Edge data frame, must contain "from", "to", "cor", and "p_adjust" columns.
- group_name: Analysis group name.
- clustersize: Maximum community size threshold (default 25).
- diffmessage: Differential analysis type identifier ("NULL", "diff", "multi").

**Usage Example:**
  
```R
# Create example network data
data("stable_subnetwork_result")

# Or read data
stable_subnetwork_result=readRDS("stable_subnetwork_result.rds")

nodes <- stable_subnetwork_result@StableNetwork@bootnet_result_filter@bootnet_node
edges <- stable_subnetwork_result@StableNetwork@bootnet_result_filter@bootnet_edge

# Perform subnetwork clustering analysis
subNetwork_results <- run_cluster(nodes = nodes, edges = edges, group_name = 'T')

# Visualization
# Overall subnetwork partitioning plot
network_show(Network = subNetwork_results, plot_type = 'overall_cluster_network', show_node_legend = TRUE, palette_style = "morandi")
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_cluster_T-vs-N.png)

```R
# Subnetwork display: using subnetwork 1 as an example
network_show(
  Network = subNetwork_results,
  plot_type = 'sub_network',
  subnetwork_name = c('subnet_1'),
  node_colortype = 'Class',
  show_node_name = TRUE,
  show_edge_legend = TRUE,
  show_node_legend = TRUE,
  centrality_scatterplot = FALSE, palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_cluster_T-vs-N_subnet_1.png)

### 2.3 Molecular Interaction Network Analysis

#### 2.3.1 Protein-Protein Interaction Network Construction - run_interaction_table Function

**Function Overview:**
  Constructs protein-protein interaction networks based on the STRING database, supporting species such as mouse (mmu) and human (hsa). The function automatically handles ID mapping and interaction information retrieval.

**Main Parameter Description:**
  
- node_list: Target node list (gene IDs or symbols).
- species: Species identifier (hsa: human, mmu: mouse, etc.).
- score_threshold: STRING interaction score threshold (0-1000, default 600).
- database_path: Local storage path for the STRING database.

**Usage Example:**
  
  ```R
# Generate overall interaction table
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

# Build differential feature interaction network
interaction_network <- run_interaction_network(
  Interaction_table = interaction_table,
  Diff_anno = Diff_anno,
  compare_group = "treatment:control"
)

# Visualization
network_show(
  Network = interaction_network,
  plot_type = 'interaction_network',
  show_node_name = TRUE, palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_interaction_network_T-vs-N.png)

#### 2.3.2 General Molecular Interaction Network Construction - run_interaction_network Function

**Function Overview:**
  Supports importing custom interaction tables to build various types of molecular interaction networks (protein-protein, transcriptional regulation, metabolic pathways, etc.).

**Main Parameter Description:**
  
- interaction_data: Interaction table, must contain feature1, feature2 columns; optionally score, relationship columns.
- Diff_anno: Annotation table with differential information, must contain a Class column.
- compare_group: Comparison group (format: "T:N").

**Usage Example:**
  
  ```R
data("protein_data")
data("multiplex_network_result")

# Or read data
protein_data=readRDS("protein_data.rds")
multiplex_network_result=readRDS("multiplex_network_result.rds")

interaction_data <- protein_data$Interaction_table
diff_anno=multiplex_network_result@Diff_anno |>
  dplyr::filter(!(State %in% c("Non-significant","Non")))

interaction_network <- run_interaction_network(
  Interaction_table = interaction_data,
  Diff_anno = diff_anno,
  compare_group = 'T:N'
)

# Visualization
network_show(
  Network = interaction_network,
  plot_type = 'interaction_network',
  show_node_name = TRUE, palette_style = "morandi"
)
```
![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_interaction_network_T-vs-N.png)

### 2.4 Condition-Specific Correlation Network with Consistent Node Layout - run_conditional_network Function

**Function Overview:**
  Condition-specific correlation network analysis is used to compare molecular network differences under different experimental conditions (e.g., disease vs. health). This function ensures the same node layout is used across different conditions, facilitating intuitive comparison of network structural changes.

**Main Parameter Description:**
  
- count_table: Molecular expression table, containing a feature_ID column.
- samplelist: Sample grouping information table, containing group and sample columns.
- compare_group: Comparison group (format: "T:N").
- Diff_anno: Differential analysis annotation table.
- node_list: List of nodes to analyze, typically significantly differential molecules.
- nBoots: Number of Bootstrap iterations (default 50, set to 5 in the example for faster demonstration).
- stability_threshold: Network stability threshold.
- bootnet_R_threshold: Correlation coefficient threshold.

**Usage Example:**
  
```R
# Create example network data
data("differential_network_result")

# Or read data
differential_network_result=readRDS("differential_network_result.rds")

count_table = differential_network_result@PreData@count_table
samplelist = differential_network_result@PreData@samplelist
Diff_anno = differential_network_result@Diff_anno
node_list = differential_network_result@node_list

# Run conditional network analysis
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

# Visualization
# Experimental group
network_show(
  Network = Conditional_network,
  plot_type = 'case_overall_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean', palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_conditional_network_T.png)

```R
# Control group
network_show(
  Network = Conditional_network,
  plot_type = 'control_overall_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean', palette_style = "morandi"
)
```
![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_conditional_network_N.png)

### 2.5 Condition-Specific Multiplex Network with Consistent Node Layout - run_conditional_multiplexnetwork Function

**Function Overview:**
  Condition-specific multiplex network analysis integrates molecular interaction networks (e.g., protein-protein interactions) with condition-specific correlation networks to construct a multi-level, multi-type multiplex network. This function ensures the same node layout is used across different conditions (e.g., experimental and control groups), facilitating comparison of multiplex network structural changes.

**Main Parameter Description:**
  
- Interaction_network: Molecular interaction network object, generated by the run_interaction_network function.
- Conditional_network: Condition-specific network object containing network data for experimental and control groups, generated by the run_conditional_network function.
- compare_group: Comparison group (format: "T:N").

**Usage Example:**
  
```R
# Create example network data
data("multiplex_network_result")

# Or read data
multiplex_network_result=readRDS("multiplex_network_result.rds")

Interaction_network = multiplex_network_result@Interaction_network
Conditional_network = multiplex_network_result@Conditional_network

# Run conditional multiplex network analysis
Conditional_multiplexnetwork <- run_conditional_multiplexnetwork(
  Interaction_network = Interaction_network,
  Conditional_network = Conditional_network,
  compare_group = 'T:N'
)

# Visualization
# Experimental group
network_show(
  Network = Conditional_multiplexnetwork,
  plot_type = 'case_multi_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean', palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_conditional_multiplexnetwork_T.png)

```R
# Control group
network_show(
  Network = Conditional_multiplexnetwork,
  plot_type = 'control_multi_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Normalized mean', palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_conditional_multiplexnetwork_N.png)

### 2.6 Differential Network Analysis - run_diff_network Function

**Function Overview:**
  Differential network analysis is used to compare network differences between two conditions (e.g., experimental and control groups), identifying nodes and edges that change significantly under different conditions. This analysis can detect dynamic changes in network structure, including added/removed edges, edges with significantly changed correlation strength, and node differential expression information.

**Main Parameter Description:**
  
- Conditional_network: Condition-specific network object containing network data for experimental and control groups, generated by the run_conditional_network function.
- Conditional_multiplexnetwork: Multiplex network object; input this when performing differential network analysis using a multiplex network object, generated by the run_conditional_multiplexnetwork function.
- edge_FC_threshold: Edge fold change threshold for filtering edges with significantly changed correlation strength, default: 1.2.
- edge_p_threshold: Edge significance p-value threshold for statistical testing, default: 0.05.
- compare_group: Comparison group (format: "T:N"), e.g., "T:N"表示 tumor group vs. normal group.

**Usage Example:**
  
#### 2.6.1 Perform Differential Network Analysis using Conditional Network Object
  
```R
# Create example network data
data("differential_network_result")

# Or read data
differential_network_result=readRDS("differential_network_result.rds")

Conditional_network = differential_network_result@Conditional_network

# Run differential network analysis
Differential_network <- run_diff_network(
  Conditional_network = Conditional_network,
  edge_FC_threshold = 1.2,
  edge_p_threshold = 0.05,
  compare_group = 'T:N'
)

# Visualization
network_show(
  Network = Differential_network,
  plot_type = 'diff_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Log2FC', palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_diff_network_T-vs-N.png)

#### 2.6.2 Perform Differential Network Analysis using Multiplex Network Object

```R
# Create example network data
data("multiplex_network_result")

# Or read data
multiplex_network_result=readRDS("multiplex_network_result.rds")

Conditional_network = multiplex_network_result@Conditional_network
Conditional_multiplexnetwork = multiplex_network_result@Conditional_multiplexnetwork

Differential_multiplexnetwork <- run_diff_network(
  Conditional_network = Conditional_network,
  Conditional_multiplexnetwork = Conditional_multiplexnetwork,
  edge_FC_threshold = 1.2,
  edge_p_threshold = 0.05,
  compare_group = 'T:N'
)

# Visualization
network_show(
  Network = Differential_multiplexnetwork,
  plot_type = 'diff_network',
  show_node_legend = TRUE,
  show_edge_legend = TRUE,
  node_colortype = 'Log2FC', palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_diff_multinetwork_T-vs-N.png)

### 2.7 Network Node Enrichment Analysis

#### 2.7.1 Single Group Network Enrichment Analysis - run_enrichment Function

**Function Overview:**
  Network node enrichment analysis performs KEGG pathway enrichment analysis on nodes within a network, identifying biological pathways significantly enriched in the network. This function supports multi-omics data integration analysis, capable of simultaneously processing metabolomics data (KEGG IDs starting with C) and genomics/proteomics data (KEGG IDs starting with K), providing comprehensive pathway functional annotation.

**Main Parameter Description:**
  
- annotation_table_select: Data frame or list containing feature annotations, must include a KEGG.ID column.
- species: Species abbreviation (e.g., "hsa" for human, "mmu" for mouse).
- nCores: Number of cores for parallel computation.
- omics_name: Omics type (e.g., c("proteomics","metabolomics")); typically used only in multi-omics analysis; when setting this parameter, the annotation data must include an omics_name column.
- group_name: Analysis group name, default "data".
- annotation_table_select: Data frame or list containing feature annotations, must include a KEGG.ID column.
- enrichment_p_threshold: Enrichment analysis p-value threshold, default 0.05.
- database_path: Path to the enrichment database.

**Analysis Features:**
  
- Multi-omics Integration: Automatically identifies and integrates enrichment results from different omics types (metabolomics, genomics, etc.).
- Parallel Computing: Supports multi-core parallel processing to improve efficiency for large-scale data analysis.
- Automatic Database Management: Automatically downloads and manages KEGG database mapping files.
- Flexible Input: Supports two input formats: single network node list and subnetwork node list.

**Usage Example:**
  
##### 2.7.1.1 Overall Network Enrichment Analysis
  
```R
# Create example network data
data("stable_subnetwork_result")

# Or read data
stable_subnetwork_result=readRDS("stable_subnetwork_result.rds")

subnetworks_nodes <- stable_subnetwork_result@StableNetwork@bootnet_result_filter@bootnet_node

# Run enrichment analysis
Enrichment <- run_enrichment(
  annotation_table_select = subnetworks_nodes,
  species = 'hsa',
  group_name = 'T'
)

# Visualization
network_show(Network = Enrichment, plot_type = 'enrichment', palette_style = "morandi")
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_enrichment_T.png)

##### 2.7.1.2 Subnetwork Enrichment Analysis

```R
# Create example network data
data("stable_subnetwork_result")

# Or read data
stable_subnetwork_result=readRDS("stable_subnetwork_result.rds")

subnetworks_nodes <- lapply(
  stable_subnetwork_result@SubNetwork@subnetworks, 
  function(subnet) {
    subnet$nodes
  }
)

# Run enrichment analysis
Subnet_Enrichment <- run_enrichment(
  annotation_table_select = subnetworks_nodes,
  species = 'hsa',
  group_name = 'T'
)

# Visualization (example: subnet_1)
network_show(
  Network = Subnet_Enrichment,
  plot_type = 'subnetwork_enrichment',
  subnetwork_name = c("subnet_1"), palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_enrichment_subnet_1_T.png)

#### 2.7.2 Differential Overall Network Enrichment Analysis (Applicable to Multiplex Overall Network) - run_diff_enrichment Function

**Function Overview:**
  Differential overall network enrichment analysis is based on differential network analysis results, performing enrichment analysis on nodes corresponding to edges grouped by correlation state. This function can identify biological pathways significantly enriched under different correlation states (e.g., only in experimental group, only in control group, enhanced in experimental group, enhanced in control group, etc.), revealing the functional significance of network differences.

**Main Parameter Description:**
  
- Differential_network: Differential network object, generated by the run_diff_network function.
- species: Species abbreviation (e.g., "hsa" for human, "mmu" for mouse).
- nCores: Number of CPU cores for parallel computation.
- omics_name: Omics type (e.g., c("proteomics","metabolomics")); typically used only in multi-omics analysis; when setting this parameter, the annotation data must include an omics_name column.
- enrichment_p_threshold: Enrichment significance p-value threshold, default 0.05.
- compare_group: Comparison group name used in enrichment analysis, default "data".
- database_path: File path to the enrichment database directory.

**Usage Example:**
  
```R
data("differential_network_result")

# Or read data
differential_network_result=readRDS("differential_network_result.rds")

Differential_network = differential_network_result@Differential_network

Enrichment <- run_diff_enrichment(
  Differential_network = Differential_network,
  species = 'hsa',
  enrichment_p_threshold = 0.05,
  compare_group = 'T:N'
)

# Visualization
network_show(Network = Enrichment, plot_type = 'enrichment', palette_style = "morandi")
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_diff_enrichment_T-vs-N.png)

#### 2.7.3 Differential Subnetwork Enrichment Analysis (Applicable to Multiplex Subnetworks) - run_diffsubnet_enrichment Function

**Function Overview:**
  Differential subnetwork enrichment analysis performs functional enrichment analysis on individual subnetworks within a differential network. This function can identify biological pathways significantly enriched for different correlation state groups (e.g., only in experimental group, only in control group, etc.) within each subnetwork, providing in-depth biological interpretation for understanding functional differences in network modules.

**Main Parameter Description:**
  
- Differential_subnetwork: Differential subnetwork object, generated by the run_diff_network function.
- species: Species abbreviation (e.g., "hsa" for human, "mmu" for mouse).
- nCores: Number of CPU cores for parallel computation.
- omics_name: Omics type (e.g., c("proteomics","metabolomics")); typically used only in multi-omics analysis; when setting this parameter, the annotation data must include an omics_name column.
- enrichment_p_threshold: Enrichment significance p-value threshold, default 0.05.
- compare_group: Comparison group name used in enrichment analysis, default "data".
- database_path: File path to the enrichment database directory.

**Usage Example:**
  
```R
Differential_subnetwork = differential_network_result@Differential_subnetwork

DiffSubnet_Enrichment <- run_diffsubnet_enrichment(
  Differential_subnetwork = Differential_subnetwork,
  species = 'hsa',
  enrichment_p_threshold = 0.05,
  compare_group = 'T:N'
)

# Visualization (example: subnet_1)
network_show(
  Network = DiffSubnet_Enrichment,
  plot_type = 'differential_subnetwork_enrichment',
  subnetwork_name = c("subnet_1"), palette_style = "morandi"
)
```

![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/run_diff_enrichment_subnet_1_T-vs-N.png)

# III. Plotting Functions

## 3.1 Overall Pipeline Result Data and Image Saving Function

### pipline_save Function, usage reference Part I

**Function Overview:**
  The pipline_save function is a comprehensive result saving function for the network analysis pipeline, capable of automatically saving all images and data files generated during the network analysis process. This function supports three main network analysis types (Stable_SubNetwork, Stable_DifferentialNetwork, Stable_MultiplexNetwork) and automatically organizes the output file structure based on the analysis type.

**Main Parameter Description:**
  
- Network: Network analysis result object (e.g., Stable_SubNetwork, StableNetwork, etc.).
- stable_num: Number of bootstrap networks to display in the stability test.
- richfactor_threshold: Enrichment analysis rich factor threshold, default 0.
- R_threshold: Default 0.3.
- StabilityTest: Logical value, whether to save stability test results, default TRUE.
- OverallNetwork: Logical value, whether to save overall network analysis results, default TRUE.
- NetworkClustering: Logical value, whether to save network clustering results, default TRUE.
- SubNetWork: Logical value, whether to save subnetwork analysis results, default TRUE.
- outdir: Output directory, default "./".

**Saved Content:**
  Depending on the network analysis type, the function automatically saves the following content:
  
**Stable_SubNetwork Type Saved Content:**
- Stability Test Results: Bootstrap network images and node/edge data files.
- Overall Network Analysis: Stable correlation network images, enrichment analysis bubble charts, and related data files.
- Network Clustering Analysis: Clustered network images and annotation data files.
- Subnetwork Analysis: Images for each subnetwork, topological analysis results (centrality analysis), and - enrichment analysis results.

**Stable_DifferentialNetwork Type Saved Content:**
- Stability Test Results: Bootstrap network images and data for experimental and control groups.
- Differential Network Analysis: Differential network images, node/edge data files, and enrichment results.
- Differential Clustering Analysis: Differential clustering network images and data files.
- Differential Subnetwork Analysis: Images and enrichment analysis results for each differential subnetwork.

**Stable_MultiplexNetwork Type Saved Content:**
- Interaction Network Analysis: Protein-protein interaction network images and data files.
- Stability Test: Bootstrap network results for experimental and control groups.
- Correlation Networks: Overall correlation networks for experimental and control groups.
- Differential Multiplex Network Analysis: Multiplex differential network images, data, and enrichment results.
- Network Clustering Analysis: Clustered network images and annotation data files.
- Differential Multiplex Subnetwork Analysis: Images and enrichment analysis results for each differential - - multiplex subnetwork.

**Usage Example:**
  
```R
# Save stable subnetwork analysis results
data("stable_subnetwork_result")
pipline_save(stable_subnetwork_result)

# Save differential network analysis results
data("differential_network_result")
pipline_save(differential_network_result)

# Save multiplex network analysis results
data("multiplex_network_result")
pipline_save(multiplex_network_result)
```

## 3.3 One-Click Re-styling

**Function Overview:**
OmicsNetwork supports multiple publication-ready color palettes, allowing users to easily customize the visual style of their plots. By setting the `palette_style` parameter in `network_show()` or `pipline_save()`, you can switch between different aesthetic themes.

**Supported Styles:**
- `"morandi"` (default): Low-saturation, muted colors, ideal for modern publications.
- `"npg"`: Nature Publishing Group style.
- `"aaas"`: Science (AAAS) style.
- `"nejm"`: New England Journal of Medicine style.
- `"lancet"`: The Lancet style.
- `"jco"`: Journal of Clinical Oncology style.
- `"jama"`: JAMA style.
- `"brewer"`: Classic high-contrast RColorBrewer style.

**Usage Example:**

```R
# Morandi Style (Default)
network_show(Network = subnet_result, plot_type = "overall_cluster_network", palette_style = "morandi")

# Lancet Style
network_show(Network = subnet_result, plot_type = "overall_cluster_network", palette_style = "lancet")
```

**Visual Comparison:**

*Morandi Style:*
![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Palette_Showcase_Cluster_morandi.png)

*Lancet Style:*
![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Palette_Showcase_Cluster_lancet.png)

*NPG Style:*
![](https://github.com/mzlab-research/OmicsNetwork/blob/main/OmicsNetwork_Plot/Palette_Showcase_Cluster_npg.png)

