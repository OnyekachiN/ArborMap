# ArborMap Workflow
Follow the steps below to deploy ArborMap on your single-cell data

## Prerequisites 
Ensure you have **Python 3.13+** installed on your system. You can check your version by running:
```bash
python --version
```
## Step 1: Clone Repository
git clone https://github.com/OnyekachiN/ArborMap.git

## Step 2: Set Up a Virtual Environment (Recommended)
```bash
# On macOS/Linux
python3 -m venv arbormap_env
source arbormap_env/bin/activate

# On Windows (Command Prompt)
python -m venv arbormap_env
arbormap_env\Scripts\activate
```
## Step 3: Install Dependencies
Install the required packages listed in the `requirements.txt` file using `pip`
```bash
pip install -r requirements.txt
```
## Step 4: Run ArborMap tool
### Workflow for unimodal single-cell data
You can run the parameter selection step to help you choose clustering parameters, this step is optional.
Run the script with help flag to see usage for scRNA-seq.
 ```bash
python ScRNA_seq_Parameter_Search.py -h
```
**Expected Output**
```text
ArborMap Parameter Search(ScRNA-seq)

positional arguments:
  file                  csv file that contains the pca or integrated reduced cell embeddings
  jobs                  The number of jobs or workers to run for the KNN and ensemble tree algorithms
  file_name_clusters    Name of output file (default: ./data/cluster_evaluation.csv)

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Path to the output folder (default: ./data) (default: ./data)
```
Next, an illustration of how to implement script is shown below.
```bash
python ScRNA_seq_Parameter_Search.py Allen_tutorial_subset.csv 12 '/RNA_Parameter_Search'  
```
**Expected Output**
```text
starting ensemble tree model
Type of A: <class 'scipy.sparse._csr.csr_matrix'>
starting KNN algorithm
Working on k 20
Building SNN graph
```
### Workflow for multimodal single-cell data
You may use our custom multimodal integration script to integrate assays available in your single-cell datasets. Run the script with help flag to see usage f

```bash
Python Make_weightedsum_matrices.py -h
```
**Expected Output**
```text
ArborMAP weighted integration module

positional arguments:
  file                  csv file with RNA reduced cell embeddings
  file                  csv file with ADT reduced cell embeddings
  weight_rna            weight assigned to multimodal assay
  weight_adt            weighte assigned to multimodal assay
  file_name_integrated  Name of output integrated file (default: ./data/multimodal_assay_integration.csv)

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Path to the output folder (default: ./data) (default: ./data)
```
Next, an illustration of how to implement script is shown below.
```bash
python Make_weightedsum_matrices.py \
Data/Multimodal/multimodal_tutorial_subset_sct.csv \
Data/Multimodal/multimodal_tutorial_subset_adt.csv \
0.1 \
0.9 \
--output_dir '/Data/Multimodal/' \
tutorial_weighted_matrix.csv  
```
Similarly for multimodal data, you can run the parameter selection step to help you choose clustering parameters, this step is optional. Use the help flag to see usage for multimodal data
```bash
python Multimodal_Parameter_Search.py -h
```
**Expected Output**
```text
ArborMap Parameter Search(CITE-seq)

positional arguments:
  file                  csv file with RNA reduced cell embeddings
  file                  csv file with ADT reduced cell embeddings
  file                  weighted integrated multimodal matrix
  jobs                  The number of jobs or workers to run for the KNN and ensemble tree algorithms
  file_name_clusters    Name of output file (default: ./data/cluster_evaluation.csv)

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Path to the output folder (default: ./data) (default: ./data)
```
Next, to run multimodal parameter script, see demonstration below.
```bash
python Codes/RF_clustering_evaluation_analysis_multimodal.py \
Data/Multimodal/multimodal_tutorial_subset_sct.csv \
Data/Multimodal/multimodal_tutorial_subset_adt.csv \
Data/Multimodal/tutorial_weighted_matrix.csv \
12 \
--output_dir '/Data/Multimodal/' \
Multimodal_parameter_search.csv
```
### Run complete ArborMap clustering pipeline 
Illustration of how to run complete ArborMap pipeline either after running step above or not. Start off by running script with help flag.
```bash
python ArborMap_v1.0.py -h  
```
**Expected Output**
```text
ArborMap clustering Tool

positional arguments:
  file                  csv file that contains the reduced cell embeddings
  file                  csv file that contains the colors to use for umap
  K                     Number of nearest neighbors
  file_name_clusters    Name of output file (default: ./data/louvain_KNN_SNN_clusters.csv)

options:
  -h, --help            show this help message and exit
  --n_workers N_WORKERS
                        How many workers or n_jobs to use for RandomTree embedding and KNN algorithmn model (default: 0)
  --resolution RESOLUTION
                        Resolution for louvain community detection (default: None)
  --SNN_prune SNN_PRUNE
                        How much to prune the SNN algorithm (default: None)
  --output_dir OUTPUT_DIR
                        Path to the output folder (default: ./data) (default: ./data)
  --save                Save Tree model sparse matrix if --save is used (default: False)
  --try_resolutions     Try multiple resolutions using one K value --save is used (default: False)
  --file_name_umap FILE_NAME_UMAP
                        Name of output UMAP embedding file (default: ./data/ArborMAP_UMAP.csv)
```
Next, an illustration of how to run complete ArborMap pipeline shown below
```bash
python ArborMap_v1.0.py \
Data/Allen_tutorial_subset.csv \
Data/35_color_set.csv \
--n_workers 12 \
100 \
--resolution 0.5 \
--output_dir 'Data/' \
tutorial_subset_clusters.csv \
--file_name_umap tutorial_subset_UmapEmbedding.csv \
--save  
```
**Expected Output**
```text
Embeddeding dimensions are (10000, 20)
starting ensemble tree model
starting KNN model
Starting SNN computation
sparse csr matrix saved to: Data/ArborMap_KNN_sparse_matrix_k100_0.5.npz
Identifying and reassigning singleton clusters
[ 0  1  2  3  4  5  6  7  8  9 10]
k=100, res=0.5: 11 clusters, singletons reassigned=0
True
Computation time for ensemble and knn algorithm is: 32.4946 seconds
File saved to: Data/tutorial_subset_clusters.csv
making umap
axist detail and returning embedding
 x min = -7.49815411567688, x max =16.999959421157836
 y min = -10.249425792694092, y max = 21.603423976898192
Computation time for UMAP is: 284.9648 seconds
File saved to: Data/tutorial_subset_UmapEmbedding.csv
```
## Step 5: (Optional) Add cluster result to Seurat object
The following steps show how you can add the cluster information and UMAP embedding into a Seurat object.
```R
library(Seurat)

# Read in object, ArborMap cluster file and ArborMap UMAP embedding file
object = readRDS('Allen_tutorial_downsample_object.Rds')
ArborMap_clusters = read.csv('Data/tutorial_subset_clusters.csv', header =T, row.names = 1)

ArborMap_UMAP = read.csv('Data/tutorial_subset_UmapEmbedding.csv', row.names = 1)
ArborMap_UMAP = as.matrix(ArborMap_UMAP)

Idents(object)

# Make show the order of cells in the  object matches the ArborMap cluster file
match( rownames(object@meta.data), rownames(ArborMap_clusters))
reorder_idx = match( rownames(object@meta.data), rownames(ArborMap_clusters))
ArborMap_clusters = ArborMap_clusters[reorder_idx,]
identical( rownames(object@meta.data), rownames(ArborMap_clusters))

#Add ArborMap cluster into objects metadata
object@meta.data$ArborMap_clusters = ArborMap_clusters$cluster

Idents(object) = object@meta.data$ArborMap_clusters
levels(object) = as.character(unique(sort(ArborMap_clusters$cluster)))
Idents(object)

# Add ArborMap UMAP embedding into objects reduction slot
umap.dr <- CreateDimReducObject(embeddings = ArborMap_UMAP, key = "ArborMapUMAP_", assay = "RNA")
object[["ArborMapumap"]] <- umap.dr
```
Next, the ArborMap results can now be visualized in Seurat using the DimPlot() visualization function, demonstrated below.
```R
library(ggplot2)

Idents(object) = object@meta.data$ArborMap_clusters
col_number <- length(table(Idents(object)))
my_colors <- read.csv("Data/35_color_set.csv")$color_codes[1:col_number]
scales::show_col(my_colors)

umap = DimPlot(object, 
               reduction = 'ArborMapumap', 
               label = FALSE, 
               label.size = 12, 
               cols = my_colors,
               group.by = 'ArborMap_clusters', 
               shuffle = TRUE,
               raster = FALSE
               
) + guides(color = guide_legend(override.aes = list(size = 5)))+
  ggtitle(NULL) +
  theme(strip.text = element_text(size = 8))
```
