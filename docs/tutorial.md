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
Next, an illustration of how to implement script is how below.
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

Run the script with help flag to see usage for multimodal data
```bash
python Multimodal_Parameter_Search.py -h
```
**Expected Output**
```text
ArborMap Parameter Search(CITE-seq)

positional arguments:
  file                  csv file with ADT reduced cell embeddings
  file                  csv file with RNA reduced cell embeddings
  file                  weighted integrated multimodal matrix
  jobs                  The number of jobs or workers to run for the KNN and ensemble tree algorithms
  file_name_clusters    Name of output file (default: ./data/cluster_evaluation.csv)

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Path to the output folder (default: ./data) (default: ./data)
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
                        How many workers or n_jobs to use for RandomTree embedding and KNN algorithmn model (default: None)
  --resolution RESOLUTION
                        Resolution for louvain community detection (default: None)
  --SNN_prune SNN_PRUNE
                        How much to prune the SNN algorithm (default: None)
  --output_dir OUTPUT_DIR
                        Path to the output folder (default: ./data) (default: ./data)
  --save                Save Tree model sparse matrix if --save is used (default: False)
  --try_resolutions     Try multiple resolutions using one K value --save is used (default: False)y multiple resolutions using one K value --save is used (default: False)

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
