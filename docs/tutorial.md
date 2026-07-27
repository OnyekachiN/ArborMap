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
Next, an illustration of how to implement script is hown below.
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
