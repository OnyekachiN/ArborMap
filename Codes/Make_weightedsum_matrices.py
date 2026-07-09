import pandas as pd
import argparse
import sys
import os

def save_dataframe(output_folder, df, file_name):
    # Ensure the directory exists
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, str(file_name))
    df.to_csv(output_path, index=True)
    print(f"File saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Custom weighted integration module',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('emb_sct', metavar='file', help='csv file with ADT reduced cell embeddings')
    parser.add_argument('emb_adt', metavar = 'file', help = 'csv file with RNA reduced cell embeddings')
    parser.add_argument('weight_rna', nargs='+',type = float, help = 'weighted integrated multimodal matrix')
    parser.add_argument('weight_adt', nargs ='+',type = float, help = 'weighted integrated multimodal matrix')
    parser.add_argument("--output_dir", type=str, default="./data",help="Path to the output folder (default: ./data)")
    parser.add_argument("file_name_clusters", type=str, default="./data/cluster_evaluation.csv",help="Name of output file (default: ./data/cluster_evaluation.csv)")

    args = parser.parse_args()

    #read in embedding for sct and adt
    Embedding_data_rna = pd.read_csv(args.emb_sct, sep = ',', index_col=0)
    Embedding_data_adt = pd.read_csv(args.emb_adt, sep = ',', index_col=0)

    #Define column names
    col_names_list = []
    for i in range(1, len(Embedding_data_rna.columns)+1):
        col_names = 'Embedding' + str(i)
        col_names_list.append(col_names)

    #Change column names for sct and adt dataframes
    Embedding_data_rna.columns = col_names_list
    Embedding_data_adt.columns = col_names_list

    ##carry out weighted sum of matrices 
    weights_rna = args.weight_rna
    weights_adt = args.weight_adt
    

    for i in range(len(weights_rna)):
        weighted_sum_distances = Embedding_data_rna *weights_rna[i] + Embedding_data_adt*weights_adt[i]
        save_dataframe(args.output_dir, weighted_sum_distances, args.file_name_clusters)

if __name__ == '__main__':
    sys.exit(main())