import numpy as np
from sklearn.ensemble import RandomTreesEmbedding
import networkx as nx
from community import community_louvain
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score
from time import perf_counter
import sys
import argparse
import os

# define all the functions functions needed for clustering and UMAP visualizations
def save_dataframe(output_folder, df, file_name):
    # Ensure the directory exists
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, str(file_name))
    df.to_csv(output_path, index=False)
    print(f"File saved to: {output_path}")

def select_K_res(sparse_mat, embedding,jobs =1, SNN_prune=None):
    print('starting KNN algorithm')
    output_df = pd.DataFrame(columns=['K', 'resolution', 'Number_of_clusters','CH_score','Silhouette'])
    K_list = [230, 240, 250, 260, 270, 280, 290, 300, 310, 320]
    
    for k in K_list:
        print(f'Working on k {k}')
        snn_prune_val = 1/15 if SNN_prune is None else SNN_prune
        
        # Fit nearest neighbors on distance matrix
        nn = NearestNeighbors(n_neighbors=k, metric='cosine', algorithm='brute', n_jobs=jobs)
        nn.fit(sparse_mat)
        _, indices = nn.kneighbors()
        
        # Build SNN graph
        n = sparse_mat.shape[0]
        G = nx.Graph()
        G.add_nodes_from(range(n))
        
        print('Building SNN graph')
        for i in range(n):
            neigh_i = indices[i]
            set_i = set(neigh_i)
            for j in neigh_i:
                if i == j:
                    continue
                neigh_j = set(indices[j])
                shared = set_i.intersection(neigh_j)
                if not shared:
                    continue
                s = len(shared)
                weight = s / (2 * k - s)  # Seurat-style SNN similarity
                if weight >= snn_prune_val:
                    G.add_edge(i, j, weight=weight)
        
        # Run Louvain clustering for multiple resolutions
        resolutions = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        for res in resolutions:
            partition = community_louvain.best_partition(G, weight='weight', resolution=res, random_state=42)
            df = pd.DataFrame.from_dict(partition, orient='index', columns=['cluster'])
            df.index.name = 'node'
            df.reset_index(inplace=True)
            
            # Map nodes to cell IDs
            df['cell_id'] = [embedding.index[i] for i in df['node']]
            df = df.set_index('cell_id').reindex(embedding.index)

            # Identify singleton clusters (clusters of size 1) 
            cluster_sizes = df['cluster'].value_counts()
            singleton_clusters = cluster_sizes[cluster_sizes == 1].index
            
            if len(singleton_clusters) > 0:
                clusters = df.loc[~df['cluster'].isin(singleton_clusters), 'cluster'].unique()
                #A = nx.to_pandas_adjacency(G, weight='weight')
                new_assignments = {}

                # Convert network to sparse array matrix (much faster indexing and memory demand)
                nodes = list(G.nodes())
                A = nx.to_scipy_sparse_array(G, nodelist = nodes,weight='weight', format='csr')
                print(type(A))

                #A_np = A.values
                #nodes = np.array(A.index)  # numeric node labels
                node_to_pos = {node: i for i, node in enumerate(nodes)}  # map node ID → row/col index

                # Precompute cluster - node indices mapping
                cluster_to_nodes = {
                    clust: df.loc[df['cluster'] == clust, 'node'].map(node_to_pos).values # type: ignore
                    for clust in clusters
                }

                # For each singleton cell
                for sing in singleton_clusters:
                    sing_cell_index = df.index[df['cluster'] == sing][0]
                    sing_cell_node  = df.loc[sing_cell_index, 'node']
                    sing_pos = node_to_pos[sing_cell_node]

                    # Compute mean connectivity in a vectorized way
                    connectivity = {
                        clust: A[sing_pos, clust_indices].mean() if len(clust_indices) > 0 else 0
                        for clust, clust_indices in cluster_to_nodes.items()
                    }
                    best_cluster = max(connectivity, key=connectivity.get) # type: ignore
                    new_assignments[sing_cell_index] = best_cluster
                
                # Reassign singleton cells
                for cell, new_cluster in new_assignments.items():
                    df.loc[cell, 'cluster'] = new_cluster

            # Convert to int for downstream metrics
            df['cluster'] = df['cluster'].astype(int)
            
            # Compute cluster metrics
            CH = calinski_harabasz_score(embedding, df['cluster'])
            silhouette = silhouette_score(embedding, df['cluster'], metric='euclidean')
            n_clusters = len(df['cluster'].unique())
            
            new_row = pd.DataFrame([{
                'K': k,
                'resolution': res,
                'Number_of_clusters': n_clusters,
                'CH_score': CH,
                'Silhouette': silhouette
            }])
            output_df = pd.concat([output_df, new_row], ignore_index=True)
            print(f"k={k}, res={res}: {n_clusters} clusters, CH={CH}, Silhouette={silhouette}, singleton clusters reassigned={len(singleton_clusters)}")
    
    return output_df

def main():
    parser = argparse.ArgumentParser(description='COMMOT pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('reduced_obj', metavar='file', help='csv file that contains the pca reduced cell embeddings')
    parser.add_argument('jobs', type=int, help='The number of jobs or workers to run for the KNN and ensemble tree algorithms')
    #parser.add_argument('min_sample_split',type =int, help = 'The minimum number of samples required to split an internal node', default=2)
    parser.add_argument("--output_dir", type=str, default="./data",help="Path to the output folder (default: ./data)")
    parser.add_argument("file_name_clusters", type=str, default="./data/cluster_evaluation.csv",help="Name of output file (default: ./data/cluster_evaluation.csv)")

    args = parser.parse_args()

    pca_emb = pd.read_csv(args.reduced_obj, sep = ',', index_col=0)
    print(pca_emb.shape)
    
    print('starting ensemble tree model')
    start_computation = perf_counter()
    rf_classifier = RandomTreesEmbedding(n_estimators= 1000, random_state = 1, n_jobs=args.jobs)
    rf_classifier.fit(pca_emb)
    Ajc_mtx = rf_classifier.transform(pca_emb)
    print(f"Type of A: {type(Ajc_mtx)}") 

    # 1. Step 1, to assist in deciding on what K and resolution is best, check various k's and resolutions 
    param_selection_pca = select_K_res(Ajc_mtx,pca_emb, jobs =args.jobs)
    end_computation = perf_counter()
    print(f"Computation time for ensemble and cluster evaluation algorithm is: {end_computation - start_computation:.4f} seconds")

    # 2. Step 2, write outputs to file 
    save_dataframe(args.output_dir, param_selection_pca, args.file_name_clusters)

if __name__ == '__main__':
    sys.exit(main())