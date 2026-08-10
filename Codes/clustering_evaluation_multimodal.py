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

def select_K_res(sparse_mat, embedding_sct, embedding_adt, jobs =1,SNN_prune=None):
    print('starting KNN algorithm')
    output_df = pd.DataFrame(columns=['K', 
                                      'resolution', 
                                      'Number_of_clusters',
                                      'CH_score_sct',
                                      'CH_score_adt',
                                      'CH_score_avg',
                                      'Silhouette_sct', 
                                      'Silhouette_adt'
                                      ])
    K_list = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110,120]
    
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
            df['cell_id'] = [embedding_sct.index[i] for i in df['node']]
            df = df.set_index('cell_id').reindex(embedding_sct.index)

            # Identify singleton clusters (clusters of size 1) 
            cluster_sizes = df['cluster'].value_counts()
            singleton_clusters = cluster_sizes[cluster_sizes == 1].index
            
            if len(singleton_clusters) > 0:
                print('reassigning singletons')
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
            unique_clusters = df['cluster'].unique()
            mapping = {old: new for new, old in enumerate(unique_clusters)}
            df['cluster'] = df['cluster'].map(mapping)
            
            # Compute cluster metrics
            CH_sct = calinski_harabasz_score(embedding_sct, df['cluster'])
            silhouette_sct = silhouette_score(embedding_sct, df['cluster'], metric='euclidean')

            CH_adt = calinski_harabasz_score(embedding_adt, df['cluster'])
            silhouette_adt = silhouette_score(embedding_adt, df['cluster'], metric='euclidean')
            n_clusters = len(df['cluster'].unique())

            #find the mean of sct and adt Calinski harabaz scores
            CH_avg = (CH_sct + CH_adt)/2
            
            new_row = pd.DataFrame([{
                'K': k,
                'resolution': res,
                'Number_of_clusters': n_clusters,
                'CH_score_sct': CH_sct,
                'CH_score_adt': CH_adt,
                'CH_score_avg': CH_avg,
                'Silhouette_sct': silhouette_sct,
                'Silhouette_adt': silhouette_adt
            }])
            output_df = pd.concat([output_df, new_row], ignore_index=True)
            print(f"k={k}, res={res}: {n_clusters} clusters, CH_sct={CH_sct}, CH_adt={CH_adt} Silhouette_sct={silhouette_sct},Silhouette_adt={silhouette_adt}, reassigned singleton cells={len(singleton_clusters)}")
    
    return output_df

def rank_similairy_score(df):
    ranks_sct = df['CH_score_sct'].rank(ascending=True)
    ranks_adt = df['CH_score_adt'].rank(ascending=True)

    similarity_01 = (ranks_sct +ranks_adt)/(max(ranks_sct)+max(ranks_adt))  *(1 - abs(ranks_sct - ranks_adt) / (ranks_sct + ranks_adt))
    

    df['similarity_score'] = similarity_01
    df['rank_sct'] = ranks_sct
    df['rank_adt'] = ranks_adt
    return df

def main():
    parser = argparse.ArgumentParser(description='Tree based cluster evaluation for multimodal data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('emb_sct', metavar='file', help='csv file with ADT reduced cell embeddings')
    parser.add_argument('emb_adt', metavar = 'file', help = 'csv file with RNA reduced cell embeddings')
    parser.add_argument('integrated_mat', metavar = 'file', help = 'weighted integrated multimodal matrix')
    parser.add_argument('jobs', type=int, help='The number of jobs or workers to run for the KNN and ensemble tree algorithms')
    parser.add_argument("--output_dir", type=str, default="./data",help="Path to the output folder (default: ./data)")
    parser.add_argument("file_name_clusters", type=str, default="./data/cluster_evaluation.csv",help="Name of output file (default: ./data/cluster_evaluation.csv)")

    args = parser.parse_args()

    emb_adt = pd.read_csv(args.emb_adt, sep = ',', index_col=0)
    emb_sct = pd.read_csv(args.emb_sct, sep = ',', index_col=0)

    integrated_mat = pd.read_csv(args.integrated_mat, sep = ',', index_col=0)

    print('starting ensemble tree model')
    start_computation = perf_counter()
    rf_classifier = RandomTreesEmbedding(n_estimators= 1000, random_state = 1, n_jobs=args.jobs)
    rf_classifier.fit(integrated_mat)
    Ajc_mtx = rf_classifier.transform(integrated_mat)
    print(f"Type of A: {type(Ajc_mtx)}") 

    # 1. Step 1, to assist in deciding on what K and resolution is best, check various k's and resolutions 
    param_selection = select_K_res(Ajc_mtx, emb_sct, emb_adt, jobs =args.jobs) # type: ignore
    end_computation = perf_counter()
    print(f"Computation time for ensemble and multimodal cluster evaluation algorithm is: {end_computation - start_computation:.4f} seconds")

    #calculate ranked correlation between the rna and adt calinski harabaz scores
    # Rank values so that largest value gets highest rank
    print('Ranking and calculating score for sct and adt Calinski Harabsz scores')
    param_selection = rank_similairy_score(param_selection)

    # 2. Step 2, write outputs to file 
    save_dataframe(args.output_dir, param_selection, args.file_name_clusters)

if __name__ == '__main__':
    sys.exit(main())
