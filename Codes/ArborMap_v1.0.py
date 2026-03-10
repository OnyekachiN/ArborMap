import pandas as pd
from sklearn.ensemble import RandomTreesEmbedding
import numpy as np
from time import perf_counter
import numpy as np
import networkx as nx
from community import community_louvain
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import umap.umap_ as umap
import seaborn as sns
from matplotlib.colors import ListedColormap
import sys
import argparse
import os

#define functions  
def create_RandomTreemodel(data,n_trees =1000, random_seed = 1, n_workers = 10):
    rf_classifier = RandomTreesEmbedding(n_estimators= n_trees, random_state = random_seed, n_jobs=n_workers)
    rf_classifier.fit(data)
    Ajc_mtx = rf_classifier.transform(data)
    return Ajc_mtx

def save_dataframe(output_folder, df, file_name):
    # Ensure the directory exists
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, str(file_name))
    df.to_csv(output_path, index=True)
    print(f"File saved to: {output_path}")

def KNN_SNN(k, res, sparse_mat, emb_indx,jobs =1, SNN_prune=None):
    print('starting KNN model')
    snn_prune_val = 1/15 if SNN_prune is None else SNN_prune
    # Compute KNN graph (on precomputed distances) 
    nn = NearestNeighbors(n_neighbors=k, metric='cosine', algorithm='brute', n_jobs=jobs)
    nn.fit(sparse_mat)
    _, indices = nn.kneighbors()

    n = sparse_mat.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(n))

    # Jaccard SNN computation (Seurat-style) 
    print('Starting SNN computation')
    for i in range(n):
        neigh_i = indices[i]
        set_i = set(neigh_i)
        for j in neigh_i:
            if i == j:
                continue
            neigh_j = indices[j]
            shared = set_i.intersection(neigh_j)
            if not shared:
                continue
            s = len(shared)
            weight = s / (2 * k - s)
            if weight >= snn_prune_val:
                G.add_edge(i, j, weight=weight)

    # Louvain clustering 
    partition = community_louvain.best_partition(G, weight='weight', resolution=res, random_state=42)
    df = pd.DataFrame.from_dict(partition, orient='index', columns=['cluster'])
    df.index.name = 'node'
    df.reset_index(inplace=True)
    print(df)

    #Map node indices to cell IDs 
    df['cell_id'] = [emb_indx.index[i] for i in df['node']]
    df = df.set_index('cell_id').reindex(emb_indx.index)

    print('Identifying and reassigning singleton clusters')
    # Identify singleton clusters (clusters of size 1) 
    cluster_sizes = df['cluster'].value_counts()
    singleton_clusters = cluster_sizes[cluster_sizes == 1].index

    if len(singleton_clusters) > 0:
        clusters = df.loc[~df['cluster'].isin(singleton_clusters), 'cluster'].unique()
        #A = nx.to_pandas_adjacency(G, weight='weight')
        new_assignments = {}

        # Convert network to sparse array matrix (much faster indexing and memory demand)
        #A_np = A.values
        nodes = list(G.nodes())
        A = nx.to_scipy_sparse_array(G, nodelist = nodes,weight='weight', format='csr')
        print(type(A))
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
    df['cluster']   = df['cluster'].astype(int)
    
    unique_clusters = df['cluster'].unique()

    # order the numbers of clusters from 0 to n by increment of 1
    mapping = {old: new for new, old in enumerate(unique_clusters)}
    df['cluster'] = df['cluster'].map(mapping)

    n_clusters      = df['cluster'].nunique()            
    df['res']       = res
    df['KNN']       = k
    print(df['cluster'].unique())

    print(f"k={k}, res={res}: {n_clusters} clusters, singletons reassigned={len(singleton_clusters)}")
    print(df.index.equals(emb_indx.index))
    return df

def make_UMAP(sparse_mat, data_df, vmax, cmap_col, k, folder, dims, spread_n = 1, n_neighbors = 30, jobs=1, res = None):
    print('making umap')
    if spread_n != 1:
        spread_n = spread_n
    if n_neighbors !=30:
        n_neighbors = n_neighbors
    if res is not None:
        res = res
    #plot UMAP
    reducer = umap.UMAP(n_neighbors=n_neighbors,        
        min_dist=0.3,          
        metric='cosine',       
        random_state=42,
        spread = spread_n,
        repulsion_strength=1,
        n_jobs= jobs)

    embedding = reducer.fit_transform(sparse_mat)
    embedding.shape # type: ignore
    
    fig, ax = plt.subplots(figsize=(10, 8)) 
    sc = ax.scatter(
        embedding[:, 0], # type: ignore
        embedding[:, 1], # type: ignore
        c=data_df.cluster, 
        cmap=cmap_col, 
        s=1,
        vmin=0, vmax=vmax-1
    )

    cbar = fig.colorbar(sc, ax=ax, boundaries=np.arange(0, vmax +1) - 0.5)  
    cbar.set_ticks(np.arange(0, vmax))
    #ax.set_title(f'UMAP projection of the {proj_name} dataset using {dims} dimensions and {weights} weights at resolution of {res}', fontsize=10)
    fig.savefig(folder + f'UMAP_k{k}_{dims}_dims_{res}_resolution.png', dpi=300, bbox_inches='tight')
    #plt.show()
    print('axist detail and returning embedding')
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    print(f' x min = {x_min}, x max ={x_max}', sep = '/')
    print(f' y min = {y_min}, y max = {y_max}', sep = '/')
    embedding = pd.DataFrame(embedding)
    embedding.columns = ['UMAP_1', 'UMAP_2']
    embedding.index = data_df.index
    return embedding

def main():
    parser = argparse.ArgumentParser(description='Ensemble tree sparse matrix and graph network construction',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('reduced_obj', metavar='file', help='csv file that contains the reduced cell embeddings')
    parser.add_argument('color_file', metavar='file', help='csv file that contains the colors to use for umap')
    #parser.add_argument("--n_trees", type=int,help="The number of trees or n_estimators to use for the RandomTree Embedding model")
    #parser.add_argument("--Randomseed", type=int,help="Random seed for RandomTree embedding model")
    #parser.add_argument("--n_workers", type=int,help="How many workers or n_jobs to use for RandomTree embedding model")
    parser.add_argument('K',type=int, help = 'Number of nearest neighbors')
    parser.add_argument('resolution',type =float, help = 'Resolution for louvain community detection')
    parser.add_argument('--SNN_prune',type =float, help = 'How much to prune the SNN algorithm')
    parser.add_argument("--output_dir", type=str, default="./data",help="Path to the output folder (default: ./data)")
    parser.add_argument("file_name_clusters", type=str, default="./data/louvain_KNN_SNN_clusters.csv",help="Name of output file (default: ./data/louvain_KNN_SNN_clusters.csv)")
    parser.add_argument("file_name_umap", type=str, default="./data/ArborMAP_UMAP.png",help="Name of output UMAP figure (default: ./data/ArborMAP_UMAP.png)")


    args = parser.parse_args()

    ## read in embedding for only rna(RNAseq) or rna and adt (CITEseq)
    Embedding_data = pd.read_csv(args.reduced_obj, sep = ',', index_col=0)
    color = pd.read_csv(args.color_file, index_col= 0)

    print(Embedding_data.shape)
    #Build ensmble tree model using RandomTreeEmbedding
    # Time the fitting process
    print('starting ensemble tree model')
    start_computation = perf_counter()
    rf_classifier = RandomTreesEmbedding(n_estimators= 1000, random_state = 1, n_jobs=12)
    rf_classifier.fit(Embedding_data)
    Ajc_mtx = rf_classifier.transform(Embedding_data)
    print(f"Type of A: {type(Ajc_mtx)}") 

    dims = Embedding_data.shape[1]
    k = args.K
    res = args.resolution

    df = KNN_SNN(k, res, Ajc_mtx, Embedding_data, jobs = 10)
    end_computation = perf_counter()
    print(f"Computation time for ensemble and knn algorithm is: {end_computation - start_computation:.4f} seconds")

    ## Save df to local folder
    save_dataframe(args.output_dir, df, args.file_name_clusters)

    number = len(df['cluster'].unique())
    color = color['color_codes'].head(number)
    my_cmap = ListedColormap(sns.color_palette(color).as_hex()) # type: ignore
    my_cmap.colors

    start_computation = perf_counter()
    embedding = make_UMAP(Ajc_mtx, df, number, my_cmap, k, args.output_dir, dims, n_neighbors=30, res = res, jobs =12)
    end_computation = perf_counter()
    print(f"Computation time for UMAP is: {end_computation - start_computation:.4f} seconds")

    save_dataframe(args.output_dir, embedding, args.file_name_umap)
if __name__ == '__main__':
    sys.exit(main())