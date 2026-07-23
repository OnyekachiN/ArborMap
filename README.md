# ArborMap
ArborMap is a bioinformatics python-based framework for unsupervised clustering of single-cell data, including scRNA-seq and CITE-seq. ArborMap takes as input dimensionality-reduced or batch-corrected cell embeddings and returns cell clusters alongside their UMAP representations. ArborMap is designed as a robust alternative to traditional PCA + graph clustering pipelines, particulary for capturing non-linear relationships in high-dimensional data. The core of the clustering engine is

**RandomTreesEmbedding → KNN/SNN Jaccard graph (Seurat-style) → Louvain community detection → UMAP projection.**

## Method Overview

### 1. Ensemble Tree Embedding
ArborMap uses an ensemble of randomized decision trees to transform cell embeddings into a sparse representation.

- Each tree partitions the data differently
- Cells that fall into the same leaf are considered similar
- This captures **non-linear structure** better than PCA alone

### 2. KNN Graph Construction
A K-nearest neighbor graph is built from the tree-derived embedding.

- Distance metric applied on transformed space

### 3. SNN Graph
The KNN graph is converted into a Shared Nearest Neighbor (SNN) graph:

- Edges weighted by shared neighbors
- Improves robustness to noise

### 4. Clustering
- Louvain clustering applied on SNN graph
- Resolution parameter controls granularity

### 5. UMAP Visualization
UMAP is used for visualization of clusters.
