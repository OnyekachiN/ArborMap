# ArborMap
ArborMap is a bioinformatics pipeline for single-cell clustering, including scRNA-seq and CITE-seq. ArborMap takes as input dimensionality-reduced or batch-corrected cell embeddings and returns cell clusters alongside their UMAP representations. The core of the clustering engine is

**RandomTreesEmbedding → KNN/SNN Jaccard graph (Seurat-style) → Louvain community detection → UMAP projection.**

