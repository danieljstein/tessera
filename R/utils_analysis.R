#' Run UMAP and save fgraph and embeddings in Seurat object
#' 
#' @param obj A Seurat object.
#' @param reduction Name of dimensional reduction in `obj` to use as input to UMAP.
#' @param dims Dimensions of `reduction` to use as input to UMAP. If NULL, use all dimensions.
#' @param fgraph_only If TRUE, only compute and store the fuzzy simplicial graph (fgraph) and skip UMAP embedding.
#' @param graph.name Name of graph to store the fgraph in `obj`. Defaults to `<assay>_fgraph`.
#' @param reduction.name Name of dimensional reduction to store the UMAP embeddings in `obj`. Defaults to "umap".
#' @param assay Assay to set as default assay for the new dimensional reduction. If NULL, use the default assay of `obj`.
#' @param key Key prefix to use for the new dimensional reduction. Defaults to "UMAP_".
#' @param n_neighbors Number of nearest neighbors to use in UMAP. See `?uwot::umap` for details.
#' @param n_components Number of UMAP dimensions to compute. Ignored if `fgraph_only` is TRUE.
#' @param metric Distance metric to use in UMAP. See `?uwot::umap` for details.
#' @param spread UMAP spread parameter. See `?uwot::umap` for details.
#' @param min_dist UMAP minimum distance parameter. See `?uwot::umap` for details.
#' @param n_threads Number of threads to use in UMAP. See `?uwot::umap` for details.
#' @param fast_sgd Whether to use the fast stochastic gradient descent optimization in UMAP. See `?uwot::umap` for details.
#' @param verbose Whether to print progress messages.
#' @param ... Additional parameters to pass to `uwot::umap`.
#' 
#' @returns The input Seurat object with the fgraph and (optionally) UMAP embeddings added.
#'
#' @export
RunUMAPCustom = function(
    obj, reduction = "pca", dims = NULL, fgraph_only = FALSE,
    graph.name = NULL, reduction.name = "umap", assay = NULL, key = "UMAP_",
    n_neighbors = 30, n_components = 2, metric = "cosine", spread = 1, min_dist = 0.3,
    n_threads = NULL, fast_sgd = TRUE, verbose = TRUE, ...
) {
    if (is.null(assay)) {
        assay = Seurat::DefaultAssay(obj)
    }
    if (is.null(graph.name)) {
        graph.name = paste0(assay, '_fgraph')
    }

    embeddings = Seurat::Embeddings(obj, reduction)
    if (!is.null(dims)) {
        embeddings = embeddings[,dims]
    }

    if (fgraph_only) {
        out = list()
        out$fgraph = uwot::similarity_graph(
            X = embeddings,
            n_neighbors = n_neighbors,
            metric = metric,
            method = "umap",
            verbose = verbose
        )
    } else {
        out = uwot::umap(
            X=embeddings, ret_extra = "fgraph",
            n_threads = n_threads, fast_sgd = fast_sgd, verbose = verbose,
            n_neighbors = n_neighbors, n_components = n_components,
            metric = metric, spread = spread, min_dist = min_dist, ...
        )
        colnames(out$embedding) = paste0("UMAP_", seq_len(n_components))
        obj[[reduction.name]] = Seurat::CreateDimReducObject(embeddings = out$embedding,
                                                     key = key, assay = assay)
    }
    rownames(out$fgraph) = colnames(obj)
    colnames(out$fgraph) = colnames(obj)
    obj[[graph.name]] = Seurat::as.Graph(out$fgraph)

    return(obj)
}

#' Merge clusters in a Seurat object
#' 
#' @param obj A Seurat object.
#' @param to_merge A vector of cluster labels to merge.
#' @param new_label The new label for the merged cluster.
#' @param clusters.name Name of the metadata column in `obj` containing the cluster labels
#'  to be merged. Defaults to 'seurat_clusters'.
#' @param new.clusters.name Name of the metadata column to store the new cluster labels.
#'  If NULL, defaults to `<clusters.name>_m<new_label>`.
#' 
#' @returns The input Seurat object with the merged cluster labels added to metadata.
#' 
#' @export
MergeClusters = function(obj, to_merge, new_label, clusters.name = 'seurat_clusters', new.clusters.name = NULL) {
    stopifnot(is.factor(obj@meta.data[[clusters.name]]))
    
    if (is.null(new.clusters.name)) {
        new.clusters.name = paste0(clusters.name, '_m', new_label)
    }

    old_levels = levels(obj@meta.data[[clusters.name]])
    new_levels = c(old_levels[!(old_levels %in% to_merge)], new_label)

    obj@meta.data[[new.clusters.name]] = as.character(obj@meta.data[[clusters.name]])
    obj@meta.data[[new.clusters.name]][obj@meta.data[[new.clusters.name]] %in% to_merge] = new_label

    obj@meta.data[[new.clusters.name]] = factor(
        obj@meta.data[[new.clusters.name]],
        levels = new_levels
    )

    return(obj)
}

#' Add sub-cluster labels to a Seurat object
#' 
#' @param obj A Seurat object.
#' @param obj_sub A Seurat object containing the sub-clustered cells.
#' @param cluster The parent cluster label in `obj` that was sub-clustered.
#' @param clusters.name Name of the metadata column in `obj` containing the parent cluster labels.
#' @param sub.clusters.name Name of the metadata column to store the sub-cluster labels.
#'  If NULL, defaults to `<clusters.name>_s<cluster>`.
#'
#' @returns The input Seurat object with sub-cluster labels added to metadata.
#' 
#' @export
AddSubClusterLabels = function(
    obj, obj_sub, cluster, clusters.name = "seurat_clusters", sub.clusters.name = NULL
) {
    stopifnot(is.factor(obj@meta.data[[clusters.name]]))
    stopifnot(sum(obj@meta.data[[clusters.name]] == cluster) == ncol(obj_sub))
    stopifnot(all(rownames(obj@meta.data %>%
                           select(!!sym(clusters.name)) %>%
                           filter(!!sym(clusters.name) == cluster)) ==
                  rownames(obj_sub@meta.data)))

    if (is.null(sub.clusters.name)) {
        sub.clusters.name = paste0(clusters.name, '_s', cluster)
    }
    
    clusters.old = obj@meta.data[[clusters.name]]
    obj@meta.data[[sub.clusters.name]] = as.character(clusters.old)
    obj@meta.data[rownames(obj_sub@meta.data), sub.clusters.name] = paste0(cluster, '_', obj_sub$seurat_clusters)
    idx = which(levels(clusters.old) == cluster)
    
    new_levels = c(
        levels(clusters.old)[seq_len(idx-1)],
        paste0(cluster, '_', levels(obj_sub$seurat_clusters)),
        levels(clusters.old)[idx + seq_len(nlevels(clusters.old)-idx)]
    )
    obj@meta.data[[sub.clusters.name]] = factor(obj@meta.data[[sub.clusters.name]], levels = new_levels)

    return(obj)
}

#' Create a Seurat object for sub-clustering a specific cluster
#' 
#' @param obj A Seurat object.
#' @param cluster The cluster label to sub-cluster.
#' @param clusters.name Name of the metadata column in `obj` containing the cluster labels.
#' @param npcs Number of principal components to use for sub-clustering.
#' @param fast_sgd Whether to use fast SGD in UMAP.
#' @param n_neighbors Number of neighbors to use in UMAP.
#' @param scale.factor Scale factor for normalization. If NULL, use median nCount_RNA.
#' @param use.existing.embeddings Name of existing dimensional reduction in `obj`
#'   to use for sub-clustering. If NULL, compute PCA on the subsetted data.
#' @param meta.vars.include Metadata variables to include in the sub-clustered object.
#' @param harmony.group.by.vars Metadata variables to use for Harmony integration.
#' @param early_stop Whether to use early stopping in Harmony.
#' 
#' @returns A Seurat object for the sub-clustered cells.
#'
#' @export
MakeSubClusterObj = function(
    obj, cluster,
    clusters.name = "seurat_clusters", npcs=30,
    fast_sgd = FALSE,
    n_neighbors = 15,
    scale.factor = NULL,
    use.existing.embeddings = NULL,
    meta.vars.include = NULL,
    harmony.group.by.vars = NULL,
    early_stop = TRUE
) {
    subset_idx = which(obj@meta.data[[clusters.name]] == cluster)
    meta.vars.include = unique(c(meta.vars.include, harmony.group.by.vars))
    # counts = obj[["RNA"]]$counts[,subset_idx]
    # colnames(counts) = NULL
    obj_sub = CreateSeuratObject(
        counts = obj[["RNA"]]$counts[,subset_idx],
        # counts = counts,
        meta.data = obj@meta.data[subset_idx, meta.vars.include]
    )

    if (is.null(use.existing.embeddings)) {
        if (is.null(scale.factor)) {
            scale.factor = median(obj_sub$nCount_RNA)
        }
        obj_sub <- NormalizeData(object = obj_sub, scale.factor = scale.factor)
        obj_sub <- FindVariableFeatures(object = obj_sub)
        obj_sub <- ScaleData(object = obj_sub)  # might get an error in this line if the cells don't have good names
        obj_sub <- RunPCA(object = obj_sub, npcs = npcs)
        reduction = "pca"
    } else {
        obj_sub[[use.existing.embeddings]] <- CreateDimReducObject(
            embeddings = Embeddings(obj, use.existing.embeddings)[subset_idx,],
            loadings = Loadings(obj, use.existing.embeddings),
            key = use.existing.embeddings
        )
        reduction = use.existing.embeddings
    }

    if (!is.null(harmony.group.by.vars)) {
        obj_sub = harmony::RunHarmony(
            obj_sub, harmony.group.by.vars,
            early_stop = early_stop,
            reduction.use = reduction,
            reduction.save = paste0(reduction, '_harmony'))
        reduction = paste0(reduction, '_harmony')
    }

    obj_sub = RunUMAPCustom(obj_sub, reduction = reduction,
                            n_neighbors = n_neighbors, fast_sgd = fast_sgd)
    return(obj_sub)
}

#' Find sub-clusters within a specific cluster of a Seurat object
#' 
#' @param obj A Seurat object.
#' @param cluster The cluster label to sub-cluster.
#' @param clusters.name Name of the metadata column in `obj` containing the cluster labels.
#' @param sub.clusters.name Name of the metadata column to store the sub-cluster labels.
#'  If NULL, defaults to `<clusters.name>_s<cluster>`.
#' @param resolution Resolution parameter for clustering.
#' @param algorithm Clustering algorithm to use. See `?Seurat::FindClusters` for details.
#' @param npcs Number of principal components to use for sub-clustering.
#' @param method Method to use for clustering. See `?Seurat::FindClusters` for details.
#' @param n_neighbors Number of neighbors to use in UMAP.
#' @param fast_sgd Whether to use fast SGD in UMAP.
#' @param scale.factor Scale factor for normalization. If NULL, use median nCount_RNA.
#' @param use.existing.embeddings Name of existing dimensional reduction in `obj`
#'   to use for sub-clustering. If NULL, compute PCA on the subsetted data.
#' @param meta.vars.include Metadata variables to include in the sub-clustered object.
#' @param harmony.group.by.vars Metadata variables to use for Harmony integration.
#' @param early_stop Whether to use early stopping in Harmony.
#' @param return_obj_sub If TRUE, return a list with the updated `obj` and the sub-clustered object.
#' 
#' @returns The input Seurat object with sub-cluster labels added to metadata.
#' 
#' @export
FindSubClusterCustom = function(
    obj, cluster,
    clusters.name = "seurat_clusters",
    sub.clusters.name = NULL,
    resolution = 0.5, algorithm = 1, npcs=30,
    method = 'igraph',
    n_neighbors = 15,
    fast_sgd = FALSE,
    scale.factor = NULL,
    use.existing.embeddings = NULL,
    meta.vars.include = NULL,
    harmony.group.by.vars = NULL,
    early_stop = TRUE,
    return_obj_sub = FALSE
) {
    obj_sub = MakeSubClusterObj(obj, cluster, clusters.name = clusters.name,
                                npcs = npcs, n_neighbors = n_neighbors,
                                scale.factor = scale.factor,
                                use.existing.embeddings = use.existing.embeddings,
                                meta.vars.include = meta.vars.include,
                                harmony.group.by.vars = harmony.group.by.vars,
                                early_stop = early_stop,
                                fast_sgd = fast_sgd)
    if (algorithm == 4) {
        suppressWarnings({    # the igraph conversion spits out a lot of warnings, which is slow if printed
            obj_sub <- FindClusters(object = obj_sub, resolution = resolution, algorithm = algorithm,
                            method = method, graph.name = 'RNA_fgraph')
        })
    } else {
        obj_sub <- FindClusters(object = obj_sub, resolution = resolution, algorithm = algorithm,
                            method = method, graph.name = 'RNA_fgraph')
    }

    obj = AddSubClusterLabels(obj, obj_sub, cluster, clusters.name = clusters.name, sub.clusters.name = sub.clusters.name)

    if (return_obj_sub) {
        list('obj' = obj, 'obj_sub' = obj_sub)
    } else {
        return(obj)
    }
}
