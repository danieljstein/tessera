#' Generate mesh data structure from coordinates, for input to DMT analysis
#' 
#' Creates a mesh from point coordinates using Delauney triangulation.
#' Stores the points, triangles, and edges for the mesh, as well
#' as a mapping from each triangle to its associated vertices.
#' A gene-by-cell matrix can also be included, as well as additional
#' metadata for each point.
#' 
#' @param X,Y A pair of numeric vectors with the coordinates for each of V points.
#' @param counts A G x V gene-by-cell matrix of transcript counts.
#' @param meta_data A data frame with additional cell metadata to include.
#' @param meta_vars_include Names of columns in meta_data to include.
#' 
#' @returns A list containing the mesh data structures:
#' * `pts` is a V-by-M data table with columns `X` and `Y` containing the
#'   coordinates of cells and additional metadata.
#' * `tris` is a F-by-4 data table containing the X,Y coordinates of each
#'    triangle's centroid in the first two columns, and area and
#'    largest height of each triangle in the last two columns.
#' * `edges` is a E-by-14 data table with columns `from_pt`, `to_pt`, `from_tri`, `to_tri`,
#'   `x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`, `x1_tri`, `y0_tri`, `y1_tri`,
#'   `length_pt`, `length_tri`. If only one triangle uses an edge, then the `from_tri`,
#'   `x0_tri`, and `y0_tri` fields will contain NaN values.
#' * `tri_to_pt` is a F-by-V sparse matrix with value 1 at (i,j) if
#'   triangle i uses point j as a vertex.
#' * `counts` is a G x V gene-by-cell matrix of transcript counts.
#' 
#' Note that all indices stored in these data structures are 1-indexed.
#' 
#' @export
init_data = function(X, Y, counts, meta_data=NULL, meta_vars_include=c()) {
    pts = cbind(X=X, Y=Y)
    
    triplets = geometry::delaunayn(pts)
    tris = init_tris_cpp(triplets, pts)
    colnames(tris) = c('X', 'Y', 'area', 'height')
    tri_to_pt = Matrix::sparseMatrix(i = rep(seq_len(nrow(triplets)), each = 3), j = c(t(triplets)), x = 1)
    
    edges = init_edges_cpp(triplets, as.matrix(pts), as.matrix(tris))
    colnames(edges) = c(
        'from_pt', 'to_pt', 'from_tri', 'to_tri',
        'x0_pt', 'x1_pt', 'y0_pt', 'y1_pt', 
        'x0_tri', 'x1_tri', 'y0_tri', 'y1_tri',
        'length_pt', 'length_tri'
    )
    data = list(
        pts = data.table(pts), 
        tris = data.table(tris), 
        edges = data.table(edges), 
        tri_to_pt = tri_to_pt,
        # udv_cells = udv_cells, 
        counts = counts
    )    
    data$pts[, ORIG_ID := 1:.N]
    if (!is.null(meta_data) & length(meta_vars_include) > 0) {
        for (colname in meta_vars_include) {
            if (colname %in% colnames(meta_data)) {
                data$pts[[colname]] = meta_data[[colname]]                
            } else {
                warning(glue('{colname} not in meta_data'))
            }
            # dmt$pts[[colname]] = meta_data[[colname]][dmt$pts$ORIG_ID]
        }
    }
    
    return(data)
}

#' Prunes mesh by eliminating long edges and small connected components
#' 
#' Pruning works in three steps:
#' 1. Any triangle with an edge longer than a threshold length is removed.
#'    Afterwards, any edge that no longer belongs to a triangle is removed.
#'    Then any point that no longer belongs to an edge is removed.
#' 2. Connected components of triangles with a shared edge are computed.
#'    If everything is connected, then nothing is pruned.
#' 3. Otherwise, components (and corresponding triangles) that contain less
#'    than `mincells` points are removed. Afterwards, any edge that no
#'    longer belongs to a triangle is removed. Then any point that no longer
#'    belongs to an edge is removed.
#' 
#' @param data A list containing the mesh data structures:
#' * `pts` is a V-by-M data table with columns `X` and `Y` containing the
#'   coordinates of cells and additional metadata.
#' * `tris` is a F-by-4 data table containing the X,Y coordinates of each
#'    triangle's centroid in the first two columns, and area and
#'    largest height of each triangle in the last two columns.
#' * `edges` is a E-by-14 data table with columns `from_pt`, `to_pt`, `from_tri`, `to_tri`,
#'   `x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`, `x1_tri`, `y0_tri`, `y1_tri`,
#'   `length_pt`, `length_tri`. If only one triangle uses an edge, then the `from_tri`,
#'   `x0_tri`, and `y0_tri` fields will contain NaN values.
#' * `tri_to_pt` is a F-by-V sparse matrix with value 1 at (i,j) if
#'   triangle i uses point j as a vertex.
#' @param thresh_quantile Floating point value between 0 and 1, inclusive.
#'   Quantile of edge length above which edges are pruned. Defaults to 0.95.
#' @param mincells Minimum number of cells required for a connected
#'   component of triangles to be kept. Defaults to 10.
#' @param thresh Edge length above which edges are pruned. If equal to NA,
#'   then this value is ignored and `thresh_quantile` is used to compute
#'   the threshold. Otherwise, if `thresh` is set, then `thresh_quantile`
#'   is ignored. Defaults to NA.
#' 
#' @returns A list containing the mesh data structures with (possibly)
#'   fewer points, edges, and triangles. Indices have been updated since
#'   some objects might have been removed.
#' 
#' @export
prune_graph = function(data, thresh_quantile = .95, mincells = 10, thresh = NA) {    
    if (is.na(thresh)) {
        ## no cutting 
        if (thresh_quantile >= 1) {
            return(data)
        } else {
            thresh = quantile(data$edges$length_pt, thresh_quantile)
        }
    }
    
    pts = data$pts
    edges = data$edges 
    tris = data$tris
    tri_to_pt = data$tri_to_pt

    edges[, qc := length_pt <= thresh]
    pts_keep = edges[qc == TRUE]
    pts_keep = with(pts_keep, union(from_pt, to_pt))
    pts$qc = FALSE
    pts$qc[pts_keep] = TRUE
    
    ## (2) Remove incomplete triangles
    ## First, find valid triangles
    i_tri_keep = edges[qc == TRUE][, .(c(to_tri, from_tri))][!is.na(V1)][, .N, by = V1][N == 3, V1]
    tris[, id := 1:.N][, qc := id %in% i_tri_keep][, id := NULL]
    ## Next, filter edges
    edges[, qc := (
        (!is.na(from_tri) & (from_tri %in% i_tri_keep | to_tri %in% i_tri_keep)) | 
        (is.na(from_tri) & to_tri %in% i_tri_keep)
    )]
    
    ## Filter points now 
    pts_keep = edges[qc == TRUE, .(union(from_pt, to_pt))]$V1
    pts$qc = FALSE
    pts$qc[pts_keep] = TRUE
    
    ## (3) Remove small components 
    ## First, get components of connected triangles 
    g = edges[
        qc == TRUE & from_tri %in% which(tris$qc == TRUE) & to_tri %in% which(tris$qc == TRUE), 
        .(from_tri, to_tri)
    ]
    g = as.matrix(g) 
    g = igraph::from_edgelist()$fun(g, directed = FALSE) 
    # add missing isolated triangles to graph
    g = igraph::add_vertices(g, nv = nrow(tris) - igraph::gorder(g))
    
    comps = igraph::components(g)$membership
    tris$comp = factor(comps) 
    
    ## if only one component, then keep everything
    if (nlevels(tris$comp) == 1) {
        return(data)
    } 
    
    ## Next, find number of (unique) points in each component
    ##       and label the good components 
    pts_by_comps = Matrix::t(tri_to_pt) %*% Matrix::sparse.model.matrix(~0+comp, tris)
    colnames(pts_by_comps) = gsub('^comp', '', colnames(pts_by_comps))
    pts_by_comps@x = rep(1, length(pts_by_comps@x)) ## count each cell once 
    pts_per_comp = Matrix::colSums(pts_by_comps)
    comps_keep = as.integer(names(which(pts_per_comp >= mincells)))
    
    ## keep triangles in good components 
    tris$qc = tris$qc & tris$comp %in% comps_keep
    
    ## keep points in good components 
    ## NOTE: 
    ##     >0 means you belong to at least one good component
    ##     ==1 means that you cannot bridge between two components 
    ##     ==1 does not prevent weakly connected holes in the same component 
    pts$qc = pts$qc & (Matrix::rowSums(pts_by_comps[, comps_keep, drop=F]) > 0)
    # pts$qc = pts$qc & (rowSums(pts_by_comps[, comps_keep]) == 1) 
    
    ## Propagate it down to edges 
    pts_keep = which(pts$qc == TRUE)
    edges[, qc := qc & (from_pt %in% pts_keep & to_pt %in% pts_keep)]
    ## remove edges that connect exterior triangles
    ## then also remove points, as needed 
    edges$qc = edges$qc & (tris$qc[edges$from_tri] | tris$qc[edges$to_tri])
    ## Check if any points are now disconnected 
    pts_disconnected = setdiff(seq_len(nrow(pts)), edges[qc == TRUE, .(from_pt, to_pt)][, .(union(from_pt, to_pt))]$V1)
    pts$qc[pts_disconnected] = FALSE
    
    ## Finally, remove all QC'd out data 
    edges = edges[qc == TRUE]
    pts[qc == TRUE, id := 1:.N]
    edges$from_pt = pts$id[edges$from_pt]
    edges$to_pt = pts$id[edges$to_pt]
    
    tris[qc == TRUE, id := 1:.N]
    edges$from_tri = tris$id[edges$from_tri]
    edges$to_tri = tris$id[edges$to_tri]
    i_remove = which(is.na(edges$to_tri))
    edges$x1_tri[i_remove] = edges$x0_tri[i_remove]
    edges$y1_tri[i_remove] = edges$y0_tri[i_remove]
    edges$to_tri[i_remove] = edges$from_tri[i_remove]
    edges$from_tri[i_remove] = NA
    edges[is.na(from_tri), `:=`(x0_tri = NA, y0_tri = NA)]
    
    ## At this point, qc flag is finalized
    ## you can remove tris and pts that did not pass QC 
    tris_qc = which(tris$qc == TRUE)
    tris = tris[tris_qc, ]
    pts_qc = which(pts$qc == TRUE)
    pts = pts[pts_qc, ]    
    data$counts = data$counts[, pts_qc]
    # data$udv_cells$embeddings = data$udv_cells$embeddings[pts_qc, ]    
    data$tri_to_pt = data$tri_to_pt[tris_qc, pts_qc]
    
    pts[, `:=`(id=NULL, qc=NULL)]
    tris[, `:=`(id=NULL, qc=NULL, comp=NULL)]
    edges[, qc := NULL]
    data$pts = pts
    data$tris = tris
    data$edges = edges
    
    return(data)
}

#' For every edge on the boundary, adds a second degenerate triangle
#' 
#' Edges at the boundary will only border a single triangle.
#' In order to deal with boundary conditions properly, this step
#' ensures that every edge has two associated triangles (which use
#' the edge) as well as two associated points (endpoints). The
#' triangles that are added are degenerate triangles centered at
#' the midpoint of each boundary edge.
#' 
#' @param data A list containing the mesh data structures:
#' * `pts` is a V-by-M data table with columns `X` and `Y` containing the
#'   coordinates of cells and additional metadata.
#' * `tris` is a F-by-4 data table containing the X,Y coordinates of each
#'    triangle's centroid in the first two columns, and area and
#'    largest height of each triangle in the last two columns.
#' * `edges` is a E-by-14 data table with columns `from_pt`, `to_pt`, `from_tri`, `to_tri`,
#'   `x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`, `x1_tri`, `y0_tri`, `y1_tri`,
#'   `length_pt`, `length_tri`. If only one triangle uses an edge, then the `from_tri`,
#'   `x0_tri`, and `y0_tri` fields will contain NaN values.
#' * `tri_to_pt` is a F-by-V sparse matrix with value 1 at (i,j) if
#'   triangle i uses point j as a vertex.
#' 
#' @returns A list containing the mesh data structures with (possibly)
#'   additional triangles. The `tris` table is updated to include the added
#'   external triangles. The `edge` table is updated so that boundary edges
#'   point to the newly added triangles. Also `pts` is unchanged, `tri_to_pt`
#'   is updated. Because new external triangles only have 2 associated points,
#'   `tri_to_pt` has rows that sum to either 2 or 3.
#' 
#' @export
add_exterior_triangles = function(data) {
    edges = data$edges
    tris = data$tris
    pts = data$pts
    tri_to_pt = data$tri_to_pt

    n_boundary_edges = sum(is.na(edges$from_tri))
    
    ## For each boundary edge, 
    ##   create a new external triangle
    ##   coordinates are mean of point edges 
    ##   area = 0
    ##   comp = NA 
    edges[, boundary := FALSE]
    edges[
        is.na(from_tri), 
        `:=`(
            from_tri = (nrow(tris)+1):(nrow(tris)+n_boundary_edges), 
            boundary = TRUE
        )
    ]
    edges[boundary==TRUE, `:=`(x0_tri = .5 * (x0_pt + x1_pt), y0_tri = .5 * (y0_pt + y1_pt))]
    tris = rbindlist(list(
        tris[, external := FALSE], 
        edges[boundary == TRUE][, .(X=x0_tri, Y=y0_tri, area=0, height=0, external=TRUE)]
    ))
    
    ## Add map from external tris to their points 
    tri_to_pt_external = Matrix::sparseMatrix(
        i = rep(1:n_boundary_edges, each = 2), 
        j = as.numeric(t(as.matrix(edges[boundary == TRUE][, .(from_pt, to_pt)]))), 
        x = 1, 
        dims = c(n_boundary_edges, nrow(pts))
    )
    tri_to_pt = Matrix::rbind2(tri_to_pt, tri_to_pt_external)    
    
    data$edges = edges
    data$tris = tris
    data$tri_to_pt = tri_to_pt
    return(data)
}




