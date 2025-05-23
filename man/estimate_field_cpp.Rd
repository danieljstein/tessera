% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{estimate_field_cpp}
\alias{estimate_field_cpp}
\title{Compute a spatial gradient field at each point (cell)}
\usage{
estimate_field_cpp(coords, embeddings, adj_i, adj_p)
}
\arguments{
\item{coords}{A \code{N} x \code{2} matrix of cell coordinates.}

\item{embeddings}{A \code{N} x \code{D} matrix of cell embeddings.}

\item{adj_i, adj_p}{A \code{N} x \code{N} sparse adjacency matrix
in dgCMatrix format.}
}
\value{
A \code{2} x \code{D} x \code{N} array in column-major ordering
containing the spatial gradient in expression for each of
\code{D} embedding dimensions at every point in space.
}
\description{
Distance between neighboring cells is normalized to unit distance
so that only the direction from each cell to its neighbors matters.
The gradient is then the average gradient in expression of each
embedding dimension between the index cell and its neighbors.
}
