# Smooth embeddings along edges

Smooth embeddings along edges

## Usage

``` r
smooth_embedding(dmt, smooth_emb = 0)
```

## Arguments

- dmt:

  A DMT object.

- smooth_emb:

  Number of smoothing iterations to perform.

## Value

The input DMT object with smoothed embeddings stored in
`dmt$udv_cells$embeddings`.
