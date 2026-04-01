# Analysis Guide

> Draft: This section is the guide to downstream analyses that start from assembled outputs or SNP-capable HDF5 inputs.

## Summary

Analysis Guide will document the methods that run on filtered SNP or sequence data inside ipyrad2 or through wrapped external tools.

## Who This Section Is For

- Users moving from assembled outputs into downstream analysis
- Users choosing between built-in analyses and exported external workflows

## Page Map

- [pca](./pca.md): PCA, t-SNE, and UMAP on SNP-capable HDF5 data
- [dapc](./dapc.md): sklearn-backed DAPC-style clustering on SNP-capable HDF5 data
- [popgen](./popgen.md): built-in population-genetic summary statistics from sequence or SNP-capable HDF5 inputs

## Recommended Reading Order

1. Start here for method selection.
2. Read [pca](./pca.md) if you want PCA-family ordinations from SNP HDF5 data.
3. Read [dapc](./dapc.md) if you want DAPC-style clustering and discriminant coordinates from SNP HDF5 data.
4. Read [popgen](./popgen.md) if you want diversity, differentiation, heterozygosity, or SFS summaries inside ipyrad2.
5. Return to Writing Outputs if you need an exported format instead.

## Related Sections

- Writing Outputs
- Files and Data Types

## Open TODOs

- Add an overview table covering PCA-family, DAPC, admixture, sNMF, popgen, and converters.
