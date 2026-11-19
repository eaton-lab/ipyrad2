"""Resolve assembly calling controls before starting native work."""
from ..utils.exceptions import IPyradError


def resolve_calling_options(model, sample_depth, locus_depth, geno_q, site_q,
                            fit_loci=1024, folds=5, max_iter=10000):
    if model not in ('native', 'sfs'):
        raise IPyradError('genotype_model must be native or sfs.')
    sample_depth = (1 if model == 'sfs' else 5) if sample_depth is None else sample_depth
    locus_depth = (5 if model == 'sfs' else sample_depth) if locus_depth is None else locus_depth
    for label, value in [('min_sample_depth', sample_depth), ('min_locus_depth', locus_depth)]:
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise IPyradError(f'{label} must be a positive integer.')
    for label, value in [('min_geno_q', geno_q), ('min_site_q', site_q)]:
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise IPyradError(f'{label} must be a nonnegative integer.')
    if model == 'sfs':
        for label, value, lower in [('sfs_fit_loci', fit_loci, 1), ('sfs_folds', folds, 2),
                                    ('sfs_max_iter', max_iter, 1)]:
            if isinstance(value, bool) or not isinstance(value, int) or value < lower:
                raise IPyradError(f'{label} must be an integer >= {lower}.')
    return sample_depth, locus_depth
