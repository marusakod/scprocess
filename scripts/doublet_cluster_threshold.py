"""Estimate a cluster-level doublet cutoff without reading expression data."""
import math
import statistics


def calculate_doublet_cluster_threshold(counts, method='adaptive', ceiling=None):
    """Estimate the adaptive split median/MAD mixture or use a manual cutoff.

    counts contains (number called doublets, total cells) pairs. The half-count
    correction is used for fitting only; exclusion compares raw proportions.
    """
    if method not in ('adaptive', 'manual'):
        raise ValueError('Doublet threshold method must be adaptive or manual')
    if ceiling is None:
        if method == 'manual':
            raise ValueError('Manual doublet threshold requires an explicit cutoff')
        ceiling = 0.5
    if not math.isfinite(ceiling) or not 0 <= ceiling <= 1:
        raise ValueError('Doublet reference ceiling must be between 0 and 1')
    counts = list(counts)
    if not counts:
        raise ValueError('No doublet clusters are available')
    for doublets, total in counts:
        if total <= 0 or not 0 <= doublets <= total:
            raise ValueError('Doublet cluster counts must satisfy 0 <= doublets <= total, total > 0')
    result = dict(method=method, threshold=ceiling, reference_ceiling=ceiling,
                  mad_multiplier=None, n_reference=0,
                  mean_logit=None, scaled_mad_logit=None, uncapped_threshold=None)
    if method == 'manual':
        return result
    return _split_mad_mixture(counts, result)


def _split_mad_mixture(counts, result):
    """Fix groups below/above 50%; select only the between-centres crossing."""
    groups = [[math.log((d + 0.5) / (n - d + 0.5))
               for d, n in counts if (d / n < 0.5 if low else d / n > 0.5)]
              for low in (True, False)]
    if any(len(group) < 2 for group in groups):
        raise ValueError('Split-MAD mixture requires at least two clusters strictly below '
                         'and strictly above 50% doublets. Inspect the data or explicitly '
                         'select manual; no automatic fallback is applied.')
    means = [statistics.median(group) for group in groups]
    scales = [1.4826 * statistics.median(abs(x - centre) for x in group)
              for group, centre in zip(groups, means)]
    if any(scale == 0 for scale in scales):
        raise ValueError('Split-MAD mixture group has zero MAD; inspect the data or '
                         'explicitly select manual.')
    sizes = [len(group) for group in groups]
    weights = [size / sum(sizes) for size in sizes]
    lo, hi = means
    sl, sh = scales
    wl, wh = weights
    a = 0.5 / sl**2 - 0.5 / sh**2
    b = hi / sh**2 - lo / sl**2
    c = math.log(wh / wl) + math.log(sl / sh) + lo**2 / (2 * sl**2) - hi**2 / (2 * sh**2)
    if abs(a) < 1e-10:
        roots = [-c / b] if abs(b) > 1e-10 else []
    else:
        disc = b*b - 4*a*c
        roots = [(-b - math.sqrt(disc)) / (2*a), (-b + math.sqrt(disc)) / (2*a)] if disc >= 0 else []
    interior = [value for value in roots if lo < value < hi]
    if len(interior) != 1:
        raise ValueError('Split-MAD mixture has no unique between-centres crossing; '
                         'inspect the data or explicitly select manual.')
    logit = interior[0]
    threshold = 1 / (1 + math.exp(-logit)) if logit >= 0 else math.exp(logit) / (1 + math.exp(logit))
    result.update(threshold=threshold, uncapped_threshold=threshold, reference_ceiling=0.5,
                  mad_multiplier=None, n_reference=sizes[0], n_low=sizes[0], n_high=sizes[1],
                  n_equal_50=len(counts)-sum(sizes), median_logit_low=lo, median_logit_high=hi,
                  scaled_mad_logit_low=sl, scaled_mad_logit_high=sh, weight_low=wl, weight_high=wh)
    return result
