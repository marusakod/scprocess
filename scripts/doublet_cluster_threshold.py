"""Estimate a cluster-level doublet cutoff without reading expression data."""
import math
import statistics


def calculate_doublet_cluster_threshold(counts, method='mean_mad', ceiling=0.5,
                                       mad_multiplier=2.0):
    """Use equally weighted clusters below ceiling to fit mean + scaled MAD.

    counts contains (number called doublets, total cells) pairs. The half-count
    correction is used for fitting only; exclusion compares raw proportions.
    """
    if method not in ('mean_mad', 'fixed'):
        raise ValueError('Doublet threshold method must be mean_mad or fixed')
    if not math.isfinite(ceiling) or not 0 <= ceiling <= 1:
        raise ValueError('Doublet reference ceiling must be between 0 and 1')
    if not math.isfinite(mad_multiplier) or mad_multiplier < 0:
        raise ValueError('Doublet MAD multiplier must be finite and non-negative')
    counts = list(counts)
    if not counts:
        raise ValueError('No doublet clusters are available')
    for doublets, total in counts:
        if total <= 0 or not 0 <= doublets <= total:
            raise ValueError('Doublet cluster counts must satisfy 0 <= doublets <= total, total > 0')
    result = dict(method=method, threshold=ceiling, reference_ceiling=ceiling,
                  mad_multiplier=mad_multiplier, n_reference=0,
                  mean_logit=None, scaled_mad_logit=None, uncapped_threshold=None)
    if method == 'fixed':
        return result
    reference = [(d, n) for d, n in counts if d / n < ceiling]
    if not reference:
        raise ValueError('No doublet clusters fall below the reference ceiling; '
                         'cannot estimate mean-MAD threshold. Inspect the data or explicitly use fixed.')
    logits = [math.log((d + 0.5) / (n - d + 0.5)) for d, n in reference]
    mean = statistics.mean(logits)
    centre = statistics.median(logits)
    mad = 1.4826 * statistics.median(abs(value - centre) for value in logits)
    if mad == 0:
        raise ValueError('Reference cluster logit rates have zero MAD; cannot estimate an adaptive threshold. '
                         'Inspect the data or explicitly use fixed.')
    cutoff = mean + mad_multiplier * mad
    # Stable inverse logit, including extreme but valid multiplier values.
    estimate = 1 / (1 + math.exp(-cutoff)) if cutoff >= 0 else math.exp(cutoff) / (1 + math.exp(cutoff))
    result.update(threshold=min(ceiling, estimate), n_reference=len(reference),
                  mean_logit=mean, scaled_mad_logit=mad, uncapped_threshold=estimate)
    return result
