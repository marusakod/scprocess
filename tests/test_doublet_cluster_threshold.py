import math
import unittest

from scripts.doublet_cluster_threshold import calculate_doublet_cluster_threshold


class TestDoubletClusterThreshold(unittest.TestCase):
    def test_reference_filter_and_equal_cluster_weights(self):
        counts = [(1, 10), (2, 10), (3000, 10000), (5, 10), (9, 10)]
        fit = calculate_doublet_cluster_threshold(counts)
        expected_logits = [math.log(1.5 / 9.5), math.log(2.5 / 8.5), math.log(3000.5 / 7000.5)]
        self.assertEqual(fit['n_reference'], 3)
        self.assertAlmostEqual(fit['mean_logit'], sum(expected_logits) / 3)
        self.assertGreater(fit['scaled_mad_logit'], 0)
        self.assertLessEqual(fit['threshold'], 0.5)

    def test_cap_never_weakens_fixed_ceiling(self):
        fit = calculate_doublet_cluster_threshold([(1, 10), (2, 10), (3, 10)], mad_multiplier=100)
        self.assertGreater(fit['uncapped_threshold'], 0.5)
        self.assertEqual(fit['threshold'], 0.5)

    def test_zero_rates_have_finite_logits(self):
        fit = calculate_doublet_cluster_threshold([(0, 100), (1, 100), (10, 100)])
        self.assertTrue(math.isfinite(fit['mean_logit']))
        self.assertTrue(math.isfinite(fit['scaled_mad_logit']))
        self.assertGreater(fit['threshold'], 0)

    def test_zero_mad_errors_without_fallback(self):
        with self.assertRaisesRegex(ValueError, 'zero MAD'):
            calculate_doublet_cluster_threshold([(1, 100)] * 10)

    def test_no_reference_clusters_errors_without_fallback(self):
        with self.assertRaisesRegex(ValueError, 'No doublet clusters'):
            calculate_doublet_cluster_threshold([(50, 100), (100, 100)])

    def test_fixed_method_does_not_fit(self):
        fit = calculate_doublet_cluster_threshold([(50, 100), (100, 100)], method='fixed')
        self.assertEqual(fit['threshold'], 0.5)
        self.assertIsNone(fit['mean_logit'])

    def test_invalid_parameters(self):
        for kwargs in ({'ceiling': 1.1}, {'ceiling': -1}, {'mad_multiplier': -1},
                       {'mad_multiplier': float('nan')}, {'method': 'unknown'}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                calculate_doublet_cluster_threshold([(1, 10)], **kwargs)
        with self.assertRaises(ValueError):
            calculate_doublet_cluster_threshold([(11, 10)])


if __name__ == '__main__':
    unittest.main()
