import math
import unittest

from scripts.doublet_cluster_threshold import calculate_doublet_cluster_threshold


class TestDoubletClusterThreshold(unittest.TestCase):
    def test_manual_cutoff_and_strict_boundary(self):
        fit = calculate_doublet_cluster_threshold([(30, 100), (31, 100)], method='manual', ceiling=.3)
        self.assertEqual(fit['method'], 'manual')
        self.assertEqual(fit['threshold'], .3)
        self.assertEqual([d/n > fit['threshold'] for d,n in [(30,100),(31,100)]], [False,True])
        self.assertIsNone(fit['mean_logit'])

    def test_reject_removed_names(self):
        for method in ('split_mad_mixture', 'mean_mad', 'fixed', 'unknown'):
            with self.subTest(method=method), self.assertRaisesRegex(ValueError, 'adaptive or manual'):
                calculate_doublet_cluster_threshold([(1,10)], method=method)

    def test_schema_matches_public_methods(self):
        import json
        from pathlib import Path
        schema = json.loads((Path(__file__).resolve().parents[1] / 'resources/schemas/config.schema.json').read_text())
        integration = schema['properties']['integration']['properties']
        self.assertEqual(integration['int_dbl_threshold_method']['enum'], ['adaptive','manual'])
        self.assertEqual(integration['int_dbl_threshold_method']['default'], 'adaptive')
        self.assertNotIn('int_dbl_mad_multiplier', integration)

    def test_schema_rejects_adaptive_with_cutoff(self):
        import json
        from pathlib import Path
        import jsonschema
        schema = json.loads((Path(__file__).resolve().parents[1] / 'resources/schemas/config.schema.json').read_text())['properties']['integration']
        for supplied in ({'int_dbl_threshold_method': 'adaptive', 'int_dbl_cl_prop': .3},
                         {'int_dbl_cl_prop': .3},
                         {'int_dbl_threshold_method': 'manual'}):
            with self.subTest(supplied=supplied), self.assertRaises(jsonschema.ValidationError):
                jsonschema.validate(dict(int_use_paga=False, **supplied), schema)
        for supplied in ({'int_dbl_threshold_method': 'adaptive'}, {},
                         {'int_dbl_threshold_method': 'manual', 'int_dbl_cl_prop': .3}):
            jsonschema.validate(dict(int_use_paga=False, **supplied), schema)

    def test_default_resolution_preserves_method_validation(self):
        import json
        from pathlib import Path
        import jsonschema
        from scripts.scprocess_utils import _get_default_config_from_schema
        from snakemake.utils import update_config
        schema = json.loads((Path(__file__).resolve().parents[1] / 'resources/schemas/config.schema.json').read_text())
        for method in (None, 'adaptive', 'manual'):
            supplied = {'integration': {}} if method is None else {'integration': {'int_dbl_threshold_method': method}}
            effective = _get_default_config_from_schema(schema, supplied)['integration']
            self.assertNotIn('int_dbl_cl_prop', effective)
            if method == 'manual':
                supplied['integration']['int_dbl_cl_prop'] = .3
            update_config(effective, supplied['integration'])
            jsonschema.validate(effective, schema['properties']['integration'])
        supplied = {'integration': {'int_dbl_threshold_method': 'adaptive', 'int_dbl_cl_prop': .3}}
        effective = _get_default_config_from_schema(schema, supplied)['integration']
        update_config(effective, supplied['integration'])
        with self.assertRaises(jsonschema.ValidationError):
            jsonschema.validate(effective, schema['properties']['integration'])

    def test_manual_requires_explicit_cutoff(self):
        with self.assertRaisesRegex(ValueError, 'explicit cutoff'):
            calculate_doublet_cluster_threshold([(1,10)], method='manual')

    def test_invalid_parameters(self):
        for kwargs in ({'ceiling': 1.1}, {'ceiling': -1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                calculate_doublet_cluster_threshold([(1,10)], **kwargs)
        for counts in ([], [(11,10)], [(0,0)]):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                calculate_doublet_cluster_threshold(counts)


class TestSplitMadMixture(unittest.TestCase):
    def test_default_robust_parameters_and_weights(self):
        counts = [(1, 100), (5, 100), (10, 100), (60, 100), (80, 100), (50, 100)]
        fit = calculate_doublet_cluster_threshold(counts)
        self.assertEqual(fit['method'], 'adaptive')
        self.assertEqual((fit['n_low'], fit['n_high'], fit['n_equal_50']), (3, 2, 1))
        self.assertAlmostEqual(fit['weight_low'], 3/5)
        self.assertAlmostEqual(fit['median_logit_low'], math.log(5.5/95.5))
        self.assertEqual(fit['threshold'], fit['uncapped_threshold'])
        x = math.log(fit['threshold']/(1-fit['threshold']))
        def density(side):
            mu = fit['median_logit_' + side]
            sd = fit['scaled_mad_logit_' + side]
            return fit['weight_' + side] / sd * math.exp(-.5*((x-mu)/sd)**2)
        self.assertAlmostEqual(density('low'), density('high'))
        self.assertLess(fit['median_logit_low'], x)
        self.assertLess(x, fit['median_logit_high'])

    def test_no_cap_and_legacy_ceiling_does_not_change_split(self):
        counts = [(20,100), (30,100), (40,100), (85,100), (90,100), (95,100)]
        fit = calculate_doublet_cluster_threshold(counts, ceiling=.2)
        self.assertGreater(fit['threshold'], .5)
        self.assertEqual(fit['reference_ceiling'], .5)

    def test_insufficient_or_zero_mad_group_errors(self):
        for counts in ([(1,100),(5,100),(90,100)],
                       [(1,100),(1,100),(70,100),(90,100)]):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                calculate_doublet_cluster_threshold(counts)


if __name__ == '__main__':
    unittest.main()
