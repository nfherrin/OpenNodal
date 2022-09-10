import unittest
import driver
import os
import numpy as np


class TestLight(unittest.TestCase):

    dirname = os.path.dirname(__file__)
    executable = os.path.join(dirname, "../OpenNodal.exe")

    k_tol = 1e-6
    err_tol = 1e-5

    def assert_relative_equal(self, a, b, *args, **kwargs):
        # tolerance is as input or sqrt(eps)
        tol = kwargs.get("tol", np.sqrt(np.finfo(float).eps))
        return self.assertTrue(((np.abs(a - b) / np.max([a, b])) <= tol))

    def assert_absolute_equal(self, a, b, *args, **kwargs):
        # tolerance is as input or sqrt(eps)
        tol = kwargs.get("tol", np.sqrt(np.finfo(float).eps))
        return self.assertTrue((np.abs(a - b) <= tol))

    def test_2d_ss_uniform(self):
        input_file = os.path.join(self.dirname, "../examples/2d_ss_uniform.inp")
        data = driver.run(self.executable, input_file)

        keff = driver.get_value(data, "XKEFF")
        kref = 0.9931877117978810
        self.assert_absolute_equal(keff, kref, tol=self.k_tol)

        # NOTE: this doesn't really make sense here.
        # The mesh is too unrefined to give any sort of meaninful data for this metric.
        #### inf = driver.get_value(data, "XFLUX_inf")
        #### infref = 2.277553e-10
        #### self.assert_relative_equal(inf, infref, tol=self.err_tol)

        ratio = driver.get_value(data, "RATIO_calc")
        ratioref = 2.229582e-01
        self.assert_relative_equal(ratio, ratioref, tol=self.err_tol)

    def test_refinement(self):

        base_input_file = os.path.join(self.dirname, "../examples/2d_ss_uniform.inp")

        nrefine = 5  # [1, 2, 4, 8, 16]

        xkeff_refine, xflux_refine = driver.refinement_order(
            self.executable, base_input_file, nrefine
        )

        self.assert_absolute_equal(xkeff_refine, 2.0, tol=0.01)
        self.assert_absolute_equal(xflux_refine, 2.0, tol=0.01)


if __name__ == "__main__":
    unittest.main()
