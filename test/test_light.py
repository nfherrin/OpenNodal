import unittest
import driver
import os
import numpy as np


class TestLight(unittest.TestCase):

    dirname = os.path.dirname(__file__)
    executable = os.path.join(dirname, "../OpenNodal.exe")

    k_tol = 1e-6
    err_tol = 1e-5

    def assertg_relative_equal(self, a, b, *args, **kwargs):
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


if __name__ == "__main__":
    unittest.main()
