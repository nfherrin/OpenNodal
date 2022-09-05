import unittest
import driver
import os
import glob
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

        with open(base_input_file, "r") as f:
            dat = f.readlines()

        # TODO work with numbers other than two...
        for i in range(len(dat)):
            line = dat[i]
            if "nsplit" in line:
                dat[i] = "nsplit 2\n"

        tmp_input_file = base_input_file + ".tmp"
        with open(tmp_input_file, "w") as f:
            f.writelines(dat)

        data = driver.run(self.executable, tmp_input_file)
        print("k2=", driver.get_value(data, "XKEFF"))

        # cleanup after ourselves
        for fl in glob.glob(tmp_input_file + "*"):
            os.remove(fl)

        # TODO
        self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
