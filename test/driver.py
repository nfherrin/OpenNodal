import subprocess
import sys
import numpy as np
import glob
import os

def run(executable, input_file, *args):
    # runs executable with argument input_file
    # expecting input_file to be a file name
    # any additional arguments may be listed as well
    # returns --- stdout & stderr as a list split into lines

    this_list = [executable, input_file]
    for arg in args:
        this_list.append(arg)

    return (
        subprocess.check_output(this_list, stderr=subprocess.STDOUT)
        .decode(sys.stdout.encoding)
        .strip()
        .splitlines()
    )


def get_value(data, value):
    # given list of string data (typically from running a code), look for a line
    # beginning with "value" and return a floating point number after this "value"
    # e.g. "klong = 1.0" return 1.0
    # skips the "=" as well
    for line in data:
        this_line = line.split()
        if len(this_line) > 0:
            if this_line[0] == value:
                return float(this_line[2])


def input_replace(base_input_file, tmp_input_file, card, value):

    with open(base_input_file, "r") as f:
        dat = f.readlines()

    found = False

    for i in range(len(dat)):
        line = dat[i]
        if card in line:
            dat[i] = card + " " + value + "\n"
            found = True
            break

    if not found:
        print(
            "Tried to modify card "
            + card
            + " in input file "
            + base_input_file
            + " but card was not found."
        )
        return

    with open(tmp_input_file, "w") as f:
        f.writelines(dat)

    return


def refinement_order(executable, base_input_file, refine_max):

    tmp_input_file = base_input_file + ".tmp"

    xkerr_arr = np.zeros(refine_max)
    xferr_arr = np.zeros(refine_max)

    for refine in range(refine_max):

        nsplit = 2 ** refine
        input_replace(base_input_file, tmp_input_file, "nsplit", str(nsplit))
        data = run(executable, tmp_input_file)

        xkerr_arr[refine] = get_value(data, "XKEFF_err")
        # finite difference convergence is second-order in inf norm
        xferr_arr[refine] = get_value(data, "XFLUX_inf")

        for fl in glob.glob(tmp_input_file + "*"):
            os.remove(fl)

    # Use the last 3 points to calculate the refinement order (e.g., Pilch formula).
    # Alternatively, could use a power fit.
    # See Eq. (3.2) in W. C. Dawn dissertation.

    refine = lambda q, r : np.log((q[-3] - q[-2])/(q[-2] - q[-1]))/np.log(r)
    xkeff_refine = refine(xkerr_arr, 2.0)
    xflux_refine = refine(xferr_arr, 2.0)

    return xkeff_refine, xflux_refine
