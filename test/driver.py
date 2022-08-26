import subprocess
import sys


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
