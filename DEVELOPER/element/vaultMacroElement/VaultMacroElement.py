import csv
import sys

import numpy as np

from Kinematics import *


acc_file_name = "acc.csv"
result_file_name = "result.csv"


def read_csv(file_name):
    with open(file_name, "r") as f:
        reader = csv.reader(f)
        data = list(reader)
    return data


def write_csv(data, file_name):
    with open(file_name, "w", newline="") as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)


dt = sys.argv[1]
theta = sys.argv[2]
acc = read_csv(acc_file_name)
model = Model_TH()


# Do something


result = [[1, 2], [3, 4], [5, 6]]
write_csv(result, result_file_name)
