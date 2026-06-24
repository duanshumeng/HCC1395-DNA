import pandas as pd
import os
import itertools
import sys

input_dir = sys.argv[1]
year = sys.argv[2]
caller = sys.argv[3]

def unique_combinations(arr):
    combos = itertools.combinations(arr, 2)
    unique_combos = set(combos)
    return unique_combos

aim_arr = []

for i in os.listdir(input_dir):
    if year in i and caller in i and i.endswith('.gz'):
        aim_arr.append(i)

print(aim_arr)
