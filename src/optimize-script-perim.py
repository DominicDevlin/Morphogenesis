import numpy as np
from skopt import Optimizer
from skopt.sampler import Hammersly
import os
import sys
import pandas as pd
import time


# def f(x):
#     return ((x[0]**2) * np.sin(x[0])  + np.random.randn() * 0.1)
# n_int = [(-3.0,3.0)]

# slurm index number to (use csv fiel to read in Jstem and Jdiff params)
# index = int(sys.argv[1])
#name of c++ exec

np.random.seed(int(time.time()))

index = 0
prepend = ""

if (len(sys.argv) > 1):
    index = int(sys.argv[1])
    print("INDEX IS: ", index)
    # prepend = "xvfb-run -a "


p_stem = 1.
p_diff = 1.

# Define the possible values for J_stem and J_diff
# J_stem_values = [1., 2., 3., 4., 5.]
J_diff_values = [0.4,0.5]

J_stem_values = []

val = 0.0
for i in range(12):
    J_stem_values.append(val)
    val += 0.02

# Total number of combinations (5x5 = 25)
n_cols = len(J_stem_values)
n_rows = len(J_diff_values)
num_combinations = len(J_stem_values) * len(J_diff_values)
# Ensure the index is within the valid range
if index >= num_combinations:
    print("Index out of range, should be between 0 and 24.")
    exit()
else:
    # Determine J_stem and J_diff based on the index
    rounder = index // n_cols  # Integer division to determine the row (J_diff)
    leftover = index % n_cols  # Modulo to determine the column (J_stem)

    # Assign the values from the sets
    p_stem = J_stem_values[leftover]
    p_diff = J_diff_values[rounder]

    print(f"Index: {index} => J_stem: {p_stem}, J_diff: {p_diff}")



file_path = 'org-data-' + str(p_stem) + '-' + str(p_diff) + '/optimize.txt'

def rounder(number, amount):
    return round(number * amount) / amount

def f(x, time=0):
    # round params
    # not rounding this x[0] = rounder(x[0], 1/0.1e-3) # not rounding this
    #x[1] = rounder(x[1], 1/0.5)
    # x[1] = rounder(x[1], 1/0.5)
    # x[1] = rounder(x[1], True)
    name = prepend + "./phase-optimize-perim "
    for var in x:
        name = name + str(var) + " "
    name = name + str(p_stem) + " " + str(p_diff) + " " + str(time)
    print(name)
    os.system(name)

    data = pd.read_csv(file_path, delimiter='\t')
    # 0 is average of three best, 1 is best. we are going to use average. 
    n = (data.iloc[-1][0])
    print(n)
    return float(n)


### specify ranges to optimize, each is a tuple with min and max

# differentiation rate, will just be the secretion constant (2.4e-3 is default, 1.5 is about minimum before 0 becomes equilibrium)
diff_rate = [1e-3,5e-3]
# J of stem to diff
J_med = [2.,8.]
# max growth rate per DTS OF stem cells. Taking this out for now.
# growth rate should depnd on J_stem
# Vmax = 1 / (1 + J_stem)
# V_smax = [0.,Vmax]
# doing addtition rate now
min_rate = 100
max_rate = 1200
V_smax = [min_rate, max_rate]


# n dimensional param space
var_list = [diff_rate, V_smax, J_med]


#number of bayesian iterations
iterations = 500


# acq_func_kwargs = {"xi: ": 10}

opt = Optimizer(var_list, n_initial_points=50) 
# n_initial_points should be AT LEAST 50 when loooking for a broad scope? More points can make it worse by dilution?
# Should be like this - I run for 20,000 to 30,000 MCS, or stop when they hit the back wall (or any wall???)
# When I analyse results, each SHOULD be able to hit the back. If it can't hit the back (within a reasonable time frame) MORPHOGENESIS HAS FAILED!!


for i in range(iterations):

    suggested = opt.ask()

    print(suggested)
    y = f(suggested, i)
    opt.tell(suggested, y)
    print('iteration:', i, suggested, y)
    
    