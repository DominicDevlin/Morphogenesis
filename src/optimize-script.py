import numpy as np
from skopt import gp_minimize
from skopt import Optimizer
import os
import sys
import pandas as pd

# def f(x):
#     return ((x[0]**2) * np.sin(x[0])  + np.random.randn() * 0.1)
# n_int = [(-3.0,3.0)]

# slurm index number to (use csv fiel to read in Jstem and Jdiff params)
# index = int(sys.argv[1])
#name of c++ exec

J_stem = 3.
J_diff = 12.

file_path = 'org-data/optimize.txt'

def rounder(number, amount):
    return round(number * amount) / amount

def f(x, time):
    # round params
    x[0] = rounder(x[0], 1/0.2e-3)
    x[1] = rounder(x[1], 1/0.5)
    x[2] = rounder(x[2], 1/0.5)
    x[3] = rounder(x[3], 1/0.5)
    x[4] = rounder(x[4], 1/0.5)
    name = "./phase-optimize "
    for var in x:
        name = name + str(var) + " "
    name = name + str(J_stem) + " " + str(J_diff) + " " + str(time)   
    print(name)
    os.system(name)

    data = pd.read_csv(file_path, delimiter='\t')
    n = (data.iloc[-1][0])
    print(n)
    return float(n)

    


### specify ranges to optimize, each is a tuple with min and max

# differentiation rate, will just be the secretion constant (2.4e-3 is default, 1.5 is about minimum before 0 becomes equilibrium)
diff_rate = (1.4e-3,1e-2)
# J of cells with medium
Jmed = [J_stem, J_diff + 3]
# J of stem to diff
Jsd = [J_stem, J_diff]
# max growth rate per DTS OF stem cells (need to sort out this implementation)
V_smax = [0.5,5]
V_dmax = [0.5,5]

# 5 dimensional param space
var_list = [diff_rate, Jmed, Jsd, V_smax, V_dmax]

#number of bayesian iterations
iterations = 200

# 6.5
# 12
# 3.5
# 1
# 4
# 12
# 19


inits = [2.4e-3, 6.5, 12, 1, 1]

opt = Optimizer(var_list)

for i in range(iterations):
    suggested = []
    if i == 0:
        suggested = inits
    else:
        suggested = opt.ask()
    print(suggested)

    y = f(suggested, i)
    opt.tell(suggested, y)
    print('iteration:', i, suggested, y)
    
    # for var in y
    
    
#     os.system("./" + name + " " + threshold + " " + mutation + " " + kp + " " + dp + " " + kq + " " + dq)
#     # dont forget to send iteration number to c++