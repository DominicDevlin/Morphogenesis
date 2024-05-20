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
    prepend = "xvfb-run -a "


J_stem = 10.
J_diff = 12.

J_stem += index

file_path = 'org-data-' + str(J_stem) + '/optimize.txt'

def rounder(number, amount):
    return round(number * amount) / amount

def f(x, time=0):
    # round params
    # x[0] = rounder(x[0], 1/0.1e-3) # not rounding this
    x[1] = rounder(x[1], 1/0.5)
    x[2] = rounder(x[2], 1/0.5)
    x[3] = rounder(x[3], True)
    name = prepend + "./phase-optimize "
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
diff_rate = [1.4e-3,0.025]
# J of cells with medium
Jmed = [0.5*J_diff, J_diff + 3]
# J of stem to diff
Jsd = [J_stem, J_diff + J_diff]
# max growth rate per DTS OF stem cells. Taking this out for now.
# V_smax = [0.,1.]
# V_dmax = [0.,1.]
# v-V threshold, usually 2
gthresh = [0.,3.]

# 5 dimensional param space
var_list = [diff_rate, Jmed, Jsd, gthresh]

# do random sampling
# hammer = Hammersly()
# n_samples = 10
# x_samples = hammer.generate(var_list, n_samples=n_samples)
# hammer_results = []
# start_time = 0 - n_samples
# for i in x_samples:
#     # print("random sample: ", i)
#     hres = f(i, start_time)
#     hammer_results.append(hres)
#     print(i, "Hammersly random sampling result is: ", hres)
#     start_time+=1

# for i in range(len(x_samples)):
#     opt.tell(x_samples[i], hammer_results[i])

#number of bayesian iterations
iterations = 800

inits=[]

## first punt
punt_sec_rate = 2.039e12*pow((J_stem+14.567),-12.1771)+0.0018588
punt_J_med = 0.5 + 0.5*J_diff
punt_J_sd = 12
# punt_vsmax = 1
# punt_vdmax = 1
punt_gthresh = 2
inits.append([punt_sec_rate, punt_J_med, punt_J_sd, punt_gthresh])

### second punt
punt_sec_rate = 128.123*pow((J_stem+3.66212),-5.64574)+0.00194831
punt_J_med = 0.5*J_diff+2
punt_J_sd = 14.5
# punt_vsmax = 1
# punt_vdmax = 1
punt_gthresh = 1
inits.append([punt_sec_rate, punt_J_med, punt_J_sd, punt_gthresh])


punt_sec_rate = 0.0223029*pow((J_stem+0.757648),-1.79529)+0.00182658
if punt_sec_rate > diff_rate[1]:
    punt_sec_rate = diff_rate[1]
punt_J_med = J_stem+1
punt_J_sd = punt_J_med*2
# punt_vsmax = 1
# punt_vdmax = 1
punt_gthresh = 0
inits.append([punt_sec_rate, punt_J_med, punt_J_sd, punt_gthresh])


# acq_func_kwargs = {"xi: ": 10}

opt = Optimizer(var_list, n_initial_points=0) 
# n_initial_points should be AT LEAST 50 when loooking for a broad scope?
# Should be like this - I run for 20,000 to 30,000 MCS, or stop when they hit the back wall (or any wall???)
# When I analyse results, each SHOULD be able to hit the back. If it can't hit the back (within a reasonable time frame) MORPHOGENESIS HAS FAILED!!


for i in range(iterations):
    suggested = []
    if i < len(inits):
        suggested = inits[i]
    else:
        suggested = opt.ask()

    print(suggested)
    y = f(suggested, i)
    opt.tell(suggested, y)
    print('iteration:', i, suggested, y)
    
    
    
