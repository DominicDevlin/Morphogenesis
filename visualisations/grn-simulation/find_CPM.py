import math
import os 


def CPM_attractions():

    directory = 'conc'

    states = {}

    count = 0

    # iterate over files in
    # that directory
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            # DO ALL THINGS NEEDED FOR SINGLE FILE HERE
            # first - sum all states. Only count staes if > 8 percent?

            insides = open(f, 'r')
            lines = insides.readlines()
            i = 120
            for i in range(120, len(lines)):
                bools = []
                strs = lines[i].split("\t")
                for n in range(9):
                    if strs[n].replace('.','',1).isdigit() or 'e' in strs[n]:
                        # print(float(str))
                        if float(strs[n]) > 0.5:
                            bools.append(1)
                        else:
                            bools.append(0)
                tbl = tuple(bools)
                if tbl not in states:
                    states[tbl] = 1
                else:
                    states[tbl] += 1
                count += 1

    min_size = count * 0.05

    # print(states)

    attractors = []

    for key in states:
        if states[key] > min_size:
            attractors.append(key)
    print("attractors: ")
    print(attractors)

    return attractors

             

