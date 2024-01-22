import math
from find_CPM import *

## weird result... evolved networks seem to have INCREASED CONNECTIVITY (all of them). This might be
## disconnected graphs are so rare 

str_m = '{ { 2, 0, 0, -1, 0, 0, 1, 2, 1 }, { 1, 1, 0, -1, 1, 0, 0, 0, 0 }, { -2, 0, -1, 0, 0, -1, 0, 1, 0 }, { 1, 0, 2, 1, 0, 1, 0, 0, -2 }, { -2, 0, 1, 0, 2, -1, -1, -1, 2 }, { 0, 0, 0, 0, 0, 1, 0, 0, 0 }, { -2, 0, -1, 0, 2, 0, 1, 0, -1 }, { 1, 0, -2, 0, -1, -2, 0, 0, 0 }, { -1, 0, 1, 0, -1, 1, 0, 0, 0 }, { 1, 0, 1, 0, 0, 1, -1, 2, 1 }, { 0, 0, 0, 0, -1, 1, 0, 0, -1 }, { 0, 2, 1, 0, -1, -2, 1, 0, 2 }, { 0, 0, 0, -1, -1, 1, 0, 0, -1 }, { -2, 0, 1, -2, 0, 1, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0, 0, 0, 0 }, { -1, 0, -1, 0, 0, 0, -1, 0, 0 }, { 0, -1, 0, 0, 1, -1, 0, 0, 0 }, { 0, -1, 0, 0, 0, 1, 0, -2, 2 }, { 0, 1, 0, 1, 0, 0, 0, 0, -1 }, { 1, -1, 0, 0, 0, 0, 1, -1, 0 }, { 0, -1, 0, 0, 0, 1, 0, 0, 1 }, { 0, 0, 0, 1, 0, 0, 0, 0, -1 }, { 0, 0, -2, 0, 1, -2, 0, 0, 0 }, { 1, 0, -1, 0, -2, -1, 1, -1, 2 }, { 0, 0, 1, 0, 0, 0, -2, 0, 1 }, { 1, 0, 0, -1, 0, 0, 0, -1, 2 }, { 0, 0, 0, 0, 0, 0, 1, 0, 1 }, }'




matrix = []


n_genes = 9
gene_list = [0.,0.,0.,1.,1.,1.,1.,1.,1.]

n_inits = int(math.pow(2, n_genes))
list_of_inits = []

bools = list
update_steps = 400

theta = - 0.3

edges = []

d_rate = math.exp(-0.29)

special_edges = []




### FUNCTIONS START

# update the network
def update():
    for i in range(n_genes):
      x_1 = 0
      for j in range(n_genes):
          x_1 += (matrix[i][j] * gene_list[j])

      x_1 += theta
      gene_list[i] = (1 / (1 + math.exp(-20 * x_1))) * 0.25 + gene_list[i] * d_rate



# turn the expression pattern into boolean
def booleanise():
    new_lst = []
    for value in gene_list:
        if value < 0.5:
            new_lst.append(0)
        else:
            new_lst.append(1)
    global bools
    bools = new_lst
    
def check_if_transition():
    for value in gene_list:
       if value > 0.35 and value < 0.65:
           return True
    return False


def print_bools():
    print(bools)

def print_glist():
    print(gene_list)



# make sure this is called after booleanise
def bool_list_to_int():
    result = 0
    for index, val in enumerate(bools):
        if val:
            result |= 1 << index
    return result


# reduce the size of the transition matrix by pruning outside nodes
def prune():
  to_remove = []
  for i in range(n_inits):
      count = 0
      has_inc = False
      for j in edges:
          if j[0] == i:
              count += 1
          elif j[1] == i:
              has_inc = True
              break
      if has_inc:
          continue
      if count > 0:
        for x in range(len(edges)):
            if edges[x][0] == i or edges[x][1] == i:
                to_remove.append(x)

  to_remove.sort(reverse=True)
  # remove pruned nodes
  for i in to_remove:
      del edges[i]  


# turn string into list of lists
def convert_matrix():
    sublist = []
    for i in range(len(str_m)):
        if str_m[i].isnumeric():
            v = int(str_m[i])
            if str_m[i-1] == '-':
                sublist.append(-v)
            else:
                sublist.append(v)
        elif str_m[i] == '}' and len(sublist) > 0:
            matrix.append(sublist)
            sublist = []


def find_attractor(states):
    
    subs = []
    i = int(update_steps * 0.8)
    while (i < update_steps):
        if states[i] not in subs:
            subs.append(states[i])
        i+=1
    global special_edges
    if len(subs) > 1:
        i = int(update_steps*0.8)
        last_s = states[i]
        i +=1
        while (i < update_steps):
            curr_s = states[i]
            if last_s != curr_s:
                tup = (last_s, curr_s)
                if tup not in special_edges:
                    special_edges.append(tup)
            i += 1
    
    subs.sort()
    return subs


 

#### FUNCTIONS END


        
convert_matrix()

for i in range(n_inits):
    binary_str = bin(i)[2:].zfill(9)
    int_bit = [int(bit) for bit in binary_str]
    list_of_inits.append(int_bit)


attractors = []

# do first loop to find attractors

bool_list = []

att_count = {}

for i in range(n_inits):
    gene_list = list(list_of_inits[i])
    list_of_states = []
    for j in range(update_steps):  
      update()
      booleanise()
      current_state = bool_list_to_int()
      list_of_states.append(current_state)
    if bools not in bool_list:
        bool_list.append(bools)
    tpbool = tuple(bools)
    if tpbool in att_count:
        att_count[tpbool] += 1
    else:
        att_count[tpbool] = 1
    attractors.append(find_attractor(list_of_states))



tupatt = []
[tupatt.append(tuple(x)) for x in attractors]

# remove minor attractors from attractor count (using this to compare against actual CPM attractors)
major_attractors = []

for key in att_count:
    if att_count[key] > 10:
        major_attractors.append(key)




## this is for numerised version of attractions
att = {}
for x in tupatt:
    if x not in att:
        att[x] = 1
    else:
        att[x] += 1
        
print("attractors are: ", att)





## NEXT LOOP TO GET EDGES
for i in range(n_inits):
    gene_list = list(list_of_inits[i])
    booleanise()
    last_state = bool_list_to_int()
    # list_of_states = []
    for j in range(update_steps):
        update()
      
        if check_if_transition():
            continue
        else:
            booleanise()
            current_state = bool_list_to_int()
            # list_of_states.append(current_state)
            if (last_state != current_state and i > 0):
                t = (last_state, current_state)
                if t not in edges:
                    edges.append(t)
            last_state = current_state
    # attractors.append(find_attractor(list_of_states))


for i in special_edges:
    if i not in edges:
        edges.append(i)



## Also want to find how commonly they appear in each. 

print("number of edges: ", len(edges))
prune()
print("first prune done, n edges: ", len(edges))
# prune()
# print("second proune done, n edges: ", len(edges))
# prune()
# print("third prune done")
# prune()
# prune()
# prune()
# prune()    

print(edges)

print(bool_list)

cpm_att = CPM_attractions()

both = 0
cpm_count = len(cpm_att)
single_count = len(major_attractors)

for i in major_attractors:
    for j in cpm_att:
        if i == j:
            both += 1

print("single cell attractors: ")
print(major_attractors)


print("CPM count: ", cpm_count)
print("single count: ", single_count)
print("both count: ", both)
