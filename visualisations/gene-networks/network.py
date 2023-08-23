# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import random


def generate_color():
  return '#%06x' % random.randint(0, 0xFFFFFF)


full_network = [ [ 0, 0, -1, 0, 0, 2, 0, 0, 0 ], [ 2, 1, 1, 0, 1, 0, 0, 0, -1 ], [ 2, 2, 0, 0, 0, 0, 0, 2, 0 ], 
 [ 0, 2, 0, 0, 0, 0, -1, -1, 1 ], [ -1, 0, 0, 0, 1, 1, 0, 0, 0 ], 
 [ 0, 2, 2, 0, 1, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 0, 0, 0, 0 ], [ -1, 0, 0, 0, 0, 0, 0, 0, 0 ], [ 1, 0, 0, 0, 0, 0, -1, 0, -1 ] ]

# can turn morph 2, morph 3 off, TF 1, TF 2 off

# must be in descending order
# to_kill = [5,4,2,1]


# for x in to_kill:
#     for i in range(len(full_network)):
#         for j in range(len(full_network[i])):
#             if j == x:
#                 full_network[i].pop(j)
#                 break

# for x in to_kill:
#     for i in range(len(full_network)):
#         if i == x:
#             full_network.pop(i)

print(full_network)

# tally
n2 = 0
n1 = 0
zero = 0
p1 = 0
p2 = 0

for i in range(len(full_network)):
    for j in range(len(full_network[i])):
        if full_network[i][j] == -2:
            n2 += 1
        elif full_network[i][j] == -1:
            n1 += 1
        elif full_network[i][j] == 0:
            zero += 1
        elif full_network[i][j] == 1:
            p1 += 1
        else:
            p2 += 1

print(n2, " ", n1, " ", zero, " ", p1, " ", p2)

nodes = [0,1,2,3,4,5,6,7,8]



net = Network(notebook=True, directed=True)


colors = []

for i in range(9):
    c = generate_color()
    while (c in colors):
        c = generate_color()
    colors.append(c)


    net.add_node(i, label=str(i), color=c)



for i in range(len(full_network)):
    for j in range(len(full_network[i])):
        if full_network[i][j] != 0:
            net.add_edge(j, i)


net.repulsion()

net.show('output.html')







# G = nx.MultiDiGraph()
# G.add_nodes_from(nodes)



# exclusion_list = []

# color_list = ["red", "green", "black", "yellow", "purple", "orange"]

# print("length of network is: ", len(full_network))

# for i in range(len(full_network)):
#     if i == 0:
#         net.add_node(0, "morph 1", color="red", size=10)
#     elif i == 1:
#         net.add_node(1, "MF1", color="blue", size=10)
#     elif i == 2:
#         net.add_node(2, "MF2", color="blue", size=10)
#     elif i == 3:
#         net.add_node(3, "TF3", color="green", size=10)
#     elif i == 4:
#         net.add_node(4, "TF4", color="green", size=10)
    # elif i == 5:
    #     net.add_node(5, "TF5", color="green", size=10)
    # elif i == 6:
    #     net.add_node(6, "TF6", color="green", size=10)
    # elif i == 7:
    #     net.add_node(7, "22", color="green", size=10)
    # elif i == 8:
    #     net.add_node(8, "33", color="green", size=10)
    #     print("node added")

    # elif i < 10:
    #     net.add_node(i, "locks", color="purple", size=10)
    # elif i < 11:
    #     net.add_node(i, "morph", color="red", size=10)
    # elif i < 12:
    #     net.add_node(i, "TF", color="green", size=10)
    # elif i < 9:
    #     net.add_node(i, "apop", color="black", size=10)
    # elif i < 11:
    #     net.add_node(i, "length", color="blue", size=10)
    # elif i < 16:
    #     net.add_node(i, "locks", color="purple", size=10)
    # else:
    #     net.add_node(i, "keys", color="orange", size=10)

for i in range(len(full_network)):
    for j in range(len(full_network[i])):
        if full_network[i][j] != 0:
            net.add_edge(j, i)

# net.add_node(nodes)
# net.add_edges(edges)

# net.from_nx(G)
net.show("network.html")
