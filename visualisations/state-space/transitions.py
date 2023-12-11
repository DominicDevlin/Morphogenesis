# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import math
from tarjan_algo import *

G = nx.MultiDiGraph()


edges = [ ('123010', '124150'), ('123010', '90242'), ('124150', '124158'), ('124150', '91382'), ('124150', '91390'), ('124158', '91390'), ('90114', '90118'), ('90114', '90242'), ('90118', '90114'), ('90150', '90118'), ('90242', '123010'), ('90242', '124150'), ('90242', '90114'), ('90242', '91382'), ('90246', '90118'), ('90246', '90242'), ('90278', '90118'), ('90278', '90150'), ('90278', '90246'), ('91270', '90242'), ('91302', '90278'), ('91366', '90242'), ('91366', '91302'), ('91382', '90242'), ('91382', '90246'), ('91382', '91270'), ('91382', '91366'), ('91382', '91390'), ('91390', '124158'), ('91390', '91382'), ]
edge_widths = [ 5.29987, 1.29861, 5.23044, 1.80944, 2.83906, 5.24986, 1.29861, 4.79512, 2.39722, 0.893147, 5.00402, 0.893147, 5.16981, 0.893147, 0.2, 0.893147, 0.893147, 0.893147, 0.2, 0.2, 1.58629, 0.2, 1.29861, 5.12725, 0.893147, 0.2, 1.80944, 2.14591, 1.58629, 5.23695, ]
nodes = [ '90114', '90118', '90150', '90242', '90246', '90278', '91270', '91302', '91366', '91382', '91390', '123010', '124150', '124158', ]
node_sizes = [ 1074.7, 76.5, 1.6, 529.2, 0.6, 6.5, 0.2, 11.7, 45.3, 166, 1546.9, 87.3, 224.4, 115.3, ]

print(len(edges))


prune_nodes(nodes, node_sizes, edges, edge_widths)
prune_edges(edges, edge_widths, 0.5)




G.add_nodes_from(nodes)

for i in range(len(edges)):
    G.add_edge(edges[i][0], edges[i][1], width = edge_widths[i])


node_d = {}
for i in range(len(nodes)):
    val = math.sqrt(node_sizes[i])
    val = val * 2
    node_d[nodes[i]] = val

nx.set_node_attributes(G, node_d, 'size')

net = Network(directed=True)
net.repulsion()
net.from_nx(G)
net.show('transitions.html')


CellStates(nodes, edges, node_sizes)



# pos = nx.nx_agraph.graphviz_layout(G)
# nx.draw_networkx_nodes(G, pos, node_size = node_sizes)
# nx.draw_networkx_edges(G, pos, edgelist=G.edges(), nodelist=node_sizes) #, width=edge_widths, connectionstyle='arc3, rad = 0.1') #, nodelist=node_sizes)
# nx.draw_networkx_labels(G,pos)
# plt.show()


