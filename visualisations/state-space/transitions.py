# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import math
from tarjan_algo import *

G = nx.MultiDiGraph()

edges = [ ('19843', '19847'), ('19843', '19875'), ('19843', '19879'), ('19843', '19895'), ('19847', '19843'), ('19847', '19895'), ('19875', '19843'), ('19875', '19879'), ('19875', '19895'), ('19879', '19843'), ('19879', '19875'), ('19879', '19895'), ('19895', '19843'), ('19895', '19847'), ('19895', '19879'), ('30736', '31760'), ('31744', '31760'), ('31760', '30736'), ('31760', '31744'), ]
edge_widths = [ 0.893147, 1.58629, 0.893147, 4.51749, 0.2, 0.893147, 0.893147, 0.2, 0.2, 0.893147, 0.2, 0.893147, 4.53073, 0.2, 0.2, 1.29861, 1.29861, 1.29861, 0.893147, ]
nodes = [ '19843', '19847', '19875', '19879', '19895', '30736', '31744', '31760', '107177', ]
node_sizes = [ 794.2, 0.4, 0.9, 1.3, 86.1, 220.6, 77.9, 34.7, 1056.4, ]

print(len(edges))


# prune_nodes(nodes, node_sizes, edges, edge_widths)
# prune_edges(edges, edge_widths)




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


