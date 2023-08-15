# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import math
from tarjan_algo import *

G = nx.MultiDiGraph()


edges = [ ('12024', '12280'), ('12160', '12024'), ('12160', '12280'), ('12160', '77696'), ('12168', '12160'), ('12200', '12160'), ('12200', '12280'), ('12280', '12024'), ('12280', '12160'), ('12280', '12168'), ('12280', '12200'), ('77696', '12160'), ]
edge_widths = [ 5.89036, 3.81092, 6.04644, 5.42575, 1.29861, 0.2, 0.2, 5.76834, 6.13489, 1.29861, 0.2, 5.19043, ]
nodes = [ '12024', '12160', '12168', '12200', '12280', '77696', ]
node_sizes = [ 384.3, 1521.2, 0.6, 0.3, 302.3, 2187.5, ]


print(len(edges))

# extra pruning algorithm. Only used for inductive organisms to get better look
for i in range(len(edge_widths) -1, -1, -1):
    if edge_widths[i] < 1:
        del edge_widths[i]
        del edges[i]





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



