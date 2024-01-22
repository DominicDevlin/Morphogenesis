# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import math
from tarjan_algo import *

G = nx.MultiDiGraph()


edges = [ ('1542', '1734'), ('1542', '1750'), ('1558', '1542'), ('1558', '1598'), ('1558', '1734'), ('1558', '1750'), ('1558', '3606'), ('1558', '3646'), ('1598', '1558'), ('1598', '3606'), ('1598', '3646'), ('1734', '1750'), ('1750', '1558'), ('1750', '1734'), ('2305', '2345'), ('2305', '2497'), ('2305', '2537'), ('2345', '2305'), ('2345', '2497'), ('2345', '2537'), ('2345', '3369'), ('2497', '2305'), ('2497', '2345'), ('2497', '2537'), ('2537', '2345'), ('2537', '2497'), ('3113', '3115'), ('3115', '3119'), ('3115', '3627'), ('3115', '3631'), ('3119', '3631'), ('3369', '3113'), ('3369', '3115'), ('3369', '3371'), ('3369', '3881'), ('3369', '3883'), ('3371', '3115'), ('3371', '3627'), ('3371', '3883'), ('3606', '1558'), ('3606', '3646'), ('3625', '3627'), ('3627', '3631'), ('3631', '3647'), ('3646', '1598'), ('3646', '3606'), ('3646', '3647'), ('3647', '3631'), ('3647', '3646'), ('3881', '3369'), ('3881', '3625'), ('3881', '3883'), ('3883', '3627'), ]
edge_widths = [ 3.14444, 0.2, 3.14444, 1.29861, 3.4581, 4.37439, 3.83759, 0.893147, 0.2, 0.2, 1.99176, 4.37439, 4.90048, 3.29104, 1.29861, 3.6012, 1.58629, 3.63399, 4.66591, 5.34166, 4.29434, 1.29861, 2.97259, 5.37615, 5.80947, 4.55671, 0.2, 1.58629, 3.49584, 0.2, 1.58629, 0.2, 1.29861, 3.93767, 2.5979, 1.80944, 3.5673, 2.27944, 1.99176, 1.58629, 3.86356, 0.2, 4.15124, 4.24305, 1.80944, 1.80944, 0.2, 0.893147, 4.11202, 0.2, 0.2, 2.50259, 3.19573, ]
nodes = [ '1542', '1558', '1598', '1734', '1750', '2305', '2345', '2497', '2537', '3113', '3115', '3119', '3369', '3371', '3606', '3625', '3627', '3631', '3646', '3647', '3881', '3883', ]
node_sizes = [ 4.5, 916.6, 4.7, 24.4, 112.7, 7.6, 1587.8, 74.2, 234.2, 0.2, 20.1, 1.4, 66.4, 12.9, 112.6, 0.1, 48.1, 71.3, 2610.5, 147.5, 4.8, 6.5, ]

print(len(edges))


# prune_nodes(nodes, node_sizes, edges, edge_widths)
# prune_edges(edges, edge_widths, 1)




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


