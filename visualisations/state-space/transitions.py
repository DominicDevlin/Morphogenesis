# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import math
from tarjan_algo import *

G = nx.MultiDiGraph()


edges = [ ('1011', '1010'), ('1011', '883'), ('1011', '995'), ('1363', '1491'), ('1363', '339'), ('1363', '467'), ('1363', '499'), ('1491', '1363'), ('1491', '339'), ('1491', '467'), ('1491', '499'), ('339', '1363'), ('339', '371'), ('339', '467'), ('339', '499'), ('371', '1011'), ('371', '339'), ('371', '499'), ('371', '883'), ('467', '339'), ('467', '371'), ('467', '499'), ('499', '1011'), ('499', '371'), ('499', '467'), ('499', '883'), ('742', '998'), ('867', '995'), ('883', '1011'), ('883', '867'), ('883', '995'), ('994', '998'), ('995', '867'), ('995', '994'), ('995', '998'), ('995', '999'), ('998', '742'), ('999', '995'), ]
edge_widths = [ 2.5979, 0.2, 2.83906, 4.7326, 0.2, 1.29861, 0.2, 4.67734, 0.2, 2.76495, 0.893147, 0.2, 0.2, 1.58629, 0.2, 0.893147, 0.2, 0.893147, 0.893147, 1.58629, 0.2, 2.90805, 3.03321, 1.80944, 0.2, 0.2, 3.09037, 4.27754, 0.893147, 0.2, 0.2, 2.5979, 4.29434, 2.5979, 0.2, 3.09037, 3.33549, 0.2, ]
nodes = [ '339', '371', '467', '499', '742', '867', '883', '994', '995', '998', '999', '1010', '1011', '1363', '1491', ]
node_sizes = [ 3.9, 2.8, 8.4, 16, 239.4, 28.9, 1.5, 24.3, 212.4, 208.1, 2412.2, 363.5, 88, 87.1, 102.3, ]

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


