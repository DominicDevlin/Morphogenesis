# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import math
from tarjan_algo import *

G = nx.MultiDiGraph()


edges = [ ('1054', '30'), ('14', '46'), ('14', '526'), ('14', '558'), ('30', '1054'), ('30', '62'), ('46', '554'), ('46', '558'), ('526', '558'), ('546', '610'), ('546', '802'), ('554', '546'), ('558', '554'), ('570', '554'), ('570', '574'), ('570', '634'), ('570', '638'), ('570', '826'), ('574', '570'), ('574', '638'), ('610', '614'), ('610', '866'), ('614', '622'), ('614', '742'), ('614', '750'), ('62', '574'), ('622', '750'), ('630', '638'), ('634', '638'), ('638', '766'), ('742', '614'), ('742', '622'), ('742', '750'), ('750', '622'), ('802', '866'), ('818', '882'), ('826', '570'), ('826', '818'), ('826', '882'), ('826', '890'), ('866', '610'), ('882', '630'), ('882', '886'), ('886', '630'), ('886', '638'), ('886', '894'), ('890', '882'), ('894', '638'), ]
edge_widths = [ 3.41888, 0.893147, 0.893147, 2.14591, 0.2, 3.37805, 0.2, 0.893147, 0.893147, 0.2, 3.6012, 3.37805, 2.5979, 2.39722, 1.29861, 1.80944, 0.893147, 1.29861, 3.41888, 1.29861, 4.22535, 0.2, 3.03321, 3.83759, 0.2, 3.41888, 4.46268, 0.2, 1.58629, 2.27944, 0.893147, 0.2, 2.27944, 4.38965, 4.00666, 0.2, 0.2, 0.2, 0.2, 0.2, 4.22535, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, ]
nodes = [ '14', '30', '46', '62', '526', '546', '554', '558', '570', '574', '610', '614', '622', '630', '634', '638', '742', '750', '766', '802', '818', '826', '866', '882', '886', '890', '894', '1054', ]
node_sizes = [ 39.8, 21.7, 0.4, 22.4, 0.3, 17, 42.3, 3.2, 133.1, 10.2, 41.3, 151.2, 444.9, 0.1, 2, 26.9, 953.4, 337.5, 225.4, 24.4, 0.1, 8.6, 104.1, 1.5, 0.3, 0.2, 0.1, 277.7, ]

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


