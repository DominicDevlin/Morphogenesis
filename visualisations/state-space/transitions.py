# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import math
from tarjan_algo import *

G = nx.MultiDiGraph()

edges = [ ('107011', '107139'), ('107011', '108163'), ('107011', '74371'), ('107075', '107011'), ('107075', '108035'), ('107075', '108099'), ('107075', '75395'), ('107091', '107219'), ('107139', '108163'), ('107139', '74371'), ('107139', '75395'), ('107155', '75411'), ('107203', '108163'), ('107203', '75395'), ('107219', '107155'), ('107219', '107203'), ('107219', '108179'), ('107219', '108243'), ('107219', '74371'), ('107219', '74387'), ('107219', '74451'), ('107219', '75395'), ('107219', '75475'), ('107227', '108243'), ('108035', '107011'), ('108035', '107075'), ('108035', '107139'), ('108035', '108099'), ('108035', '108163'), ('108035', '112131'), ('108035', '75267'), ('108035', '75331'), ('108035', '75347'), ('108035', '75395'), ('108051', '108035'), ('108051', '108059'), ('108051', '108115'), ('108051', '108163'), ('108051', '108179'), ('108059', '108035'), ('108059', '108051'), ('108059', '75323'), ('108083', '108051'), ('108091', '108051'), ('108099', '108035'), ('108099', '75267'), ('108099', '75329'), ('108099', '75331'), ('108099', '75347'), ('108099', '75457'), ('108099', '75459'), ('108099', '75473'), ('108115', '108035'), ('108115', '108051'), ('108115', '108163'), ('108115', '108179'), ('108115', '108243'), ('108115', '75347'), ('108123', '108051'), ('108123', '108059'), ('108123', '108115'), ('108123', '108243'), ('108123', '108251'), ('108123', '75355'), ('108163', '107139'), ('108163', '108035'), ('108163', '75267'), ('108163', '75395'), ('108179', '108035'), ('108179', '108051'), ('108179', '108163'), ('108179', '112259'), ('108179', '75395'), ('108243', '108051'), ('108243', '108163'), ('108243', '108179'), ('108243', '75395'), ('108243', '75411'), ('108243', '75475'), ('108251', '108243'), ('112131', '108035'), ('112131', '128515'), ('112259', '112131'), ('128515', '112131'), ('74323', '74451'), ('74331', '107091'), ('74331', '107219'), ('74331', '108123'), ('74331', '108243'), ('74331', '108251'), ('74331', '74451'), ('74331', '74459'), ('74331', '75355'), ('74331', '75483'), ('74371', '75395'), ('74387', '74371'), ('74387', '75395'), ('74451', '107219'), ('74451', '74371'), ('74451', '74387'), ('74451', '75395'), ('74451', '75411'), ('74451', '75475'), ('74459', '107219'), ('74459', '107227'), ('74459', '108243'), ('74459', '108251'), ('74459', '74451'), ('74459', '75475'), ('74459', '75483'), ('75267', '75275'), ('75267', '75299'), ('75267', '75323'), ('75267', '75393'), ('75267', '75409'), ('75275', '75353'), ('75281', '75345'), ('75283', '108051'), ('75283', '75409'), ('75291', '108051'), ('75291', '75385'), ('75299', '75315'), ('75313', '75321'), ('75313', '75385'), ('75315', '108083'), ('75315', '75313'), ('75315', '75323'), ('75321', '75313'), ('75321', '75323'), ('75321', '75385'), ('75323', '108083'), ('75323', '75291'), ('75323', '75315'), ('75329', '75345'), ('75329', '75473'), ('75331', '75267'), ('75331', '75283'), ('75331', '75347'), ('75331', '75473'), ('75345', '75353'), ('75345', '75385'), ('75347', '75283'), ('75347', '75345'), ('75347', '75355'), ('75347', '75473'), ('75347', '75475'), ('75353', '75385'), ('75355', '108115'), ('75355', '108123'), ('75355', '108243'), ('75355', '108251'), ('75355', '74323'), ('75355', '74331'), ('75355', '75347'), ('75355', '75353'), ('75355', '75385'), ('75355', '75387'), ('75355', '75475'), ('75355', '75483'), ('75377', '75385'), ('75385', '75321'), ('75385', '75387'), ('75385', '75513'), ('75387', '108091'), ('75387', '108123'), ('75387', '75355'), ('75387', '75385'), ('75393', '75401'), ('75393', '75409'), ('75393', '75417'), ('75393', '75425'), ('75393', '75441'), ('75393', '75449'), ('75393', '75473'), ('75393', '75481'), ('75393', '91777'), ('75393', '91809'), ('75395', '108163'), ('75395', '75267'), ('75395', '75291'), ('75395', '75353'), ('75395', '75393'), ('75395', '75409'), ('75395', '75411'), ('75395', '75417'), ('75395', '75425'), ('75395', '75441'), ('75395', '75443'), ('75395', '75449'), ('75395', '75457'), ('75395', '75459'), ('75395', '75473'), ('75395', '75475'), ('75395', '75483'), ('75395', '91777'), ('75395', '91779'), ('75401', '75449'), ('75409', '75313'), ('75409', '75345'), ('75409', '75377'), ('75409', '75385'), ('75409', '75441'), ('75409', '75449'), ('75409', '75473'), ('75409', '75481'), ('75409', '75513'), ('75411', '75355'), ('75411', '75395'), ('75411', '75409'), ('75411', '75441'), ('75411', '75443'), ('75411', '75473'), ('75411', '75475'), ('75411', '75483'), ('75417', '75449'), ('75417', '75481'), ('75417', '75513'), ('75425', '75441'), ('75425', '75449'), ('75441', '75385'), ('75441', '75449'), ('75441', '75513'), ('75443', '75315'), ('75443', '75449'), ('75449', '75321'), ('75449', '75385'), ('75449', '75513'), ('75457', '75409'), ('75457', '75473'), ('75457', '75481'), ('75459', '75409'), ('75459', '75473'), ('75473', '75281'), ('75473', '75345'), ('75473', '75385'), ('75473', '75409'), ('75473', '75481'), ('75475', '108243'), ('75475', '75347'), ('75475', '75355'), ('75475', '75395'), ('75475', '75409'), ('75475', '75411'), ('75475', '75481'), ('75481', '75353'), ('75481', '75385'), ('75481', '75513'), ('75483', '108243'), ('75483', '108251'), ('75483', '75353'), ('75483', '75355'), ('75483', '75475'), ('75513', '75385'), ('91777', '75393'), ('91777', '75425'), ('91777', '91825'), ('91779', '91809'), ('91809', '75441'), ('91825', '75441'), ('91825', '75449'), ('91825', '91833'), ('91833', '75449'), ('91833', '75513'), ]
edge_widths = [ 1.29861, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.29861, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.80944, 0.2, 0.893147, 0.2, 0.2, 0.893147, 0.2, 1.58629, 1.29861, 0.2, 1.99176, 3.69651, 3.6012, 0.2, 0.2, 0.2, 1.58629, 3.6012, 0.2, 0.2, 0.2, 0.893147, 0.2, 1.99176, 0.2, 0.893147, 0.2, 0.2, 0.2, 0.2, 0.893147, 0.2, 0.2, 0.2, 0.2, 1.99176, 2.68491, 0.893147, 0.893147, 0.2, 0.2, 2.5979, 1.99176, 3.14444, 1.29861, 0.893147, 0.2, 0.2, 3.41888, 0.2, 4.17029, 0.2, 0.893147, 3.5322, 0.2, 1.29861, 0.2, 1.99176, 3.6012, 0.893147, 0.2, 0.2, 2.83906, 3.09037, 2.83906, 0.2, 2.76495, 0.2, 0.2, 1.29861, 0.2, 0.2, 0.2, 0.893147, 3.33549, 1.80944, 1.58629, 1.58629, 0.2, 0.893147, 1.29861, 0.2, 0.2, 0.2, 0.2, 1.29861, 2.27944, 0.2, 1.58629, 0.893147, 1.80944, 0.2, 1.29861, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.29861, 0.893147, 0.2, 0.2, 0.893147, 0.2, 0.2, 2.50259, 1.29861, 3.86356, 0.2, 1.58629, 0.2, 3.86356, 0.2, 0.2, 0.2, 1.80944, 0.893147, 2.5979, 0.2, 0.2, 4.82497, 0.893147, 0.2, 0.2, 4.78497, 1.99176, 0.2, 3.14444, 2.27944, 3.19573, 3.14444, 2.76495, 0.893147, 1.99176, 0.893147, 0.2, 2.14591, 0.2, 0.2, 0.2, 4.54381, 0.893147, 1.58629, 0.2, 0.2, 0.2, 0.2, 0.2, 0.893147, 0.2, 0.2, 0.2, 0.2, 1.58629, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.29861, 2.27944, 1.29861, 1.99176, 1.80944, 0.2, 1.99176, 0.2, 0.2, 0.2, 0.893147, 0.2, 0.2, 0.2, 1.58629, 1.80944, 2.50259, 2.76495, 0.2, 3.66574, 1.29861, 0.2, 0.2, 0.2, 0.2, 4.43411, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.893147, 0.2, 0.2, 2.27944, 1.29861, 0.2, 0.2, 1.58629, 0.2, 1.80944, 0.2, 1.80944, 2.14591, 2.5979, 2.14591, 2.14591, 0.2, 0.2, 1.58629, 4.74329, 0.2, 0.2, 1.80944, 0.2, 0.893147, 0.893147, 0.2, 1.29861, 0.893147, 0.2, ]
nodes = [ '74323', '74331', '74371', '74387', '74451', '74459', '75267', '75275', '75281', '75283', '75291', '75299', '75313', '75315', '75321', '75323', '75329', '75331', '75345', '75347', '75353', '75355', '75377', '75385', '75387', '75393', '75395', '75401', '75409', '75411', '75417', '75425', '75441', '75443', '75449', '75457', '75459', '75473', '75475', '75481', '75483', '75513', '91777', '91779', '91809', '91825', '91833', '107011', '107075', '107091', '107139', '107155', '107203', '107219', '107227', '108035', '108051', '108059', '108083', '108091', '108099', '108115', '108123', '108163', '108179', '108243', '108251', '112131', '112259', '128515', ]
node_sizes = [ 0.1, 15.6, 0.8, 0.5, 2.1, 6.9, 0.9, 0.1, 0.1, 0.2, 0.2, 0.1, 0.3, 2.2, 1, 1.1, 0.2, 0.6, 1.5, 1.1, 4.5, 131.9, 0.1, 673.2, 119.7, 32, 107.4, 0.1, 8.4, 4, 1.9, 4.5, 7.9, 0.6, 18.6, 0.4, 0.3, 2.6, 4.1, 12.3, 5.4, 56.5, 1.3, 0.1, 0.4, 1.8, 1.2, 1.2, 0.9, 0.1, 0.9, 0.2, 0.2, 3.9, 0.2, 174, 14.2, 1.5, 0.6, 0.1, 1.4, 6.8, 19.2, 55.7, 12.3, 13.7, 3, 1072.5, 0.2, 357.6, ]
print(len(edges))


# prune_nodes(nodes, node_sizes, edges, edge_widths)
prune_edges(edges, edge_widths, 3)




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


