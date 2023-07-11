# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network


net = Network(notebook=True, directed=True)

net.add_node(1, "1", color="orange", size=30)
net.add_node(2, "2", color="orange", size=45)
net.add_node(3, "3", color="orange", size=10)
net.add_node(4, "4", color="orange", size=10)
net.add_node(5, "5", color="red", size=35)
net.add_node(6, "6", color="red", size=10)
net.add_node(7, "7", color="orange", size=10)


net.add_edge(1, 2)
net.add_edge(2, 1)
net.add_edge(2, 3)
net.add_edge(3, 2)
net.add_edge(4, 5)
net.add_edge(5, 4)
net.add_edge(5, 6)
net.add_edge(6, 5)
net.add_edge(6, 7)



# net.add_node(nodes)
# net.add_edges(edges)

# net.from_nx(G)
net.show("phenotypes.html")
