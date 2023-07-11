# import
import networkx as nx
# import matplotlib.pyplot as plt
from pyvis.network import Network


net = Network(notebook=True, directed=True)

net.add_node(1923, "1923", color="red", size=30)

net.add_node(2003, "2003 up", color="green", size=15)
net.add_node(20031, "2003 down", color="green", size=15)

net.add_node(2007, "2007 up", color="green", size=20)
net.add_node(20071, "2007 down", color="green", size=20)

net.add_node(4055, "4055", color="green", size=30)


net.add_node(899, "899", color="orange", size=40)
net.add_node(903, "903", color="blue", size=10)
net.add_node(919, "919", color="blue", size=30)


net.add_edge(1923, 2003)
net.add_edge(20031, 1923)
net.add_edge(2003, 2007)
net.add_edge(20071, 20031)
net.add_edge(2007, 4055)
net.add_edge(4055, 20071)
net.add_edge(20071, 1923)


net.add_edge(1923, 899)
net.add_edge(899, 1923)
net.add_edge(899, 903)
net.add_edge(903, 919)



# net.add_node(nodes)
# net.add_edges(edges)

# net.from_nx(G)
net.show("phenotypes.html")
