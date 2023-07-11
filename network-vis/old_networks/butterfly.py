# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network


net = Network(notebook=True, directed=True)

net.add_node(1410, "1410", color="blue", size=25)
net.add_node(1922, "1922", color="blue", size=15)


net.add_node(2006, "2006", color="orange", size=15)
net.add_node(982, "982", color="orange", size=30)

net.add_node(4055, "4055", color="orange", size=5)

net.add_node(1923, "1923", color="red", size=5)

net.add_node(3971, "3971", color="orange", size=5)

net.add_node(3975, "3975", color="orange", size=5)


net.add_edge(1923, 3971)
net.add_edge(3971, 3975)
net.add_edge(3975, 4055)
net.add_edge(4055, 982)
net.add_edge(982, 2006)
net.add_edge(2006, 1923)


net.add_edge(1923, 1922)
net.add_edge(1922, 1410)
net.add_edge(1410, 1922)




# net.add_node(nodes)
# net.add_edges(edges)

# net.from_nx(G)
net.show("phenotypes.html")
