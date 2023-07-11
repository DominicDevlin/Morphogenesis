# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network


net = Network(notebook=True, directed=True)

net.add_node(1923, "1923", color="orange", size=30)
net.add_node(3971, "3971", color="orange", size=10)
net.add_node(1411, "1411", color="green", size=30)


net.add_node(1995, "1995", color="blue", size=10)
net.add_node(2027, "2027", color="blue", size=10)
net.add_node(2031, "2031", color="green", size=20)
net.add_node(879, "879", color="green", size=30)

net.add_edge(1923, 3971)
net.add_edge(3971, 1923)
net.add_edge(1923, 1411)
net.add_edge(1411, 1923)
net.add_edge(1995, 1923)
net.add_edge(2027, 1995)
net.add_edge(2031, 2027)
net.add_edge(2031, 879)
net.add_edge(879, 2031)


# net.add_node(nodes)
# net.add_edges(edges)

# net.from_nx(G)
net.show("phenotypes.html")
