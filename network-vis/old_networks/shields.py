# import
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network


net = Network(notebook=True, directed=True)

net.add_node(3803, "3803", color="orange", size=30)
net.add_node(3675, "3675", color="orange", size=45)
net.add_node(3659, "3659", color="orange", size=10)
net.add_node(3587, "3587", color="orange", size=10)
net.add_node(3715, "3715", color="red", size=35)
net.add_node(3971, "3971", color="red", size=10)
net.add_node(3987, "3987", color="orange", size=10)
net.add_node(4059, "4059", color="orange", size=10)

net.add_node(3203, "3203", color="blue", size=10)
net.add_node(3459, "3459", color="blue", size=30)
net.add_node(3458, "3458", color="blue", size=10)
net.add_node(1410, "1410", color="blue", size=35)


net.add_edge(3803, 3675)
net.add_edge(3675, 3659)
net.add_edge(3659, 3587)
net.add_edge(3587, 3715)
net.add_edge(3715, 3971)
net.add_edge(3971, 3987)
net.add_edge(3987, 4059)
net.add_edge(4059, 3803)


net.add_edge(3715, 3203)
net.add_edge(3203, 3459)
net.add_edge(3459, 3458)
net.add_edge(3458, 3459)
net.add_edge(3458, 1410)


# net.add_node(nodes)
# net.add_edges(edges)

# net.from_nx(G)
net.show("phenotypes.html")
