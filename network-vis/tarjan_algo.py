# Python program to find strongly connected components in a given
# directed graph using Tarjan's algorithm (single DFS)
#Complexity : O(V+E)

from collections import defaultdict
import sys

# This class represents an directed graph
# using adjacency list representation


class Graph:

	def __init__(self, vertices):
		# No. of vertices
		self.V = vertices

		# default dictionary to store graph
		self.graph = defaultdict(list)

		self.Time = 0
		self.components = []

	# function to add an edge to graph
	def addEdge(self, u, v):
		self.graph[u].append(v)

	'''A recursive function that find finds and prints strongly connected
	components using DFS traversal
	u --> The vertex to be visited next
	disc[] --> Stores discovery times of visited vertices
	low[] -- >> earliest visited vertex (the vertex with minimum
				discovery time) that can be reached from subtree
				rooted with current vertex
	st -- >> To store all the connected ancestors (could be part
		of SCC)
	stackMember[] --> bit/index array for faster check whether
				a node is in stack
	'''

	def SCCUtil(self, u, low, disc, stackMember, st):

		# Initialize discovery time and low value
		disc[u] = self.Time
		low[u] = self.Time
		self.Time += 1
		stackMember[u] = True
		st.append(u)

		# Go through all vertices adjacent to this
		for v in self.graph[u]:

			# If v is not visited yet, then recur for it
			if disc[v] == -1:

				self.SCCUtil(v, low, disc, stackMember, st)

				# Check if the subtree rooted with v has a connection to
				# one of the ancestors of u
				# Case 1 (per above discussion on Disc and Low value)
				low[u] = min(low[u], low[v])

			elif stackMember[v] == True:

				'''Update low value of 'u' only if 'v' is still in stack
				(i.e. it's a back edge, not cross edge).
				Case 2 (per above discussion on Disc and Low value) '''
				low[u] = min(low[u], disc[v])

		# head node found, pop the stack and print an SCC
		w = -1 # To store stack extracted vertices
		comp = []
		if low[u] == disc[u]:
			while w != u:
				w = st.pop()
				# print(w, end=" ")
				stackMember[w] = False
				comp.append(w)
				
			# print()
			self.components.append(comp)
	# The function to do DFS traversal.
	# It uses recursive SCCUtil()

	def SCC(self):

		# Mark all the vertices as not visited
		# and Initialize parent and visited,
		# and ap(articulation point) arrays
		disc = [-1] * (self.V)
		low = [-1] * (self.V)
		stackMember = [False] * (self.V)
		st = []

		# Call the recursive helper function
		# to find articulation points
		# in DFS tree rooted with vertex 'i'
		for i in range(self.V):
			if disc[i] == -1:
				self.SCCUtil(i, low, disc, stackMember, st)

	# function to check if a path exists between two vertices. 
	def isReachable(self, s, d):
			# Mark all the vertices as not visited
			visited =[False]*(self.V)

			# Create a queue for BFS
			queue=[]

			# Mark the source node as visited and enqueue it
			queue.append(s)
			visited[s] = True

			while queue:

					#Dequeue a vertex from queue
					n = queue.pop(0)
						
					# If this adjacent node is the destination node,
					# then return true
					if n == d:
									return True

					#  Else, continue to do BFS
					for i in self.graph[n]:
							if visited[i] == False:
									queue.append(i)
									visited[i] = True
				# If BFS is complete without visited d
			return False


	def ReturnList(self):
		return self.components



def CellStates(nodes, edges, node_sizes):
	## ALGORITHM TO FIND COMPONENTS AND CHECK FOR STEM AND DIFFERENTIATED CELLS
	# we need to turn the nodes into numbers
	new_nodes = []
	for i in range(len(nodes)):
		new_nodes.append(i)

	# make graph and run function
	g1 = Graph(len(nodes))
	for i in range(len(edges)):
			j1 = nodes.index(edges[i][0])
			j2 = nodes.index(edges[i][1])
			g1.addEdge(j1, j2)
			
	g1.SCC()
	comps = g1.ReturnList()
	# print(comps)

	#turn the components into the actual values, not index
	for i in range(len(comps)):
			for j in range(len(comps[i])):
					comps[i][j] = nodes[comps[i][j]]


	## Now we prune the nodes so that only the major states remain (over 200)
	pruned_nodes = []
	for i in range(len(nodes)):
		if node_sizes[i] > 50:
			pruned_nodes.append(nodes[i])


	# prune the component list so we only get components with major nodes
	pruned_components = []
	for i in range(len(comps)):
		pruned_list = []
		for j in range(len(comps[i])):
			if pruned_nodes.count(comps[i][j]):
				pruned_list.append(comps[i][j])
		if len(pruned_list) > 0:
			pruned_components.append(pruned_list)

	# print(pruned_nodes)
	print(pruned_components)


	#turned the pruned nodes back into numbers
	big_nodenumbers = []
	for i in range(len(pruned_nodes)):
		big_nodenumbers.append(nodes.index(pruned_nodes[i]))



	# iteratore through major nodes and see whether they are stem (can reach all other nodes)
	full_diff_list = []
	is_stem = []
	for i in big_nodenumbers:
		check_list = []
		check = True
		for j in big_nodenumbers:
			if i != j:
				if g1.isReachable(i, j) == False:
					check = False
					check_list.append(False)
				else:
					check_list.append(True)
			else:
				check_list.append(True)
		is_stem.append(check)
		full_diff_list.append(check_list)
	# print(is_stem)
	# print(full_diff_list)




	# Need to ascertain whether components can access other components. 
	# This will let us know which have partial stemness and which don't.
	# if len(pruned_components > 2):
	# 	for i in range(len(big_nodenumbers)):
	# 		for j in range(len(big_nodenumbers)):
	# 			if i != j and check_list[i][j] == True:
	# 				# Need to check whether i and j are in different components. 
	# 				n1 = pruned_nodes[i]
	# 				n2 = pruned_nodes[j]
	# 				for x in range(pruned_components):
	# 					if pruned_components[i].count(n1) and pruned_components[x].count(n2) == 0:


	original_stdout = sys.stdout

	with open('differentiation.txt', 'w') as f:
		sys.stdout = f
		#Re align everything so that we can determine if each component is stem or not. 
		for i in range(len(pruned_components)):
			print("Component #", i, " contains ", len(pruned_components[i]), " big states.")
			stemness = is_stem[pruned_nodes.index(pruned_components[i][0])]
			if stemness:
				print("The states ", pruned_components[i], " are stem.")
			else:
				print("The states ", pruned_components[i], " are differentiated.")

	sys.stdout = original_stdout


# Need edges, nodes, and size of nodes for pruning. 
# edges = [ ('51327', '63539'), ('63539', '59395'), ('60931', '60419'), ('63615', '63539'), ('63491', '65027'), ('63547', '59395'), ('59455', '59451'), ('65027', '60931'), ('60419', '65027'), ('59447', '59519'), ('59451', '59443'), ('59451', '59455'), ('63491', '59451'), ('63551', '63491'), ('59519', '63551'), ('59395', '59447'), ('65027', '60419'), ('55423', '63607'), ('63491', '59511'), ('63491', '59443'), ('51327', '63607'), ('59519', '63491'), ('59519', '59455'), ('59447', '63491'), ('59395', '59455'), ('63615', '59395'), ('64515', '59395'), ('59511', '59519'), ('59519', '63607'), ('63539', '63491'), ('59395', '59511'), ('51327', '63615'), ('59455', '63491'), ('49277', '49279'), ('63615', '63491'), ('55423', '59395'), ('59395', '51327'), ('51327', '49279'), ('59443', '59519'), ('63615', '63607'), ('63491', '59519'), ('63615', '59511'), ('51327', '59395'), ('59395', '64515'), ('49279', '51327'), ('51327', '63491'), ('63607', '63491'), ('59455', '51327'), ('59519', '59443'), ('60419', '59519'), ('63491', '64515'), ('59447', '51327'), ('59455', '59519'), ('63607', '59395'), ('60419', '64515'), ('51327', '59511'), ('59519', '59451'), ('59443', '51327'), ('59519', '59447'), ('60419', '59395'), ('59519', '63547'), ('59451', '51327'), ('59447', '59395'), ('51327', '55423'), ('59511', '63491'), ('55423', '63491'), ('65027', '64515'), ('59519', '63539'), ('63547', '63491'), ('59395', '59443'), ('59519', '59395'), ('63615', '59519'), ('59395', '63491'), ('59395', '59451'), ('63491', '59395'), ('51327', '59519'), ('59443', '63491'), ('55423', '59519'), ('59511', '51327'), ('59451', '59519'), ('59519', '51327'), ('59519', '59511'), ('59395', '59519'), ('59511', '59395'), ('60419', '51327'), ('64515', '63491'), ('59443', '59395'), ('64515', '65027'), ('59451', '59395'), ('55423', '63615'), ('64515', '60419'), ('60419', '63491'), ('59451', '63491'), ('59455', '59395'), ('59519', '63615'), ]
# edge_widths = [ 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.29861, 0.2, 0.893147, 0.2, 0.893147, 5.41494, 1.29861, 0.2, 0.893147, 2.5979, 4.43411, 2.14591, 2.39722, 2.83906, 3.09037, 4.71086, 0.893147, 0.893147, 4.96217, 0.2, 5.19721, 2.5979, 1.58629, 0.2, 1.80944, 0.2, 0.2, 0.2, 3.83759, 1.29861, 2.5979, 0.2, 1.80944, 0.893147, 6.45575, 0.2, 0.2, 0.2, 1.29861, 0.2, 3.14444, 0.893147, 0.2, 4.4485, 3.24452, 0.893147, 0.2, 2.90805, 2.50259, 2.39722, 5.37048, 2.83906, 3.24452, 1.99176, 4.75388, 0.2, 6.51355, 1.80944, 6.77368, 6.09164, 1.29861, 0.2, 2.5979, 1.58629, 5.74126, 3.03321, 5.81677, 2.39722, 0.2, 5.8058, 1.29861, 5.73733, 1.99176, 1.80944, 4.58203, 2.39722, 2.76495, 0.893147, 3.72636, ]
# nodes = [ '60931', '59451', '63539', '63551', '60419', '59511', '59395', '63615', '49277', '59443', '59455', '55423', '65027', '63607', '63491', '59519', '59447', '51327', '64515', '63547', '49279', ]
# node_sizes = [ 0.1, 9.59, 5.6, 0.24, 38.23, 14.09, 863.53, 49.96, 1, 3.97, 1.56, 5.8, 5453.58, 4.68, 2154.76, 395.27, 0.16, 1155.88, 2790.6, 8.4, 50.13, ]







			


# print(nodes.index('91136'), "    ", nodes.index('91136'))

# u = nodes.index('91136')
# v = nodes.index('28162')

# if g1.isReachable(u, v):
#     print("There is a path from %d to %d" % (u,v))
# else :
#     print("There is no path from %d to %d" % (u,v))

# if g1.isReachable(v, u):
#     print("There is a path from %d to %d" % (v, u))
# else :
#     print("There is no path from %d to %d" % (v, u))

