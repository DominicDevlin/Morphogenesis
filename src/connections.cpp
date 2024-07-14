#include "connections.h"
#include "parameter.h"
#include <math.h>
#include <stdio.h>
#include <bits/stdc++.h>

extern Parameter par;

using namespace std;

#define NIL -1

Graph::Graph(int V)
{
	this->V = V;
	adj = new list<int>[V];
}

Graph::~Graph()
{
	delete[] adj;
	components.clear();
	UnComponents.clear();
	OneComp.clear();
}

void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w);
}

void Graph::addUnEdge(int v, int w)
{
	adj[v].push_back(w);
	adj[w].push_back(v);
}

// A recursive function that finds and prints strongly
// connected components using DFS traversal u --> The vertex
// to be visited next disc[] --> Stores discovery times of
// visited vertices low[] -- >> earliest visited vertex (the
// vertex with minimum
//			 discovery time) that can be reached from
//			 subtree rooted with current vertex
// *st -- >> To store all the connected ancestors (could be
// part
//		 of SCC)
// stackMember[] --> bit/index array for faster check
// whether
//				 a node is in stack
void Graph::SCCUtil(int u, int disc[], int low[], stack<int> *st, bool stackMember[])
{
	// A static variable is used for simplicity, we can
	// avoid use of static variable by passing a pointer.
	static int time = 0;

	// Initialize discovery time and low value
	disc[u] = low[u] = ++time;
	st->push(u);
	stackMember[u] = true;

	// Go through all vertices adjacent to this
	list<int>::iterator i;
	for (i = adj[u].begin(); i != adj[u].end(); ++i)
	{
		int v = *i; // v is current adjacent of 'u'

		// If v is not visited yet, then recur for it
		if (disc[v] == -1)
		{
			SCCUtil(v, disc, low, st, stackMember);

			// Check if the subtree rooted with 'v' has a
			// connection to one of the ancestors of 'u'
			// Case 1 (per above discussion on Disc and Low
			// value)
			low[u] = min(low[u], low[v]);
		}

		// Update low value of 'u' only of 'v' is still in
		// stack (i.e. it's a back edge, not cross edge).
		// Case 2 (per above discussion on Disc and Low
		// value)
		else if (stackMember[v] == true)
			low[u] = min(low[u], disc[v]);
	}

	// head node found, pop the stack and print an SCC
	int w = 0; // To store stack extracted vertices
	vector<int> new_comp;
	if (low[u] == disc[u])
	{
		while (st->top() != u)
		{
			w = (int)st->top();
			// cout << w << " ";
			new_comp.push_back(w);
			stackMember[w] = false;
			st->pop();
		}
		w = (int)st->top();
		// cout << w << "\n";
		new_comp.push_back(w);
		components.push_back(new_comp);
		new_comp.clear();
		stackMember[w] = false;
		st->pop();
	}
}

// The function to do DFS traversal. It uses SCCUtil()
void Graph::SCC()
{
	int *disc = new int[V];
	int *low = new int[V];
	bool *stackMember = new bool[V];
	stack<int> *st = new stack<int>();

	// Initialize disc and low, and stackMember arrays
	for (int i = 0; i < V; i++)
	{
		disc[i] = NIL;
		low[i] = NIL;
		stackMember[i] = false;
	}

	// Call the recursive helper function to find strongly
	// connected components in DFS tree with vertex 'i'
	for (int i = 0; i < V; i++)
		if (disc[i] == NIL)
			SCCUtil(i, disc, low, st, stackMember);
}

// Method to print connected components in an
// undirected graph. Needed to test whether there are stems.
// void Graph::connectedComponents()
// {
//     // Mark all the vertices as not visited
//     bool* visited = new bool[V];
//     for (int v = 0; v < V; v++)
//         visited[v] = false;

//     for (int v = 0; v < V; v++) {
//         if (visited[v] == false) {
//             // print all reachable vertices
//             // from v
//             DFSUtil(v, visited);

// 						UnComponents.push_back(OneComp);
// 						OneComp.clear();
//             cout << "\n";
//         }
//     }
//     delete[] visited;
// }

// void Graph::DFSUtil(int v, bool visited[])
// {
//     // Mark the current node as visited and print it
//     visited[v] = true;
//     cout << v << " ";
// 		OneComp.push_back(v);

//     // Recur for all the vertices
//     // adjacent to this vertex
//     list<int>::iterator i;
//     for (i = adj[v].begin(); i != adj[v].end(); ++i)
//         if (!visited[*i])
//             DFSUtil(*i, visited);
// }

// Function to return the number of
// connected components in an undirected graph
int Graph::NumberOfconnectedComponents(vector<vector<int>> &disconGroups)
{

	// Mark all the vertices as not visited
	bool *visited = new bool[V];

	// To store the number of connected components
	int count = 0;
	for (int v = 0; v < V; v++)
		visited[v] = false;

	for (int v = 0; v < V; v++)
	{
		if (visited[v] == false)
		{
			vector<int> group;
			DFSUtil(v, visited, group);
			count += 1;
			disconGroups.push_back(group);
		}
	}

	return count;
}

void Graph::DFSUtil(int v, bool visited[], vector<int> &group)
{

	// Mark the current node as visited
	if (find(group.begin(), group.end(), v) == group.end())
	{
		group.push_back(v);
	}

	visited[v] = true;
	vector<vector<int>> disconGroups;

	// Recur for all the vertices
	// adjacent to this vertex
	list<int>::iterator i;

	for (i = adj[v].begin(); i != adj[v].end(); ++i)
		if (!visited[*i])
			DFSUtil(*i, visited, group);
}

// Depracated
int CalculateThreshold(vector<int> &sizes, double scale = 1.)
{
	double threshold{};
	for (int i : sizes)
	{
		threshold += i;
	}
	cout << "threshold before scaling: " << threshold << endl;

	threshold = threshold * par.node_percent * scale;
	threshold = floor(threshold);

	// new threshold is going to be flat number
	// 10 cells * (12100 - 8000)/40 *2 = 2040

	return (int)(threshold);
}

vector<vector<int>> Graph::CreateDiGraph(map<int, int> &nodes, map<int, int> &types, vector<int> start, vector<int> end, double scale, bool cycling)
{
	vector<int> state_sizes{};
	vector<int> state_values{};

	// state values to numbers.
	if (par.potency_edges)
	{
		for (auto i : types)
		{
			state_sizes.push_back(types[i.first]);
			state_values.push_back(i.first);
		}
	}
	else
	{
		for (auto i : nodes)
		{
			if (types.find(i.first) == types.end())
			{
				state_sizes.push_back(1);
			}
			else
			{
				state_sizes.push_back(types[i.first]);
			}
			state_values.push_back(i.first);
		}
	}

	// int threshold = CalculateThreshold(state_sizes, scale);
	// cout << "THRESHOLD IS: " << threshold << endl;

	// replace state values with the new set of numbers.

	for (int i = 0; i < start.size(); ++i)
	{
		auto f1 = find(state_values.begin(), state_values.end(), start[i]);
		auto f2 = find(state_values.begin(), state_values.end(), end[i]);
		if (f1 != state_values.end() && f2 != state_values.end())
		{
			int ind1 = f1 - state_values.begin();
			int ind2 = f2 - state_values.begin();
			this->addEdge(ind1, ind2);
			// cout << "Edge: " << ind1 << "  " << ind2 << endl;
			// cout << "Working Edge: " << start[i] << "  " << end[i] << endl;
		}
		else
		{
			cerr << "Error creating graph!\n";
			// cout << "Not Working Edge: " << start[i] << "  " << end[i] << endl;
		}
	}
	if (pruned.size())
		cout << "error in function call" << endl;

	// create SCC. This creates all the components stored in "components"
	this->SCC();

	// prune SCC to exclude transients
	for (int i = 0; i < components.size(); ++i)
	{
		vector<int> new_list{};
		for (int j = 0; j < components[i].size(); ++j)
		{
			// check if node is big enough i.e. not a transient
			if (cycling && state_sizes[components[i][j]] > par.node_threshold * 3)
			{
				new_list.push_back(components[i][j]);
			}
			else if (state_sizes[components[i][j]] > par.node_threshold)
			{
				new_list.push_back(components[i][j]);
			}
		}
		if (new_list.size() > 0)
		{
			pruned.push_back(new_list);
		}
	}

	vector<vector<int>> SCC_types{};

	if (par.print_fitness && pruned.size() > 0)
	{
		// cout << "new component...." << endl;
		for (int i = 0; i < pruned.size(); ++i)
		{
			vector<int> new_scc{};
			for (int j = 0; j < pruned[i].size(); ++j)
			{
				int el = pruned[i][j];
				auto it = types.begin();
				std::advance(it, el);
				// cout << it->first << "  ";
				new_scc.push_back(it->first);
			}
			// cout << endl;
			SCC_types.push_back(new_scc);
		}
	}
	if (pruned.size() == 1)
	{
		strongly_connected = true;
	}
	else
	{
		strongly_connected = false;
	}

	return SCC_types;
}

bool Graph::StronglyConnected()
{
	return strongly_connected;
}

double Graph::rescale(map<int, int> &types, map<int, int> &new_types)
{
	// calculate total time
	int ttime = 0;
	for (auto &i : types)
	{
		ttime += i.second;
	}
	int subtime = 0;
	for (auto &i : new_types)
	{
		subtime += i.second;
	}
	double ratio = double(ttime) / double(subtime);
	return ratio;
}

// removes bidirectional edges for undirected graph. Function is not needed.
void PruneBiEdges(vector<int> &start, vector<int> &end)
{
	// get a pair
	for (int i = 0; i < start.size(); ++i)
	{
		int ender = start[i];
		int starter = end[i];

		for (int j = i + 1; j < start.size(); ++j)
		{
			if (starter == start[j] && ender == end[j])
			{
				auto it1 = start.begin() + j;
				auto it2 = end.begin() + j;

				*it1 = move(start.back());
				*it2 = move(end.back());
				start.pop_back();
				end.pop_back();
			}
		}
	}
}

void PruneEdges(map<pair<int, int>, int> &tally, int n_orgs)
{
	for (auto it = tally.begin(); it != tally.end();)
	{
		if (it->second < n_orgs * par.prune_amount)
		{
			// cout << "Erased edge: " << it->first.first << " " << it->first.second << " " << it->second << endl;
			it = tally.erase(it);
		}
		else
		{
			++it;
		}
	}
}

map<int, int> Graph::CreateUnGraph(map<int, int> &nodes, map<int, int> &types, map<pair<int, int>, int> &tally, int n_orgs, bool cycling)
{
	vector<int> state_sizes{};
	vector<int> state_values{};

	if (par.prune_edges)
	{
		PruneEdges(tally, n_orgs);
	}

	// cout << "Node list: " << endl;
	// for (auto &i : types)
	// {
	// 	cout << i.first << endl;
	// }

	// state values to numbers.
	if (par.potency_edges)
	{
		for (auto i : types)
		{
			state_sizes.push_back(types[i.first]);
			state_values.push_back(i.first);
			// cout << i.first << " " << i.second << endl;
		}
	}
	else
	{
		for (auto i : nodes)
		{
			if (types.find(i.first) == types.end())
			{
				state_sizes.push_back(0);
			}
			else
			{
				state_sizes.push_back(types[i.first]);
			}
			state_values.push_back(i.first);
		}
	}

	// replace state values with the new set of numbers.
	for (auto &i : tally)
	{
		auto f1 = find(state_values.begin(), state_values.end(), i.first.first);
		auto f2 = find(state_values.begin(), state_values.end(), i.first.second);
		if (f1 != state_values.end() && f2 != state_values.end())
		{
			int ind1 = f1 - state_values.begin();
			int ind2 = f2 - state_values.begin();
			this->addUnEdge(ind1, ind2);
			// cout << "UNCON Edge: " << ind1 << "  " << ind2 << endl;
		}
		else
		{
			cerr << "Error creating graph!\n";
			cout << "Missing edge is: " << i.first.first << "  " << i.first.second << "  " << i.second << endl;
		}
	}

	vector<vector<int>> disconGroups{};

	int val = this->NumberOfconnectedComponents(disconGroups);

	if (par.print_fitness)
	{
		cout << "Number of connected components: " << val << endl;
	}

	// now we find if there are weakly connected components within each group:

	// we need to create subsets of nodes and types for each group of discongroup.
	// Then have to find the edges that are specifically involved with those nodes...
	map<int, int> toreturn{};

	int counter = 0;
	for (vector<int> i : disconGroups)
	{
		map<int, int> new_nodes{};
		map<int, int> new_types{};
		for (int j : i)
		{
			int c = 0;
			for (auto it = types.begin(); it != types.end(); ++it)
			{
				if (c == j)
				{
					new_types[it->first] = it->second;
					// cout << "printing addition: " << it->first << "   " << new_types[it->first] << endl;
				}
				++c;
			}
		}
		double resc = rescale(types, new_types);

		vector<int> new_start{};
		vector<int> new_end{};
		for (auto key : new_types)
		{
			for (auto &i : tally)
			{
				if (key.first == i.first.first)
				{
					new_start.push_back(i.first.first);
					new_end.push_back(i.first.second);
				}
			}
		}
		Graph subgraph(new_types.size());
		vector<vector<int>> subcomps = subgraph.CreateDiGraph(new_nodes, new_types, new_start, new_end, resc, cycling);
		for (auto i : subcomps)
		{
			pruned.push_back(i);
		}

		// for (auto i : subcomps)
		// {
		// 	for (auto j : i)
		// 	{
		// 		cout << j << "  ";
		// 	}
		// 	cout << endl;
		// }

		// cout << "subcount done with size " << subcomps.size() << endl;

		if (subcomps.size() > 0)
			toreturn[counter] = subcomps.size();

		++counter;
	}

	// we make a map with each component, and value number of weakly connected components (0 = SC, 1 + = weakly connected)
	return toreturn;
}

vector<vector<int>> Graph::GetComps(map<int, int> types, int threshold)
{

	vector<vector<int>> comps;
	for (auto i : pruned)
	{
		int sumt{};
		for (auto j : i)
		{
			sumt += types[j];
		}
		if (sumt > threshold)
		{
			comps.push_back(i);
		}
	}
	return comps;
}

void Graph::PrintComponents()
{
	for (auto i : components)
	{
		for (int j : i)
		{
			cout << j << "  ";
		}
		cout << std::endl;
	}
}
