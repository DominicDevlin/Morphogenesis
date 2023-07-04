#ifndef _COPY_H_
#define _COPY_H_


#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <bits/stdc++.h>


using namespace std;


// A class that represents an directed graph
class Graph {
	int V; // No. of vertices
	list<int>* adj; // A dynamic array of adjacency lists

  vector<vector<int>> components{};

	vector<vector<int>> UnComponents{};

	vector<int> OneComp{};

	bool strongly_connected{};

	vector<vector<int>> pruned{};



	// A Recursive DFS based function used by SCC()
	void SCCUtil(int u, int disc[], int low[],
				stack<int>* st, bool stackMember[]);

public:
	Graph(int V); // Constructor
	void addEdge(int v, int w); // function to add an edge to graph

	void addUnEdge(int v, int w);


	void SCC(); // prints strongly connected components

	// turns states into graph, checks for SCC
	vector<vector<int>> CreateDiGraph(map<int, int>& nodes, map<int, int>& diffs, vector<int> start, vector<int> end, double scale = 1., bool cycling = false);

	void PrintComponents();
	
  ~Graph();

	// void connectedComponents();

	void DFSUtil(int v, bool visited[], vector<int> &group);

	int NumberOfconnectedComponents(vector<vector<int>> &disconGroups);

	// void connectedComponents();

	double rescale(map<int,int> &types, map<int,int> &new_types);

	bool StronglyConnected();

	map<int,int> CreateUnGraph(map<int, int>& nodes, map<int, int>& diffs, map<pair<int,int>,int>& tally, int n_orgs=1, bool cycling = false);


	inline vector<vector<int>> ReturnPrunedComponents()
	{
		return pruned;
	}


};

#endif