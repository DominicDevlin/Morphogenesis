#include "storage.h"
#include "parameter.h"
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>



extern Parameter par;

using namespace std;



storage::~storage()
{
  switch_tally.clear();
  phenotype_tally.clear();
  edge_tally.clear();
}



void storage::add_to_switches(unordered_map<string, int> switches)
{
  for (auto i : switches)
  {
    switch_tally[i.first] += i.second;
  }
}

void storage::add_to_time(map<int, int>& phen, map<int,int>& adult)
{
  for (auto i : phen)
  {
    phenotype_tally[i.first] += i.second;
  }

  for (auto i : adult)
  {
    adult_tally[i.first] += i.second;
  }
}

void storage::do_averaging()
{
  for (auto i : switch_tally)
  {
    int val = ceil((double)i.second / (double)par.n_orgs);
    switch_tally[i.first] = val;
  }

  for (auto i : adult_tally)
  {
    int val = ceil((double)i.second / (double)par.n_orgs);
    adult_tally[i.first] = val;
  }


  for (auto i : phenotype_tally)
  {
    if (!par.potency_edges)
    {
      if (adult_tally.find(i.first) == adult_tally.end())
      {
        adult_tally[i.first] = 10;
      }
      int val = ceil((double)i.second / (double)par.n_orgs);
      phenotype_tally[i.first] = val;
    }

  }


}

map<string, int>& storage::get_switch_tally()
{
  return switch_tally;
}

map<int, int>& storage::get_phenotype_tally()
{
  return phenotype_tally;
}

map<int, int>& storage::get_adult_tally()
{
  return adult_tally;
}



void storage::add_to_edges(map<pair<int,int>,int>& tally)
{

  for (auto &i : tally)
  {
    edge_tally[i.first] += i.second;
  }


}

map<pair<int,int>,int>& storage::get_edges()
{
  return edge_tally;
}



void storage::write_to_file(bool cycles)
{
/* Files for creating graph in python
  Need: edges, edge width, nodes, node size. */

  if (mkdir("transition-data", 0777) == -1)
    cerr << "Error : transitions directory already created." << endl;
  else
    cout << "Directory created." << endl;

  ofstream outfile;
  string python = "transition-data/python_graph.txt";
  outfile.open(python, ios::app);
  outfile << "edges = [ ";
  for (auto i : switch_tally)
  {
    outfile << i.first << ", ";
  }
  outfile << "]" << endl;

  outfile << "edge_widths = [ ";
  for (auto i : switch_tally)
  {
    double logged = log(i.second) + 0.2;
    outfile << logged << ", ";
  }
  outfile << "]" << endl;

  if (cycles)
  {
    outfile << "nodes = [ ";
    for (auto i : phenotype_tally)
    {
      outfile << "'" << i.first << "', ";
    }
    outfile << "]" << endl;

    outfile << "node_sizes = [ ";
    for (auto i : phenotype_tally)
    {
      outfile << (i.second) / 10. << ", ";
    }
    outfile << "]" << endl;
  }
  else
  {
    outfile << "nodes = [ ";
    for (auto i : adult_tally)
    {
      outfile << "'" << i.first << "', ";
    }
    outfile << "]" << endl;

    outfile << "node_sizes = [ ";
    for (auto i : adult_tally)
    {
      outfile << (i.second) / 10. << ", ";
    }
    outfile << "]" << endl;
  }


}