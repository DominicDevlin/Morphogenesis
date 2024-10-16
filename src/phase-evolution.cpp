#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <math.h>
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "info.h"
#include "parameter.h"
#include "sqr.h"
#include "omp.h"
#include <chrono>
#include <random>
#include "fft.h"
#include <sys/stat.h>
#include "connections.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

using namespace std;

auto start = chrono::steady_clock::now();
//rng for making random networks
auto mseed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne( static_cast<mt19937::result_type>(mseed) );
std::uniform_real_distribution<double> double_num(0.0, 1.0);
std::uniform_int_distribution<> genes_dist(0, par.n_genes-1);
std::uniform_int_distribution<> activ_dist(0, par.n_activators-1);
std::uniform_int_distribution<> TF_dist(0, par.n_TF-1);
std::normal_distribution<> morph_mut_dist(0,0.1);


int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

INIT {

  try 
  {
    CPM->set_seed();
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);
    CPM->ConstructInitCells(*this);
    CPM->SetRandomTypes();
    
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);

  }
}

TIMESTEP 
{
  cerr << "Error" << endl;
}

vector<double> random_gene(int length)
{
  vector<double> g{};

  for (int k=0; k < length; ++k)
  {
    double val = double_num(mersenne);
    int connection=0;
    if (val < 0.05)
    {
      connection = -2;
    }
    else if (val < 0.20)
    {
      connection = -1;
    }
    else if (val < 0.74)
    {
      connection = 0;
    }
    else if (val < 0.93)
    {
      connection = 1;
    }
    else
    {
      connection = 2;
    }
    g.push_back(connection);
  }
  return g; 
}

void ConstructNetwork(bool randomise, vector<vector<vector<double>>>& nets, vector<vector<vector<double>>>& morphs)
{
  int length = par.n_diffusers + par.n_MF + par.n_TF;
  int total = par.n_diffusers + par.n_MF + par.n_TF + 1;
  // we can just alter these appropriately.
  int old_d = 1;
  int old_mf = 2; // includes tf + mf


  for (auto &net : nets)
  {
    for (int i = 0; i < total; ++i)
    {
      if (i < old_d) // old gene
      {
        vector<double>& g = net[i];
        for (int j=1;j<length;++j)
        {
          if (j >= old_d && j < par.n_diffusers)
          {
            g.insert(g.begin() +j, 0.);
          }
          else if (j < length)
          {
            g.push_back(0.);
          }
        }
      }
      else if (i < par.n_diffusers) // new gene
      {
        if (randomise)
        {
          vector<double> new_g = random_gene(length);
          net.insert(net.begin()+i,new_g);
        }
        else
        {
          vector<double> new_g(length,0.);
          net.insert(net.begin()+i,new_g);
        }

      }
      else if (i < par.n_diffusers+old_mf) // old gene
      {
        vector<double>& g = net[i];
        for (int j=1;j<length;++j)
        {
          if (j >= old_d && j < par.n_diffusers)
          {
            g.insert(g.begin() +j, 0.);
          }
          else if (j < length)
          {
            g.push_back(0.);
          }
        }    
      }
      else if (i < length) // new gene
      {
        if (randomise)
        {
          vector<double> new_g = random_gene(length);
          net.insert(net.begin()+i,new_g);
        }
        else
        {
          vector<double> new_g(length,0.);
          net.insert(net.begin()+i,new_g);
        }    
      }
      else
      {
        vector<double>& g = net[i];
        for (int j=1;j<length;++j)
        {
          if (j >= old_d && j < par.n_diffusers)
          {
            g.insert(g.begin() +j, 0.);
          }
          else if (j < length)
          {
            g.push_back(0.);
          }
        }     
      }
    }
  }
  if (randomise)
  {
    double loglow = log10(5e-4);
    double logup = log10(5e-2);
    uniform_real_distribution<> sec_dis(loglow, logup);

    loglow = log10(5e-9);
    logup = log10(5e-6);
    uniform_real_distribution<> dif_dis(loglow, logup);

    for (auto &morph : morphs)
    {
      for (int i = 0; i < par.n_diffusers; ++i)
      {
        if (i >= old_d)
        {
          double secr = sec_dis(mersenne);
          secr = pow(10, secr);
          morph[i][0] = secr;

          double dif = dif_dis(mersenne);
          dif = pow(10, dif);
          morph[i][2] = dif;
        }
      }
    }
  }

}

void swap(double *xp, double *yp)
{
  double temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void swapv(vector<vector<double>> *xp, vector<vector<double>> *yp)
{
  vector<vector<double>> temp = *xp;
  *xp = *yp;
  *yp = temp;

}

void swapb(vector<bool> *xp, vector<bool> *yp)
{
  vector<bool> temp = *xp;
  *xp = *yp;
  *yp = temp;
}


void swapd(Dish *d, int max_idx, int i)
{
  Dish tmp = d[max_idx];
  d[max_idx] = d[i];
  d[i] = tmp;
}

void sorter(vector<vector<vector<double>>> &networks, vector<vector<vector<double>>> &morphogens, vector<double> &fitlist, Dish *dishes)
{
  int i, j, max_idx;
  int n = par.n_orgs;

  for (i = 0; i < n-1; i++)
  {
    max_idx = i;
    for (j=i+1;j<n;j++)
    {
      if (fitlist.at(j) > fitlist.at(max_idx))
        max_idx = j;
    }
    // swap largest element with first element
    swap(&fitlist.at(max_idx), &fitlist.at(i));
    swapv(&networks.at(max_idx), &networks.at(i));
    swapv(&morphogens.at(max_idx), &morphogens.at(i));
    // std::swap(dishes[max_idx], dishes[i]);

    // swapd(&dishes[max_idx], &dishes[i]);
  }
}

void addEdge(map<int, set<int>>& graph, int u, int v) {
    graph[u].insert(v);
    // graph[v].insert(u); // Remove this line if the graph is directed
}

bool dfs(map<int, set<int>>& graph, int current, int target, set<int>& visited) {
    if (current == target) return true; // target node found
    visited.insert(current); // mark the current node as visited
    for (int neighbor : graph[current]) {
        if (visited.find(neighbor) == visited.end()) { // if neighbor hasn't been visited
            if (dfs(graph, neighbor, target, visited)) return true;
        }
    }
    return false;
}

bool isPathExists(vector<pair<int, int>>& edges, vector<int>& nodes, int start, int end) {
    // Build the graph
    map<int, set<int>> graph;
    for (auto edge : edges) {
        addEdge(graph, edge.first, edge.second);
    }
    
    // Set of visited nodes
    set<int> visited;
    
    // Start DFS from the start node
    return dfs(graph, start, end, visited);
}


// randomise a new network.  
vector<vector<double>> get_random_network()
{
  vector<vector<double>> matrix;
  matrix.resize(par.n_genes);
  for (int i=0; i < par.n_genes;++i)
  {
    matrix.at(i).resize(par.n_activators);
  }

  for (int i = 0; i < par.n_genes; ++i)
  {
    for (int j = 0; j < par.n_activators; ++j)
    {
      double val = double_num(mersenne);
      if (val < 0.07)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.20)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.7)
      {
        matrix[i][j] = 0;
      }
      else if (val < 0.93)
      {
        matrix[i][j] = 1;
      }
      else
      {
        matrix[i][j] = 2;
      }
    }
  }
  return matrix;
}



// mutate a network. Currently no bias towards ON when mutating networks. 
void mutate(vector<vector<double>> &network)
{
  for (int k=0; k < par.n_mutations; ++k)
  {
    double val = double_num(mersenne);
    int i = genes_dist(mersenne);
    int j = activ_dist(mersenne);
    if (val < 0.5)
    {
      network[i][j] -= 0.5;
      if (network[i][j] < -2)
        network[i][j] = -2;
    }
    else
    {
      network[i][j] += 0.5;
      if (network[i][j] > 2)
        network[i][j] = 2;
    }
  }
}

void mutate_morphogens(vector<vector<double>> &morph)
{
  
  double val = double_num(mersenne);
  val *= par.n_diffusers;
  int to_mutate = floor(val);
  if (to_mutate>= par.n_diffusers)
  {
    cerr << "mutation error\n";
    to_mutate = par.n_diffusers-1;
  }
  
  double mfac1 = morph_mut_dist(mersenne);
  double mfac2 = morph_mut_dist(mersenne);

  // we only mutate (0)=secretion rate and (2)=diffusion coefficient
  morph[to_mutate][0]= morph[to_mutate][0] * exp(-mfac1);
  // 5e-2 to 5e-4
  if (morph[to_mutate][0] > 5e-2)
    morph[to_mutate][0] = 5e-2;
  else if (morph[to_mutate][0] < 5e-4)
    morph[to_mutate][0] = 5e-4;

  morph[to_mutate][2]= morph[to_mutate][2] * exp(-mfac2);
  // 5e6 to 5e9
  if (morph[to_mutate][2] > 5e-6)
    morph[to_mutate][2] = 5e-6;
  else if (morph[to_mutate][2] < 5e-9)
    morph[to_mutate][2] = 5e-9;

}


void output_networks(vector<vector<vector<double>>>& netw)
{
  for (int org=0;org<par.n_orgs;++org)
    for (int i=0;i<par.n_genes;++i)
    {
      if (i == 0)
        cout << "{ ";
      for (int j=0;j<par.n_activators;++j)
      {
        if (j==0)
          cout << "{ " << netw[org][i][j] << ", ";
        else if (j==par.n_activators-1)
          cout << netw[org][i][j] << " }, ";
        else 
          cout << netw[org][i][j] << ", ";
      }
      if (i == par.n_genes -1)
        cout << "}" << endl;
    }
}

void record_networks(vector<vector<vector<double>>>& netw, vector<vector<vector<double>>>& morphs, string oname)
{
  string nname = oname + "/" + "genomes.txt";
  std::ofstream outfile;
  outfile.open(nname, ios::app);
  for (int org=0;org<par.n_orgs;++org)
  {
    for (int i=0;i<par.n_genes;++i)
    {
      if (i == 0)
        outfile << "{ ";
      for (int j=0;j<par.n_activators;++j)
      {
        if (j==0)
          outfile << "{ " << netw[org][i][j] << ", ";
        else if (j==par.n_activators-1)
          outfile << netw[org][i][j] << " }, ";
        else 
          outfile << netw[org][i][j] << ", ";
      }
      if (i == par.n_genes -1)
        outfile << "}" << endl;
    }
    for (int i=0;i<par.n_diffusers;++i)
    {
      outfile << "secr_rate[" << i << "] = " << morphs[org][i][0] << "; decay_rate[" << i << "] = " << morphs[org][i][1] << "; diff_coeff[" << i << "] = " << morphs[org][i][2] << "; ";
    }
    outfile << endl; 
  }
  outfile.close();
}



void printn(vector<vector<double>> netw, vector<vector<double>> morph, vector<double> fitn, string oname)
{
  // create and open file
  std::string var_name = oname + "/gene_networks.txt";
  std::ofstream outfile;
  outfile.open(var_name, ios::app);

  for (int i=0;i<par.n_genes;++i)
  {
    if (i == 0)
      outfile << "start_matrix = { ";
    for (int j=0;j<par.n_activators;++j)
    {
      if (j==0)
        outfile << "{ " << netw[i][j] << ", ";
      else if (j==par.n_activators-1)
        outfile << netw[i][j] << " }, ";
      else 
        outfile << netw[i][j] << ", ";
    }
    if (i == par.n_genes -1)
      outfile << "};" << endl;
  }
  for (int i=0;i<par.n_diffusers;++i)
  {
    outfile << "secr_rate[" << i << "] = " << morph[i][0] << "; decay_rate[" << i << "] = " << morph[i][1] << "; diff_coeff[" << i << "] = " << morph[i][2] << "; ";
  } 
  outfile << endl;

  outfile.close();

  // max fitness 
  double max_fit = fitn.front();

  //average fitness
  double avgfit = 0;
  for (double i : fitn)
  {
    avgfit += i;
  }
  avgfit = avgfit / par.n_orgs;


  if (par.asym_only && par.asymmetry_selection && avgfit > par.swap_selection)
  {
    par.asym_only = false;
  }

  //calculate time since begin
  auto end = chrono::steady_clock::now();
  auto diff = end - start;

  //output fitness and time since beginning simulation. 
  var_name = oname + "/fitness.txt";
  outfile.open(var_name, ios::app);
  outfile << max_fit << '\t' << avgfit << '\t' << chrono::duration <double, milli> (diff).count() << endl;

}




// function that simulates a population for a single evolutionary step. 
vector<double> process_population(vector<vector<vector<double>>>& network_list, vector<vector<vector<double>>>& morphogens, int time)
{
  vector<double> inter_org_fitness{};
  inter_org_fitness.resize(par.n_orgs);

  // create memory for dishes. 
  Dish* dishes = new Dish[par.n_orgs];

  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].SetMorphogens(morphogens[i]);
  }

  // run organisms in parallel. 
  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for 
  for (int i=0; i < par.n_orgs; ++i)  
  {
    dishes[i].CPM->set_num(i+1);
    // does init block above.
    dishes[i].Init();

    int t;

    dishes[i].CPM->start_network(network_list.at(i));
    dishes[i].CPM->Set_evoJ(par.J_stem_diff);
    


    // make temperature lower for division section
    dishes[i].CPM->CopyProb(par.T);

    // run simulation for single organism for mcs montecarlo steps.
    for (t=0;t<par.mcs;t++) 
    {
      // PROGRAMMED CELL DIVISION SECTION
      if (t < par.end_program)
      {

        //programmed divisions
        if (t % par.div_freq == 0 && t <= par.div_end)
        {
          dishes[i].CPM->Programmed_Division(par.phase_evolution); // need to get the number of divisions right. 
        }

        
        if (t >= par.begin_network && t % par.update_freq == 0)
        {
          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell();
          for (int r=0;r<par.program_its;r++) 
          {
            dishes[i].PDEfield->Secrete(dishes[i].CPM);
            dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
          } 
        }
      }
      else 
      {
        // Normal division stage

        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell();
        }
        // diffusion stuff
        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); 
        }  
        dishes[i].CPM->ConstrainedGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);
    
      bool check_end = dishes[i].CPM->EndOptimizer();
      if (check_end)
        t = par.mcs-1;
 
      // ensure all cells are connected for shape calculations. 
      if (t % 1000 == 0 || (t<1000 && t>par.end_program && t%100 == 0))
      {
        bool check_shape = dishes[i].CPM->CheckShape();
        if (check_shape == false)
        {
          inter_org_fitness[i] = 0;
          t = par.mcs;
          // cout << "Org number: " << i << " has bad shape. " << endl;
        }
      }
      // get fitness at end of development
      if (t == par.mcs-1)
      {
        // need fitness function here
        inter_org_fitness[i] = dishes[i].CPM->DistanceTravelled();

        if (par.type_selection)
        {
          map<int, int> phens = dishes[i].CPM->get_phenotype_time();
          map<int, int> types = dishes[i].CPM->get_AdultTypes();  

          map<pair<int,int>,int> edge_tally{};
          dishes[i].CPM->set_switches(edge_tally);

          for (auto it = edge_tally.begin(); it != edge_tally.end();)
          {
            if (it->second < par.prune_amount)
            {
              it = edge_tally.erase(it);
            }
            else
            {
              ++it;
            }
          }
          vector<vector<int>> scc;
          map<int,int>subcomps{};
          Graph ungraph(types.size());
          subcomps = ungraph.CreateUnGraph(phens, types, edge_tally);
          scc = ungraph.GetComps(types, 8000);
          for (auto n : scc)
          {
              cout << "component: ";
              for (int j : n)
              {
                  cout << j << "  ";
              }
              cout << std::endl;
          }

          vector<vector<int>> temp_dendrogram(scc.size());
            
          vector<pair<int,int>> temp_edges{};
          vector<int> temp_nodes{};
          for (auto n : edge_tally)
          {
            if (n.second >= par.prune_amount)
            {
              temp_edges.push_back(n.first);
              if (std::find(temp_nodes.begin(), temp_nodes.end(), n.first.first) == temp_nodes.end())
              {
                temp_nodes.push_back(n.first.first);
              }
              if (std::find(temp_nodes.begin(), temp_nodes.end(), n.first.second) == temp_nodes.end())
              {
                temp_nodes.push_back(n.first.second);
              }
            }
          }

          for (int n = 0; n < scc.size(); n++)
          {
            for (int j = 0; j < scc.size(); j++)
            {
              if (j != n)
              {
                for (int x=0; x <scc[n].size();++x)
                {
                  for (int y = 0; y < scc[j].size();++y)
                  {
                    // cout << i << '\t' << j << endl;
                    int start = scc[n][x];
                    int end = scc[j][y];
                    // cout << start << '\t' << end << endl;
                    // Check if there is a path from start to end
                    if (isPathExists(temp_edges, temp_nodes, start, end)) 
                    {
                      temp_dendrogram[n].push_back(j);
                      y=scc[j].size();
                      x=scc[n].size();
                      // std::cout << "There is a path between " << start << " and " << end << std::endl;
                    } 
                    else 
                    {
                      // std::cout << "No path exists between " << start << " and " << end << std::endl;
                    }
                  }
                }
              }
            }
          }

          vector<int> dead_ends{};
          for (int i = 0; i < temp_dendrogram.size(); i++)
          {
            if (temp_dendrogram[i].size() == 0)
            {
              cout << "Dead end at " << scc[i][0] << endl;
              dead_ends.push_back(i);
            }
          }
          set<int> dead_ends_set(dead_ends.begin(), dead_ends.end());


          vector<vector<int>> dendrogram(scc.size());
          vector<pair<int,int>> edges{};
          vector<int> nodes{};
          for (auto n : edge_tally)
          {
            if (n.second > 20)
            {
              edges.push_back(n.first);
              if (std::find(nodes.begin(), nodes.end(), n.first.first) == nodes.end())
              {
                nodes.push_back(n.first.first);
              }
              if (std::find(nodes.begin(), nodes.end(), n.first.second) == nodes.end())
              {
                nodes.push_back(n.first.second);
              }
            }
          }

          for (int n = 0; n < scc.size(); n++)
          {
            for (int j = 0; j < scc.size(); j++)
            {
              if (j != n)
              {
                for (int x=0; x <scc[n].size();++x)
                {
                  for (int y = 0; y < scc[j].size();++y)
                  {
                    // cout << i << '\t' << j << endl;
                    int start = scc[n][x];
                    int end = scc[j][y];
                    // cout << start << '\t' << end << endl;
                    // Check if there is a path from start to end
                    if (isPathExists(edges, nodes, start, end)) 
                    {
                      dendrogram[n].push_back(j);
                      y=scc[j].size();
                      x=scc[n].size();
                      std::cout << "There is a path between " << start << " and " << end << std::endl;
                    } 
                    else 
                    {
                      std::cout << "No path exists between " << start << " and " << end << std::endl;
                    }
                  }
                }
              }
            }
          }

          // Filter dendrogram
          for (auto& vec : dendrogram) 
          {
              // Remove elements not in dead_end_set
            vec.erase(std::remove_if(vec.begin(), vec.end(),[&dead_ends_set](int x) { return dead_ends_set.find(x) == dead_ends_set.end(); }), vec.end());
          }

          int max_diffs=0;

          for (int n = 0; n < dendrogram.size(); ++n)
          {
            if (dendrogram[n].size() > max_diffs)
              max_diffs = dendrogram[n].size();
          }
          cout << "differentiates into: " << max_diffs << endl;
          inter_org_fitness[i] += inter_org_fitness[i] * (0.2*max_diffs);     
        }
      }        
    }
        
    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }

  if (par.evo_pics && time % par.pic_gen_interval == 0)
  {
    string dirn = par.data_file + "/" + to_string(time+1);
    if (mkdir(dirn.c_str(), 0777) != -1)
      cout << "Directory created." << endl;
    record_networks(network_list, morphogens, dirn);
    for (int i=0; i < par.n_orgs; ++i)
    {
      // change this function for max state space
      dishes[i].CPM->ColourCells(par.phase_evolution);
      fft new_org(par.sizex,par.sizey);
      new_org.ImportCPM(dishes[i].get_cpm());
      string f2 = "org-";
      string n2 = to_string(i);
      string ftype = ".png";
      string foutput = dirn + "/" + f2 + n2 + ftype;
      new_org.cpmOutput(foutput);

    }
  }

  delete[] dishes;

  // do sorting algorithm and return fitness
  sorter(network_list, morphogens, inter_org_fitness, dishes);

  //output to standard output
  output_networks(network_list);

  // output to file
  printn(network_list.front(), morphogens.front(), inter_org_fitness, par.data_file);

  // output_Js(Js, par.data_file);


  vector<vector<vector<double>>> nextgen{};
  vector<vector<vector<double>>> nextgen_morphs{};
  int j = 0;
  for (int i=0; i < par.n_orgs;++i)
  {
    nextgen.push_back(network_list[j]);
    nextgen_morphs.push_back(morphogens[j]);
    //mutate network with probability = mut_rate
    double mu = double_num(mersenne);
    if (mu < par.mut_rate)
    {
      mutate(nextgen.back());
    }
    double mu2 = double_num(mersenne);
    if (mu2 < par.mut_rate)
    {
      mutate_morphogens(nextgen_morphs.back());
    }
    ++j; 
    if (j >= par.n_orgs / 4)
      j=0;
  }

  // set nextgen = current pop. 
  for (int i=0;i<par.n_orgs;++i)
  {
    network_list.at(i) = nextgen.at(i);
    morphogens[i] = nextgen_morphs[i];
  }
  return inter_org_fitness;
}




// Main function
int main(int argc, char *argv[]) 
{
  vector<vector<double>> start_genome;
  if (par.starter)
    start_genome = par.start_matrix;
  else
    start_genome = par.start_n;
  vector<vector<double>> start_morph;
  for (int j = 0; j < par.n_diffusers; ++j)
  {
    vector<double> values{par.secr_rate[j], par.decay_rate[j], par.diff_coeff[j]};
    // cout << par.secr_rate[j] << '\t' <<  par.decay_rate[j] << '\t' << par.diff_coeff[j] << endl;
    start_morph.push_back(values);
  }
  string dirn = "evo-data";

  if (par.file_genomes && argc > 1)
  {
    dirn = dirn + "-" + argv[1];
    ifstream file("genomes.txt");
    vector<vector<vector<double>>> genomes;
    string line;
    while (getline(file, line)) 
    {
      vector<vector<double>> genome;

      vector<double> row{};
      stringstream ss(line);
      string value;
      while (ss >> value)
      {
        for (int i=0;i < value.size();++i)
        {
          
          if (value[i] == '-')
          {
            string ns{ value[i] };
            ++i;
            bool check = true;
            while (check)
            {
              if (isdigit(value[i]) || value[i] == '.' )
              {
                ns = ns + value[i];
                string nsn = ns;
                ++i;
              }
              else
              {
                check = false;
              }
            }
            row.push_back(stod(ns));
          }
          else if (isdigit(value[i]))
          {
            string ns{ value[i] };
            ++i;
            bool check = true;
            while (check)
            {
              if (isdigit(value[i]) || value[i] == '.' )
              {
                ns = ns + value[i];
                string nsn = ns;
                ++i;
              }
              else
              {
                check = false;
              }
            }
            row.push_back(stod(ns));
          }
        }      

        if (row.size() == par.n_activators)
        {
          genome.push_back(row);
          row.clear();
        }          

      }
      genomes.push_back(genome);
    }

    int itera = stoi(argv[1]);
    start_genome = genomes[itera];

    ifstream nextfile("params.txt");
    vector<vector<vector<double>>> morphs;
    while (getline(nextfile, line)) 
    {
      vector<vector<double>> mph;

      vector<double> row{};
      stringstream ss(line);
      string value;
      while (ss >> value)
      {
        for (int i=0;i < value.size();++i)
        {
          
          if (value[i] == '.')
          {
            string ns{ value[i-1] };
            ns += value[i];
            ++i;
            bool check = true;
            while (check)
            {
              if (isdigit(value[i]) || value[i] == 'e' || value[i] == '+' || value[i] == '-')
              {
                ns = ns + value[i];
                string nsn = ns;
                ++i;
              }
              else
              {
                check = false;
              }
            }
            row.push_back(stod(ns));
          }
        }      

        if (row.size() == 3)
        {
          mph.push_back(row);
          row.clear();
        }          

      }
      morphs.push_back(mph);
    }

    for (auto i : morphs)
    {
      for (auto j : i)
        for (double k : j)
        {
          cout << k << '\t';
        }
      cout << endl;
    }
    start_morph = morphs[itera];
  }

#ifdef QTGRAPHICS
  if (par.evo_pics)
  {
    QApplication* a = new QApplication(argc, argv);
    // par.data_file = "evo-data";
    // if (mkdir(par.data_file.c_str(), 0777) != -1)
    //   cout << "Directory created." << endl;
  }
  
#endif

  
  if (mkdir(dirn.c_str(), 0777) != -1)
    cout << "Directory created." << endl;
  
  par.data_file = dirn;

  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_output=false;
  par.store = false;
  par.velocities=false;
  par.output_gamma=false;
  par.pickseed = 0;
  par.umap = false;
  par.output_sizes=false;

  par.mcs = 20000;
  par.end_program=100;
  par.phase_evolution = true;
  par.type_selection=true;
  par.max_statespace=true;
  par.sizex=150;
  par.sizey=250;
  
  par.offset=75;
  
  par.mut_rate=0.25;
  par.n_orgs = 60;

  par.prune_edges=true;
  par.prune_amount=5;

  // this is true for type selection
  par.gene_record = true;
  Parameter();


  // make initial random networks. 
  vector<vector<vector<double>>> networks{};
  for (int i=0;i<par.n_orgs;++i)
  {
    networks.push_back(start_genome);
  }
  vector<vector<vector<double>>> morphogens;
  for (int i = 0; i < par.n_orgs; ++i)
  {
    morphogens.push_back(start_morph);
  }
  if (!par.starter && !par.file_genomes)
  {
    bool randomise = true;
    ConstructNetwork(randomise, networks, morphogens);
  }



  for (int t=0;t<par.evs;++t)
  {
    cout << "current ev timestep is: " << t+1 << endl;
    // process population. 
    vector<double> fit = process_population(networks, morphogens, t);

    // output every x evolution steps. 
    // if (t%1==0)
    // {
    //   printn(networks.front(), polarities.front(), fit);
    // }
  }
  // finished
  par.CleanUp();

  return 0;
}
