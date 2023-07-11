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
#include "storage.h"
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


int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

INIT {

  try 
  {
    CPM->set_seed();
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
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

// randomise a new network.  
vector<vector<int>> get_random_network()
{
  vector<vector<int>> matrix;
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
      // slight ON bias for random networks. This is due to theta = -0.3. 
      if (val < 0.05)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.2)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.74)
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



// function that simulates a population for a single evolutionary step. 
vector<int> process_population(vector<vector<vector<int>>>& network_list, vector<vector<bool>> &pols, 
                                        vector<int>& potency, vector<int>& nodes, vector<double>& org_sizes, vector<int>& weak_cons, vector<int>& s_count)
{
  vector<double> inter_org_fitness{};
  inter_org_fitness.resize(par.n_orgs);

  // create memory for dishes. 
  Dish* dishes = new Dish[par.n_orgs];

  // run organisms in parallel. 
  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for 
  for (int i=0; i < par.n_orgs; ++i)  
  {
    dishes[i].CPM->set_num(i+1);
    // does init block above.
    dishes[i].Init();

    int t{};

    dishes[i].CPM->start_network(network_list[i], pols.at(i));




    // make temperature lower for division section
    dishes[i].CPM->CopyProb(par.eT);
    par.T = par.eT;

    // run simulation for single organism for mcs montecarlo steps.
    for (t=0;t<par.mcs;t++) 
    {
      if (t == 100)
        dishes[i].CPM->record_GRN();


      if (t == par.end_program)
      {
        dishes[i].CPM->CopyProb(par.lT); // normal temperature for normal development timing. 
        par.T = par.lT;
      } 
      // PROGRAMMED CELL DIVISION SECTION
      if (t < par.end_program)
      {

        //programmed divisions
        if (t % par.div_freq == 0 && t <= par.div_end)
        {
          dishes[i].CPM->Programmed_Division(); // need to get the number of divisions right. 
        }

        
        if (t >= par.begin_network && t % par.update_freq == 0)
        {
          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
          for (int r=0;r<par.program_its;r++) 
          {
            dishes[i].PDEfield->Secrete(dishes[i].CPM);
            dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
          } 
          dishes[i].CPM->record_GRN();  
        }
      }
      else 
      {
        // Normal division stage

        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
          dishes[i].CPM->record_GRN(); 

        }

        // diffusion stuff
        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); 
        }  
        dishes[i].CPM->CellGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);
    

      // calculate the diversity over last 20% of time steps. 
      if (t > par.mcs * par.fitness_begin && t % par.fitness_typerate == 0)
      {
        // am now doing for curvature as well (taking mean)
        dishes[i].CPM->update_fitness();
      }
 
      // ensure all cells are connected for shape calculations. 
      if (t > 0 && t % 1000 == 0)
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
        inter_org_fitness[i] = dishes[i].CPM->get_fitness();
      }

      // potency test!!
      if (t >= 6000 && t < 8000 && t % 20 == 0 && par.scramble)
        dishes[i].CPM->swap_cells();        
    }
        
    // cout << "Sim complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }
  // n_orgs must divide n_replicates
  int n_threads = par.n_orgs / par.n_replicates;
  

  //  dont parallelise
  for (int i=0; i < n_threads; ++i)
  {
    double org_size{};

    nodes[i] = 0;
    int shape_count = 0;
    storage stores{};

    bool cycles=false;
    for (int i = 0; i < par.n_orgs; ++i)
    {
      if (dishes[i].CPM->CycleCheck())
      {
        cout << "FOUND A LONG CYCLE" << endl;
        cycles = true;
        break;
      }
    }
    // need to use i * par.n_replicates
    int j = par.n_replicates * i;
    int maxj = j + par.n_replicates;

    for (; j < maxj; ++j)  
    {
      if (dishes[j].CPM->CheckShape() == true)
        ++shape_count;

      map<int,int> phentime = dishes[j].CPM->get_phenotype_time();
      map<int,int> difftime = dishes[j].CPM->get_AdultTypes();
      stores.add_to_time(phentime, difftime);
      map<pair<int,int>,int> tally{};
      if (cycles)
      {
        dishes[j].CPM->set_long_switches(tally);
        stores.add_to_edges(tally);
      }
      else
      {
        dishes[j].CPM->set_switches(tally);
        stores.add_to_edges(tally);        
      }

      org_size += dishes[j].CPM->CountCells();
    }
    org_size = org_size / par.n_replicates;
    org_sizes[i] = org_size;


    map<int, int> phentally = stores.get_phenotype_tally();
    map<int, int> adulttally = stores.get_adult_tally();
    
    nodes[i] = adulttally.size();    

    if (shape_count > (par.n_replicates / 2))
    {
      if (adulttally.size() > 1)
      {
        map<pair<int,int>,int> full_tally{};
        full_tally = stores.get_edges();
        // Graph LargeGraph(adulttally.size());
        // LargeGraph.CreateDiGraph(phentally, adulttally, starts, ends);
        // bool toti = LargeGraph.StronglyConnected();
        // bool connected = true;

        // NEED TO ACCOUNT FOR WEIRD CYCLE CASE JUST IN CASE!!
        map<int,int> subcomps{};
        if (cycles)
        {
          Graph UnGraph(phentally.size());
          subcomps = UnGraph.CreateUnGraph(phentally, phentally, full_tally, par.n_replicates, cycles);
        }
        else
        {
          Graph UnGraph(adulttally.size());
          subcomps = UnGraph.CreateUnGraph(phentally, adulttally, full_tally, par.n_replicates);
        }


        // connected vs unconnected
        potency[i] = subcomps.size();

        // number of weakly connected subgraphs
        int weak=0;
        for (auto &k : subcomps)
        {
          weak += k.second - 1;
        }
        // if (weak > 0)
        // {
        //   weak = 0;
        //   par.prune_edges = false;
        //   Graph UnGraph(adulttally.size());
        //   subcomps = UnGraph.CreateUnGraph(phentally, adulttally, full_tally, par.n_replicates);
        //   for (auto &k : subcomps)
        //   {
        //     weak += k.second - 1;
        //   }
        //   par.prune_edges = true;          
        // }
        weak_cons[i] = weak;
        s_count[i] = shape_count;
      }
      else
      {
        potency[i] = -1;
        weak_cons[i] = 0;
        s_count[i] = shape_count;
      }
    }
    else 
    {
      potency[i] = -2;
      s_count[i] = shape_count;
      weak_cons[i] = 0;
    }
  }


  delete[] dishes;
  // do sorting algorithm and return fitness


  return potency;
}


void OutputPotency(vector<vector<vector<int>>>& genomes, vector<int>& potency, vector<int>& nodes, vector<double>& org_sizes, vector<int>& weak_cons, vector<int>& s_count)
{
  string var_name = "p_networks.txt";
  ofstream outfile;
  outfile.open(var_name, ios::app);
  for (auto netw : genomes)
  {
    for (int i=0;i<par.n_genes;++i)
    {
      if (i == 0)
        outfile << "{ ";
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
        outfile << "}" << endl;
    }
  }
  outfile.close();

  var_name = "p_potencies.txt";
  outfile.open(var_name, ios::app);
  for (int i=0; i < potency.size(); ++i)
  {
    outfile << potency[i] << '\t' << weak_cons[i] << '\t' << nodes[i] << '\t' << org_sizes[i] << '\t' << s_count[i] << endl;
  }
  outfile.close();


}


// Main function
int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_record=true;
  par.gene_output=false;
  par.potency_edges=true;
  par.velocities = false;
  Parameter();

  par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10 * par.n_replicates);
  // This is currently depracated. 
  vector<bool> start_p = { 0, 0, 0, 0 };

  vector<vector<vector<int>>> networks{};
  vector<vector<bool>> polarities{};
  vector<vector<int>> new_net{};

  vector<vector<vector<int>>> current_genomes{};
  int s = par.n_orgs / par.n_replicates;


  vector<int> potency(s, 0);
  vector<int> weak_cons(s,0);
  vector<int> nodes(s, 0);
  vector<double> org_sizes(s,0);
  vector<int> s_count(s,0);

  for (int i=0;i<par.n_orgs;++i)
  {
      polarities.push_back(start_p);
  }

  int n_steps = 5000;
  // we have to create organisms in sets of 10

  for (int i = 0; i < n_steps; ++i)
  {
    int count = 0;
    new_net = get_random_network();
    // new_net = { { 1, 0, 0, 1, -1, 1, -1, 0, -1 }, { 1, -1, 0, 1, 0, 0, 0, 0, 0 }, { 0, 0, 1, 0, 1, 1, 0, 0, 1 }, { 0, 0, 0, 1, 0, -1, 0, 0, 2 }, { 0, 1, 0, -2, -1, 1, -1, 1, 1 }, { 0, -1, -1, 0, 0, 1, -1, 0, 1 }, { 0, 0, -1, -1, -1, 0, -1, 0, 2 }, { -1, 0, 0, -1, 1, 1, 2, 0, 0 }, { 0, -1, -1, 0, 2, 0, -2, -1, 0 }, { 0, 2, 0, -1, 1, 0, 0, 0, 0 }, { 0, -1, -2, -1, 0, 1, 1, 0, 0 }, { 0, 0, -2, 1, 0, 0, 0, 1, -2 }, { 0, -2, 1, 0, 2, -1, 0, -1, 0 }, { 1, 1, 0, 0, 1, -1, 0, 0, 0 }, { 2, 0, 0, 0, 0, -1, 0, 1, -2 }, { -1, 0, 0, 1, 0, 0, 1, 1, 1 }, { 1, 0, -1, 0, 1, 0, -1, -1, 1 }, { 2, 0, 2, 0, 0, 1, 0, 0, 1 }, { 0, 0, 0, 0, 0, 0, -2, 1, 0 }, { 0, 0, 2, 1, 0, -1, -1, 1, 0 }, { 0, -1, -1, 1, 0, -1, 0, 1, 0 }, { 1, 0, -1, 1, 0, -1, 1, -1, 0 }, { -1, 0, 0, 1, 1, 0, 0, 0, 0 }, { 0, 0, 2, 1, 1, 0, 0, 0, 0 }, { 1, -1, 0, 0, 0, 0, -2, 0, 0 }, { -2, 0, 0, 0, 2, -1, -1, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, -1, -1 }, };
    
    current_genomes.push_back(new_net);

    for (int i=0;i<par.n_orgs;++i)
    {
      if (count % par.n_replicates == 0 && count)
      {
        new_net = get_random_network();
        // new_net = { { 1, 0, 0, 1, -1, 1, -1, 0, -1 }, { 1, -1, 0, 1, 0, 0, 0, 0, 0 }, { 0, 0, 1, 0, 1, 1, 0, 0, 1 }, { 0, 0, 0, 1, 0, -1, 0, 0, 2 }, { 0, 1, 0, -2, -1, 1, -1, 1, 1 }, { 0, -1, -1, 0, 0, 1, -1, 0, 1 }, { 0, 0, -1, -1, -1, 0, -1, 0, 2 }, { -1, 0, 0, -1, 1, 1, 2, 0, 0 }, { 0, -1, -1, 0, 2, 0, -2, -1, 0 }, { 0, 2, 0, -1, 1, 0, 0, 0, 0 }, { 0, -1, -2, -1, 0, 1, 1, 0, 0 }, { 0, 0, -2, 1, 0, 0, 0, 1, -2 }, { 0, -2, 1, 0, 2, -1, 0, -1, 0 }, { 1, 1, 0, 0, 1, -1, 0, 0, 0 }, { 2, 0, 0, 0, 0, -1, 0, 1, -2 }, { -1, 0, 0, 1, 0, 0, 1, 1, 1 }, { 1, 0, -1, 0, 1, 0, -1, -1, 1 }, { 2, 0, 2, 0, 0, 1, 0, 0, 1 }, { 0, 0, 0, 0, 0, 0, -2, 1, 0 }, { 0, 0, 2, 1, 0, -1, -1, 1, 0 }, { 0, -1, -1, 1, 0, -1, 0, 1, 0 }, { 1, 0, -1, 1, 0, -1, 1, -1, 0 }, { -1, 0, 0, 1, 1, 0, 0, 0, 0 }, { 0, 0, 2, 1, 1, 0, 0, 0, 0 }, { 1, -1, 0, 0, 0, 0, -2, 0, 0 }, { -2, 0, 0, 0, 2, -1, -1, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, -1, -1 }, };
        current_genomes.push_back(new_net);
      }
      networks.push_back(new_net);
      polarities.push_back(start_p);
      ++count;
    }
    process_population(networks, polarities, potency, nodes, org_sizes, weak_cons, s_count);
    OutputPotency(current_genomes, potency, nodes, org_sizes, weak_cons, s_count);
    current_genomes.clear();
    networks.clear();
  }

  // finished
  par.CleanUp();

  return 0;
}
