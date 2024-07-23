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
#include <sys/stat.h>

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

using namespace std;

auto start = chrono::steady_clock::now();
// rng for making random networks
auto mseed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne(static_cast<mt19937::result_type>(mseed));
std::uniform_real_distribution<double> double_num(0.0, 1.0);
std::uniform_int_distribution<> genes_dist(0, par.n_genes - 1);
std::uniform_int_distribution<> activ_dist(0, par.n_activators - 1);
std::uniform_int_distribution<> TF_dist(0, par.n_TF - 1);

int PDE::MapColour(double val)
{

  return (((int)((val / ((val) + 1.)) * 100)) % 100) + 155;
}

INIT
{

  try
  {
    CPM->set_seed();
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);
    CPM->ConstructInitCells(*this);
    CPM->SetRandomTypes();
  }
  catch (const char *error)
  {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
}

TIMESTEP
{
  cerr << "Error" << endl;
}

// function that simulates a population for a single evolutionary step.
vector<double> process_population(vector<vector<vector<double>>> &network_list)
{
  vector<double> inter_org_fitness{};
  inter_org_fitness.resize(par.n_orgs);
  int shape_checking{};

  // create memory for dishes.
  Dish *dishes = new Dish[par.n_orgs];
  

  // run organisms in parallel.
  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();

    int t;

    dishes[i].CPM->start_network(network_list.at(i));
    dishes[i].CPM->Set_evoJ(par.J_stem_diff);


    // make temperature lower for division section
    dishes[i].CPM->CopyProb(par.T);

    // run simulation for single organism for mcs montecarlo steps.
    for (t = 0; t < par.mcs; t++)
    {
      // PROGRAMMED CELL DIVISION SECTION
      if (t < par.end_program)
      {

        // programmed divisions
        if (t % par.div_freq == 0 && t <= par.div_end)
        {
          dishes[i].CPM->Programmed_Division(); // need to get the number of divisions right.
        }

        if (t >= par.begin_network && t % par.update_freq == 0)
        {
          if (par.phase_evolution)
            dishes[i].CPM->update_phase_network(t);
          else
            dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
          for (int r = 0; r < par.program_its; r++)
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
          if (par.phase_evolution)
            dishes[i].CPM->update_phase_network(t);
          else
            dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
        }

        // diffusion stuff
        for (int r = 0; r < par.pde_its; r++)
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1);
        }
        if (par.phase_evolution)
          dishes[i].CPM->ConstrainedGrowthAndDivision(t);
        else
          dishes[i].CPM->CellGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);

      bool check_end = dishes[i].CPM->EndOptimizer();
      if (check_end)
      {
        cout << "made to end" << endl;
        t = par.mcs-1;
      }

      if (t >= 6000 && t < 8000 && t % 40 == 0 && par.scramble)
        dishes[i].CPM->swap_cells();     

      // ensure all cells are connected for shape calculations.
      if (t > 0 && t % 1000 == 0)
      {
        bool check_shape = dishes[i].CPM->CheckShape();
        if (check_shape == false)
        {
          shape_checking += 1;
          //inter_org_fitness[i] = 0;
          // t = par.mcs;
        }
        // cout << "Org number: " << i << " has bad shape. N-cells: " << dishes[i].CPM->CountCells() << endl;
      }
      // get fitness at end of development
      // if (t == par.mcs - 1)
      // {
      //   inter_org_fitness[i] = dishes[i].CPM->get_fitness();
      // }
    }

    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;
  }

  bool cycles=false;
  for (int i = 0; i < par.n_orgs; ++i)
  {
    if (dishes[i].CPM->CycleCheck())
    {
      cycles = true;
    }
  }
  if (cycles)
    cout << "THERE IS CYCLING" << endl;

  if (shape_checking)
    cout << shape_checking << " times were organisms not connected in space." << endl;

  
  storage stores;
  map<pair<int,int>,int> edge_tally{};

  for (int i = 0; i < par.n_orgs; ++i)
  {
    if (par.count_bud_cells)
      dishes[i].CPM->CheckCellsInBud();

      
    map<int, int> phentime = dishes[i].CPM->get_phenotype_time();
    map<int, int> difftime = dishes[i].CPM->get_AdultTypes();
    stores.add_to_time(phentime, difftime);
    stores.add_to_switches(dishes[i].CPM->transitions(cycles));
    

    if (cycles)
      dishes[i].CPM->set_long_switches(edge_tally);
    else
      dishes[i].CPM->set_switches(edge_tally);
  }


  map<int, int> phentally = stores.get_phenotype_tally();
  map<int, int> adulttally = stores.get_adult_tally();

  map<int, int> subcomps{};

  if (cycles)
  {
    Graph UnGraph(phentally.size());
    subcomps = UnGraph.CreateUnGraph(phentally, phentally, edge_tally, par.n_orgs);
  }
  else
  {
    Graph UnGraph(adulttally.size());
    subcomps = UnGraph.CreateUnGraph(phentally, adulttally, edge_tally, par.n_orgs);
  }



  if (mkdir("transition-data", 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;

  ofstream outfile;
  string switch_out = "transition-data/potency.dat";
  outfile.open(switch_out, ios::app);

  for (auto kv : subcomps)
  {
    outfile << "Component number: " << kv.first;
    if (kv.second > 1)
      outfile << " is weakly connected." << endl;
    else
      outfile << " is strongly connected." << endl;
  }
  outfile.close();

  stores.do_averaging();
  stores.write_to_file(cycles);

  delete[] dishes;
  // do sorting algorithm and return fitness

  return inter_org_fitness;
}

// Main function
int main(int argc, char *argv[])
{

  par.graphics = false;
  par.contours = false;
  par.print_fitness = false;
  par.randomise = false;
  par.gene_record = true;
  par.gene_output = false;
  par.velocities = false;
  Parameter();
  par.n_orgs = 4;
  par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 20 * par.n_orgs);

  par.sizex=150;
  par.sizey=250;
  par.end_program = 100;
  par.adult_begins = 1000;
  par.offset=75;



  // make initial random networks.
  vector<vector<vector<double>>> networks{};
  for (int i = 0; i < par.n_orgs; ++i)
  {
    networks.push_back(par.start_matrix);
  }

  process_population(networks);

  // finished
  par.CleanUp();

  return 0;
}
