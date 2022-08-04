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



void swap(double *xp, double *yp)
{
  double temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void swapv(vector<vector<int>> *xp, vector<vector<int>> *yp)
{
  vector<vector<int>> temp = *xp;
  *xp = *yp;
  *yp = temp;

}

void swapb(vector<bool> *xp, vector<bool> *yp)
{
  vector<bool> temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void sorter(vector<vector<vector<int>>> &networks, vector<vector<bool>> &pol, vector<double> &fitlist)
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
    swapb(&pol.at(max_idx), &pol.at(i));
  }
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
      if (val < 0.01)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.18)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.8)
      {
        matrix[i][j] = 0;
      }
      else if (val < 0.98)
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

vector<bool> get_random_pol()
{
  vector<bool> pol;
  pol.resize(par.n_TF);
  for (int i=0;i<par.n_TF;++i)
  {
    double val = double_num(mersenne);
    if (val < 0.7)
      pol[i]=0;
    else
      pol[i]=1;
  }
  return pol;
}


// mutate a network. Currently no bias towards ON when mutating networks. 
void mutate(vector<vector<int>> &network)
{
  for (int k=0; k < par.n_mutations; ++k)
  {
    double val = double_num(mersenne);
    int i = genes_dist(mersenne);
    int j = activ_dist(mersenne);
    if (val < 0.05)
    {
      network[i][j] = -2;
    }
    else if (val < 0.2)
    {
      network[i][j] = -1;
    }
    else if (val < 0.8)
    {
      network[i][j] = 0;
    }
    else if (val < 0.95)
    {
      network[i][j] = 1;
    }
    else
    {
      network[i][j] = 2;
    }
  }
}

// mutate the TF polarities (whether each TF is passed onto daughter upon cell reproduction)
void polmutate(vector<bool> &pol)
{
  int i = TF_dist(mersenne);
  double val = double_num(mersenne);
  if (val > 0.7)
    pol[i] = true;
  else
    pol[i] = false;
}



// function that simulates a population for a single evolutionary step. 
vector<double> process_population(vector<vector<vector<int>>>& network_list, vector<vector<bool>> &pols)
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

    int t;

    dishes[i].CPM->start_network(network_list.at(i), pols.at(i));

    // make temperature lower for division section
    dishes[i].CPM->CopyProb(par.eT);
    par.T = par.eT;

    // run simulation for single organism for mcs montecarlo steps.
    for (t=0;t<par.mcs;t++) 
    {
      if (t == par.end_program)
      {
        dishes[t].CPM->CopyProb(par.lT); // normal temperature for normal development timing. 
        par.T = par.lT;
      } 
      // PROGRAMMED CELL DIVISION SECTION
      if (t < par.end_program)
      {

        if (t%par.div_freq==0 && t > 0 && t <= par.div_end)
          dishes[i].CPM->Programmed_Division();
        
        if (t >= par.begin_network && t % par.update_freq == 0)
        {
          if (t == par.begin_network && par.morphogen)
            dishes[i].CPM->morphogenWave();

          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
        }
      }
      else 
      {
        // CELLS CAN NOW ADHESE ACCORDING TO PROTEINS + HIGHER TEMP. 
        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
        }

        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); 
        }  
        

        dishes[i].CPM->CellGrowthAndDivision();
      }
      dishes[i].CPM->AmoebaeMove(t);

      // calculate the diversity over last 20% of time steps. 
      if (t > par.mcs * par.fitness_begin && t % par.fitness_typerate)
      {
        dishes[i].CPM->type_fitness();
        dishes[i].CPM->ShapeFitness();
      }

      // calculate rate of somatic cell production
      if (t % par.fitness_somrate == 0 && t>0)
        dishes[i].CPM->som_fitness();


      if (t == par.mcs-1)
      {
        inter_org_fitness.at(i) = dishes[i].CPM->get_fitness();
        cout << i << " inter org fitness: " << inter_org_fitness.at(i) << endl;
      }        
    }
        
    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }
  delete[] dishes;
  // do sorting algorithm and return fitness
  sorter(network_list, pols, inter_org_fitness);
  vector<vector<vector<int>>> nextgen{};
  vector<vector<bool>> nextgenpol{};
  int j = 0;
  for (int i=0; i < par.n_orgs;++i)
  {
    // Currently no random networks are added if largest fitness > 5
    if (inter_org_fitness.front() > 6)
    {
      nextgen.push_back(network_list.at(j));
      nextgenpol.push_back(pols.at(j));

      //mutate network with probability = mut_rate=.5
      double mu = double_num(mersenne);
      if (mu < par.mut_rate)
      {
        mutate(nextgen.back());
      }
      double mu2 = double_num(mersenne);
      if  (mu2 < par.polm_rate)
      {
        polmutate(nextgenpol.back());
      }

      if (j > par.n_orgs / 4)
        j=0;
      else
        ++j;
    }
    else
    {
      // the last 1/4 are random networks
      if (i > (par.n_orgs * 3)/4)
      {
        nextgen.push_back(get_random_network());
        nextgenpol.push_back(get_random_pol());
      }
      else 
      {
        nextgen.push_back(network_list.at(j));
        nextgenpol.push_back(pols.at(j));
      }

      //mutate network with probability = 0.5
      double mu = double_num(mersenne);
      if (mu > par.mut_rate)
      {
        mutate(nextgen.back());
      }
      double mu2 = double_num(mersenne);
      if  (mu2 > 0.5)
      {
        polmutate(nextgenpol.back());
      }

      if (j > par.n_orgs / 4)
        j=0;
      else
        ++j;
    }
  }
  // set nextgen = current pop. 
  for (int i=0;i<par.n_orgs;++i)
  {
    network_list.at(i) = nextgen.at(i);
    pols.at(i) = nextgenpol.at(i);
  }
  return inter_org_fitness;
}


void printn(vector<vector<int>> netw, vector<bool> pol, vector<double> fitn)
{
  // create and open file
  std::string var_name = "gene_networks.txt";
  std::ofstream outfile;
  outfile.open(var_name, ios::app);

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
  outfile << "{ ";
  for (int i=0;i<par.n_TF;++i)
  {
    if (i < par.n_TF - 1)
      outfile << pol[i] << ", ";
    else 
      outfile << pol[i] << " }";
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

  //calculate time since begin
  auto end = chrono::steady_clock::now();
  auto diff = end - start;

  //output fitness and time since beginning simulation. 
  var_name = "fitness.txt";
  outfile.open(var_name, ios::app);
  outfile << max_fit << '\t' << avgfit << '\t' << chrono::duration <double, milli> (diff).count() << endl;

}


// Main function
int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_output=false;
  Parameter();

  // If starting from a specific organism. 
  bool starter = false;
  vector<vector<int>> start_n = { { 0, 0, 1, -2, 0, -1, 0, 0, 0 }, { -1, 0, 0, 0, 0, 0, 1, 0, 0 }, { 1, 0, 0, -1, 1, 0, 1, 0, 0 }, { 1, 0, -1, 1, 0, 0, 1, 0, 1 }, { -1, 0, 0, 0, -1, 0, 1, 1, 1 }, { 0, 0, 0, 0, 1, 2, 1, 0, 0 }, { 0, -1, 1, 0, 1, 0, 0, -1, 1 }, { -1, -1, 0, 0, 1, -1, 1, 0, -1 }, { 0, 1, 1, 0, 0, -2, 1, 0, 0 }, { 2, 0, 0, 0, -1, 1, 0, -1, 0 }, { -1, -1, 0, 1, 1, 0, 1, 0, 0 }, { 0, 0, 0, 0, 0, 0, 1, -1, 1 }, { 0, 2, 0, 0, 0, -1, 0, 0, 1 }, { 0, 0, 0, 0, -1, 1, 0, 0, 0 }, { 0, 0, 1, 1, 0, 0, 1, 0, 0 }, { 0, -2, 1, 1, 0, 0, 0, 0, 0 }, { 0, 0, 0, -1, 0, 1, 0, -1, 0 }, { -1, -1, 0, 1, 1, 0, 0, 2, 0 }, { -1, 0, 0, 1, 0, -1, 1, 1, 1 }, { 1, 1, -2, 1, 0, 0, -1, 0, 0 }, };


  vector<bool> start_p = { 0, 0, 1, 0 };

  // make initial random networks. 
  vector<vector<vector<int>>> networks{};
  vector<vector<bool>> polarities{};
  for (int i=0;i<par.n_orgs;++i)
  {
    if (starter)
    {
      networks.push_back(start_n);
      polarities.push_back(start_p);
    }
    else
    {
      networks.push_back(get_random_network());
      polarities.push_back(get_random_pol());
    }
  }


  for (int t=0;t<par.evs;++t)
  {
    cout << "current ev timestep is: " << t+1 << endl;
    // process population. 
    vector<double> fit = process_population(networks, polarities);

    // output every x evolution steps. 
    if (t%1==0)
    {
      printn(networks.front(), polarities.front(), fit);
    }
  }
  // finished
  par.CleanUp();

  return 0;
}
