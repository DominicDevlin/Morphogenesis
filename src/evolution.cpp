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
      if (val < 0.05)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.20)
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
    else if (val < 0.20)
    {
      network[i][j] = -1;
    }
    else if (val < 0.74)
    {
      network[i][j] = 0;
    }
    else if (val < 0.93)
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


void output_networks(vector<vector<vector<int>>>& netw)
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
  // outfile << "{ ";
  // for (int i=0;i<par.n_TF;++i)
  // {
  //   if (i < par.n_TF - 1)
  //     outfile << pol[i] << ", ";
  //   else 
  //     outfile << pol[i] << " }";
  // }
  // outfile << endl;
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
  var_name = "fitness.txt";
  outfile.open(var_name, ios::app);
  outfile << max_fit << '\t' << avgfit << '\t' << chrono::duration <double, milli> (diff).count() << endl;

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
        }
      }
      else 
      {
        // Normal division stage

        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
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
    }
        
    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }
  delete[] dishes;


  // do sorting algorithm and return fitness
  sorter(network_list, pols, inter_org_fitness);

  //output to standard output
  output_networks(network_list);

  // output to file
  printn(network_list.front(), pols.front(), inter_org_fitness);

  vector<vector<vector<int>>> nextgen{};
  vector<vector<bool>> nextgenpol{};
  int j = 0;
  for (int i=0; i < par.n_orgs;++i)
  {
    // Currently no random networks are added if largest fitness > this
    if (inter_org_fitness.front() > 30 || !par.insert_randoms)
    {
      nextgen.push_back(network_list.at(j));
      nextgenpol.push_back(pols.at(j));

      //mutate network with probability = mut_rate
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
      ++j; 
      if (j >= par.n_orgs / 4)
        j=0;
    }
    else
    {
      // the last 1/4 are random networks
      if (i >= (par.n_orgs * 3)/4)
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

      ++j;
      if (j >= par.n_orgs / 4)
        j=0;
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




// Main function
int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record = false;
  par.store = false;
  par.pickseed = 0;
  par.mcs = 12100;
  Parameter();

  // This is currently depracated. 
  vector<bool> start_p = { 0, 0, 0, 0 };

  // make initial random networks. 
  vector<vector<vector<int>>> networks{};
  vector<vector<bool>> polarities{};
  for (int i=0;i<par.n_orgs;++i)
  {
    if (par.starter)
    {
      networks.push_back(par.start_n);
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
    // if (t%1==0)
    // {
    //   printn(networks.front(), polarities.front(), fit);
    // }
  }
  // finished
  par.CleanUp();

  return 0;
}
