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
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <map>


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



// function that simulates a population for a single evolutionary step. 
void process_population(vector<vector<vector<double>>>& network_list, vector<vector<bool>> &pols)
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
    dishes[i].CPM->CopyProb(par.T);

    // run simulation for single organism for mcs montecarlo steps.
    for (t=0;t<par.mcs;t++) 
    {
      if (t == 100)
        dishes[i].CPM->record_GRN();


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

          if (par.noise && t > par.noise_start)
            dishes[i].CPM->add_noise();
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
    }
      
    // cout << "Sim complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;



  }
  // dishes[0].CPM->PrintColourList();

  vector<vector<double>> sum_history{};
  vector<int> sum_types{};


  for (int i = 0; i < par.n_orgs;++i)
  {
    vector<vector<double>> history = dishes[i].CPM->OrganismGenes(153);
    vector<int> types = dishes[i].CPM->OrganismTypes(153);

    cout << history.size() << endl;
    cout << types.size() << endl;
    if (history.size() != types.size())
    {
      cout << "error with expression data!" << endl;
    }
    for (int k=0;k<history.size();++k)
    {
      sum_history.push_back(history[k]);
      sum_types.push_back(types[k]);
    }  
    
  }


  

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;  



  string fnamen = par.data_file + "/cdata.dat";
  ofstream outfile;
  outfile.open(fnamen, ios::app);

  for (int i = 0; i < par.n_genes;++i)
  {
    string n = "c" + to_string(i+1);
    outfile << n << '\t';
  }
  outfile << "class" << '\t' << "cell" << endl;

  // vector<int> col_index{};
  int count=1;
  int f = sum_history.size();
  for (int i = 0; i < f;++i)
  {
    for (auto j : sum_history[i])
    {
      outfile << j << '\t';
    }
    int p = sum_types[i];
    if (par.colour_index.find(p) == par.colour_index.end())
    {
      // we will just look for the closest one
      int min_dist=100000;
      int n=0;
      for (auto i : par.colour_index)
      {
        int dist = abs(p - i.first);
        if (dist < min_dist)
        {
          min_dist = dist;
          n = i.second;
        }
      }
      outfile << n;

    }
    else 
    {
      outfile << par.colour_index[p];
    }

    if (((i+1) % 150) == 0)
    {
      count+=1;
    }
    outfile << '\t' << count << endl;

  }


  outfile.close();


  // fnamen = par.data_file + "/index.dat";
  // outfile.open(fnamen, ios::app);
  // for (int i = 0; i < col_index.size(); ++i)
  // {
  //   outfile << i << "\t" << col_index[i] << endl;
  // }
  // outfile.close();

  delete[] dishes;



}




// Main function
int main(int argc, char *argv[]) {

  
  par.graphics = false;
  par.contours = false;
  par.print_fitness = false;
  par.randomise = false;
  par.gene_record = true;
  par.gene_output = false;
  par.velocities = false;
  par.umap = true;
  par.noise = true;

  Parameter();
  par.n_orgs = 4;
  // This is currently depracated.
  vector<bool> start_p = {0, 0, 0, 0};

  // make initial random networks.
  vector<vector<vector<double>>> networks{};
  vector<vector<bool>> polarities{};
  for (int i = 0; i < par.n_orgs; ++i)
  {
    networks.push_back(par.start_matrix);
    polarities.push_back(start_p);
  }

  process_population(networks, polarities);

  // finished
  par.CleanUp();

  return 0;
}
