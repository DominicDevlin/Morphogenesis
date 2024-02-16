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

#include <sys/stat.h>
#include <cstring>

using namespace std;


int PDE::MapColour(double val)
{

  return (((int)((val / ((val) + 1.)) * 100)) % 100) + 155;
}

INIT 
{
  try 
  {
    CPM->set_seed();
    CPM->set_datafile(par.data_file);
    // Define initial distribution of cells
    // CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);

    CPM->FillGrid();
    CPM->ConstructInitCells(*this);


    if (par.velocities)
      par.output_sizes = true;
    
    // for (int i=0;i<par.divisions;i++) {
    //   CPM->DivideCells();
    // }

    CPM->FractureSheet();
    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();

    
  } 
  catch(const char* error) 
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

void process_population(vector<vector<vector<int>>> &network_list, vector<vector<bool>> &pols)
{

  vector<vector<double>> cell_displacements;

  Dish *dishes = new Dish[par.n_orgs];

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();
    

    int t;
    dishes[i].CPM->start_network(network_list.at(i), pols.at(i));


    for (t = 0; t < par.mcs; t++)
    {
      if (t == par.end_program)
      {
        dishes[i].CPM->CopyProb(par.lT); // normal temperature for normal development timing.
        par.T = par.lT;
      }

      // record initial expression state. This occurs before any time step updates. 
      if (t == 0)
      {
        if (par.flush_cells)
        {
          dishes[i].CPM->SetAllStates();
          dishes[i].PDEfield->FlushGrid();
        }

      }
        
      
      {
        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_network(t);
          if (par.noise && t > par.noise_start)
            dishes[i].CPM->add_noise();
            
          dishes[i].AverageChemCell(); 

 

        }
        
        if (par.velocities)
        {
          dishes[i].CPM->RecordMasses();
        }

        if (par.output_sizes)
        {
          dishes[i].CPM->RecordSizes();
        }

        // for (int r=0;r<par.pde_its;r++) 
        // {
        //   dish->PDEfield->Secrete(dish->CPM);
        //   dish->PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        // }


        // dish->CPM->CellGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);

    }
  }

  if (par.output_sizes)
  {
    for (int i = 0; i < par.n_orgs; ++i)
    {
      vector<vector<double>> displc = dishes[i].CPM->ReturnMSD();
      for (auto i : displc)
        cell_displacements.push_back(i);
    }
  }
  string var_name = par.data_file + "/meansqauredisplacement.dat"; 
  ofstream outfile;
  outfile.open(var_name, ios::app);  

  int timer = 0;
  for (auto j=0;j<par.mcs-par.equilibriate;++j)
  {
    double msd=0;
    int n = 0;

    for (auto i=0;i<cell_displacements.size();++i)
    {
      msd += cell_displacements[i][j];
      ++n;
    }
    outfile << timer << "\t" << msd/((double)n) << endl;
      
    ++timer;
  }


  outfile.close();

}




int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record=false;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  par.output_sizes = true;
  par.mcs=40000;
  par.sizex=150;
  par.sizey=150;
  par.end_program=0;
  Parameter();



  par.n_orgs = 10;

  vector<bool> start_p = {0, 0, 0, 0};

  // make initial random networks.
  vector<vector<vector<int>>> networks{};
  vector<vector<bool>> polarities{};
  for (int i = 0; i < par.n_orgs; ++i)
  {
    networks.push_back(par.start_matrix);
    polarities.push_back(start_p);
  }

  process_population(networks, polarities);


  

  // make initial random networks. 
  // vector<vector<vector<int>>> networks{};
  // vector<vector<bool>> polarities{};
  // for (int i=0;i<par.n_orgs;++i)
  // {
  //     networks.push_back(par.start_matrix);
  //     polarities.push_back(start_p);
  // }

  // process_population(networks, polarities);

  // finished
  par.CleanUp();

  return 0;
}
