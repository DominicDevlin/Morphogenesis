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

    par.sheet=true;
    par.periodic_boundaries=true;

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

void process_population()
{

  vector<vector<double>> shape_index(par.n_orgs);
  vector<vector<double>> cell_displacements;

  Dish *dishes = new Dish[par.n_orgs];

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();

    // equilibriate cells with high T
    if (par.highT)
    {
      dishes[i].CPM->CopyProb(par.highT_temp);
    }
    else
    {
      dishes[i].CPM->CopyProb(par.T);
    }

    int t;

    for (t = 0; t < par.mcs; t++)
    {              
      if (t==par.highT_time)
        dishes[i].CPM->CopyProb(par.T);

      if (par.velocities)
      {
        dishes[i].CPM->RecordMasses();
      }
      // if (par.output_sizes)
      // {
      //   dishes[i].CPM->RecordSizes();
      // }

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

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;

  // string sTemp = to_string(par.T);
  // string sJ = to_string(par.minJ);

  std::stringstream stream;
  stream << std::fixed << std::setprecision(1) << par.T << "-" << par.minJ;
  string s = stream.str();

  string var_name = par.data_file + "/msd" + s + ".dat"; 
  ofstream outfile;
  outfile.open(var_name, ios::app);  

  int timer = 1;
  for (auto j=0;j<par.mcs-par.equilibriate-1;++j)
  {
    double msd=0;
    int n = 0;

    for (auto i=0;i<cell_displacements.size();++i)
    {
      msd += cell_displacements[i][j];
      ++n;
    }
    // outfile << timer << "\t" << msd/((double)n) << endl;
    outfile << (msd)/((double)n) << endl;

    ++timer;
  }

  outfile.close();
  delete[] dishes;

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
  par.mcs=100000 + par.equilibriate;
  par.sizex=150;
  par.sizey=150;
  par.end_program=0;
  Parameter();

  par.periodic_boundaries = true;
  par.flush_cells = true;


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
  
  vector<double> Jlist;
  int n_trials = ceil((par.sheet_maxJ-par.sheet_minJ)/par.J_width);;
  for (int i = 0; i < n_trials; ++i)
  {
    double J = par.sheet_minJ + par.J_width*i;
    Jlist.push_back(J);
  }
  
  for (int i = 0; i < n_trials; ++i)
  {
    par.sheet_J = Jlist[i];
    process_population();
  }

  // finished
  par.CleanUp();

  return 0;
}
