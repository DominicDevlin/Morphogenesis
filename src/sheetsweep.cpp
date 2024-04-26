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
  vector<vector<double>> true_adhesion(par.n_orgs);
  vector<vector<double>> cell_displacements;

  Dish *dishes = new Dish[par.n_orgs];

  vector<double> Js;
  for (int i = 0;i<par.n_orgs;++i)
  {
    Js.push_back(par.sheet_minJ+i*par.J_width);
  }

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
    dishes[i].CPM->Set_J(Js[i]);

    int t;

    for (t = 0; t < par.mcs; t++)
    {              
      if (t==par.highT_time)
        dishes[i].CPM->CopyProb(par.T);

      if (par.velocities)
      {
        dishes[i].CPM->RecordMasses();
      }

      dishes[i].CPM->AmoebaeMove(t);

      if (t % 500 == 0 && t > 0)
      {
        dishes[i].CPM->initVolume();
        dishes[i].CPM->adjustPerimeters();
        vector<double> tperims = dishes[i].CPM->PerimitersRadiusN(sqrt(13));
        vector<double> volumes = dishes[i].CPM->GetVolumes();

        vector<double> tadhesion = dishes[i].CPM->TrueAdhesion();

        double avg{};
        for (int j = 0; j < tperims.size(); ++j)
        {
          double sindex = tperims[j] / sqrt(volumes[j]);
          // cout << i << '\t';
          avg+=sindex;
          shape_index[i].push_back(sindex);
        }
        avg/=tperims.size();

        double avg_adh{};
        for (int j = 0; j < tadhesion.size();++j)
        {
          true_adhesion[i].push_back(tadhesion[j]);
        }
        avg/=tadhesion.size();

      }
    }
  }


  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;


  string var_name = par.data_file + "/shapes.dat";
  ofstream outfile;
  outfile.open(var_name, ios::app);  

  size_t length=0;
  for (const auto& inner_vec : shape_index) 
  {
      if (inner_vec.size() > length) 
      {
          length = inner_vec.size();
      }
  }
  // Print the vector of vectors as columns
  for (size_t i = 0; i < length; i++) {
      for (size_t j = 0; j < shape_index.size(); j++) 
      {
          if (i < shape_index[j].size()) 
          {
              outfile << shape_index[j][i] << "\t";
          } 
      }
      outfile << std::endl; // Newline after each column is printed
  }
  outfile.close();

  // now do same for adhesion
  var_name = par.data_file + "/adhesions.dat";
  outfile.open(var_name, ios::app);  
  length=0;
  for (const auto& inner_vec : true_adhesion) 
  {
      if (inner_vec.size() > length) 
      {
          length = inner_vec.size();
      }
  }
  // Print the vector of vectors as columns
  for (size_t i = 0; i < length; i++) {
      for (size_t j = 0; j < true_adhesion.size(); j++) 
      {
          if (i < true_adhesion[j].size()) 
          {
              outfile << true_adhesion[j][i] << "\t";
          } 
      }
      outfile << std::endl; // Newline after each column is printed
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
  par.mcs=10000 + par.equilibriate;
  par.sizex=200;
  par.sizey=200;
  par.end_program=0;
  Parameter();

  par.periodic_boundaries = true;
  par.flush_cells = true;

  par.n_orgs = ceil((par.sheet_maxJ-par.sheet_minJ)/par.J_width);
  par.sheet_J=par.sheet_minJ;
  
  process_population();

  // finished
  par.CleanUp();

  return 0;
}
