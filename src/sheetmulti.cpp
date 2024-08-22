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

  vector<vector<double>> hex_order(par.n_orgs);
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

    dishes[i].CPM->CopyProb(par.T);
    dishes[i].CPM->Set_J(par.sheet_J);
    dishes[i].CPM->set_mixJ(par.sheetmixJ);
    if (par.sheetmix)
    {
      dishes[i].CPM->StartSheetTypes();

    }
    int t;

    for (t = 0; t < par.mcs; t++)
    {              
      if (par.highT && t < par.highT_time)
      {
        // double diff = par.highT_temp - par.T;
        // double toT = par.highT_temp - diff * (double(t) / double(par.highT_time));
        double toT = par.highT_temp;
        dishes[i].CPM->CopyProb(toT);
      }
      if (t==par.highT_time)
        dishes[i].CPM->CopyProb(par.T);

      if (par.velocities)
      {
        dishes[i].CPM->RecordMasses();
      }

      if (par.sheetmix)
      {
        dishes[i].CPM->RandomSheetType();
      }

      if (t % 10 == 0 && t >= par.start_sheet_measure && t<= par.end_sheet_measure && par.sheet_hex)
      {
        dishes[i].CPM->initVolume();
        dishes[i].CPM->adjustPerimeters();
        // vector<double> tperims = dishes[i].CPM->PerimitersRadiusN(sqrt(13));
        vector<double> tperims = dishes[i].CPM->TruePerimeters();
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

        // now do hexatic order
        vector<double> hexes = dishes[i].CPM->GetHexes();
        for (auto &j : hexes)
        {
          hex_order[i].push_back(j);
        }
      }



      // if (par.output_sizes)
      // {
      //   dishes[i].CPM->RecordSizes();
      // }

      dishes[i].CPM->AmoebaeMove(t);
    }
  }

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;




  if (par.output_sizes)
  {
    for (int i = 0; i < par.n_orgs; ++i)
    {
      vector<vector<double>> displc = dishes[i].CPM->ReturnMSD();
      for (auto i : displc)
        cell_displacements.push_back(i);
    }
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << par.T << "-" << par.sheet_J;
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
  }

  if (par.sheet_hex)
  {

    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << par.sheet_J;
    string s = stream.str();
    string var_name = par.data_file + "/shapes-" + s + ".dat"; 
    ofstream outfile;
    outfile.open(var_name, ios::app);  

    // Print the vector of vectors as columns
    for (size_t i = 0; i < shape_index.size(); i++) 
    {
      for (size_t j = 0; j < shape_index[i].size(); j++) 
      {

        outfile << shape_index[i][j] << endl;
      }
    }
    outfile.close();  


    // now do same for hexatic order
    var_name = par.data_file + "/hex-order-" + s + ".dat";
    outfile.open(var_name, ios::app);  
    // Print the vector of vectors as columns
    for (size_t i = 0; i < hex_order.size(); i++) 
    {
      for (size_t j = 0; j < hex_order[i].size(); j++) 
      {

        outfile << hex_order[i][j] << endl;
      }
    }
    outfile.close();  
  }



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
  par.mcs=100000 + par.equilibriate;
  par.sizex=200;
  par.sizey=200;
  par.end_program=0;
  par.sheet=true;
  par.periodic_boundaries=true;

  if (par.velocities)
    par.output_sizes = true;
  else
    par.output_sizes = false;
  Parameter();


  if (par.sheet_hex)
  {
    par.mcs = par.end_sheet_measure + par.equilibriate;
  }

  par.periodic_boundaries = true;
  par.flush_cells = true;


  par.n_orgs = 60;

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
  int n_trials = ceil((par.sheet_maxJ-par.sheet_minJ)/par.J_width)+1;
  for (int i = 0; i < n_trials; ++i)
  {
    double J = par.sheet_minJ + par.J_width*i;
    Jlist.push_back(J);
  }

  vector<vector<double>> hex(par.n_orgs);
  vector<vector<double>> shape(par.n_orgs);

  vector<vector<vector<double>>> hex_order{};
  vector<vector<vector<double>>> shape_index{};

  for (int i = 0; i < n_trials; ++i)
  {
    hex_order.push_back(hex);
    shape_index.push_back(shape);
    par.sheet_J = Jlist[i];
    if (par.sheetmix)
      par.sheetmixJ = 2*par.sheet_J;
    cout << par.sheet_J << endl;
    process_population();
  }


  // if (par.sheet_hex)
  // {
  //   // collapse into one vector
  //   vector<vector<double>> shape_final(n_trials);
  //   vector<vector<double>> hex_final(n_trials);

  //   for (int i = 0; i < n_trials; ++i)
  //   {
  //     for (auto &j : shape_index[i])
  //       for (double &k : j)
  //       {
  //         shape_final[i].push_back(k);
  //       }
  //     for (auto &j : hex_order[i])
  //       for (double &k : j)
  //       {
  //         hex_final[i].push_back(k);
  //       }
  //   }
  //   string var_name = par.data_file + "/shapes.dat";
  //   ofstream outfile;
  //   outfile.open(var_name, ios::app);  
  //   int length=0;
  //   for (const auto& inner_vec : shape_final) 
  //   {
  //       if (inner_vec.size() > length) 
  //       {
  //           length = inner_vec.size();
  //       }
  //   }
  //   // Print the vector of vectors as columns
  //   for (size_t i = 0; i < length; i++) 
  //   {
  //     for (size_t j = 0; j < shape_final.size(); j++) 
  //     {
  //       int k = 0;
  //       if (i < shape_final[j].size()) 
  //       {
  //           outfile << shape_final[j][i] << "\t";
  //       } 
  //     }
  //     outfile << std::endl; // Newline after each column is printed
  //   }
  //   outfile.close();  
  //   // now do same for hexatic order
  //   var_name = par.data_file + "/hex-order.dat";
  //   outfile.open(var_name, ios::app);  
  //   length=0;
  //   for (const auto& inner_vec : hex_final) 
  //   {
  //       if (inner_vec.size() > length) 
  //       {
  //           length = inner_vec.size();
  //       }
  //   }
  //   // Print the vector of vectors as columns
  //   for (size_t i = 0; i < length; i++) 
  //   {
  //     for (size_t j = 0; j < hex_final.size(); j++) 
  //     {
  //       int k = 0;
  //       if (i < hex_final[j].size()) 
  //       {
  //           outfile << hex_final[j][i] << "\t";
  //       } 
  //     }
  //     outfile << std::endl; // Newline after each column is printed
  //   }
  //   outfile.close();  
  // }
  
  // finished
  par.CleanUp();

  return 0;
}
