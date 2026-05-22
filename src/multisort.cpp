/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <math.h>
#include "model/dish.h"
#include "model/random.h"
#include "model/cell.h"
#include "model/info.h"
#include "model/parameter.h"
#include "model/sqr.h"
#include "model/storage.h"
#include "model/connections.h"
#include "model/fft.h"
#include <sys/stat.h>
#include <chrono>
#include "omp.h"

#ifdef QTGRAPHICS
#include "model/qtgraph.h"
#else
#include "model/x11graph.h"
#endif

#include <sys/stat.h>
#include <cstring>

using namespace std;




INIT 
{
  try 
  {
    CPM->set_seed();
    CPM->set_datafile(par.data_file);
    // Define initial distribution of cells
    // CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);

    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);
    CPM->ConstructInitCells(*this);



    // if (par.do_voronoi)
    // {
    //   par.highT=false;
    //   CPM->Voronoi(par.sizex, par.sizey);
    // }
    // else
    // {
    //   CPM->FractureSheet();
    // }

    CPM->ClearGrid();
    if (par.make_sparse_cells)
    {
      CPM->PopulateSparseCells(0.8, 80, 0, 0);
    }
    else if (par.do_voronoi)
    {
      CPM->Voronoi(par.sizex, par.sizey);
    }


    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();


    if (par.set_colours)
    {
      CPM->SetColours();
    }
    par.end_program=0;
    if (par.store)
    {
      if (mkdir("data_film", 0777) == -1)
        cerr << "Data film 2 " << strerror(errno) << endl;
      else
        cout << "data_film created." << endl;       
    }

    if (mkdir(par.data_file.c_str(), 0777) == -1)
      cerr << "Error : " << strerror(errno) << endl;
    else
      cout << "Directory created." << endl;

  } 
  catch(const char* error) 
  {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);

  }

}

int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

TIMESTEP
{
  cerr << "Error\n";
}


void process_population()
{
  Dish *dishes = new Dish[par.n_orgs];
  omp_set_num_threads(par.n_orgs);

  vector<vector<double>> boundary_lengths(par.n_orgs);

  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  { 
    vector<double> bound_lengths{};
    dishes[i].CPM->set_num(i+1);
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
    dishes[i].CPM->SetAreas(par.cell_target_area);
    if (par.H_perim == true)
    {
      par.H_perim = true; 
      dishes[i].CPM->MeasureCellPerimeters();
      dishes[i].CPM->SetPerims(par.ptarget_perimeter);
    }
    // dish->CPM->ColourCellsByIndex();'
    dishes[i].CPM->StartDynamicAdhesion();
    dishes[i].CPM->SetSortingTypesRandomly();
      
    int t;
    for (t = 0; t < par.mcs; t++)
    {
      if (t==par.highT_time)
      {
        dishes[i].CPM->CopyProb(par.T);
      }

      if (t % 10 == 0)
      {
        dishes[i].CPM->UpdateDynamicAdhesion();
      }

      dishes[i].CPM->AmoebaeMove(t);
      if (par.active_motion)
      {
        dishes[i].CPM->update_cell_velocities_MCS();
      }
      // if (t % 100 == 0)
      // {
      //   cout << dish->CPM->CalculateABBoundaryLength() << endl;
      // }
      if (t%10==0)
      {
        bound_lengths.push_back(dishes[i].CPM->CalculateABBoundaryLength());
      }
      if (t%1000==0)
        cout << i << "  reached time step: " << t << endl;
    }
    boundary_lengths[i]=bound_lengths;
  }
  // do stuff here


  string oname = par.data_file + "/boundary-length.dat";
  ofstream outfile;
  outfile.open(oname, ios::app);  // Append mode
  outfile << fixed << setprecision(3);
  int timelength = boundary_lengths[0].size();
  for (int j = 0; j < timelength; ++j)
  {
    double avg_bound_length{};
    for (int i = 0; i < par.n_orgs; ++i)
    {
      avg_bound_length += boundary_lengths[i][j];
    }
    avg_bound_length /= par.n_orgs;
    outfile << j*10 << '\t' << avg_bound_length << endl;
  }
  outfile.close();

}



int main(int argc, char *argv[]) 
{

  par.pics_for_opt = false;
  #ifdef QTGRAPHICS
    {
      if (par.pics_for_opt)
      {
        QApplication* a = new QApplication(argc, argv);
        // if (mkdir(par.pic_dir.c_str(), 0777) != -1)
        //   cout << "Directory created." << endl;
      }
    }
  #endif


        
    Parameter();
    par.dynamic_sorting=true;
    par.graphics=false;
    par.contours=false;
    process_population();
    

}
