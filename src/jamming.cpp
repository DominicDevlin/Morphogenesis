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
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "info.h"
#include "parameter.h"
#include "sqr.h"
#include "storage.h"
#include "connections.h"
#include "fft.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

#include <sys/stat.h>
#include <cstring>

using namespace std;


void WriteData(const map<int, vector<pair<int, double>>>& shapedata, const string& oname)
{
  ofstream outfile;
  outfile.open(oname, ios::app);  // Append mode

  // First, find the maximum number of rows required
  int max_rows = 0;
  vector<int> rows{};
  for (const auto& [key, vec] : shapedata) {
    for (const auto& [index, value] : vec) 
    {
      if (index + 1 > max_rows) 
      {
          max_rows = index + 1;
          rows.push_back(index);
      }
    }
  }

  // Write the header
  outfile << fixed << setprecision(6);
  
  for (int &row : rows) 
  {
    bool first_col = true;
    
    // Iterate over the map entries
    for (const auto& [key, vec] : shapedata) 
    {
      // Output the first column (the integer index)
      if (!first_col) 
      {
        outfile << "\t";  // Separate columns with a tab
      }

      outfile << row;

      // Calculate the average for this row if there are matching pairs
      double sum = 0.0;
      int count = 0;
      for (const auto& [index, value] : vec) 
      {
        if (index == row) 
        {
          sum += value;
          ++count;
        }
      }

      if (count > 0) 
      {
        double average = sum / count;
        outfile << "\t" << average << '\t' << count;  // Output the average in the second column
      } 
      else 
      {
        // cout << "Error in time output" << endl;
        outfile << "\t" << 0.0 << '\t' << count;  // No data for this row, leave empty
      }

      first_col = false;  // Set this to false after the first column
    }

    outfile << endl;  // Newline after each row
  }

  outfile.close();  
}




INIT 
{
  try 
  {
    
    CPM->set_seed();
    CPM->set_datafile(par.data_file);
    // Define initial distribution of cells

    if (par.make_sheet)
    {
      CPM->ConstructSheet(par.sheetx,par.sheety);
      par.divisions = 6;
    }
    else
      CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);

    

    CPM->ConstructInitCells(*this);
    if (par.velocities)
      par.output_sizes = true;

    // par.divisions = 6;
    if (par.do_voronoi)
    {
      par.highT=false;
      int xtoshift = par.sizex/2 - par.dewet_length/2;
      int ytoshift = par.sizey/2 - par.L2/2;
      // cout << "dewet length: " << par.dewet_length << "  .vertical length: " << par.L2 << endl;
      CPM->Voronoi(par.dewet_length,round(par.L2+5), ytoshift, xtoshift);
    }
    else
    {
      for (int i=0;i<par.divisions;i++) 
      {
        CPM->DivideCells();
      }
    }
    
    // If we have only one big cell and divide it a few times
    // we start with a nice initial clump of cells. 
    // 
    // The behavior can be changed in the parameter file using 
    // parameters n_init_cells, size_init_cells and divisions

    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();

    CPM->start_network(par.start_matrix, par.start_polarity);

    CPM->Set_evoJ(par.J_SL);


    par.print_fitness = true;

    if (par.set_colours)
    {
      CPM->SetColours();
    }


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

TIMESTEP { 
 
  try {
    static int t=0;
 
    static Dish *dish=new Dish();
    
    if (t < 1)
    { 
      cout << "calling init" << endl;
      dish->Init();
      dish->CPM->CopyProb(par.T);

      dish->CPM->SetAreas(par.cell_areas);
      
    }

    // bool GRN = true;

    static Info *info=new Info(*dish, *this);
    // record initial expression state. This occurs before any time step updates. 
    if (t == 100)
    {
      if (par.flush_cells)
      {
        dish->CPM->SetAllStates();
        dish->PDEfield->FlushGrid();
      }
    }

    if (par.record_pressure & t > 10)
    {
      dish->CPM->RecordStress();
    }
    // if (t % 100 == 0)
    // {
    //   dish->CPM->CheckIfCellTouchingMedium();
    // }
    // if (t == par.mcs - 1)
    // {
    //   map<int, bool> medtouches = dish->CPM->ReturnMediumTouching();
    //   ofstream outfile;
    //   string fnamee = par.data_file + "/cell-medium-touchlist.dat";
    //   outfile.open(fnamee, ios::app);  // Append mode
    //   outfile << fixed << setprecision(3);
    //   for (auto &pr : medtouches)
    //   {
    //     outfile << pr.first << '\t' << pr.second << endl;
    //   }
    //   outfile.close();
    // }

    if (par.record_pressure && t % par.pressure_time_length == 0 && t > 100)
    {
      vector<double> pressures = dish->CPM->HydrostaticPressure();
      vector<double> adh_stress = dish->CPM->AdhesionStress();
      double len = double(pressures.size());
      double pressure_var{};
      double pressure_avg = accumulate(pressures.begin(), pressures.end(), 0.0);
      pressure_avg = pressure_avg / len;
      for (double &i : pressures)
      {
        double diff = i - pressure_avg;
        pressure_var += diff * diff;
      }
      pressure_var = pressure_var / len;


      // cout << pressure_avg << '\t' << pressure_var << endl;
      string fnamee = par.data_file + "/pressures.dat";
      ofstream outfile;
      outfile.open(fnamee, ios::app);  // Append mode
      outfile << fixed << setprecision(3);
      for (auto &pr : pressures)
      {
        outfile << pr << '\t';
      }
      outfile << endl;
      outfile.close();
      
      fnamee = par.data_file + "/total_stress.dat";
      outfile.open(fnamee, ios::app);  // Append mode
      outfile << fixed << setprecision(3);
      for (double i = 0; i < len; ++i)
      {
        outfile << pressures[i] + adh_stress[i] << '\t';
      }
      outfile << endl;
      outfile.close();
    }

    if (par.velocities)
    {
      dish->CPM->RecordMasses();
    }
    if (t>10 && t % 10 == 0)
    {
      vector<vector<double>> shared_centres = dish->CPM->find_shared_centres();
      if (shared_centres.size() > 0)
      {
        string fnamee = par.data_file + "/transitions.dat";
        ofstream outfile;
        outfile.open(fnamee, ios::app);  // Append mode
        outfile << fixed << setprecision(3);
        for (auto &vv : shared_centres)
        {
          outfile << t << '\t' << vv[4] << '\t' << vv[5] << '\t' << vv[0] << '\t' << vv[1] <<'\t' << vv[2] <<'\t' << vv[3] <<endl;
        }
        outfile.close();
      }
    }
      

    // static vector<double> cooperativities;

    if (t==0)
    {
      dish->CPM->WetAllCells();
    }

    if (par.measure_time_order_params && t > 1000)
    {
      dish->CPM->PhaseHexaticOrder(t);
      dish->CPM->PhaseShapeIndex(t, true);
    
      // if (t > par.coop_start)
      // {
      //   double coop = dish->CPM->Cooperativity(200);
      //   cooperativities.push_back(coop);
      // }
      
    }    
    dish->CPM->ColourCells(true);
    dish->CPM->AmoebaeMove(t);


    if (t % par.cell_addition_rate == 0 && t > 200 && par.add_cells)
    {
      int cnum = dish->CPM->FindHighestCell();
      int mnum = dish->CPM->TopStalk();
      int count_phase = dish->CPM->CountPhaseOnCells();

      bool set=false;
      if (count_phase == 0)
      {
        set = true;
      }
      while (!set)
      {
        // pair<int,int> val = dish->CPM->ChooseAddPoint(mnum);
        pair<int,int> val = dish->CPM->ChooseAddPoint();
        set = dish->CPM->SpawnCell(val.first, val.second, cnum, t);
      }
        
    }


    if (t == par.mcs - 1)
    {



      if (par.measure_time_order_params)
      {
        map<int, vector<pair<int,double>>> shapedata = dish->CPM->Get_time_shape_index();
        map<int, vector<pair<int,double>>> hexdata = dish->CPM->Get_time_hexatic_order();

        string oname = par.data_file + "/hex_time.dat";
        WriteData(hexdata, oname);

        oname = par.data_file + "/shape_time.dat";
        WriteData(shapedata, oname);

      }

      if (par.output_gamma)
        dish->CPM->OutputGamma();

      if (par.output_sizes)
      {
        dish->CPM->OutputSizes();
        dish->CPM->Vectorfield();
      }

      if (par.velocities)
      {
        dish->CPM->CellVelocities();
      }

   
    }

    //printing every 1000 steps. Do other debugging things here as well. 
    if (t % 1000 == 0)
    {

      cout << "Number of cell types: " << dish->CPM->get_ntypes() << endl;
      cout << t << " TIME STEPS HAVE PASSED." << endl;
    }

    static bool c1 = false;
    static bool c2 = false;
    static bool c3 = false;

    //cerr << "Done\n";
    if (par.graphics && t%5==0)// !(t%par.screen_freq)) 
    {
      
      BeginScene();
      ClearImage();

      // Plot the dish. 
      dish->Plot(this);
      


    
      char title[400];
      snprintf(title,399,"CellularPotts: %.2f hr",dish->PDEfield->TheTime()/3600);      



      //char title[400];
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      EndScene();
      info->Menu();
     
    }
  
    // storage function. 
    if (par.store && !(t%par.storage_stride))//  || t == 3041) 
    {
      char fname[200];
      sprintf(fname,"%s/extend%07d.png",par.datadir,t);
    
      BeginScene();
      ClearImage();    
      dish->Plot(this);

      
      EndScene();
    
      Write(fname);
        
    }

    t++;
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


int main(int argc, char *argv[]) {
  
	try 
  {

#ifdef QTGRAPHICS
    QApplication a(argc, argv);
#endif
    Parameter();
    par.phase_evolution = true;
    // Read parameters
    bool read = false;
    if (read)
      par.Read(argv[1]);
    // Seed(par.rseed);
    
    //QMainWindow mainwindow w;
#ifdef QTGRAPHICS
    QtGraphics g(par.sizex*2,par.sizey*2);
    //a.setMainWidget( &g );
    a.connect(&g, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );

    if (par.graphics)
      g.show();
    
    a.exec();
#else
    X11Graphics g(par.sizex*2,par.sizey*2);
    int t;

    for (t=0;t<par.mcs;t++) {

      g.TimeStep();
    
    }
#endif
    
  } catch(const char* error) {
    std::cerr << error << "\n";
    exit(1);
  }
  catch(...) {
    std::cerr << "An unknown exception was caught\n";
  }
  return 0;
}
