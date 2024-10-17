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


void Outputter(map<int, vector<double>> data2, vector<vector<int>> scc, string switch_out)
{
  ofstream outfile;
  outfile.open(switch_out, ios::app);
  vector<vector<double>> vec2(scc.size());


  for (const auto& pair : data2 )
  {
    int count = 0;
    for (auto &j : scc)
    {
      int type = pair.first;
      auto it = std::find(j.begin(), j.end(), type);
      if (it != j.end())
      {
        for (double val : pair.second)
        {
          vec2[count].push_back(val);
        }
      }

      ++count;
    }
  }

  size_t max_size = 0;
  for (const auto& inner_vec2 : vec2) {
      if (inner_vec2.size() > max_size) 
      {
          max_size = inner_vec2.size();
      }
  }
  for (size_t i=0;i<vec2.size();++i)
  {
    outfile << "SCC-"+to_string(i+1) << '\t';
  }
  outfile << endl;

  for (size_t i = 0; i < max_size; i++) 
  {
      for (size_t j = 0; j < vec2.size(); j++) 
      {
          if (i < vec2[j].size()) 
          {
              outfile << vec2[j][i];
          } 
          else 
          {
              outfile << "NaN"; // Print extra spaces for alignment if no element exists
          }
          outfile << "\t";
      }
      outfile << std::endl; // Newline after each column is printed
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
      CPM->Voronoi(par.sizex,par.sheet_depth, par.sheet_shift);
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

    CPM->Set_evoJ(par.J_stem_diff);


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

    bool GRN = true;

    static Info *info=new Info(*dish, *this);
    // record initial expression state. This occurs before any time step updates. 
    if (t == 100)
    {
      if (par.flush_cells)
      {
        dish->CPM->SetAllStates();
        dish->PDEfield->FlushGrid();
      }
      if (par.gene_output)
        dish->CPM->record_GRN();

      if (par.output_init_concs)
        dish->CPM->OutputInitConcs();

    }

    // static vector<double> cooperativities;

    if (t == 1000)
    {
      dish->CPM->WetTopCells(par.dewet_length, par.dewet_cell_depth);
      // dish->CPM->WetRandomCells();
    } 
    if (par.measure_time_order_params && t > 1000)
    {
      dish->CPM->RecordMasses();
      dish->CPM->PhaseHexaticOrder(t);
      dish->CPM->PhaseShapeIndex(t);
    
      // if (t > par.coop_start)
      // {
      //   double coop = dish->CPM->Cooperativity(200);
      //   cooperativities.push_back(coop);
      // }
      if (t % 1000 == 0)
      {
        cout << "Wetting length is... " << dish->CPM->WettingLength() << endl;
      }
      

    }

    
    // if (t % 1000 == 0 && t > 0)
    // {
    //   dish->CPM->RemoveUnconnectedCells();
    //   // dish->CPM->MeasureCellPerimeters();
    // }
    
    if ((t == 0) && par.lambda_perimeter > 0)
    {
      par.H_perim = true; 
      dish->CPM->SetPerims(par.ptarget_perimeter);
      dish->CPM->MeasureCellPerimeters();
    }

    if (par.velocities)
    {
      dish->CPM->RecordMasses();
    }



    // if (t == 22000)
    // {
    //   par.J_stem=1.5; 
    //   par.J_diff = 6;
    //   par.J_med = 4.25;
    //   par.J_med2 = 4.25;
    //   par.J_stem_diff=5.6;
    //   par.gthresh=3;
    //   par.secr_rate[0]=0.00274;
    //   dish->CPM->Set_evoJ(par.J_stem_diff);
    // }
    // if (t>22000)
    // {
    //   dish->CPM->CellGrowthAndDivision(t);
    // }

    if (GRN && t >= 20000)
    {
      if (t==20000)
      {
        par.secr_rate[0]=0.0023;
        dish->CPM->StartWettingNetwork();
      }
        
      if (t % par.update_freq == 0)
      {
        dish->CPM->update_phase_network(t);
        // if (par.noise && t > par.noise_start)
          // dish->CPM->add_noise();
          
        dish->AverageChemCell(); 
        if (par.gene_output)
        {
          dish->CPM->record_GRN();
          dish->CPM->CountTypesTime();
        }
      }
      for (int r=0;r<par.pde_its;r++) 
      {
        if (!par.hold_morph_constant)
        {
          dish->PDEfield->Secrete(dish->CPM);
          dish->PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        }

      }
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
    if (t % 4000 == 0)
    {
      bool check_shape = dish->CPM->CheckAllConnected(0.9);
    }


    if (t == par.mcs - 1)
    {

      if (mkdir(par.data_file.c_str(), 0777) == -1)
        cerr << "Error : " << strerror(errno) << endl;
      else
        cout << "Directory created." << endl;

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
        // dish->CPM->scc_momenta(scc);
        // dish->CPM->momenta();
        // dish->CPM->diff_anisotropy(scc);
      }

   
    }

    //printing every 1000 steps. Do other debugging things here as well. 
    if (t % 1000 == 0)
    {

      cout << "Number of cell types: " << dish->CPM->get_ntypes() << endl;
      cout << t << " TIME STEPS HAVE PASSED." << endl;
        // dish->CPM->PrintPhenotypes();
        // // dish->CPM->WhiteSpace();
        // // dish->CPM->DeviationFromCircle();
        // cout << "ORG MASS IS: " << dish->CPM->Mass() << endl;
        // double center[] = {0.0,0.0};
        // dish->CPM->get_center(center);
        // cout << "x center: " << center[0] << "   y center: " << center[1] << endl;

        // dish->CPM->PrintColours();

        // cout << "DISTANCE TO TOP: " << dish->CPM->Optimizer() << endl;

      // dish->CPM->ShapeAlignmentByPhase();
    }


    if (t >= 6000 && t < 8000 && t % 40 == 0 && par.scramble)
      dish->CPM->swap_cells();


    // for spawning a lot of morphogen at a specific point, or changing cell types + morphogen

    if (par.convert_cells && par.convert_time == t)
    {
      dish->CPM->ConvertToStem(par.convert_x,par.convert_y,par.convert_size,par.convert_to_type, dish->PDEfield, true, par.clear_radius);
      dish->CPM->ConvertToStem(100,230,par.convert_size,par.convert_to_type, dish->PDEfield, true, par.clear_radius);  
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
      
      // static vector<array<int,2>> perim;
      // static vector<int> pcells;
      // if (t > 1000 && t % 500 == 0)
      // {
      //   // perim = dish->CPM->PerimeterCAC();
      //   // pcells = dish->CPM->CellsFromCAC(perim);
      //   // pcells = dish->CPM->LinkPerimeter();
      //   // dish->CPM->AngleCurvature();
      // }
      
      // if (t > 1500)
      // {
      //   // dish->CPM->DrawListofCAC(this, perim);
      //   // dish->CPM->DrawPerimeter(this, pcells);
      // }


      // static bool c4 = false;

      if (t>0 && t % par.begin_network == 0)
      {
        c1 = dish->PDEfield->CheckSecreting(0);
        if (par.n_diffusers > 1)
          c2 = dish->PDEfield->CheckSecreting(1);
        if (par.n_diffusers > 2)
        {
          c3 = dish->PDEfield->CheckSecreting(2);
          // c4 = dish->PDEfield->CheckSecreting(3);
        }

      }

      if (t>par.end_program && par.contours)
      {
        if (c1)
          dish->PDEfield->ContourPlot(this,0,5);
        if (c2)
          dish->PDEfield->ContourPlot(this,1,7);
        if (par.n_diffusers > 2)
        {
          if (c3)
            dish->PDEfield->ContourPlot(this,2,14);
          // if (c4)
          //   dish->PDEfield->ContourPlot(this,3,37);
        }


        // this function plots a shade for the PDE field and is very computationally costly. Isn't needed. 
        // dish->PDEfield->Plot(this,0);
        // dish->PDEfield->Plot(this,1);
        // dish->PDEfield->Plot(this,2);
      }


    
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

      if (t>par.end_program && par.contours)
      {
        c1 = dish->PDEfield->CheckSecreting(0);
        if (par.n_diffusers > 1)
          c2 = dish->PDEfield->CheckSecreting(1);
        if (par.n_diffusers > 2)
        {
          c3 = dish->PDEfield->CheckSecreting(2);
          // c4 = dish->PDEfield->CheckSecreting(3);
        }
        if (c1)
          dish->PDEfield->ContourPlot(this,0,293);
        if (c2)
          dish->PDEfield->ContourPlot(this,1,291);
        if (c3)
          dish->PDEfield->ContourPlot(this,2,292);
      }
      
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
