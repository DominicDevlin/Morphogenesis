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
      // CPM->VoronoiSeparated(par.dewet_length,round(par.L2+5), ytoshift, xtoshift);
      CPM->PopulateSparseCells(0.6, 300, 0, 0);
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

    par.end_program=0;

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
      dish->CPM->SetAreas(par.cell_target_area);
      dish->CPM->MeasureCellPerimeters();
      dish->CPM->WetAllCells();
      cout << "Number of cells: " << dish->CPM->CountCells() << endl; // 1200

    }

    // bool GRN = true;

    static Info *info=new Info(*dish, *this);

    if (par.velocities)
    {
      dish->CPM->RecordMasses();
    }
  
    dish->CPM->ColourCells(true);
    dish->CPM->AmoebaeMove(t);
    if (par.active_motion)
    {
      dish->CPM->update_cell_velocities_MCS();
    }

    // we will use this later
    // pair<int,int> val = dish->CPM->ChooseAddPoint();
    // dish->CPM->SpawnCell(val.first, val.second, cnum, t);    

    if (t == par.mcs - 1)
    {
      if (par.velocities)
      {
        dish->CPM->CellVelocities();
      }
    }

    // if (t==1000)
    // {
    //   par.medium_area_constraint=true;
    //   dish->CPM->SetMediumArea();
    // }

    if (t%1==0)
    {
      double hh = dish->CPM->SumEnergy();
      ofstream outfile;
      string oname = par.data_file + "/counter.dat";
      outfile.open(oname, ios::app);  // Append mode
      outfile << double(par.tmpcounter) / double(par.tmpcountertotal) << endl;


    }

    if (t%100==0)
    {
      cout << "Mean cell Perimenter: " << dish->CPM->MeanCellPerimeter() << endl;
      cout << "Mean cell Area: " << dish->CPM->MeanCellArea() << endl;
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



int main(int argc, char *argv[]) {
  
	try 
  {

#ifdef QTGRAPHICS
    QApplication a(argc, argv);
#endif
    Parameter();
    par.phase_evolution = true;

    //active term params
    par.active_motion = true;
    par.motility_strength = 1.5;
    par.persistence_time = 100.;
    if (par.active_motion)
    {
      par.dewet_length=450;
      par.dewet_cell_depth=47; // 47 is max
      par.cell_target_area = 100;
      par.L2 = sqrt((sqrt(3.)/2) * par.cell_target_area ) * par.dewet_cell_depth;   
      par.H_perim = true;
      par.ptarget_perimeter = 120;
      par.J_L = 0;
      par.lambda2 = 0;
      par.lambda_perimeter_phase = 0.2;
      par.J_med = 0;
      par.J_med2 = par.J_med;
      par.tmpcounter=0;
      par.tmpcountertotal=0;
      par.medium_area_constraint = false;
    }

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
