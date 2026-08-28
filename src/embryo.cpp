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


    par.make_synthetic=true;
    par.phase_evolution = false;

    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);


    CPM->ConstructInitCells(*this);
    if (par.velocities)
      par.output_sizes = true;

    par.highT=false;
    // cout << "dewet length: " << par.dewet_length << "  .vertical length: " << par.L2 << endl;

    // Note - this function will need to have a center of mass somewhere.
    
    CPM->ClearGrid();

    if (par.make_zona_pellucida)
    {
      CPM->MakeZonaPellucida(par.sizex/2, par.sizey/2, 40, 40, 2);
    }

    CPM->PopulateDenseCellsInZonaRadius(par.start_density, par.start_radius, 0, -70, par.sizex/2, par.sizey/2, 40, 40, 2);

    CPM->DifferentiateZonaPellucida();


    // Assign a random type to each of the cells
    // CPM->SetRandomTypes();

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
      dish->CPM->MeasureCellPerimeters();
      dish->CPM->SetAreas(par.cell_target_area);
      // Sox values must be assigned before SetPerims, which sizes the
      // perimeter constraint per lineage (epiblast/hypoblast).
      dish->CPM->SetPerims(par.ptarget_perimeter);
      cout << "Number of cells: " << dish->CPM->CountCells() << endl; // 1200
      dish->CPM->DrawDivisionTimes();
      dish->CPM->SetSoxColours(0);
    }

    if (t==par.initialise_sox_time)
    {
      
      dish->CPM->InitialiseRandomSoxValues();
      dish->CPM->SetMotilityStrengths();
      dish->CPM->SetPerims(par.ptarget_perimeter);


    }
    if (t>par.initialise_sox_time)
    {
      double tfrac = min(1., double(t-par.expression_starts)/double(par.time_till_full_expression));
      // double multiplier = (par.sox2binding - par.loser_sox2_adhesion) * 0.5;
      if (tfrac < 0)
        tfrac=0;
      dish->CPM->SetLoserPerimIncrease( par.loser_perim_increase * tfrac );
      // cout << par.loser_perim_increase * tfrac << endl;

      double tfrac2 = min(1., double(t-par.sox17bleb_slowdown_start)/double(par.bleb_end));
      if (tfrac2 < 0)
        tfrac2=0;
      tfrac2=1-tfrac2;
      dish->CPM->SetSox17PerimIncrease(par.hypoblast_perim_increase * tfrac2 );
      // if (t%100==0)
      // {
      //   dish->CPM->ToxictoLonelyCells();
      // }
      if (t%10==0)
      {
        dish->CPM->CheckIfDivisionHit(t);
        dish->CPM->LoserActiveMotion(tfrac);
        /*
        The thing is... I actually KNOW what the shape index should be 
        roughly for loser cells. Because of that, what i can do is
        get results of the shape index for loser cells for different adhesion /target perim values.
        From that, we can then deduce what kinds of perimeter/adhesion are roughly correct?
        */
        // dish->CPM->ShapeIndex();
        dish->CPM->SetPerims();
        dish->CPM->SetSoxColours(tfrac);
      }
      if (t%250==0)
      {
        dish->CPM->NeighbourBasedApoptosis();
      }

      // if (t%500==0)
      // {
      //   cout << "loser boundary: " << dish->CPM->LoserWinnerBoundaryLength() << endl;
      //   cout << "sox boundary: " << dish->CPM->Sox2Sox17BoundaryLength() << endl;
      // }
    }

    // if (t==1000)
    // {
    //   dish->CPM->MeanCellArea();
    //   dish->CPM->MeanCellPerimeter();
    // }

    static Info *info=new Info(*dish, *this);
  
    // dish->CPM->ColourCells(true);
    dish->CPM->AmoebaeMove(t);
    if (par.active_motion)
    {
      dish->CPM->update_cell_velocities_MCS();
    }




    // std::this_thread::sleep_for(std::chrono::milliseconds(1));
    // if (t==200)
    // {
    //   std::cout << "Press Enter to continue..."; // Nice to have a prompt
    //   std::cin.get();                            // The actual pause
    // }


    // if (t%100==0)
    // {
    //   cout << "Mean cell Perimenter: " << dish->CPM->MeanCellPerimeter() << endl;
    //   cout << "Mean cell Area: " << dish->CPM->MeanCellArea() << endl;
    // }

    //printing every 1000 steps. Do other debugging things here as well. 
    if (t % 2000 == 0)
    {
      cout << t << " TIME STEPS HAVE PASSED." << endl;
    }

    //cerr << "Done\n";
    if (par.graphics && t%5==0)// !(t%par.screen_freq)) 
    {
      
      BeginScene();
      ClearImage();

      // Plot the dish.
      dish->Plot(this);

#ifdef QTGRAPHICS
      DrawLegend();
#endif

      char title[400];

      //char title[400];
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      if (par.contours && t > 1000)
        dish->PDEfield->ContourPlot(this,0,5);

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

#ifdef QTGRAPHICS
      DrawLegend();
#endif

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
    cout << "here" << endl;
    Parameter();
    // Read parameters
    bool read = false;
    if (read)
      par.Read(argv[1]);
    par.periodic_boundaries=false;
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
