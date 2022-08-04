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

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif


using namespace std;

INIT 
{
  try 
  {
    CPM->set_seed();
    // Define initial distribution of cells
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
    CPM->ConstructInitCells(*this);
    
    // If we have only one big cell and divide it a few times
    // we start with a nice initial clump of cells. 
    // 
    // The behavior can be changed in the parameter file using 
    // parameters n_init_cells, size_init_cells and divisions
    for (int i=0;i<par.divisions;i++) {
      CPM->DivideCells();
    }
    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();

    CPM->start_network(par.start_matrix, par.start_polarity);

    par.print_fitness = true;
    
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
      dish->CPM->CopyProb(par.eT);
      par.T = par.eT;
      
    }
    if (t == par.end_program)
    {
      dish->CPM->CopyProb(par.lT);
      par.T = par.lT;
        
    }

    static Info *info=new Info(*dish, *this);
    
    if (par.gene_output && t == 100)
      dish->CPM->record_GRN();
    
    // programmed cell division section
    if (t < par.end_program)
    {
      if (t%par.div_freq==0 && t > 0 && t <= par.div_end)
        dish->CPM->Programmed_Division();
      
      if (t >= par.begin_network && t % par.update_freq == 0)
      {
        if (t == par.begin_network && par.morphogen)
          dish->CPM->morphogenWave();

        dish->CPM->update_network(t);
        dish->AverageChemCell(); 
        
        if (par.gene_output)
          dish->CPM->record_GRN();
        
        for (int r=0;r<par.pde_its;r++) 
        {
          dish->PDEfield->Secrete(dish->CPM);
          dish->PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
          // cout << i << endl;
        }
        
      }
    }
    else
    {
      // turn this on for a dance show.
      // dish->CPM->colourbymorph();
      if (t % par.update_freq == 0)
      {
        dish->CPM->update_network(t);
        dish->AverageChemCell(); 
        if (par.gene_output)
          dish->CPM->record_GRN();

      }

      for (int r=0;r<par.pde_its;r++) 
      {
        dish->PDEfield->Secrete(dish->CPM);
        dish->PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        // cout << i << endl;
      }

      // printing individual chemical concentrations. 
      // if (t % 100 == 0 && t > par.end_program)
      // {
      //   dish->PDEfield->print_concentrations(dish->CPM);
      // }

      dish->CPM->CellGrowthAndDivision();
    }
    dish->CPM->AmoebaeMove(t);
    
    if (t == par.mcs-1 && par.gene_output)
    {
      dish->CPM->print_cell_GRN();
    }


    // fitness stuff
    if (t > par.mcs * par.fitness_begin && t % par.fitness_typerate == 0)
    {
      dish->CPM->type_fitness();
      dish->CPM->ShapeFitness();
    }
    
    if (t % par.fitness_somrate == 0 && t > 0)
    {
      dish->CPM->som_fitness();
    }

    if (t == par.mcs - 1)
      dish->CPM->get_fitness();


    //printing every 1000 steps. Do other debugging things here as well. 
    if (t % 1000 == 0)
    {
      bool val = dish->CPM->CheckAllConnected();
      cout << "Are all cells connected?  " << val << endl;
      if (t>0)
        dish->CPM->CalculateShape();

      dish->CPM->print_random_cell();
      cout << "Number of cell types: " << dish->CPM->get_ntypes() << endl;
      cout << t << " TIME STEPS HAVE PASSED." << endl;

    }
  
    //cerr << "Done\n";
    if (par.graphics && !(t%par.screen_freq)) 
    {
      
      BeginScene();
      ClearImage();

      // Plot the dish. 
      dish->Plot(this);

      if (t>par.end_program && par.contours)
      {
        dish->PDEfield->ContourPlot(this,1,7);
        dish->PDEfield->ContourPlot(this,0,5);
        // this function plots a shade for the PDE field and is very computationally costly. Isn't needed. 
        //   dish->PDEfield->Plot(this,0);
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
    if (par.store && !(t%par.storage_stride)) {
      char fname[200];
      sprintf(fname,"%s/extend%07d.png",par.datadir,t);
    
      BeginScene();
      ClearImage();    
      dish->Plot(this);
      if (t>par.end_program && par.contours)
      {
        dish->PDEfield->ContourPlot(this,0,7);
        dish->PDEfield->ContourPlot(this,1,5);
      }
      
      EndScene();
    
      Write(fname);
        
    }

    t++;
  } catch(const char* error) {
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
