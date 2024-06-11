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

    if (par.velocities)
      par.output_sizes = true;
    
    // for (int i=0;i<par.divisions;i++) {
    //   CPM->DivideCells();
    // }


    CPM->FractureSheet();
    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();




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

double findMedian(std::vector<double>& nums) 
{
    int n = nums.size();
    std::nth_element(nums.begin(), nums.begin() + n / 2, nums.end());
    if (n % 2 != 0) {
        return static_cast<double>(nums[n / 2]);
    } 
    else 
    {
        double median1 = nums[n / 2];
        std::nth_element(nums.begin(), nums.begin() + n / 2 - 1, nums.end());
        double median2 = nums[n / 2 - 1];
        return (static_cast<double>(median1) + static_cast<double>(median2)) / 2.0;
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
      
      // equilibriate cells with high T
      if (par.highT)
      {
        dish->CPM->CopyProb(par.highT_temp);
      }
      else
      {
        dish->CPM->CopyProb(par.T);
      }
      dish->CPM->Set_J(par.sheet_J);
        
    }
    
    if (t==par.highT_time)
      dish->CPM->CopyProb(par.T);


    static Info *info=new Info(*dish, *this);
    
    // static vector<double> N_index;
    static vector<double> shape_index;
    if (t % 1000 == 0 && t > 0)
    {
      // N_index.clear();
      // shape_index.clear();
      cout << (par.periodic_boundaries) << endl;
      dish->CPM->initVolume();
      dish->CPM->adjustPerimeters();
      vector<double> tperims = dish->CPM->TruePerimeters();
      vector<double> volumes = dish->CPM->GetVolumes();
      // vector<double> Nperims = dish->CPM->PerimitersRadiusN(sqrt(5), 11);
      // cout << tperims[1] << '\t' << Nperims[1] << '\t' << volumes[1] << endl;

      double avg{};
      // double n_avg{};
      for (int i = 0; i < tperims.size(); ++i)
      {
        // cout << tperims[i] << '\t' << Nperims[i] << '\t' << volumes[i] << endl;
        double sindex = tperims[i] / sqrt(volumes[i]);
        // cout << i << '\t';
        avg+=sindex;
        shape_index.push_back(sindex);

        // sindex = (Nperims[i]) / (sqrt(volumes[i]));
        // cout << i << '\t';
        // n_avg+=sindex;
        // N_index.push_back(sindex);

      }
      double median = findMedian(shape_index);
      avg/=tperims.size();
      // n_avg/=tperims.size();
      cout << endl << avg << '\t' << median << endl;

  
      cout << t << " TIME STEPS HAVE PASSED." << endl;
      cout << "ORG MASS IS: " << dish->CPM->Mass() << endl;
    }
      

    if (par.velocities)
    {
      dish->CPM->RecordMasses();
    }

    if (par.output_sizes)
    {
      dish->CPM->RecordSizes();
    }
    dish->CPM->AmoebaeMove(t);

    if (t == par.mcs-1 && par.gene_output)
    {

      if (mkdir(par.data_file.c_str(), 0777) == -1)
        cerr << "Error : " << strerror(errno) << endl;
      else
        cout << "Directory created." << endl;  
    }
    
    if (t == par.mcs - 1)
    {

      
      
      if (par.output_gamma)
        dish->CPM->OutputGamma();

      if (par.output_sizes)
      {
        dish->CPM->OutputSizes();
        dish->CPM->MeanSquareDisplacement();
      }
        
      
      if (par.velocities)
      {
        // dish->CPM->CellVelocities();
        // dish->CPM->scc_momenta(scc);
        // dish->CPM->momenta();
        // dish->CPM->diff_anisotropy(scc);
      }
        
      // dish->CPM->SpecialVelocity();
      if (par.record_directions)
      {
        // dish->CPM->Directionality();
        // dish->CPM->SingleCellDirection();
      }
   

    }
   

    int freq = 10;
    // par.n_screen_freq=1000;
    if (par.draw_paths)
      freq = 100;
    
    //cerr << "Done\n";
    if (par.graphics && t%freq==0)// !(t%par.screen_freq)) 
    {
      dish->CPM->ColourCellsByIndex();
      
      BeginScene();
      ClearImage();

      // Plot the dish. 
      if (par.draw_paths && t > 600)
        dish->CPM->DrawDisplacement(this);
      else 
        dish->Plot(this);
        
    
    
      static bool c1 = false;
      static bool c2 = false;
      static bool c3 = false;
      // static bool c4 = false;


    
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

      // Plot the dish. 
      if (par.draw_paths && t > 600)
        dish->CPM->DrawDisplacement(this);
      else 
        dish->Plot(this);

      if (t>par.end_program && par.contours)
      {
        
        dish->PDEfield->ContourPlot(this,0,293);
        dish->PDEfield->ContourPlot(this,2,292);
        dish->PDEfield->ContourPlot(this,1,291);
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
    par.sizex=150;
    par.sizey=150;
    par.end_program=0;
    par.periodic_boundaries = true;
    par.flush_cells = true;

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
