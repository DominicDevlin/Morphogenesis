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
#include "omp.h"
#include <chrono>
#include <random>

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
      CPM->MakeZonaPellucida(par.sizex/2, par.sizey/2, 65, 85, 2);
    }

    CPM->PopulateDenseCellsInZonaRadius(par.start_density, par.start_radius, 0, -110, par.sizex/2, par.sizey/2, 65, 85, 2);

    CPM->DifferentiateZonaPellucida();


    // Assign a random type to each of the cells
    // CPM->SetRandomTypes();

    par.end_program=0;

    
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
  cerr << "Error" << endl;
}

void process_population()
{

  

  Dish *dishes = new Dish[par.n_orgs];

  ostringstream makefll;


  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();
    dishes[i].CPM->CopyProb(par.T);
    dishes[i].CPM->MeasureCellPerimeters();
    dishes[i].CPM->SetAreas(par.cell_target_area);

    dishes[i].CPM->SetPerims(par.ptarget_perimeter);
    cout << "Number of cells: " << dishes[i].CPM->CountCells() << endl; // 1200
    dishes[i].CPM->DrawDivisionTimes();
    dishes[i].CPM->SetColours(); 


    int t;

    for (t = 0; t < par.mcs; t++)
    {              


      if (t==par.initialise_sox_time)
      {
        
        dishes[i].CPM->InitialiseRandomSoxValues();
        dishes[i].CPM->SetMotilityStrengths();
        dishes[i].CPM->SetPerims(par.ptarget_perimeter);

      }
      if (t>par.initialise_sox_time)
      {
        double tfrac = min(1., double(t-par.initialise_sox_time)/double(par.time_till_full_expression));
        par.loser_perim_increase = 0.5 * tfrac;

        if (t%100==0)
        {
          dishes[i].CPM->ToxictoLonelyCells();
        }
        if (t%10==0)
        {
          dishes[i].CPM->CheckIfDivisionHit(t);
          dishes[i].CPM->NeighbourBasedActiveMotion(tfrac);
          dishes[i].CPM->SetPerims();
          dishes[i].CPM->SetSoxColours(tfrac);
        }
      }
      dishes[i].CPM->AmoebaeMove(t);

      if (par.active_motion)
      {
        dishes[i].CPM->update_cell_velocities_MCS();
      }

      if (par.pics_for_opt && t % 1000 == 0)
      {
        string dirn = par.pic_dir;
        if (mkdir(dirn.c_str(), 0777) != -1)
        {
          cout << "Directory created." << endl;
        }

        for (int org=0; org < par.n_orgs; ++org)
        {
          // dishes[i].CPM->ColourCells(par.phase_evolution);
          fft new_org(par.sizex,par.sizey);
          new_org.ImportCPM(dishes[org].get_cpm());
          string f2 = "org-";
          string n2 = to_string(org);
          string ftype = ".png";
          string foutput = dirn + "/" + f2 + n2 + "-" + to_string(t) + ftype;
          new_org.cpmOutput(foutput);
        }
      }

    }
  }

    
  //   ostringstream stream;
  //   stream << fixed << setprecision(2) << par.sheet_J; // Setting precision to 2 decimal points
  //   string formatted_value = stream.str();
  //   string oname = par.data_file + "/hex_time-" + formatted_value + ".dat";
  //   ofstream outfile;
  //   outfile.open(oname, ios::app);  // Append mode
  //   outfile << fixed << setprecision(3);
  //   for (int i = 0; i < container_size; ++i)
  //   {
  //     outfile << i*par.struct_avg_interval + 2*par.struct_avg_interval << '\t' << hex_order_output[i] << endl;
  //   }
  //   outfile.close();

  // } 


  delete[] dishes;

}




int main(int argc, char *argv[]) 
{
  par.pics_for_opt = true;
  par.pic_dir = "photos";

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
  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;


  par.end_program=0;
  par.n_orgs = 5;
  par.make_synthetic=true;
  par.phase_evolution = false;

  // par.velocities = true;
  //   par.output_sizes = true;
  // else
  //   par.output_sizes = false;
  
  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;



  process_population();
  
  // finished
  par.CleanUp();

  return 0;
}
