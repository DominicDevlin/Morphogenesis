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




    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);


    CPM->ConstructInitCells(*this);
    if (par.velocities)
      par.output_sizes = true;

    // cout << "dewet length: " << par.dewet_length << "  .vertical length: " << par.L2 << endl;

    // Note - this function will need to have a center of mass somewhere.
    
    CPM->ClearGrid();

    if (par.make_zona_pellucida)
    {
      CPM->MakeZonaPellucida(par.sizex/2, par.sizey/2, 55, 65, 2);
    }

    CPM->PopulateDenseCellsInZonaRadius(par.start_density, par.start_radius, 0, -83, par.sizex/2, par.sizey/2, 55, 65, 2);

    CPM->DifferentiateZonaPellucida();


    // Assign a random type to each of the cells
    // CPM->SetRandomTypes();
    
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

// Draws one progress bar per organism and updates them in place (rather than
// reprinting new lines), the same way tools like pacman/apt render parallel
// progress. `n_orgs` lines are reserved up front by init_progress_bars();
// each call here rewrites only organism org_index's line via ANSI cursor
// movement, so all bars stay pinned to their own row while ticking forward.
void init_progress_bars(int n_orgs, int total)
{
  for (int org_index = 0; org_index < n_orgs; ++org_index)
    cout << "Org " << (org_index + 1) << " [" << string(40, '-') << "] 0% (0/" << total << ")\n";
  cout << flush;
}

void update_progress_bar(int org_index, int step, int total, int n_orgs)
{
  const int bar_width = 40;
  double fraction = double(step) / double(total);
  int filled = static_cast<int>(bar_width * fraction);

  ostringstream bar;
  bar << "Org " << (org_index + 1) << " [";
  for (int i = 0; i < bar_width; ++i)
    bar << (i < filled ? '#' : '-');
  bar << "] " << static_cast<int>(fraction * 100) << "% (" << step << "/" << total << ")";

  int lines_up = n_orgs - org_index;

  #pragma omp critical(progress_print)
  cout << "\033[" << lines_up << "A\r\033[K" << bar.str() << "\033[" << lines_up << "B\r" << flush;
}

void process_population()
{

  

  Dish *dishes = new Dish[par.n_orgs];
  vector<vector<double>> sox2_start(par.n_orgs);
  vector<vector<double>> sox17_start(par.n_orgs);

  ostringstream makefll;

  // init_progress_bars(par.n_orgs, par.mcs);

  // One OS thread per organism oversubscribes the machine once n_orgs
  // exceeds the core count (context-switch/cache-thrashing overhead can
  // dwarf the actual simulation work). Cap at the hardware thread count;
  // OpenMP's parallel-for then queues the remaining organisms onto
  // whichever thread frees up, instead of running them all at once.
  omp_set_num_threads(min(par.n_orgs, omp_get_num_procs()));
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
    int initial_cell_count = dishes[i].CPM->CountCells();
    cout << "Number of cells: " << initial_cell_count << endl; // 1200
    dishes[i].CPM->DrawDivisionTimes();
    dishes[i].CPM->SetColours();

    // No apoptosis mechanism runs before t==initialise_sox_time, so this is
    // also the population cell sorting starts from - lets downstream analysis
    // compute what fraction of the starting population was eliminated.
    string celltype_fname = par.data_file + "/celltypes-org-" + to_string(i + 1) + ".dat";
    ofstream celltype_file(celltype_fname);
    celltype_file << "time\tsox2_high\tsox17_high\tundifferentiated\ttotal\tinitial_count" << endl;

    string death_cause_fname = par.data_file + "/death_causes-org-" + to_string(i + 1) + ".dat";
    ofstream death_cause_file(death_cause_fname);
    death_cause_file << "time"
                      << "\tsox2_high_lonely\tsox2_high_signal"
                      << "\tsox17_high_lonely\tsox17_high_signal"
                      << "\tundifferentiated_lonely\tundifferentiated_signal"
                      << "\ttotal_lonely\ttotal_signal" << endl;
    DeathCounts cumulative_deaths;

    string sorting_fname = par.data_file + "/boundary_lengths-org-" + to_string(i + 1) + ".dat";
    ofstream sorting_file(sorting_fname);
    sorting_file << "time" << "\tsox2sox17\tloserwinner" << endl;

    int t;

    for (t = 0; t < par.mcs; t++)
    {              


      if (t==par.initialise_sox_time)
      {
        
        dishes[i].CPM->InitialiseRandomSoxValues();
        sox2_start[i] = dishes[i].CPM->sox2_values();
        sox17_start[i] = dishes[i].CPM->sox17_values();
        dishes[i].CPM->SetMotilityStrengths();
        dishes[i].CPM->SetPerims(par.ptarget_perimeter);

      }
      if (t>par.initialise_sox_time)
      {
        double tfrac = min(1., double(t-par.initialise_sox_time)/double(par.time_till_full_expression));
        double multiplier = par.sox2binding - par.loser_sox2_adhesion * 0.5;
        // dishes[i].CPM->SetLoserPerimIncrease( multiplier * tfrac );

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
        if (t%250==0)
        {
          dishes[i].CPM->NeighbourBasedApoptosis(i + 1);
        }
      }
      if (t%100==0)
      {
        CellTypeCounts type_counts = dishes[i].CPM->CountCellTypes();
        celltype_file << t << '\t' << type_counts.sox2_high
                      << '\t' << type_counts.sox17_high << '\t'
                      << type_counts.undifferentiated << '\t' << type_counts.total()
                      << '\t' << initial_cell_count << endl;

        DeathCounts death_counts = dishes[i].CPM->CountAndClearDeathEvents();
        cumulative_deaths.Accumulate(death_counts);
        death_cause_file << t
                          << '\t' << cumulative_deaths.sox2_high.lonely << '\t' << cumulative_deaths.sox2_high.signal
                          << '\t' << cumulative_deaths.sox17_high.lonely << '\t' << cumulative_deaths.sox17_high.signal
                          << '\t' << cumulative_deaths.undifferentiated.lonely << '\t' << cumulative_deaths.undifferentiated.signal
                          << '\t' << cumulative_deaths.total_lonely() << '\t' << cumulative_deaths.total_signal() << endl;
        
        double loser_boundary= dishes[i].CPM->LoserWinnerBoundaryLength();
        double sox_boundary = dishes[i].CPM->Sox2Sox17BoundaryLength();
        sorting_file << t << '\t' << loser_boundary << '\t' << sox_boundary << endl;
      }

      dishes[i].CPM->AmoebaeMove(t);

      if (par.active_motion)
      {
        dishes[i].CPM->update_cell_velocities_MCS();
      }

      if (t % 5000 == 0)
      {
        // update_progress_bar(i, t, par.mcs, par.n_orgs);

        if (par.pics_for_opt)
        {
          fft new_org(par.sizex, par.sizey);
          new_org.ImportCPM(dishes[i].get_cpm());
          string foutput = par.pic_dir + "/org-" + to_string(i) + "-" + to_string(t) + ".png";
          new_org.cpmOutput(foutput);
        }
      }

    }

    celltype_file.close();
  }



  // string oname = par.data_file + "/sox_start_values.dat";
  // ofstream outfile;
  // outfile.open(oname);
  // for (int i = 0; i < par.n_orgs; ++i)
  // {    
  //   for (int s=0; s<sox2_start.size(); ++s)
  //   {
  //     outfile << sox2_start[i][s] << '\t' << sox17_start[i][s] << endl;
  //   }
  // }
  // outfile.close();

  // oname = par.data_file + "/sox_end_values.dat";
  // outfile.open(oname);
  // for (int i = 0; i < par.n_orgs; ++i)
  // {
  //   vector<double> sox2 = dishes[i].CPM->sox2_values();
  //   vector<double> sox17 = dishes[i].CPM->sox17_values();
  //   for (int s=0; s<sox2.size(); ++s)
  //   {
  //     outfile << sox2[s] << '\t' << sox17[s] << endl;
  //   }
  // }
  // outfile.close();

  delete[] dishes;

}




int main(int argc, char *argv[]) 
{
  par.pics_for_opt = false;
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
  

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;
  if (argc > 1)
  {
    vector<double> params;
    for (int i = 1; i < argc; ++i)
    {
      params.push_back(stod(argv[i]));
      cout << stod(argv[i]) << " ";
    }
    cout << endl;

    par.data_file = par.data_file + "/" + argv[1] + "-" + argv[2];
    par.apop_threshold=stod(argv[1]);
    par.loser_sox2_adhesion=stod(argv[2]);
    // sx2 L min= -0.1, max=0.7
    // LL  min = -0.7, max=0.7
    // sx17 L min = 0. max = 0.6
    // and we do intervals of 9x9 matrix (or 8 dividers) so..
    double LSX2min=-0.1;
    double LSX2max=0.7;
    double frac = (par.loser_sox2_adhesion - LSX2min) / ( LSX2max - LSX2min);

    double LLmin=-0.7;
    double LLmax=0.7;
    double to_add = (LLmax-LLmin)*frac;
    par.loser_loser_adhesion=LLmin+to_add;
    
    double LSX17min=-0.1;
    double LSX17max=0.6;
    to_add = (LSX17max-LSX17min) * frac;
    par.loser_sox17_adhesion=LSX17min + to_add;

    double ZLmin=0;
    double ZLmax=1.4;
    to_add = (ZLmax-ZLmin) * frac;
    par.Jzona_sticky_loser=ZLmin+to_add;

    cout << "Loser sox2 adhesion: " << par.loser_sox2_adhesion << "\nloser loser adhesion: " << par.loser_loser_adhesion << "\nloser sox17 adhesion: " << par.loser_sox17_adhesion << endl;

    par.loser_sox2_adhesion=par.loser_sox2_adhesion*par.adhesion_multiplier;
    par.loser_loser_adhesion=par.loser_loser_adhesion*par.adhesion_multiplier;
    par.loser_sox17_adhesion=par.loser_sox17_adhesion*par.adhesion_multiplier;
    par.Jzona_sticky_loser=par.Jzona_sticky_loser*par.adhesion_multiplier;

  }



  Parameter();
  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;
  par.highT=false;


  par.end_program=0;
  par.n_orgs = 100;
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

  if (par.pics_for_opt)
    mkdir(par.pic_dir.c_str(), 0777);

  process_population();

  
  // finished
  par.CleanUp();

  return 0;
}
