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
#include "fft.h"
#include <sys/stat.h>

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

#include <sys/stat.h>
#include <cstring>

using namespace std;

void OutputColumnData(vector<vector<double>> &odata, string fname)
{
    // Open file for writing
    ofstream outputFile;
    outputFile.open(fname, ios::app);

    size_t max_inner_size = 0;
    for (const auto& vec : odata) {
        if (vec.size() > max_inner_size) 
        {
            max_inner_size = vec.size();
            // cout << "m_inner size: " << max_inner_size << endl;
        }
    }


    // Output data as columns where each inner vector corresponds to a column
    for (size_t i = 0; i < max_inner_size; ++i) 
    {
        // Write the row index as the first column
        outputFile << i;

        // Write the corresponding element from each inner vector
        for (size_t j = 0; j < odata.size(); ++j) {
            if (i < odata[j].size()) {
                outputFile << "\t" << odata[j][i];
            } else {
                outputFile << "\t" << 0;  // If the inner vector is shorter, leave an empty space
            }
        }

        // Newline at the end of the row
        outputFile << "\n";
    }

    outputFile.close();      
}

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

    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);

    CPM->ConstructInitCells(*this);
    if (par.velocities)
      par.output_sizes = true;

    if (par.do_voronoi)
    {
      par.highT=false;
      CPM->Voronoi(par.sizex, par.sheet_depth, par.sheet_shift);
    }
    
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

void process_population(vector<vector<vector<int>>>& network_list, int argn=0)
{
  if (argn > 0)
  {
    par.data_file = "org-data-" + to_string(argn);
    par.pic_dir = "images-" + to_string(argn);
  }

  Dish *dishes = new Dish[par.n_orgs];

  int n_times_apart{};

  vector<vector<double>> contact_angles(par.n_orgs);
  vector<int> starting_heights(par.n_orgs);

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();
    dishes[i].PDEfield->SetSecretion(par.secr_rate);
    dishes[i].CPM->Set_evoJ(par.J_stem_diff);
    dishes[i].CPM->SetAreas(par.cell_areas);
    dishes[i].CPM->start_network(network_list[i]);

    bool stayed_together=true;

    // equilibriate cells with high T
    dishes[i].CPM->CopyProb(par.T);
    int t=0;

    for (; t < par.mcs; t++)
    {  
      if (t==0 && (par.lambda_perimeter > 0 || par.lambda_perimeter_phase>0))
      {
        // cout << par.cell_addition_rate << '\t' << par.J_med << '\t' << par.lambda_perimeter << endl;
        par.H_perim = true;
        dishes[i].CPM->SetPerims(par.ptarget_perimeter);
        dishes[i].CPM->MeasureCellPerimeters();
      }

      if (t == par.start_topping)
      {
        dishes[i].CPM->ToppingVoronoi(); 
        if (par.MakeEpithelia)
        {
          dishes[i].CPM->AddEpithelialLayer();
        }
      }      
      if (t == 1000)
      {
        starting_heights[i] = dishes[i].CPM->ReturnHeight(); 
      }

      if (par.linear_increase)
      {
        dishes[i].PDEfield->increase_secretion(t);
      }

      if (t>= par.begin_network)
      {
        if (t == par.begin_network)
        {
          dishes[i].CPM->StartWettingNetwork();
        }

        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell();  
        }
        // speed up initial PDE diffusion
        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        } 
      }

      if (par.velocities && t % 1 == 0)
      {
        dishes[i].CPM->RecordMasses(true);
      }

      if (t > par.start_topping)
      {
        dishes[i].CPM->SetCellCenters();
        int check{};
        if (par.MakeEpithelia)
          check = dishes[i].CPM->EpiContactAngle();
        else
          check = dishes[i].CPM->ContactAngle();
        if (check < 0)
        {
          bool check_shape = dishes[i].CPM->CheckShape();
          if (check_shape == false)
          {
            ++n_times_apart;
            stayed_together=false;
            t = par.mcs;
          }
        }
      }
        

      if (t > par.start_topping && t % par.measure_interval == 0)
      {
        double contacta = dishes[i].CPM->GetContactAngles();
        contact_angles[i].push_back(contacta);
      }

      

      // dishes[i].CPM->DiscreteGrowthAndDivision(t);
      if (t % par.cell_addition_rate == 0 && t > 200 && par.add_cells)
      {
        int cnum = dishes[i].CPM->FindHighestCell();
        int mnum = dishes[i].CPM->TopStalk();
        int check_n_points = dishes[i].CPM->CheckAddPoints();
        int counter = 0;
        bool set=false;
        if (check_n_points < 60)
        {
          t = par.mcs;            
        }
        else
        {
          while (!set)
          {
            pair<int,int> val = dishes[i].CPM->ChooseAddPoint(mnum);
            set = dishes[i].CPM->SpawnCell(val.first, val.second, cnum, t);
            ++counter;
            if (counter > 4)
            {
              set = true;
              cerr << "Error in spawn cell with phase count " << dishes[i].CPM->CountPhaseOnCells() << endl;
              t=par.mcs;
            }
          }
        }
      }       

      dishes[i].CPM->AmoebaeMove(t);


      // ensure all cells are connected for shape calculations. 
      if (t > 0 && t % 5000 == 0 && stayed_together==true)
      {
        bool check_shape = dishes[i].CPM->CheckShape();
        if (check_shape == false)
        {
          ++n_times_apart;
          stayed_together=false;
          t = par.mcs;
        }
        
      }

      if (t % 10 && t > 1000)
      {
        int n_phase = dishes[i].CPM->CountPhaseOnCells();
        if (n_phase < 5)
          t = par.mcs;
      }
      


      if (par.pics_for_opt && t % 500 == 0)
      {
        string dirn = par.pic_dir;
        if (mkdir(dirn.c_str(), 0777) != -1)
        {
          cout << "Directory created." << endl;
        }

        for (int org=0; org < par.n_orgs; ++org)
        {
          dishes[i].CPM->ColourCells(par.phase_evolution);
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

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;

  ostringstream stream;
  stream << fixed << setprecision(2) << par.gamma_hm << '-' << par.gamma_hl; // Setting precision to 2 decimal points
  string formatted_value = stream.str();
  double sum_heights;
  double square_heights;
  vector<double> heights(par.n_orgs);
  double avg_phase_remained = 0;
  for (int i=0; i < par.n_orgs;++i)
  {
    int n_phase = dishes[i].CPM->CountPhaseOnCells();
    avg_phase_remained += n_phase;
    heights[i] = double(starting_heights[i] - dishes[i].CPM->ReturnHeight());
    // cout << "start: " << starting_heights[i] << "\tend: " << dishes[i].CPM->ReturnHeight() <<endl;
  }
  sort(heights.begin(), heights.end(), std::greater<double>());
  vector<double> largest_heights(heights.begin(), heights.begin() + heights.size() / 2);
  int lsize = largest_heights.size();
  for (int i=0; i < lsize;++i)
  {
    sum_heights += largest_heights[i];
    square_heights += largest_heights[i] * largest_heights[i];

  }
  double mean_height = sum_heights / lsize;
  double variance = (square_heights / lsize ) - (mean_height * mean_height);
  avg_phase_remained = avg_phase_remained / lsize;
  ofstream outfile;
  string infoname = par.data_file + "/info.txt";
  outfile.open(infoname, ios::app);  // Append mode
  outfile << par.gamma_hm << '\t' << par.gamma_hl << '\t' << mean_height << '\t' << variance << '\t' << double(n_times_apart) / double(par.n_orgs) << '\t' << avg_phase_remained << endl;
  outfile.close();

  
  // double avg_empty_space = std::accumulate(empty_spaces.begin() + start, empty_spaces.begin() + half, 0.0) / half;


  string fname = par.data_file + "/contact.angle-" + formatted_value + ".dat";
  OutputColumnData(contact_angles, fname);


  delete[] dishes;

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

  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record=false;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  par.velocities=false;
  par.output_sizes = false;
  Parameter();
  par.measure_time_order_params=false;
  
  par.begin_network=2000;
  par.phase_evolution = true;
  par.min_phase_cells=4;
  par.mcs = 40000;
  par.sheet_hex=false;
  par.n_orgs = 120;
  par.do_voronoi = true;
  par.add_cells = false;

  par.measure_interval = 20;


  par.sizex=200;
  par.sizey=300;

  // typical wetting parameters used:
  par.sheet_depth=95;
  par.sheet_shift=10;
  
  bool perimeter_model = false;

  vector<vector<vector<int>>> networks{};
  for (int i = 0; i < par.n_orgs; ++i)
  {
    networks.push_back(par.start_matrix);
  }

  par.gamma_hm = 7;
  par.gamma_hl = 7;

  while (par.gamma_hm < 14.1)
  {
    while (par.gamma_hl < 14.1)
    {
      par.J_stem = 2;
      par.J_med = par.gamma_hm + 1;
      par.J_med2 = par.J_med;
      par.J_stem_diff = 1.75 + par.gamma_hm + par.gamma_hl;
      par.J_diff = 2 * par.gamma_hm + 1.5;
      if (par.MakeEpithelia)
      {
        if (par.gamma_hm < 7 || par.gamma_hl < 7)
        {
          cerr << "wrong gamma vals for epi\n" << endl;
          break;
        }
        par.epiJ=2;
        par.epiJelse = par.gamma_hm - 2;
        par.J_stem_diff = (par.gamma_hm/2.) + par.gamma_hl + 1 - (3.25/2);
        par.J_diff = par.gamma_hm - 3.25;
      }


      process_population(networks);
      par.gamma_hl += 0.5;

    }
    par.gamma_hl = 1.;
    if (par.MakeEpithelia)
      par.gamma_hl = 7.;
    par.gamma_hm += 0.5;
  }

  
  // finished
  par.CleanUp();

  return 0;
}
