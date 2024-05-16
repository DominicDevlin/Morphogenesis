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
#include "fft.h"
#include <sys/stat.h>


#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

using namespace std;

auto start = chrono::steady_clock::now();
//rng for making random networks
auto mseed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne( static_cast<mt19937::result_type>(mseed) );
std::uniform_real_distribution<double> double_num(0.0, 1.0);
std::uniform_int_distribution<> genes_dist(0, par.n_genes-1);
std::uniform_int_distribution<> activ_dist(0, par.n_activators-1);
std::uniform_int_distribution<> TF_dist(0, par.n_TF-1);
std::uniform_int_distribution<> J_dist(0, par.J_diff);

int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

INIT {

  try 
  {
    CPM->set_seed();
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);
    CPM->ConstructInitCells(*this);
    CPM->SetRandomTypes();
    
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);

  }
}

TIMESTEP 
{
  cerr << "Error" << endl;
}



void printn(vector<double> &fitn, string &oname, vector<double> &params)
{
  // create and open file
  std::string var_name = oname + "/optimize.txt";
  std::ofstream outfile;
  outfile.open(var_name, ios::app);

  // max fitness 
  double min_fit = SIZE_MAX;
  double max_fit = 0;

  //average fitness
  double avgfit = 0;
  for (double i : fitn)
  {
    avgfit += i;
    cout << i << endl;
    if (i > max_fit)
      max_fit = i;
    if (i < min_fit)
      min_fit = i;
  }
  avgfit = avgfit / par.optimization_replicates;

  //output fitness 
  outfile << min_fit << '\t' << max_fit << '\t' << avgfit << endl;

  outfile.close();

  var_name = oname + "/params.txt";
  outfile.open(var_name, ios::app);

  for (double i : params)
    outfile << i << '\t';
  outfile << endl;
  outfile.close();

}



// function that simulates a population for a single evolutionary step. 
vector<double> process_population(vector<vector<vector<int>>>& network_list, vector<vector<bool>> &pols, vector<double> params)
{
  vector<double> opt_out{};
  vector<int> cell_counter{par.optimization_replicates};
  opt_out.resize(par.optimization_replicates);

  // create memory for dishes. 
  Dish* dishes = new Dish[par.optimization_replicates];
  int time{};
  time = int(params[8]);

  // run organisms in parallel. 
  omp_set_num_threads(par.optimization_replicates);
  #pragma omp parallel for 
  for (int i=0; i < par.optimization_replicates; ++i)  
  {
    dishes[i].CPM->set_num(i+1);
    // does init block above.
    dishes[i].Init();

    int t;

    dishes[i].CPM->start_network(network_list.at(i), pols.at(i));
    // setting optimization params
    // 0 = secretion rate
    par.secr_rate[0] = params[0];
    par.J_med = params[1];
    par.J_stem_diff = params[2];
    dishes[i].CPM->Set_evoJ(par.J_stem_diff);
    par.Vs_max = params[3];
    par.Vd_max = params[4];
    par.gthresh = params[5];
    // constant params
    par.J_stem = params[6];
    par.J_diff = params[7];

    // if (i=0)
    // {
    //   cout << "params are: " << par.secr_rate[0] << '\t' << par.J_med << '\t' << par.J_stem_diff << '\t' << par.Vs_max << '\t' <<
    //   par.Vd_max << '\t' << par.secr_rate[0] << '\t' << par.secr_rate[0] << endl;
    // }

    dishes[i].CPM->CopyProb(par.T);

    cell_counter[i] = dishes[i].CPM->CountCells();

    // run simulation for single organism for mcs montecarlo steps.
    for (t=0;t<par.mcs;t++) 
    {
      // PROGRAMMED CELL DIVISION SECTION
      if (t < par.end_program)
      {

        //programmed divisions
        if (t % par.div_freq == 0 && t <= par.div_end)
        {
          dishes[i].CPM->Programmed_Division(par.phase_evolution); // need to get the number of divisions right. 
        }

        
        if (t >= par.begin_network && t % par.update_freq == 0)
        {
          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell();
          for (int r=0;r<par.program_its;r++) 
          {
            dishes[i].PDEfield->Secrete(dishes[i].CPM);
            dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
          } 
        }
      }
      else 
      {
        // Normal division stage

        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell();
        }
        // diffusion stuff
        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); 
        }  
        dishes[i].CPM->ConstainedGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);
    
      bool check_end = dishes[i].CPM->EndOptimizer();
      if (check_end)
        t = par.mcs;

      // finish simulation if organism is not growing.
      if (t%1000==0 && t > 0)
      {
        int n_cells = dishes[i].CPM->CountCells();
        if (n_cells == cell_counter[i])
        {
          t = par.mcs;
        }
        else
        {
          cell_counter[i] = n_cells;
        }
      }

      // // ensure all cells are connected for shape calculations. 
      // if (t > 0 && t % 1000 == 0)
      // {
      //   bool check_shape = dishes[i].CPM->CheckShape();
      //   if (check_shape == false)
      //   {
      //     opt_out[i] = 100;
      //     t = par.mcs;
      //     // cout << "Org number: " << i << " has bad shape. " << endl;
      //   }
      // }

      // get fitness at end of development
      
    }
    opt_out[i] = dishes[i].CPM->Optimizer();       
          
    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }

  if (par.pics_for_opt)
  {
    string dirn = par.pic_dir + "/" + to_string(time+1);
    if (mkdir(dirn.c_str(), 0777) != -1)
    {
      cout << "Directory created." << endl;
    }
      
    for (int i=0; i < par.optimization_replicates; ++i)
    {
      dishes[i].CPM->ColourCells(par.phase_evolution);
      fft new_org(par.sizex,par.sizey);
      new_org.ImportCPM(dishes[i].get_cpm());
      string f2 = "org-";
      string n2 = to_string(i);
      string ftype = ".png";
      string foutput = dirn + "/" + f2 + n2 + ftype;
      new_org.cpmOutput(foutput);
    }
  }

  delete[] dishes;
  // output to file
  printn(opt_out, par.data_file, params);

  return opt_out;
}




// Main function
int main(int argc, char *argv[]) {

#ifdef QTGRAPHICS
  if (par.evo_pics)
  {
    QApplication* a = new QApplication(argc, argv);
    if (mkdir(par.pic_dir.c_str(), 0777) != -1)
      cout << "Directory created." << endl;
  }
  
#endif


  vector<double> params;
  for (int i = 1; i < argc; ++i)
  {
    params.push_back(stod(argv[i]));
    cout << stod(argv[i]) << endl;
  }


  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record = false;
  par.store = false;
  par.velocities=false;
  par.output_gamma=false;
  par.pickseed = 0;
  par.umap = false;
  par.output_sizes=false;
  par.mcs = 50000;
  par.phase_evolution = true;
  

  Parameter();

  string dirn = par.data_file;
  if (mkdir(dirn.c_str(), 0777) != -1)
  {
    cout << "Directory created." << endl;
    std::string var_name = par.data_file + "/optimize.txt";
    std::ofstream outfile;
    outfile.open(var_name, ios::app);
    outfile << "min\tmax\tavg" << endl;
    outfile.close();
  }
   

  // This is currently depracated. 
  vector<bool> start_p = { 0, 0, 0, 0 };

  vector<vector<vector<int>>> networks{};
  vector<vector<bool>> polarities{};
  for (int i=0;i<par.optimization_replicates;++i)
  {
    networks.push_back(par.start_matrix);
    polarities.push_back(start_p);
  }

  
  // process population. 
  vector<double> fit = process_population(networks, polarities, params);

  // finished
  par.CleanUp();

  return 0;
}
