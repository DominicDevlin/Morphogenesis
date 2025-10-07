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
std::uniform_int_distribution<> J_dist(0, par.J_S);

int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

INIT {

  try 
  {
    CPM->set_seed();


    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.sizex/2, par.sizey/2,0,par.offset);
    
    CPM->ConstructInitCells(*this);

    if (par.do_voronoi)
    {
      par.highT=false;
      CPM->Voronoi(par.sizex, par.sheet_depth, par.sheet_shift);
    }

  } 
  catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);

  }
}

TIMESTEP 
{
  cerr << "Error" << endl;
}



void swap(double *xp, double *yp)
{
  double temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void swapv(vector<vector<int>> *xp, vector<vector<int>> *yp)
{
  vector<vector<int>> temp = *xp;
  *xp = *yp;
  *yp = temp;

}

void swapb(vector<bool> *xp, vector<bool> *yp)
{
  vector<bool> temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void swapJ(double *xp, double *yp)
{
  double temp = *xp;
  *xp = *yp;
  *yp = temp;  
}

void swapd(Dish *d, int max_idx, int i)
{
  Dish tmp = d[max_idx];
  d[max_idx] = d[i];
  d[i] = tmp;
}

void sorter(vector<vector<vector<int>>> &networks, vector<double> &fitlist, Dish *dishes)
{
  int i, j, max_idx;
  int n = par.n_orgs;

  for (i = 0; i < n-1; i++)
  {
    max_idx = i;
    for (j=i+1;j<n;j++)
    {
      if (fitlist.at(j) > fitlist.at(max_idx))
        max_idx = j;
    }
    // swap largest element with first element
    swap(&fitlist.at(max_idx), &fitlist.at(i));
    swapv(&networks.at(max_idx), &networks.at(i));
    // std::swap(dishes[max_idx], dishes[i]);

    // swapd(&dishes[max_idx], &dishes[i]);
  }
}


// randomise a new network.  
vector<vector<int>> get_random_network()
{
  vector<vector<int>> matrix;
  matrix.resize(par.n_genes);
  for (int i=0; i < par.n_genes;++i)
  {
    matrix.at(i).resize(par.n_activators);
  }

  for (int i = 0; i < par.n_genes; ++i)
  {
    for (int j = 0; j < par.n_activators; ++j)
    {
      double val = double_num(mersenne);
      if (val < 0.07)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.20)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.7)
      {
        matrix[i][j] = 0;
      }
      else if (val < 0.93)
      {
        matrix[i][j] = 1;
      }
      else
      {
        matrix[i][j] = 2;
      }
    }
  }
  return matrix;
}


double get_random_J()
{
  return double(J_dist(mersenne));
}


// mutate a network. Currently no bias towards ON when mutating networks. 
void mutate(vector<vector<int>> &network)
{
  for (int k=0; k < par.n_mutations; ++k)
  {
    double val = double_num(mersenne);
    int i = genes_dist(mersenne);
    int j = activ_dist(mersenne);
    if (val < 0.05)
    {
      network[i][j] = -2;
    }
    else if (val < 0.20)
    {
      network[i][j] = -1;
    }
    else if (val < 0.7)
    {
      network[i][j] = 0;
    }
    else if (val < 0.93)
    {
      network[i][j] = 1;
    }
    else
    {
      network[i][j] = 2;
    }
  }
}

void mutate_J(double &J)
{
  double val = double_num(mersenne);
  if (val < 0.5)
    J -= 0.5;
  else
    J += 0.5;

  if (J > par.J_S)
    J = par.J_S;
  else if (J < par.J_L)
    J = par.J_L;
}


void output_networks(vector<vector<vector<int>>>& netw)
{
  for (int org=0;org<par.n_orgs;++org)
    for (int i=0;i<par.n_genes;++i)
    {
      if (i == 0)
        cout << "{ ";
      for (int j=0;j<par.n_activators;++j)
      {
        if (j==0)
          cout << "{ " << netw[org][i][j] << ", ";
        else if (j==par.n_activators-1)
          cout << netw[org][i][j] << " }, ";
        else 
          cout << netw[org][i][j] << ", ";
      }
      if (i == par.n_genes -1)
        cout << "}" << endl;
    }
}

void record_networks(vector<vector<vector<int>>>& netw, string oname)
{
  string nname = oname + "/" + "genomes.txt";
  std::ofstream outfile;
  outfile.open(nname, ios::app);
  for (int org=0;org<par.n_orgs;++org)
    for (int i=0;i<par.n_genes;++i)
    {
      if (i == 0)
        outfile << "{ ";
      for (int j=0;j<par.n_activators;++j)
      {
        if (j==0)
          outfile << "{ " << netw[org][i][j] << ", ";
        else if (j==par.n_activators-1)
          outfile << netw[org][i][j] << " }, ";
        else 
          outfile << netw[org][i][j] << ", ";
      }
      if (i == par.n_genes -1)
        outfile << "}" << endl;
    }
  outfile.close();
}



void printn(vector<vector<int>> netw, vector<double> fitn, string oname)
{
  // create and open file
  std::string var_name = oname + "/gene_networks.txt";
  std::ofstream outfile;
  outfile.open(var_name, ios::app);

  for (int i=0;i<par.n_genes;++i)
  {
    if (i == 0)
      outfile << "{ ";
    for (int j=0;j<par.n_activators;++j)
    {
      if (j==0)
        outfile << "{ " << netw[i][j] << ", ";
      else if (j==par.n_activators-1)
        outfile << netw[i][j] << " }, ";
      else 
        outfile << netw[i][j] << ", ";
    }
    if (i == par.n_genes -1)
      outfile << "}" << endl;
  }

  outfile.close();

  // max fitness 
  double max_fit = fitn.front();

  //average fitness
  double avgfit = 0;
  for (double i : fitn)
  {
    avgfit += i;
  }
  avgfit = avgfit / par.n_orgs;


  if (par.asym_only && par.asymmetry_selection && avgfit > par.swap_selection)
  {
    par.asym_only = false;
  }

  //calculate time since begin
  auto end = chrono::steady_clock::now();
  auto diff = end - start;

  //output fitness and time since beginning simulation. 
  var_name = oname + "/fitness.txt";
  outfile.open(var_name, ios::app);
  outfile << max_fit << '\t' << avgfit << '\t' << chrono::duration <double, milli> (diff).count() << endl;

}








// function that simulates a population for a single evolutionary step. 
vector<double> process_population(vector<vector<vector<int>>>& network_list, int time)
{
  vector<double> inter_org_fitness{};
  inter_org_fitness.resize(par.n_orgs);

  // create memory for dishes. 
  Dish* dishes = new Dish[par.n_orgs];

  // run organisms in parallel. 
  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for 
  for (int i=0; i < par.n_orgs; ++i)  
  {
    dishes[i].CPM->set_num(i+1);
    // does init block above.
    dishes[i].Init();
    
    int t;

    dishes[i].PDEfield->SetSecretion(par.secr_rate);
    dishes[i].CPM->Set_evoJ(par.J_SL);
    dishes[i].CPM->SetAreas(par.cell_areas);
    dishes[i].CPM->start_network(network_list[i]);

    bool stayed_together=true;

    // make temperature lower for division section
    dishes[i].CPM->CopyProb(par.T);

    // run simulation for single organism for mcs montecarlo steps.
    for (t=0;t<par.mcs;t++) 
    {
      // PROGRAMMED CELL DIVISION SECTION

      if (t==0 && (par.lambda_perimeter > 0 || par.lambda_perimeter_phase>0))
      {
        // cout << par.cell_addition_rate << '\t' << par.J_med << '\t' << par.lambda_perimeter << endl;
        par.H_perim = true;
        dishes[i].CPM->SetPerims(par.ptarget_perimeter);
        dishes[i].CPM->MeasureCellPerimeters();
      }

      if (t == par.start_topping)
      {
        cout << "IN HERE WITH START TOPPING: " << par.start_topping << endl;
        dishes[i].CPM->ToppingVoronoi(); 
        if (par.MakeEpithelia)
        {
          dishes[i].CPM->AddEpithelialLayer();
        }
        if (t > par.begin_network)
          dishes[i].CPM->StartWettingNetwork();

      }      
      cout << t << endl;

      if (t == par.start_topping + 40)
      {
        dishes[i].CPM->ApoptoseDeadCells();
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
      dishes[i].CPM->AmoebaeMove(t);

      // ensure all cells are connected for shape calculations. 
      if (t > 1 && t % 5000 == 0 && stayed_together==true)
      {
        bool check_shape = dishes[i].CPM->CheckShape();
        if (check_shape == false)
        {
          stayed_together=false;
          t = par.mcs;
        }
        
      }

      if (t % 100 && t > par.begin_network + 1000)
      {
        int n_phase = dishes[i].CPM->CountPhaseOnCells();
        if (n_phase < 3)
        {
          t = par.mcs-1;
          cout << "got here " << endl;
        } 



      }
      if (t >= par.mcs-1)
      {
        int returned_height = dishes[i].CPM->ReturnHeight();
        if (returned_height < 5)
          inter_org_fitness[i] = 0;
        else
          inter_org_fitness[i] = 300 - returned_height;
      }

      // if (par.evo_pics && t % 500 == 0 && t > 0)
      // {
      //   string dirn = par.pic_dir;
      //   if (mkdir(dirn.c_str(), 0777) != -1)
      //   {
      //     cout << "Directory created." << endl;
      //   }
      //   dishes[i].CPM->ColourCells(par.phase_evolution);
      //   fft new_org(par.sizex,par.sizey);
      //   new_org.ImportCPM(dishes[i].get_cpm());
      //   string f2 = "org-";
      //   string n2 = to_string(i);
      //   string ftype = ".png";
      //   string foutput = dirn + "/" + f2 + n2 + "-" + to_string(t) + ftype;
      //   new_org.cpmOutput(foutput);
      // }

    }
        
    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;


  }

  if (par.evo_pics && time % par.pic_gen_interval == 0)
  {
    string dirn = par.pic_dir + "/" + to_string(time+1);
    if (mkdir(dirn.c_str(), 0777) != -1)
      cout << "Directory created." << endl;
    record_networks(network_list, dirn);
    for (int i=0; i < par.n_orgs; ++i)
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

  // do sorting algorithm and return fitness
  sorter(network_list, inter_org_fitness, dishes);

  //output to standard output
  output_networks(network_list);

  // output to file
  printn(network_list.front(), inter_org_fitness, par.data_file);



  vector<vector<vector<int>>> nextgen{};
  int j = 0;
  for (int i=0; i < par.n_orgs;++i)
  {
    // Currently no random networks are added if largest fitness > this
    if (inter_org_fitness.front() > 30 || !par.insert_randoms)
    {
      nextgen.push_back(network_list[j]);

      //mutate network with probability = mut_rate
      double mu = double_num(mersenne);
      if (mu < par.mut_rate)
      {
        mutate(nextgen.back());
      }
      ++j; 
      if (j >= par.n_orgs / 4)
        j=0;
    }
    else
    {
      // the last 1/4 are random networks
      if (i >= (par.n_orgs * 3)/4)
      {
        nextgen.push_back(get_random_network());
      }
      else 
      {
        nextgen.push_back(network_list.at(j));
      }

      //mutate network with probability = 0.5
      double mu = double_num(mersenne);
      if (mu > par.mut_rate)
      {
        mutate(nextgen.back());
      }
      ++j;
      if (j >= par.n_orgs / 4)
        j=0;
    }
  }

  // set nextgen = current pop. 
  for (int i=0;i<par.n_orgs;++i)
  {
    network_list.at(i) = nextgen.at(i);
  }
  return inter_org_fitness;
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
  Parameter();
  par.measure_time_order_params=false;
  par.phase_evolution = true;
  par.MakeEpithelia = false;
  par.begin_network=800;
  par.start_topping=1000;
  par.mcs = 40000;
  par.n_orgs = 60;
  par.do_voronoi = true;
  par.add_cells = false;
  par.ball_radius=54;
  par.sizex=200;
  par.sizey=300;
  par.sheet_depth=95;
  par.sheet_shift=10;
  par.min_phase_cells=4;
  par.evs = 2000;

  string dirn = par.data_file;
  if (mkdir(dirn.c_str(), 0777) != -1)
    cout << "Directory created." << endl;


  // make initial random networks. 
  vector<vector<vector<int>>> networks{};
  for (int i=0;i<par.n_orgs;++i)
  {
    if (par.starter)
    {
      networks.push_back(par.start_n);
    }
    else
    {
      networks.push_back(get_random_network());
    }
  }


  for (int t=0;t<par.evs;++t)
  {
    cout << "current ev timestep is: " << t+1 << endl;
    // process population. 
    vector<double> fit = process_population(networks, t);

    // output every x evolution steps. 
    // if (t%1==0)
    // {
    //   printn(networks.front(), fit);
    // }
  }
  // finished
  par.CleanUp();

  return 0;
}
