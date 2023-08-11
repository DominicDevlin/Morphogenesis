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
#include <sstream>
#include <sys/stat.h>
#include "fft.h"


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


int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

INIT {

  try 
  {
    CPM->set_seed();
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
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




// function that simulates a population for a single evolutionary step. 
vector<double> process_population(vector<vector<vector<int>>>& network_list, vector<vector<bool>> &pols)
{
  vector<double> inter_org_fitness{};
  inter_org_fitness.resize(par.n_orgs);

  // create memory for dishes. 
  Dish* dishes = new Dish[par.n_orgs];

  for (int i=0; i < par.n_orgs; ++i)  
  {
    dishes[i].CPM->set_num(i+1);
    // does init block above.
    dishes[i].Init();
  }


  // run organisms in parallel. 
  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for 
  for (int i=0; i < par.n_orgs; ++i)  
  {
    int t;

    dishes[i].CPM->start_network(network_list.at(i), pols.at(i));


    // make temperature lower for division section
    dishes[i].CPM->CopyProb(par.eT);
    par.T = par.eT;

    // run simulation for single organism for mcs montecarlo steps.
    for (t=0;t<par.mcs;t++) 
    {
      if (t == par.end_program)
      {
        dishes[i].CPM->CopyProb(par.lT); // normal temperature for normal development timing. 
        par.T = par.lT;
      } 
      // PROGRAMMED CELL DIVISION SECTION
      if (t < par.end_program)
      {

        //programmed divisions
        if (t % par.div_freq == 0 && t <= par.div_end)
        {
          dishes[i].CPM->Programmed_Division(); // need to get the number of divisions right. 
        }

        
        if (t >= par.begin_network && t % par.update_freq == 0)
        {
          dishes[i].CPM->update_network(t);
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
          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell();
        }
        // diffusion stuff
        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); 
        }  
        dishes[i].CPM->CellGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);
    

      // calculate the diversity over last 20% of time steps. 
      if (t > par.mcs * par.fitness_begin && t % par.fitness_typerate == 0)
      {
        // am now doing for curvature as well (taking mean)
        dishes[i].CPM->update_fitness();
      }
 
      // ensure all cells are connected for shape calculations. 
      if (t > 0 && t % 1000 == 0)
      {
        bool check_shape = dishes[i].CPM->CheckShape();
        if (check_shape == false)
        {
          inter_org_fitness[i] = 0;
          t = par.mcs;
          // cout << "Org number: " << i << " has bad shape. " << endl;
        }
      }
      // get fitness at end of development
      if (t == par.mcs-1)
      {
        inter_org_fitness[i] = dishes[i].CPM->get_fitness();
      }        
    }
        
    if (i == 1)
      cout << "Sim complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }


  

  vector<double> proportions{};
  // Compare here
  for (int i = 0;i<par.n_orgs;i++)
  {
    for (int j = i+1; j < par.n_orgs;j++)
    {
        double prp = dishes[i].CPM->CompareGrid(dishes[j].CPM->ReturnGrid());
        // cout << "PROPORTION IS: " << prp << endl;
        proportions.push_back(prp);
    }
  }
  
  double avgp{};
  double variance{};

  for (double i : proportions)
  {
    avgp += i;
  }
  avgp = avgp / proportions.size();

  for (double i : proportions)
  {
    double val = pow(i - avgp, 2);
    variance += val;
  }
  variance = variance / proportions.size();


  ofstream outfile;
  string var_name = "overlap.txt";
  outfile.open(var_name, ios::app);

  cout << "average proportion is: " << avgp << endl;
  cout << "variance is: " << variance << endl;

  outfile << avgp << '\t' << variance << endl;
  outfile.close();



  // Now we are going to add rotationally invariant version.

  vector<double> invariant_p{};
  invariant_p.resize(((par.n_orgs-1)*par.n_orgs)/2);

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for 
  for (int i = 0;i<par.n_orgs;i++)
  {
    fft org1;
    org1.AllocateGrid(par.sizex, par.sizey);
    org1.ImportGrid(dishes[i].CPM->ReturnGrid());
    org1.PolarTransform();
    // string name = "replicate" + to_string(i) + ".png";
    // org1.PolarToOutput(name);

    for (int j = i+1; j < par.n_orgs;j++)
    {
      fft org2;
      org2.AllocateGrid(par.sizex, par.sizey);
      org2.ImportGrid(dishes[j].CPM->ReturnGrid());
      org2.PolarTransform();
      // do comparison
      double inp = org1.PolarComparison(org2.GetPolar());

      // get the index
      int num{};
      for (int z=1;z<=i;++z)
        num += par.n_orgs - 1 - z;
      num+=j-1;

      // cout << "before check - " << i << " " << j << "  " << num << "  " << inp << endl;
      invariant_p.at(num)=inp;
    }
  }

  double inp_avg{};
  double inp_var{};

  for (double i : invariant_p)
  {
    inp_avg += i;
  }

  inp_avg = inp_avg / invariant_p.size();

  for (double i : invariant_p)
  {
    double val = pow(i - inp_avg, 2);
    inp_var += val;
  }
  inp_var = inp_var / proportions.size();


  var_name = "overlap_invariant.txt";
  outfile.open(var_name, ios::app);

  outfile << inp_avg << '\t' << inp_var << endl;
  outfile.close();





  delete[] dishes;

  return inter_org_fitness;
}


// Main function
int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record=false;
  par.velocities=false;
  
  Parameter();
  par.n_orgs = 60;
  if (par.file_genomes)
  {
    ifstream file("genomes.txt");
    vector<vector<vector<int>>> genomes;
    string line;
    while (getline(file, line)) 
    {
      vector<vector<int>> genome;

      vector<int> row{};
      stringstream ss(line);
      string value;
      while (ss >> value)
      {


        for (int i=0;i < value.size();++i)
        {
          if (value[i] == '-')
          {

            string ns{value[i]};
            string ns2{value[i+1]};
            string nsn = ns + ns2;
            row.push_back(stoi(nsn));
            break;
          }
          else if (isdigit(value[i]))
          {
            row.push_back(value[i] - '0');
          }

        }      

        if (row.size() == 9)
        {
          genome.push_back(row);
          row.clear();
        }          

      }
      genomes.push_back(genome);

    }

    file.close();

    for (vector<vector<int>> i : genomes)
    {

      // This is currently depracated. 
      vector<bool> start_p = { 0, 0, 0, 0 };
      vector<vector<vector<int>>> networks{};
      vector<vector<bool>> polarities{};
      for (int j=0;j<par.n_orgs;++j)
      {
          networks.push_back(i);
          polarities.push_back(start_p);
      }
      process_population(networks, polarities);
    }
  }
  else // just use start matrix in parameter.cpp
  {
      // This is currently depracated. 
      vector<bool> start_p = { 0, 0, 0, 0 };
      vector<vector<vector<int>>> networks{};
      vector<vector<bool>> polarities{};
      for (int j=0;j<par.n_orgs;++j)
      {
          networks.push_back(par.start_matrix);
          polarities.push_back(start_p);
      }
      process_population(networks, polarities);  
  }






  

  // make initial random networks. 
  // vector<vector<vector<int>>> networks{};
  // vector<vector<bool>> polarities{};
  // for (int i=0;i<par.n_orgs;++i)
  // {
  //     networks.push_back(par.start_matrix);
  //     polarities.push_back(start_p);
  // }

  // process_population(networks, polarities);

  // finished
  par.CleanUp();

  return 0;
}