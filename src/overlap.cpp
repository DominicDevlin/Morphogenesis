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
          // t = par.mcs;
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



  // Rotation and reflection invariant
  vector<double> invariant_p{};
  invariant_p.resize(((par.n_orgs-1)*par.n_orgs)/2);
  int counter = 1;

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for if(!par.overlap_images)
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
      double inp=0;
      if (!par.translation)
      {
        inp = org1.PolarComparison(org2.GetPolar(), true);
      }  
      else
      {
        // we dont need to move much about the origin because there is no "force".
        // Organisms will always be roughly translation invariant about the center.
        // we are going to translate at evenly distributed points to make a square around the center.
        // This will give a close enough approximation for translation, although it isn't really needed
        // because translation about the origin will almost always give the best result.
        // minimum value that we are going to iterate from
        int min_tx = -(floor(par.nt_intervals/2) * par.t_interval);
        int min_ty = -(floor(par.nt_intervals/2) * par.t_interval);

        for (int x=0;x<par.nt_intervals;++x)
        {
          for (int y=0;y<par.nt_intervals;++y)
          {
            int tx = min_tx + x*par.t_interval;
            int ty = min_ty + y*par.t_interval;
            org2.PolarTransform(tx, ty);
            double jcind = org1.PolarComparison(org2.GetPolar(), true);
            if (jcind > inp)
              inp = jcind;
          }
        }
      }

      // get the index
      int num{};
      for (int z=1;z<=i;++z)
        num += par.n_orgs - 1 - z;
      num+=j-1;

      // cout << "before check - " << i << " " << j << "  " << num << "  " << inp << endl;
      invariant_p.at(num)=inp;

      // used for sanity check for comparison.
      if (par.overlap_images)
      {
        string t1 = "org-i-" + to_string(counter) + ".png";
        string t2 = "org-i-shift-" + to_string(counter) + ".png";
        string t3 = "org-j-" + to_string(counter) + ".png";
        org1.PolarToOutput(t1);
        org1.ShowOptimal(t2);
        org2.PolarToOutput(t3);
        ++counter;

        // org1.OutputLoss();
        // string g1 = "grid-" + to_string(i) + ".png";
        // string g2 = "grid-" + to_string(j) + ".png";
        // org1.GridToOutput(g1);
        // org2.GridToOutput(g2);
      }
    }
    // string t1 = "org-" + to_string(i+1) + ".png";
    // org1.GridToOutput(t1);
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

  if (par.between_org_overlap)
  {
    var_name = "overlap-data.txt";
    outfile.open(var_name, ios::app);
    for (double i : invariant_p)
    {
      outfile << i << endl;
    }
    outfile.close();
  }





  delete[] dishes;

  return inter_org_fitness;
}





void print_fitness(vector<double>& fitn)
{
  // max fitness 
  double max_fit = 0;
  double stdev=0;

  //average fitness
  double avgfit = 0;
  double n_alive=0;
  for (double i : fitn)
  {
    if (i > 0)
    {
      ++n_alive;
      avgfit += i;
    }
    if (i > max_fit)
      max_fit = i;
  }
  avgfit = avgfit / n_alive;
  // calculate stdev
  for (double i : fitn)
  {
    if (i > 0)
    {
      stdev += pow(i-avgfit,2);

    }
  }
  stdev = sqrt(stdev/n_alive);


  //output fitness  
  string var_name = "fitness.txt";
  ofstream outfile;
  outfile.open(var_name, ios::app);
  outfile << max_fit << '\t' << avgfit << '\t' << stdev << endl;
  outfile.close();

}





// Main function
int main(int argc, char *argv[]) 
{



#ifdef QTGRAPHICS
  if (par.overlap_images)
    QApplication* a = new QApplication(argc, argv);
#endif



  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record=false;
  par.velocities=false;
  
  Parameter();
  par.n_orgs = par.overlap_orgs;
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

    if (!par.between_org_overlap) // compare same genomes
    {
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
        vector<double> fitness = process_population(networks, polarities);
        print_fitness(fitness);
      }
    }
    else // compare different genomes 
    {
      par.n_orgs = genomes.size();
      vector<bool> start_p = { 0, 0, 0, 0 };
      vector<vector<vector<int>>> networks{};
      vector<vector<bool>> polarities{};      
      for (vector<vector<int>> i : genomes)
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
