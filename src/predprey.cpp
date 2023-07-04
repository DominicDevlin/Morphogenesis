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

void sorter(vector<vector<vector<int>>> &networks, vector<double> &fitlist)
{
  int i{};
  int j{};
  int max_idx{};
  int n = par.n_pred;

  for (i = 0; i < n-1; i++)
  {
    max_idx = i;
    for (j=i+1;j<n;j++)
    {
      if (fitlist[j] > fitlist[max_idx])
        max_idx = j;
    }
    // swap largest element with first element
    swap(&fitlist[max_idx], &fitlist[i]);
    swapv(&networks.at(max_idx), &networks.at(i));
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
      // slight ON bias for random networks. This is due to theta = -0.3. 
      if (val < 0.01)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.18)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.8)
      {
        matrix[i][j] = 0;
      }
      else if (val < 0.98)
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
    else if (val < 0.2)
    {
      network[i][j] = -1;
    }
    else if (val < 0.8)
    {
      network[i][j] = 0;
    }
    else if (val < 0.95)
    {
      network[i][j] = 1;
    }
    else
    {
      network[i][j] = 2;
    }
  }
}


//function to compare lattices
array<long,2> CompareLattice(vector<vector<bool>> &pred, vector<vector<bool>> &prey)
{
  // first is pred, second is prey.
  array<long, 2> points = {0,0};

  // going to start by not using center (might need to change to this later? unsure)
  for (int x=1;x<par.sizex-1;++x)
    for (int y=1;y<par.sizey-1;++y)
    {
      if (prey[x][y] && pred[x][y])
      {
        ++points[0];
        --points[1];
      }
      else if (pred[x][y])
      {
        --points[0];
        ++points[1];
      }
      else if (prey[x][y])
      {
        ++points[1];
      }      
    }
  return points;
}







// function that simulates a population for a single evolutionary step. 
vector<double> process_population(vector<vector<vector<int>>>& pred,
vector<vector<vector<int>>>& prey)
{
  // create memory for dishes. 
  Dish* dishes = new Dish[par.n_orgs];

  // grids for comparison
  vector<vector<vector<bool>>> pred_grids;
  pred_grids.resize(par.n_pred);

  vector<vector<vector<bool>>> prey_grids;
  prey_grids.resize(par.n_pred);

  // check that all the cells are connected. 
  vector<bool> pred_connect_check;
  pred_connect_check.resize(par.n_pred);
  vector<bool> prey_connect_check;
  prey_connect_check.resize(par.n_pred);


  // run organisms in parallel. 
  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for 
  for (int i=0; i < par.n_orgs; ++i)  
  {
    if (i < par.n_pred)
    {
      dishes[i].CPM->set_num(i+1);
      // does init block above.
      dishes[i].Init();

      int t;

      dishes[i].CPM->start_network(pred.at(i));


      // run simulation for single organism for mcs montecarlo steps.
      for (t=0;t<par.mcs;t++) 
      {
        // PROGRAMMED CELL DIVISION SECTION
        if (t < par.end_program)
        {

          if (t%par.div_freq==0 && t > 0 && t <= par.div_end)
            dishes[i].CPM->Programmed_Division();
          
          if (t >= par.begin_network && t % par.update_freq == 0)
          {
            if (t == par.begin_network && par.morphogen)
              dishes[i].CPM->morphogenWave();

            dishes[i].CPM->update_network(t);
            dishes[i].AverageChemCell();
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
          dishes[i].CPM->CellGrowthAndDivision();
        }
        dishes[i].CPM->AmoebaeMove(t);

        // ensure all cells are connected for shape and curvature calculations. 
        if (t > 0 && t % 2000 == 0)
        {
          dishes[i].CPM->CheckShape();
        }    
      }
      pred_grids[i] = dishes[i].CPM->ReturnGrid();

      if (dishes[i].CPM->MaintainedShape())
      {
        pred_connect_check[i] = true;
      }
      else
      {
        pred_connect_check[i] = false;
      }

      if (i == 1)
        cout << "Pred #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;
    }
    else
    {
      int j = i - (par.n_pred);

      dishes[i].CPM->set_num(i+1);
      // does init block above.
      dishes[i].Init();

      int t;

      dishes[i].CPM->start_network(prey.at(j));


      // run simulation for single organism for mcs montecarlo steps.
      for (t=0;t<par.mcs;t++) 
      {
        // PROGRAMMED CELL DIVISION SECTION
        if (t < par.end_program)
        {

          if (t%par.div_freq==0 && t > 0 && t <= par.div_end)
            dishes[i].CPM->Programmed_Division();
          
          if (t >= par.begin_network && t % par.update_freq == 0)
          {
            if (t == par.begin_network && par.morphogen)
              dishes[i].CPM->morphogenWave();

            dishes[i].CPM->update_network(t);
            dishes[i].AverageChemCell();
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
          dishes[i].CPM->CellGrowthAndDivision();
        }
        dishes[i].CPM->AmoebaeMove(t);

        // ensure all cells are connected for shape and curvature calculations. 
        if (t > 0 && t % 2000 == 0)
        {
          dishes[i].CPM->CheckShape();
        }    
      }
      prey_grids[j] = dishes[i].CPM->ReturnGrid();

      if (dishes[i].CPM->MaintainedShape())
      {
        prey_connect_check[j] = true;
      }
      else
      {
        prey_connect_check[j] = false;
      }

      if (i == 1)
        cout << "Prey #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

    }
  }
  delete[] dishes;


  vector<double> pred_fitness{};
  vector<double> prey_fitness{};
  for (int i=0;i<par.n_pred;++i)
  {
    pred_fitness.push_back(0.0);
    prey_fitness.push_back(0.0);
  }
  vector<int> pred_count(par.n_pred, 0);
  vector<int> prey_count(par.n_pred, 0);
  // Separate predator and prey, do sorting
  // Need to be careful about which organisms are not all connected. 
  for (int i=0; i<par.n_pred;++i)
    for (int j=0;j<par.n_pred;++j)
    {
      // cout << pred_connect_check[i] << "   " << prey_connect_check[j] << endl;
      if (!pred_connect_check[i] && !prey_connect_check[j])
      {
        pred_fitness[i] = -10000000.;
        prey_fitness[j] = -10000000.;
      }
      else if (!pred_connect_check[i])
      {
        pred_fitness[i] = -10000000.;
      } 
      else if (!prey_connect_check[j])
      {
        prey_fitness[j] = -10000000.;
      }
      else
      {
        array<long,2> fight = CompareLattice(pred_grids[i], prey_grids[j]);
        pred_fitness[i] += fight[0];
        prey_fitness[j] += fight[1];
        ++pred_count[i];
        ++prey_count[j];
        // cout << "Fight between " << i << " and " << j << ": " << fight[0] << "  " << fight[1] << endl;
      }
    }

  for (int i=0;i<par.n_pred;++i)
  {
    if (pred_count[i] > 0)
    {
      pred_fitness[i] = pred_fitness[i] / pred_count[i];
    }
    if (prey_count[i] > 0)
    {
      prey_fitness[i] = prey_fitness[i] / prey_count[i];
    }
  } 


  // Sort pred and prey vectors separately 
  sorter(pred, pred_fitness);
  sorter(prey, prey_fitness);



  vector<vector<vector<int>>> nextgen_pred{};
  int j = 0;
  for (int i=0; i < par.n_pred;++i)
  {
    nextgen_pred.push_back(pred.at(j));

    //mutate network with probability = mut_rate=.5
    double mu = double_num(mersenne);
    if (mu < par.mut_rate)
    {
      mutate(nextgen_pred.back());
    }
    if (j > par.n_orgs / 8)
      j=0;
    else
      ++j;
  }
  vector<vector<vector<int>>> nextgen_prey{};
  j=0;
  for (int i=0; i < par.n_pred;++i)
  {
    nextgen_prey.push_back(prey.at(j));

    //mutate network with probability = mut_rate=.5
    double mu = double_num(mersenne);
    if (mu < par.mut_rate)
    {
      mutate(nextgen_prey.back());
    }
    if (j > par.n_orgs / 8)
      j=0;
    else
      ++j;
  }
  // set nextgen = current pop. 
  for (int i=0;i<par.n_pred;++i)
  {
    pred[i] = nextgen_pred[i];
    prey[i] = nextgen_prey[i];
  }

  // push all prey fitness onto pred fitness and return it.
  int p;
  for (p=0; p<par.n_pred;++p)
  {
    pred_fitness.push_back(prey_fitness[p]);
  }
  return pred_fitness;
}


void printn(vector<vector<int>> pred, vector<vector<int>> prey, vector<double> fitn)
{
  // create and open file
  std::string var_name = "gene_networks.txt";
  std::ofstream outfile;
  outfile.open(var_name, ios::app);
  outfile << "Fittest Predator: " << endl;
  for (int i=0;i<par.n_genes;++i)
  {
    if (i == 0)
      outfile << "{ ";
    for (int j=0;j<par.n_activators;++j)
    {
      if (j==0)
        outfile << "{ " << pred[i][j] << ", ";
      else if (j==par.n_activators-1)
        outfile << pred[i][j] << " }, ";
      else 
        outfile << pred[i][j] << ", ";
    }
    if (i == par.n_genes -1)
      outfile << "}" << endl;
  }
  outfile << endl;

  // Report best gene network for prey as well. 
  outfile << "Fittest Prey: " << endl;
  for (int i=0;i<par.n_genes;++i)
  {
    if (i == 0)
      outfile << "{ ";
    for (int j=0;j<par.n_activators;++j)
    {
      if (j==0)
        outfile << "{ " << prey[i][j] << ", ";
      else if (j==par.n_activators-1)
        outfile << prey[i][j] << " }, ";
      else 
        outfile << prey[i][j] << ", ";
    }
    if (i == par.n_genes -1)
      outfile << "}" << endl;
  }
  outfile << endl;
  outfile.close();

  // max fitness 
  double max_predfit = fitn.front();
  double max_preyfit = fitn[par.n_pred];

  //average fitness
  double pred_avgfit = 0;
  double prey_avgfit = 0;
  int i = 0;

  int count1 = 0;
  int count2 = 0;
  while (i < par.n_pred)
  {
    if (fitn[i] > -10000)
    {
      pred_avgfit += fitn[i];
      ++count1;
    }

    ++i;
  }
  while (i < par.n_orgs)
  {
    if (fitn[i] > -10000)
    {
      prey_avgfit += fitn[i];
      ++count2;
    }
    ++i;
  }
  if (count1 > 0)
    pred_avgfit = pred_avgfit / count1;
  if (count2 > 0)
    prey_avgfit = prey_avgfit / count2;
  
  //calculate time since begin
  auto end = chrono::steady_clock::now();
  auto diff = end - start;

  //output fitness and time since beginning simulation. 
  var_name = "fitness.txt";
  outfile.open(var_name, ios::app);
  outfile << max_predfit << '\t' << pred_avgfit << '\t' << max_preyfit << '\t' << prey_avgfit << '\t' << chrono::duration <double, milli> (diff).count() << endl;

}


// Main function
int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=false;
  par.randomise=false;
  par.gene_output=false;
  Parameter();

  // If starting from a specific organism. 
  bool starter = false;
  vector<vector<int>> start_n = { { 0, 0, 2, 0, 0, 0, -1, 2 }, { 0, 0, 1, 0, 0, 0, 0, 0 }, { 0, -1, 0, 0, 1, 1, 0, -1 }, { -1, 0, 1, 0, 0, 0, 1, 0 }, { 1, 0, 1, 1, 0, 1, 0, -1 }, { -1, 0, 0, 2, -1, 0, -1, 0 }, { 1, -1, 1, -1, 1, 0, -1, 0 }, { 0, 1, 0, 0, -1, 0, 1, 0 }, { 0, 0, -1, 0, 1, 0, 1, 2 }, { 0, 0, 1, 0, 0, 2, 1, 1 }, { 1, 0, 0, -1, 0, 0, 0, 1 }, { 0, 1, 0, 1, 0, 0, -1, 0 }, { 0, -1, 0, 0, -1, 1, -1, 1 }, { 0, 0, 2, -1, 0, 2, -1, 1 }, { 0, 0, -1, 2, 2, 1, -1, 0 }, { 0, 1, 0, -1, 1, 1, 0, -1 }, { 0, -1, -1, 0, 0, -1, 0, 0 }, { -1, 0, 0, 1, 0, 0, 2, 0 }, { 0, 0, 1, 0, 0, 0, -1, 0 }, { 0, 0, -1, 0, 0, 0, 0, 1 }, };


  // make initial random networks. 
  vector<vector<vector<int>>> pred{};
  vector<vector<vector<int>>> prey{};

  for (int i=0;i<par.n_pred;++i)
  {
    if (starter)
    {
      pred.push_back(start_n);
      prey.push_back(start_n);
    }
    else
    {
      pred.push_back(get_random_network());
      prey.push_back(get_random_network());
    }
  }



  for (int t=0;t<par.evs;++t)
  {
    cout << "current ev timestep is: " << t+1 << endl;
    // process population. 
    vector<double> fit = process_population(pred, prey);

    // output every x evolution steps. 
    if (t%1==0)
    {
      printn(pred.front(), prey.front(), fit);
    }
  }
  // finished
  par.CleanUp();

  return 0;
}
