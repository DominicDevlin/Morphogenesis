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
#include "storage.h"
#include "connections.h"

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



  // run organisms in parallel. 
  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for 
  for (int i=0; i < par.n_orgs; ++i)  
  {


    dishes[i].CPM->set_num(i+1);
    // does init block above.
    dishes[i].Init();


    dishes[i].CPM->start_network(network_list.at(i), pols.at(i));

    for (int t=0;t<par.mcs;++t)
    {

      if (t == par.end_program)
      {
        dishes[i].CPM->CopyProb(par.lT);
        par.T = par.lT;
          
      }

      
      // record initial expression state. This occurs before any time step updates. 
      if (par.gene_output && t == 100)
        dishes[i].CPM->record_GRN();
      
      // programmed cell division section
      if (t < par.end_program)
      {


        if (t % par.div_freq == 0 && t <= par.div_end)
        {
          dishes[i].CPM->Programmed_Division(); // need to get the number of divisions right. 
        }

        
      
        if (t >= par.begin_network && t % par.update_freq == 0)
        {

          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell(); 
          
          if (par.gene_output)
            dishes[i].CPM->record_GRN();    

          // speed up initial PDE diffusion
          for (int r=0;r<par.program_its;r++) 
          {
            dishes[i].PDEfield->Secrete(dishes[i].CPM);
            dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
          } 
    
        }
      }
      else
      {
        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_network(t);
          dishes[i].AverageChemCell(); 
          if (par.gene_output)
          {
            dishes[i].CPM->record_GRN();
            dishes[i].CPM->CountTypesTime();
          }

          if (par.output_gamma)
          {
            dishes[i].CPM->RecordGamma();
          }

        }

        if (par.velocities)
        {
          dishes[i].CPM->RecordMasses();
        }
        if (par.output_sizes)
        {
          dishes[i].CPM->RecordSizes();
        }


        // if (par.gene_record)
        // {
        //   dishes[i].CPM->RecordTypes();
        // }

        // if (par.gene_record && t > par.adult_begins)
        // {
        //   dishes[i].CPM->RecordAdultTypes();
        // }

        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        }

        // print individual chemical concentrations. 
        // if (t % 100 == 0 && t > par.end_program)
        // {
        //   dishes[i].PDEfield->print_concentrations(dishes[i].CPM);
        // }

        dishes[i].CPM->CellGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);




      if (t >= 6000 && t < 8000 && t % 40 == 0 && par.scramble)
        dishes[i].CPM->swap_cells();


      //printing every 1000 steps. Do other debugging things here as well. 
      // if (t % 1000 == 0)
      // {


      //   dishes[i].CPM->print_random_cell();
      //   cout << "Number of cell types: " << dishes[i].CPM->get_ntypes() << endl;
      //   cout << t << " TIME STEPS HAVE PASSED." << endl;

      //   dishes[i].CPM->PrintPhenotypes();
      //   // dishes[i].CPM->WhiteSpace();
      //   // dishes[i].CPM->DeviationFromCircle();
      //   cout << "AVG BINDING: " << dishes[i].CPM->AverageBinding() << endl;
      //   cout << "NUMBER OF MEDIUM PROTEINS ON AVG: " << dishes[i].CPM->AvgMedsOn() << endl;
      //   cout << "ORG MASS IS: " << dishes[i].CPM->Mass() << endl;
      //   double center[] = {0.0,0.0};
      //   dishes[i].CPM->get_center(center);
      //   cout << "x center: " << center[0] << "   y center: " << center[1] << endl;

      //   dishes[i].CPM->PrintColours();
      // }
    }
  }


  int count = 0;
  for (int i=0; i < par.n_orgs; ++i)  
  {

    ++count;

    string fname = "org" + to_string(count);


    dishes[i].CPM->set_datafile(fname);
    par.data_file = fname;

    // string command = "mkdir " + fname;
    // system(command.c_str());
    if (mkdir(fname.c_str(), 0777) == -1)
      cerr << "Error : " << strerror(errno) << endl;
    else
      cout << "Directory created." << endl;   

  
    if (par.output_gamma)
      dishes[i].CPM->OutputGamma();



    // if (par.umap)

    dishes[i].CPM->get_fitness();

    map<int, int> phens = dishes[i].CPM->get_phenotype_time();
    map<int, int> types = dishes[i].CPM->get_AdultTypes();  


    
    map<pair<int,int>,int> edge_tally{};
    
    // check if there are super long cycles. Need to account for this tiny edge case where there is a >3000 mcs cycle (very annoying)
    bool cycling = dishes[i].CPM->CycleCheck();
    if (cycling)
    {
      cout << "There is cycling!!" << endl;
      dishes[i].CPM->set_long_switches(edge_tally);
    }
    else
    {
      dishes[i].CPM->set_switches(edge_tally);
    }

    vector<vector<int>> scc;

    if (par.potency_edges)
    {
      // entire program is run from ungraph now
      map<int,int>subcomps{};
      if (cycling)
      {
        Graph ungraph(phens.size());
        subcomps = ungraph.CreateUnGraph(phens, phens, edge_tally, 1, true);          
      }
      else
      {
        Graph ungraph(types.size());
        subcomps = ungraph.CreateUnGraph(phens, types, edge_tally);
        scc = ungraph.GetComps(types, 500);
        for (auto i : scc)
        {
            cout << "component: ";
            for (int j : i)
            {
                cout << j << "  ";
            }
            cout << std::endl;
        }

      }

      if (par.gene_output)
      {
        ofstream outfile;
        string switch_out = par.data_file + "/potency.dat";
        outfile.open(switch_out, ios::app);

        for (auto kv : subcomps)
        {
          outfile << "Component number: " << kv.first;
          if (kv.second > 1)
            outfile << " is weakly connected." << endl;
          else
            outfile << " is strongly connected." << endl;
        }
        outfile.close();
      }
      
    }
    else
    {

      // Graph newgraph(phens.size());
      // newgraph.CreateDiGraph(phens, types, edge_start, edge_end);

      // cout << "Having a look at the undirected graph...." << endl;
      Graph ungraph(phens.size());
      map<int,int> subcomps = ungraph.CreateUnGraph(phens, types, edge_tally);
    }

    if (par.gene_output)
    {
      ofstream outnet;
      string netw = par.data_file + "/network.txt";
      outnet.open(netw, ios::app);
      for (int i=0;i<par.n_genes;++i)
      {
        if (i == 0)
          outnet << "{ ";
        for (int j=0;j<par.n_activators;++j)
        {
          if (j==0)
            outnet << "{ " << par.start_matrix[i][j] << ", ";
          else if (j==par.n_activators-1)
            outnet << par.start_matrix[i][j] << " }, ";
          else 
            outnet << par.start_matrix[i][j] << ", ";
        }
        if (i == par.n_genes -1)
          outnet << "}" << endl;
      }

      outnet << endl << "Seed is: " << endl << par.pickseed;
      outnet.close();
    }



    if (par.velocities)
    {
      // dishes[i].CPM->CellVelocities();
      dishes[i].CPM->momenta(scc);
    }
      
    // dishes[i].CPM->SpecialVelocity();
    if (par.record_directions)
    {
      dishes[i].CPM->Directionality();
      // dishes[i].CPM->SingleCellDirection();
    }

  }


}


// Main function
int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;
  par.gene_output=true;
  par.gene_record=true;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  if (par.velocities)
  par.output_sizes = true;
  
  Parameter();


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

  vector<bool> start_p = { 0, 0, 0, 0 };
  vector<vector<bool>> polarities{};
  for (vector<vector<int>> i : genomes)
  {
    polarities.push_back(start_p);
  }

  par.n_orgs = 10;
  for (vector<vector<int>> i : genomes)
  {
    vector<vector<vector<int>>> replays(par.n_orgs, i);
    process_population(replays, polarities);

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
