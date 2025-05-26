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
#include <sstream>



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
std::uniform_int_distribution<> choose_diffuser(0, par.n_diffusers-1);

std::normal_distribution<double> diff_mut_curve(0.0, 0.25);



int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

INIT {

  try 
  {
    CPM->set_seed();
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
    CPM->ConstructInitCells(*this);

    if (par.make_rectangle)
    {
      CPM->Voronoi(0);
    }

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

void swapvv(vector<double> *xp, vector<double> *yp)
{
  vector<double> temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void swapb(vector<bool> *xp, vector<bool> *yp)
{
  vector<bool> temp = *xp;
  *xp = *yp;
  *yp = temp;
}

void swapd(Dish *d, int max_idx, int i)
{
  Dish tmp = d[max_idx];
  d[max_idx] = d[i];
  d[i] = tmp;
}

void sorter(vector<vector<vector<int>>> &networks, vector<double> &fitlist, Dish *dishes, vector<vector<double>>& coeffs)
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
    swapvv(&coeffs.at(max_idx), &coeffs.at(i));
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
      // slight ON bias for random networks. This is due to theta = -0.3. 
      if (val < 0.05)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.20)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.74)
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
    else if (val < 0.74)
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


vector<double> get_random_diff_coeffs()
{
  vector<double> rand_diff_coeffs;
  for (int i = 0; i < par.n_diffusers; ++i)
  {
    double val = double_num(mersenne);
    val = pow(val, 3.);
    double diff_co = par.min_diff_coeff + val * (par.max_diff_coeff - par.min_diff_coeff);
    rand_diff_coeffs.push_back(diff_co);
  }
  return rand_diff_coeffs;
}

void diff_mutate(vector<double> &coeffs)
{
  int tm = choose_diffuser(mersenne);
  coeffs[tm] = coeffs[tm] * exp( -diff_mut_curve(mersenne));
  if (coeffs[tm] > par.max_diff_coeff)
    coeffs[tm] = par.max_diff_coeff;
  else if (coeffs[tm] < par.min_diff_coeff)
    coeffs[tm] = par.min_diff_coeff;
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

void record_networks(vector<vector<vector<int>>>& netw, string oname, string time="")
{
  string nname = oname + "/" + "genomes" + "-" + time + ".txt";
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





void printn(vector<vector<int>> netw, vector<double> fitn, vector<double> coeffs)
{
  // create and open file
  std::string var_name = par.sim_file + "/gene_networks.txt";
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

  if (par.do_morphogen_evolution)
  {
    std::string coeff_file = par.sim_file + "/diff_coeffs.txt";
    outfile.open(coeff_file, ios::app);
    for (int i = 0; i < par.n_diffusers; ++i)
    {
      if (i == 0)
        outfile << "{ ";

      outfile << coeffs[i] << ", ";

      if (i == par.n_diffusers-1)
        outfile << "}" << endl;
    }
    outfile.close();
  }


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
  var_name = par.sim_file + "/fitness.txt";
  outfile.open(var_name, ios::app);
  outfile << max_fit << '\t' << avgfit << '\t' << chrono::duration <double, milli> (diff).count() << endl;

}








// function that simulates a population for a single evolutionary step. 
vector<double> process_population(vector<vector<vector<int>>>& network_list, int time, vector<vector<double>>& org_diff_coeffs)
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

    dishes[i].CPM->start_network(network_list.at(i));
    dishes[i].PDEfield->SetParameters(org_diff_coeffs[i]);


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
        if (t==100 && par.make_rectangle)
        {
          dishes[i].CPM->SetRectangularMF();
        }
        else if (!par.make_rectangle)
        {
          //programmed divisions
          if (t % par.div_freq == 0 && t <= par.div_end)
          {
            dishes[i].CPM->Programmed_Division(); // need to get the number of divisions right. 
          }
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
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }

  if (par.evo_pics && time % par.pic_gen_interval == 0)
  {
    string dirn = par.data_file + "/" + to_string(time+1);
    if (mkdir(dirn.c_str(), 0777) != -1)
      cout << "Directory created." << endl;
    for (int i=0; i < par.n_orgs; ++i)
    {
      dishes[i].CPM->ColourCells();
      fft new_org(par.sizex,par.sizey);
      new_org.ImportCPM(dishes[i].get_cpm());
      string f2 = "org-";
      string n2 = to_string(i);
      string ftype = ".png";
      string foutput = dirn + "/" + f2 + n2 + ftype;
      new_org.cpmOutput(foutput);
    }
  }
  string dirn = par.sim_file;
  if (time % 100 == 0)
    record_networks(network_list, dirn, to_string(time));

  delete[] dishes;

  // do sorting algorithm and return fitness
  sorter(network_list, inter_org_fitness, dishes, org_diff_coeffs);

  //output to standard output
  // output_networks(network_list);

  // output to file
  printn(network_list.front(), inter_org_fitness, org_diff_coeffs.front());


  vector<vector<vector<int>>> nextgen{};

  vector<vector<double>> nextgen_diffs{};
  int j = 0;
  for (int i=0; i < par.n_orgs;++i)
  {
    // Currently no random networks are added if largest fitness > this
    if (inter_org_fitness.front() > 35 || !par.insert_randoms)
    {
      nextgen.push_back(network_list.at(j));
      nextgen_diffs.push_back(org_diff_coeffs.at(j));

      //mutate network with probability = mut_rate
      double mu = double_num(mersenne);
      if (mu < par.mut_rate)
      {
        mutate(nextgen.back());
      }
      double mu2 = double_num(mersenne);
      if (mu2 < par.diff_mut_rate)
      {
        diff_mutate(nextgen_diffs.back());
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
        nextgen_diffs.push_back(get_random_diff_coeffs());
      }
      else 
      {
        nextgen.push_back(network_list.at(j));
        nextgen_diffs.push_back(org_diff_coeffs.at(j));
      }

      //mutate network with probability = 0.5
      double mu = double_num(mersenne);
      if (mu > par.mut_rate)
      {
        
        mutate(nextgen.back());
      }
      double mu2 = double_num(mersenne);
      if (mu2 < par.diff_mut_rate)
      {
        diff_mutate(nextgen_diffs.back());
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
    org_diff_coeffs.at(i) = nextgen_diffs.at(i);
  }
  return inter_org_fitness;
}




// Main function
int main(int argc, char *argv[]) {

#ifdef QTGRAPHICS
  if (par.evo_pics)
    QApplication* a = new QApplication(argc, argv);
  par.data_file = "images";
  if (mkdir(par.data_file.c_str(), 0777) != -1)
    cout << "Directory created." << endl;
  
#endif

  par.sim_file = "sim_data";
  if (mkdir(par.sim_file.c_str(), 0777) != -1)
    cout << "Directory created." << endl;

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
  par.mcs = 12100;
  par.n_orgs = 60;
  
  par.pic_gen_interval = 250;
  

  Parameter();



  // make initial random networks. 
  vector<vector<vector<int>>> networks{};

  vector<vector<double>> org_diff_coeffs{};

  if (par.read_in_evolution)
  {
    ifstream file("genomes.txt");
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
  
        if (row.size() == par.n_activators)
        {
          genome.push_back(row);
          row.clear();
        }          
      }
      networks.push_back(genome);
    }
    file.close();
    for (int i=0;i<par.n_orgs;++i)
    {
      org_diff_coeffs.push_back(par.evostart_diff_coeffs);
    }    
  }
  else
    for (int i=0;i<par.n_orgs;++i)
    {
      if (par.starter)
      {
        networks.push_back(par.start_n);
        org_diff_coeffs.push_back(par.evostart_diff_coeffs);
      }
      else
      {
        networks.push_back(get_random_network());
        org_diff_coeffs.push_back(get_random_diff_coeffs());
      }
    }


  for (int t=0;t<par.evs;++t)
  {
    cout << "current ev timestep is: " << t+1 << endl;
    // process population. 
    vector<double> fit = process_population(networks, t, org_diff_coeffs);

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
