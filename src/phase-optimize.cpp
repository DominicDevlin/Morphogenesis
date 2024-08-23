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


// Function to calculate the mean
double calculateMean(const std::vector<double>& data) {
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

// Function to calculate the median
double calculateMedian(std::vector<double>& data) {
    std::sort(data.begin(), data.end());
    size_t size = data.size();
    if (size % 2 == 0) {
        return (data[size / 2 - 1] + data[size / 2]) / 2;
    } else {
        return data[size / 2];
    }
}

// Function to calculate the IQR
std::pair<double, double> calculateIQR(std::vector<double>& data) {
    std::sort(data.begin(), data.end());
    size_t size = data.size();
    double q1, q3;

    if (size % 2 == 0) {
        std::vector<double> lower(data.begin(), data.begin() + size / 2);
        std::vector<double> upper(data.begin() + size / 2, data.end());
        q1 = calculateMedian(lower);
        q3 = calculateMedian(upper);
    } else {
        std::vector<double> lower(data.begin(), data.begin() + size / 2);
        std::vector<double> upper(data.begin() + size / 2 + 1, data.end());
        q1 = calculateMedian(lower);
        q3 = calculateMedian(upper);
    }

    return std::make_pair(q1, q3);
}

// Function to calculate the range
double calculateRange(const std::vector<double>& data) {
    auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());
    return *maxIt - *minIt;
}


void OutputShapes(map<int, vector<double>> data2, string &oname, int time)
{
  ofstream outfile;
  outfile.open(oname, ios::app);

  outfile << time << '\t';
  for (const auto& [key, values] : data2)
  {
    std::vector<double> tempValues = values;
    double mean = calculateMean(values);
    double median = calculateMedian(tempValues);
    auto [q1, q3] = calculateIQR(tempValues);
    double range = calculateRange(values);
    int observations = tempValues.size();

    outfile << key << '\t' << mean << '\t' << median << '\t' << q1 << '\t' << q3 << '\t' << range << '\t' << observations << '\t'; 
  }
  outfile << endl;

  outfile.close();  
}

void OutputResults(vector<double>& lengths, vector<double>& variances, vector<int>& phase_remained, string& oname, int time)
{

  //average fitness
  double avg_length = 0;
  double avg_variance = 0;
  double avg_phase_remained = 0;
  vector<double> vec;
  for (int i = 0; i < lengths.size(); ++i)
  {
    double fitness = pow(lengths[i],2) / variances[i];
    vec.push_back(fitness);

    avg_length += lengths[i];
    avg_variance += variances[i];
    avg_phase_remained += phase_remained[i];
  }
  avg_length = avg_length / lengths.size();
  avg_variance = avg_variance / lengths.size();
  avg_phase_remained = avg_phase_remained / lengths.size();
  
  int start = 1;
  int half = par.optimization_replicates / 2;

  std::sort(vec.begin(), vec.end(), std::greater<int>());
  double avgfit = std::accumulate(vec.begin() + start, vec.begin() + start + half, 0.0) / half;


  std::string var_name = oname + "/results.txt";
  std::ofstream outfile;
  outfile.open(var_name, ios::app);

  outfile << time << '\t' << avgfit << '\t' << avg_length << '\t' << avg_variance << '\t' << avg_phase_remained << endl;
  outfile.close();

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

  // we are going to skip the highest fitness one
  int start = 1;

  // use the next 3 highest fitness
  int half = par.optimization_replicates / 2;
  // need three lowest
  std::sort(fitn.begin(), fitn.end());
  double avgfit = std::accumulate(fitn.begin() + start, fitn.begin() + half + start, 0.0) / half;

  for (double i : fitn)
  {
    if (i > max_fit)
      max_fit = i;
    if (i < min_fit)
      min_fit = i;
  }
  // avgfit = avgfit / par.optimization_replicates;

  //output fitness 
  outfile << avgfit << '\t' << min_fit << '\t' << max_fit << endl;

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
  vector<int> cell_counter(par.optimization_replicates, 0);
  opt_out.resize(par.optimization_replicates);

  // create memory for dishes. 
  Dish* dishes = new Dish[par.optimization_replicates];
  int time{};
  time = int(params[5]);

  par.secr_rate[0] = params[0];
  par.Vs_max = params[1];
  par.J_stem_diff = params[2];
  // constant params
  par.J_stem = params[3];
  par.mcs= 40000 + int(par.J_stem)*15000;
  par.J_diff = params[4];

  // if (par.J_stem > par.J_diff)
  //   par.J_stem_diff = par.J_stem;
  // else
  //   par.J_stem_diff = par.J_diff;

  if (par.J_diff > par.J_stem)
  {
    par.J_med = 0.5*par.J_diff+0.25;
    par.J_med2 = 0.5*par.J_diff+0.25;
  }
  if (par.J_med < par.J_stem)
  {
    par.J_med = par.J_stem;
    par.J_med2 = par.J_stem;
  }


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

    // par.J_stem_diff = params[2];
    dishes[i].CPM->Set_evoJ(par.J_stem_diff);
    // par.Vs_max = params[3];
    // par.Vd_max = params[4];


    // if (i=0)
    // {
    //   cout << "params are: " << par.secr_rate[0] << '\t' << par.J_med << '\t' << par.J_stem_diff << '\t' << par.Vs_max << '\t' <<
    //   par.Vd_max << '\t' << par.secr_rate[0] << '\t' << par.secr_rate[0] << endl;
    // }

    dishes[i].CPM->CopyProb(par.T);

    

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
        // WE ARE GOING TO CHANGE THIS SO THAT IT JUST RANDOMLY ADDS MASS TO ONE OF THE STEM CELLS!! (can also do sigmoidal function?)
        dishes[i].CPM->DiscreteGrowthAndDivision(t);
        // dishes[i].CPM->CellGrowthAndDivision(t);
      }
      dishes[i].CPM->AmoebaeMove(t);
    
      if (t % 250 == 0 && par.insitu_shapes)
      {
        // cout << 'here' << endl;
        dishes[i].CPM->PhaseShapeIndex();
        // dish->CPM->AdhesionByState();
        dishes[i].CPM->HexaticOrder();
      }

      bool check_end = dishes[i].CPM->EndOptimizer(t);
      if (check_end)
      {
        cout << "here@@" << endl;
        t = par.mcs;
      }

      
      if (t == par.end_program)
        cell_counter[i] = dishes[i].CPM->CountCells();
      // finish simulation if organism is not growing.
      if (t%par.max_div_time==0 && t > 0)
      {
        int n_cells = dishes[i].CPM->CountCells();
        if (n_cells <= cell_counter[i])
        {
          cout << "here!!" << endl;
          t = par.mcs;
        }
        else
        {
          cell_counter[i] = n_cells;
        }
      }

      // ensure all cells are connected for shape calculations. 
      if (t > 0 && t % 5000 == 0)
      {
        bool check_shape = dishes[i].CPM->CheckShape();
        if (check_shape == false)
        {
          opt_out[i] = par.sizey;
          t = par.mcs;
          // cout << "Org number: " << i << " has bad shape. " << endl;
        }
      }

      // get fitness at end of development
      
    }
    opt_out[i] = dishes[i].CPM->Optimizer();       
    cout << "finished" << endl;
    if (i == 1)
      cout << "Sim #1 complete. The number of cells is: " << dishes[i].CPM->CountCells() << endl;

  }

  
  // combine maps together:
  map<int, vector<double>> data = dishes[0].CPM->Get_state_shape_index();

  for (int i = 1; i < par.optimization_replicates;++i)
  {
    map<int, vector<double>> next = dishes[i].CPM->Get_state_shape_index();
    for (auto&kv : next)
    {
      int key = kv.first;
      vector<double>& vec = kv.second;
      if (data.find(key) != data.end()) {
          data[key].insert(data[key].end(), vec.begin(), vec.end());
      } 
      else 
      {
          // If key does not exist, insert the new key-value pair
          data[key] = vec;
      }
    }
  }
  string oname = par.data_file + "/shape-data.txt";
  OutputShapes(data, oname, time);

  // combine maps together:
  map<int, vector<double>> hexdata = dishes[0].CPM->GetHexaticOrderList();

  for (int i = 1; i < par.optimization_replicates;++i)
  {
    map<int, vector<double>> next = dishes[i].CPM->GetHexaticOrderList();
    for (auto&kv : next)
    {
      int key = kv.first;
      vector<double>& vec = kv.second;
      if (hexdata.find(key) != hexdata.end()) {
          hexdata[key].insert(hexdata[key].end(), vec.begin(), vec.end());
      } 
      else 
      {
          // If key does not exist, insert the new key-value pair
          hexdata[key] = vec;
      }
    }
  }
  string hname = par.data_file + "/hexatic-data.txt";
  OutputShapes(hexdata, hname, time);


  vector<double> lengths;
  vector<double> variances;
  vector<int> phase_cells;

  for (int i = 0; i < par.optimization_replicates;++i)
  {
    pair<double,double> lw = dishes[i].CPM->LengthWidth();
    lengths.push_back(lw.first);
    variances.push_back(lw.second);
    
    int n_phase = dishes[i].CPM->CountPhaseOnCells();
    phase_cells.push_back(n_phase);

  }
  OutputResults(lengths, variances, phase_cells, par.data_file, time);
  
  // output to file
  printn(opt_out, par.data_file, params);

  if (par.pics_for_opt && time % par.pics_for_opt_interval == 0)
  {
    string dirn = par.pic_dir + "/" + to_string(time+1);
    if (mkdir(dirn.c_str(), 0777) != -1)
    {
      cout << "Directory created." << endl;
    }

    if (par.evo_pics)  
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

  return opt_out;
}




// Main function
int main(int argc, char *argv[]) 
{
  par.sizex = 200;
  par.sizey = 250;

  vector<double> params;
  for (int i = 1; i < argc; ++i)
  {
    params.push_back(stod(argv[i]));
    cout << stod(argv[i]) << " ";
  }
  cout << endl;

  par.pic_dir = par.pic_dir + "-" + argv[4] + "-" + argv[5];
  par.data_file = par.data_file + "-" + argv[4] + "-" + argv[5];

#ifdef QTGRAPHICS
  if (par.pics_for_opt)
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
  par.mcs = 100000;
  par.phase_evolution = true;
  par.insitu_shapes = true;
  par.measure_time_order_params=false;
  

  Parameter();

  string dirn = par.data_file;
  if (mkdir(dirn.c_str(), 0777) != -1)
  {
    cout << "Directory created." << endl;
    std::string var_name = par.data_file + "/optimize.txt";
    std::ofstream outfile;
    outfile.open(var_name, ios::app);
    outfile << "avg\tmin\tmax" << endl;
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
