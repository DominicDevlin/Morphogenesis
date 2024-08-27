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


int WriteData(const map<int, vector<pair<int, double>>>& shapedata, const string& oname)
{
  ofstream outfile;
  outfile.open(oname, ios::app);  // Append mode

  int phase_counts{};

  // First, find the maximum number of rows required
  int max_rows = 0;
  vector<int> rows{};
  for (const auto& [key, vec] : shapedata) {
    for (const auto& [index, value] : vec) 
    {
      if (index + 1 > max_rows) 
      {
          max_rows = index + 1;
          rows.push_back(index);
      }
      if (key==1) // phase on
      {
        ++phase_counts;
      }      
    }
  }

  // Write the header
  outfile << fixed << setprecision(6);
  
  for (int &row : rows) 
  {
    bool first_col = true;
    
    // Iterate over the map entries
    for (const auto& [key, vec] : shapedata) 
    {
      // Output the first column (the integer index)
      if (!first_col) 
      {
        outfile << "\t";  // Separate columns with a tab
      }

      outfile << row;

      // Calculate the average for this row if there are matching pairs
      vector<double> values;
      for (const auto& [index, value] : vec) 
      {
        if (index == row) 
        {
          values.push_back(value);
        }
      }

      if (!values.empty()) 
      {
        sort(values.begin(), values.end());
        double median;
        int size = values.size();
        if (size % 2 == 0)
            median = (values[size / 2 - 1] + values[size / 2]) / 2.0;
        else
            median = values[size / 2];
        outfile << "\t" << median;  // Output the average in the second column
      } 
      else 
      {
        cout << "Error in time output" << endl;
        outfile << "\t";  // No data for this row, leave empty
      }

      first_col = false;  // Set this to false after the first column
    }

    outfile << endl;  // Newline after each row
  }

  outfile.close(); 
  return phase_counts; 
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
    CPM->SetRandomTypes();    
    if (par.set_colours)
    {
      CPM->SetColours();
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
  int total_steps{};

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {

    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();
    dishes[i].CPM->start_network(network_list.at(i));
    dishes[i].CPM->Set_evoJ(par.J_stem_diff);

    bool stayed_together=true;

    // equilibriate cells with high T
    dishes[i].CPM->CopyProb(par.T);

    int t=0;

    for (; t < par.mcs; t++)
    {  
      // if (t % 1000 == 0)
      //   cout << t << endl;        
      if (par.gene_record && t == 100)
        dishes[i].CPM->record_GRN();
      if (t < par.end_program)
      {
        if (t % par.div_freq == 0 && t <= par.div_end)
        {
          
          dishes[i].CPM->Programmed_Division(par.phase_evolution); // need to get the number of divisions right. 
        }
       
      
        if (t >= par.begin_network && t % par.update_freq == 0)
        {

          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell(); 
          
          if (par.gene_record)
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
          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell(); 
          if (par.gene_record)
          {
            dishes[i].CPM->record_GRN();
            dishes[i].CPM->CountTypesTime();
          }

          if (par.output_gamma)
          {
            dishes[i].CPM->RecordGamma();
          }
        }

        if (t > 200 && par.measure_time_order_params && t % 1 == 0)
        {
          dishes[i].CPM->PhaseShapeIndex(t);
          dishes[i].CPM->HexaticOrder(t);
        }

        if (par.velocities)
        {
          dishes[i].CPM->RecordMasses();
        }
        if (par.output_sizes)
        {
          dishes[i].CPM->RecordSizes();
        }

        for (int r=0;r<par.pde_its;r++) 
        {
          dishes[i].PDEfield->Secrete(dishes[i].CPM);
          dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        }

        dishes[i].CPM->DiscreteGrowthAndDivision(t);
      }
      // ensure all cells are connected for shape calculations. 
      if (t > 0 && t % 5000 == 0 && stayed_together==true)
      {
        bool check_shape = dishes[i].CPM->CheckShape();
        if (check_shape == false)
        {
          ++n_times_apart;
          stayed_together=false;
        }
      }
      // if (par.output_sizes)
      // {
      //   dishes[i].CPM->RecordSizes();
      // }

      dishes[i].CPM->AmoebaeMove(t);
      bool if_end = dishes[i].CPM->EndOptimizer(t);
      if (if_end == true)
      {
        t = par.mcs;
      }
    }
    total_steps += t;
  }

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;


  int t_shape_count{};
  int t_hex_count{};

  if (par.measure_time_order_params)
  {
    map<int, vector<pair<int,double>>> hexdata = dishes[0].CPM->Get_time_hexatic_order();
    map<int, vector<pair<int,double>>> shapedata = dishes[0].CPM->Get_time_shape_index();

    for (int i = 1; i < par.n_orgs;++i)
    {
      map<int, vector<pair<int,double>>> next = dishes[i].CPM->Get_time_hexatic_order();
      for (auto&kv : next)
      {
        int key = kv.first;
        vector<pair<int,double>>& vec = kv.second;
        if (hexdata.find(key) != hexdata.end()) 
        {
            hexdata[key].insert(hexdata[key].end(), vec.begin(), vec.end());
        } 
        else 
        {
          // If key does not exist, insert the new key-value pair
          hexdata[key] = vec;
        }
      }
      map<int, vector<pair<int,double>>> shape_next = dishes[i].CPM->Get_time_shape_index();
      for (auto&kv : shape_next)
      {
        int key = kv.first;
        vector<pair<int,double>>& vec = kv.second;
        if (shapedata.find(key) != shapedata.end()) 
        {
            shapedata[key].insert(shapedata[key].end(), vec.begin(), vec.end());
        } 
        else 
        {
          // If key does not exist, insert the new key-value pair
          shapedata[key] = vec;
        }
      }
    }
    string oname = par.data_file + "/hex_time.dat";
    t_hex_count = WriteData(hexdata, oname);

    oname = par.data_file + "/shape_time.dat";
    t_shape_count = WriteData(shapedata, oname);    
  }
  /*
    NOTES:


    NEED TO OUTPUT:
    Z-value (time based)
    shape alignment (time based)
    average fitness (length2 / variance)
    
  
  
  
  */

  ofstream outfile;
  string infoname = par.data_file + "/info.txt";
  outfile.open(infoname, ios::app);  // Append mode
  outfile << total_steps << '\t' << double(total_steps) / 60 << '\t' << t_hex_count << '\t' << t_shape_count 
  << '\t' << double(t_hex_count) / total_steps * par.measure_interval << '\t' << double(t_shape_count) / total_steps * par.measure_interval << '\t' << n_times_apart << endl;



  if (par.pics_for_opt)
  {
    string dirn = par.pic_dir;
    if (mkdir(dirn.c_str(), 0777) != -1)
    {
      cout << "Directory created." << endl;
    }

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

}




int main(int argc, char *argv[]) 
{

  bool read_in = true;
  vector<vector<double>> params_data;
  if (read_in)
  {
    std::string filename = "params-list.txt";  // the file containing the params data
    std::ifstream infile(filename);

    // Check if file is open
    if (!infile.is_open()) 
    {
      std::cerr << "Could not open the file!" << std::endl;
      return 1;
    }

    std::string line;
    // Read the file line by line
    while (std::getline(infile, line)) 
    {
      std::vector<double> row;
      std::stringstream ss(line);
      std::string value;

      // Split the line by tab characters and convert to double
      while (std::getline(ss, value, '\t')) 
      {
        try 
        {
          double num = std::stod(value);  // Convert string to double
          row.push_back(num);  // Add to the current row
        } 
        catch (const std::invalid_argument& e) 
        {
          std::cerr << "Invalid number found in the file: " << value << std::endl;
        }
      }
      
      // Add the row to the params_data
      params_data.push_back(row);
    }
    infile.close();    
  }


#ifdef QTGRAPHICS
  {
    QApplication* a = new QApplication(argc, argv);
    // if (mkdir(par.pic_dir.c_str(), 0777) != -1)
    //   cout << "Directory created." << endl;
  }
#endif

  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record=true;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  par.velocities=false;
  if (par.velocities)
    par.output_sizes = true;
  else
    par.output_sizes = false;
  Parameter();

  par.phase_evolution = true;
  par.min_phase_cells=10;


  par.n_orgs = 60;
  vector<vector<vector<int>>> networks{};
  for (int i = 0; i < par.n_orgs; ++i)
  {
    networks.push_back(par.start_matrix);
  }
  if (!read_in)
  {
    process_population(networks);
  }
  else
  {
    int argnumber=1;
    for (vector<double> &params: params_data)
    {
      par.secr_rate[0] = params[0];
      par.Vs_max = params[1];
      par.J_stem_diff = params[2];
      par.J_stem = params[3];
      par.J_diff = params[4];
      par.J_med=par.J_diff/2 + 0.25;
      if (par.J_stem > par.J_med)
        par.J_med = par.J_stem;
      par.J_med2 = par.J_med;
      cout << par.J_stem << '\t' << par.J_diff << '\t' << par.J_med << endl;
      process_population(networks, argnumber);
      ++argnumber;
    }
  }




  // if (par.sheet_hex)
  // {
  //   // collapse into one vector
  //   vector<vector<double>> shape_final(n_trials);
  //   vector<vector<double>> hex_final(n_trials);

  //   for (int i = 0; i < n_trials; ++i)
  //   {
  //     for (auto &j : shape_index[i])
  //       for (double &k : j)
  //       {
  //         shape_final[i].push_back(k);
  //       }
  //     for (auto &j : hex_order[i])
  //       for (double &k : j)
  //       {
  //         hex_final[i].push_back(k);
  //       }
  //   }
  //   string var_name = par.data_file + "/shapes.dat";
  //   ofstream outfile;
  //   outfile.open(var_name, ios::app);  
  //   int length=0;
  //   for (const auto& inner_vec : shape_final) 
  //   {
  //       if (inner_vec.size() > length) 
  //       {
  //           length = inner_vec.size();
  //       }
  //   }
  //   // Print the vector of vectors as columns
  //   for (size_t i = 0; i < length; i++) 
  //   {
  //     for (size_t j = 0; j < shape_final.size(); j++) 
  //     {
  //       int k = 0;
  //       if (i < shape_final[j].size()) 
  //       {
  //           outfile << shape_final[j][i] << "\t";
  //       } 
  //     }
  //     outfile << std::endl; // Newline after each column is printed
  //   }
  //   outfile.close();  
  //   // now do same for hexatic order
  //   var_name = par.data_file + "/hex-order.dat";
  //   outfile.open(var_name, ios::app);  
  //   length=0;
  //   for (const auto& inner_vec : hex_final) 
  //   {
  //       if (inner_vec.size() > length) 
  //       {
  //           length = inner_vec.size();
  //       }
  //   }
  //   // Print the vector of vectors as columns
  //   for (size_t i = 0; i < length; i++) 
  //   {
  //     for (size_t j = 0; j < hex_final.size(); j++) 
  //     {
  //       int k = 0;
  //       if (i < hex_final[j].size()) 
  //       {
  //           outfile << hex_final[j][i] << "\t";
  //       } 
  //     }
  //     outfile << std::endl; // Newline after each column is printed
  //   }
  //   outfile.close();  
  // }
  
  // finished
  par.CleanUp();

  return 0;
}
