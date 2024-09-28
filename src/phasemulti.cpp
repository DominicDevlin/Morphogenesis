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
        // sort(values.begin(), values.end());
        // double median;
        // int size = values.size();
        // if (size % 2 == 0)
        //     median = (values[size / 2 - 1] + values[size / 2]) / 2.0;
        // else
        //     median = values[size / 2];
        // outfile << "\t" << median;  // Output the average in the second column
        double sum = accumulate(values.begin(), values.end(), 0.0);
        double mean = sum / values.size();
        outfile << "\t" << mean;  // Output the mean in the second column
      } 
      else 
      {
        // cout << "Error in time output" << endl;
        outfile << "\t";  // No data for this row, leave empty
      }

      first_col = false;  // Set this to false after the first column
    }

    outfile << endl;  // Newline after each row
  }

  outfile.close(); 
  return phase_counts; 
}

void OutputOrder(vector<vector<pair<double,double>>> &shape_alignments, string oname)
{
    // Find the maximum length of the inner vectors
  int max_size = 0;
  for (const auto& vec : shape_alignments) 
  {
    if (vec.size() > max_size) 
    {
        max_size = vec.size();
    }
  }
  // Initialize vectors to store the sums and counts for each index
  vector<double> sum_first(max_size, 0.0);
  vector<double> sum_second(max_size, 0.0);
  vector<int> count_first(max_size, 0);
  vector<int> count_second(max_size, 0);
  // Iterate through all shape_alignments and sum up the values at each index
  for (const auto& vec : shape_alignments) 
  {
    for (int i = 0; i < vec.size(); ++i) {
        sum_first[i] += vec[i].first;
        sum_second[i] += vec[i].second;
        count_first[i]++;
        count_second[i]++;
    }
  }
  // Calculate the averages
  vector<pair<double, double>> averages(max_size);
  for (int i = 0; i < max_size; ++i) {
      if (count_first[i] > 0) {
          averages[i].first = sum_first[i] / count_first[i];
      }
      if (count_second[i] > 0) {
          averages[i].second = sum_second[i] / count_second[i];
      }
  }
  ofstream outfile;
  outfile.open(oname, ios::app);  // Append mode
  
  for (unsigned i = 0; i < averages.size(); ++i)
  {
    outfile << 300 + i * 100 << '\t' << averages[i].first << '\t' << averages[i].second << '\t' << endl;
  }

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

  vector<vector<pair<double,double>>> shape_alignments(par.n_orgs);
  vector<vector<pair<double,double>>> Z_values(par.n_orgs);
  vector<vector<pair<double,double>>> volumes(par.n_orgs);
  vector<vector<double>> cooperativities(par.n_orgs);
  vector<int> total_steps(par.n_orgs);

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    vector<pair<double,double>> org_zvals{};
    vector<pair<double,double>> org_alignments{};
    vector<pair<double,double>> org_volume{};

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
          dishes[i].CPM->PhaseHexaticOrder(t);
        }
        if (t > 200 && t % 500 == 0)
        {
          pair<double,double> zvals = dishes[i].CPM->PhaseZValues();
          org_zvals.push_back(zvals);
          pair<double,double> alignments = dishes[i].CPM->ShapeAlignmentByPhase();
          org_alignments.push_back(alignments);
          pair<double,double> vs = dishes[i].CPM->GetTargetandVolume();
          org_volume.push_back(vs);

        }

        // if (t > par.coop_start && t % 200 == 0)
        // {
        //   double coop = dishes[i].CPM->Cooperativity();
        //   cout << coop << endl;
        //   cooperativities[i].push_back(coop);
        // }

        if (par.velocities && t % 200 == 0)
        {
          dishes[i].CPM->RecordMasses();
          if (t > par.coop_start)
          {
            double coop = dishes[i].CPM->Cooperativity(200);
            cooperativities[i].push_back(coop);
          }
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

        // dishes[i].CPM->DiscreteGrowthAndDivision(t);
        if (t % par.cell_addition_rate == 0 && t > 200)
        {
          int cnum = dishes[i].CPM->FindHighestCell();
          int mnum = dishes[i].CPM->TopStalk();
          int check_n_points = dishes[i].CPM->CheckAddPoints();
          if (check_n_points < 60)
          {
            total_steps[i] = t;
            t = par.mcs;            
          }
          else
          {
            bool set=false;
            while (!set)
            {
              pair<int,int> val = dishes[i].CPM->ChooseAddPoint(mnum);
              set = dishes[i].CPM->SpawnCell(val.first, val.second, cnum, t);
            }
          }

        }        
      }
      dishes[i].CPM->AmoebaeMove(t);

      if (t % 1000 == 0 && t > 0)
      {
        dishes[i].CPM->RemoveUnconnectedCells();
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
      
      bool if_end = dishes[i].CPM->EndOptimizer(t);
      if (if_end == true)
      {
        total_steps[i] = t;
        t = par.mcs;
      }
      else if (t == par.mcs -1)
      {
        total_steps[i] = par.mcs;
      }
    }
    // this is not thread safe but changes of it happening at same time are very very rare.
    shape_alignments[i] = org_alignments;
    Z_values[i] = org_zvals;
    volumes[i] = org_volume;
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

  vector<double> coop_averages(par.n_orgs);


  double avg_phase_remained = 0;

  vector<double> fitnesses(par.n_orgs);
  vector<double> lengths(par.n_orgs);
  vector<double> variances(par.n_orgs);
  // vector<double> empty_spaces(par.n_orgs);

  for (int i=0; i < par.n_orgs;++i)
  {
    pair<double,double> lw = dishes[i].CPM->LengthWidth();
    lengths[i] = lw.first;
    variances[i] = lw.second;
    fitnesses[i] = pow(lw.first, 2) / lw.second;
   
    int n_phase = dishes[i].CPM->CountPhaseOnCells();
    avg_phase_remained += n_phase;

    // int empty_amount = dishes[i].CPM->EmptySpace();
    // empty_spaces[i] = empty_amount;

    double cop = std::accumulate(cooperativities[i].begin(), cooperativities[i].end(), 0.0);
    cop /= double(cooperativities[i].size());
    coop_averages[i] = cop;
    cout << cooperativities[i][0] << '\t' << coop_averages[i] << endl; 
  }

  vector<int> indices(par.n_orgs);
  std::iota(indices.begin(), indices.end(), 0); // Fill with 0, 1, 2, ..., n_orgs-1

  // Step 2: Sort the indices based on the values in 'fitnesses'
  std::sort(indices.begin(), indices.end(),
            [&fitnesses](int i1, int i2) { return fitnesses[i1] > fitnesses[i2]; });

  // std::sort(indices.begin(), indices.end(),
  //           [&empty_spaces](int i1, int i2) { return empty_spaces[i1] > empty_spaces[i2]; });       

  // Step 3: Reorder fitnesses, lengths, and variances based on the sorted indices
  vector<double> sorted_fitnesses(par.n_orgs);
  vector<double> sorted_lengths(par.n_orgs);
  vector<double> sorted_variances(par.n_orgs);
  vector<double> sorted_coops(par.n_orgs);
  // vector<double> sorted_empty_spaces(par.n_orgs);

  for (int i = 0; i < par.n_orgs; ++i) {
      sorted_fitnesses[i] = fitnesses[indices[i]];
      sorted_lengths[i] = lengths[indices[i]];
      sorted_variances[i] = variances[indices[i]];
      sorted_coops[i] = coop_averages[indices[i]];
      // sorted_empty_spaces[i] = empty_spaces[indices[i]];
  }
  fitnesses = sorted_fitnesses;
  lengths = sorted_lengths;
  variances = sorted_variances;
  coop_averages = sorted_coops;
  // empty_spaces = sorted_empty_spaces;

  // for (auto i : fitnesses)
  //   cout << i << '\t';
  // cout << endl;

  // lets output top half
  int start = 0;
  int half = par.n_orgs / 2;
  double avg_fitness = std::accumulate(fitnesses.begin() + start, fitnesses.begin() + half, 0.0) / half;
  double avg_length = std::accumulate(lengths.begin() + start, lengths.begin() + half, 0.0) / half;
  double avg_variance = std::accumulate(variances.begin() + start, variances.begin() + half, 0.0) / half;
  double avg_coop = std::accumulate(coop_averages.begin() + start, coop_averages.begin() + half, 0.0) / half;
  // double avg_empty_space = std::accumulate(empty_spaces.begin() + start, empty_spaces.begin() + half, 0.0) / half;
  avg_phase_remained = avg_phase_remained / par.n_orgs;

  ofstream outfile;
  string fname = par.data_file + "/fitness.txt";
  outfile.open(fname, ios::app);
  outfile << avg_fitness << '\t' << avg_length << '\t' << avg_variance << '\t' << avg_coop << '\t' <<  avg_phase_remained << endl;
  outfile.close();

  
  /*
    NOTES:
    NEED TO OUTPUT:
    Z-value (time based)
    shape alignment (time based)
    average fitness (length2 / variance) 
  */
  string o_align = par.data_file + "/alignments.txt";
  OutputOrder(shape_alignments, o_align);

  string o_Z = par.data_file + "/Z-order.txt";
  OutputOrder(Z_values, o_Z);

  string o_volume = par.data_file + "/volume-data.txt";
  OutputOrder(volumes, o_volume);


  int sum_steps = accumulate(total_steps.begin(), total_steps.end(), 0);

  string infoname = par.data_file + "/info.txt";
  outfile.open(infoname, ios::app);  // Append mode
  outfile << "average steps:\t" << double(sum_steps) / double(par.n_orgs) << '\n' 
  << "hex counts:\t" << t_hex_count << '\n' 
  << "shape counts:\t" << t_shape_count << '\n'
  << "hex per step:\t" << double(t_hex_count) / sum_steps * par.measure_interval << '\n' 
  << "shape per step:\t" << double(t_shape_count) / sum_steps * par.measure_interval << '\n'
  << "average breaks:\t" << double(n_times_apart) / double(par.n_orgs) << '\n'
  << "Jstem:\t" << par.J_stem << '\n'
  << "Jdiff:\t" << par.J_diff << '\n'
  << "differentiation rate:\t" << par.secr_rate[0] << '\n'
  << "addition rate:\t" << par.cell_addition_rate << '\n'
  << "Jsd:\t" << par.J_stem_diff << endl;
  outfile.close();


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

  bool read_in = false;
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
  par.gene_record=true;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  par.velocities=true;
  par.output_sizes = false;
  par.measure_time_order_params=true;
  Parameter();

  par.phase_evolution = true;
  par.min_phase_cells=4;
  par.mcs = 6000;
  par.sheet_hex=false;


  par.n_orgs = 2;
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
      // we are going to do three iterations, one with data, one control with no differentiation, one control with no differentiation and addtion
      par.secr_rate[0] = params[0];
      // par.Vs_max = params[1];
      par.cell_addition_rate = params[1];
      par.J_stem_diff = params[2];
      par.J_stem = params[3];
      par.J_diff = params[4];
      par.J_med=par.J_diff/2 + 0.25;
      par.mcs= 40000 + int(par.J_stem)*15000;
      if (par.J_stem > par.J_med)
        par.J_med = par.J_stem;
      par.J_med2 = par.J_med;
      cout << par.J_stem << '\t' << par.J_diff << '\t' << par.J_med << endl;
      process_population(networks, argnumber);
      ++argnumber;
      // control1 - no differentiation:
      par.secr_rate[0] = 0.0001;
      process_population(networks, argnumber);
      ++argnumber;
      //control 2 - no cell addition:
      par.cell_addition_rate = 100000;
      process_population(networks, argnumber);
      ++argnumber;
    }
  }

  
  // finished
  par.CleanUp();

  return 0;
}
