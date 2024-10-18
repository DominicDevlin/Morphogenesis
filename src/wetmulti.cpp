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


int WriteData(const vector<map<int, vector<pair<int, double>>>>& shapedata, const string& oname)
{
    ofstream outfile;
    outfile.open(oname, ios::app);  // Append mode

    int phase_counts{};

    // First, find the maximum number of rows required across all maps
    int max_rows = 0;
    vector<int>rows{};
    
    for (const auto& column : shapedata) {
        for (const auto& [key, vec] : column) {
            for (const auto& [index, value] : vec) {
                if (index + 1 > max_rows) {
                    max_rows = index + 1;
                    rows.push_back(index);
                }
                if (key == 1) {  // phase on
                    ++phase_counts;
                }
            }
        }
    }

    // Write the header
    outfile << fixed << setprecision(6);

    // Iterate over each row (index from 0 to max_rows - 1)
    for (int &row : rows) 
    {
        outfile << row;


        // Iterate over the vector of maps (each map represents a column)
        for (const auto& column : shapedata) 
        {


            // For each map in the vector, find the corresponding row's value (if it exists)
            for (const auto& [key, vec] : column) 
            {
              // only doing phase on for now
              if (key==1)
              {
                double sum = 0.0;
                int count = 0;

                for (const auto& [index, value] : vec) 
                {
                  if (index == row) {
                    sum += value;
                    ++count;
                  }
                }

                // Output the average for this row in the current column
                if (count > 0) 
                {
                  double average = sum / count;
                  outfile << '\t' << average;  // Output the average
                } 
                else 
                {
                  outfile << '\t' << 0.0;  // No data for this row, leave as 0.0
                }
              }

            }

        }

        outfile << endl;  // Newline after each row
    }

    outfile.close();
    return phase_counts;
}

void OutputCooperativities(vector<vector<double>> &cooperativities, string oname)
{
    std::vector<std::vector<double>> result;
    int intervalSize = par.measure_interval;
  
    for (const auto& vec : cooperativities) 
    {
        std::vector<double> averagedVec;

        // Ensure that the vector size is divisible by intervalSize
        int numIntervals = ceil(double(vec.size()) / double(intervalSize));


        // Loop over the intervals
        for (int i = 0; i < numIntervals; ++i) {
            // Compute the start and end of the current interval
            int start = i * intervalSize;
            int end = start + intervalSize;
            if (end > vec.size())
              end = vec.size();

            // Calculate the average over the interval
            double sum = std::accumulate(vec.begin() + start, vec.begin() + end, 0.0);
            double average = sum / intervalSize;

            // Add the average to the result vector
            averagedVec.push_back(average);
        }

        // Add the averaged vector to the result
        result.push_back(averagedVec);
    }

    // Open file for writing
    ofstream outputFile;
    outputFile.open(oname, ios::app);

    size_t max_inner_size = 0;
    for (const auto& vec : result) 
    {
        if (vec.size() > max_inner_size) 
        {
            max_inner_size = vec.size();
            // cout << "m_inner size: " << max_inner_size << endl;
        }
    }


    // Output data as columns where each inner vector corresponds to a column
    for (size_t i = 0; i < max_inner_size; ++i) 
    {
        // Write the row index as the first column
        outputFile << intervalSize*i;

        // Write the corresponding element from each inner vector
        for (size_t j = 0; j < result.size(); ++j) {
            if (i < result[j].size()) {
                outputFile << "\t" << result[j][i];
            } else {
                outputFile << "\t" << 0;  // If the inner vector is shorter, leave an empty space
            }
        }

        // Newline at the end of the row
        outputFile << "\n";
    }

    outputFile.close();  
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


void OutputColumnData(vector<vector<double>> &odata, string fname)
{
    // Open file for writing
    ofstream outputFile;
    outputFile.open(fname, ios::app);

    size_t max_inner_size = 0;
    for (const auto& vec : odata) {
        if (vec.size() > max_inner_size) 
        {
            max_inner_size = vec.size();
            // cout << "m_inner size: " << max_inner_size << endl;
        }
    }


    // Output data as columns where each inner vector corresponds to a column
    for (size_t i = 0; i < max_inner_size; ++i) 
    {
        // Write the row index as the first column
        outputFile << i;

        // Write the corresponding element from each inner vector
        for (size_t j = 0; j < odata.size(); ++j) {
            if (i < odata[j].size()) {
                outputFile << "\t" << odata[j][i];
            } else {
                outputFile << "\t" << 0;  // If the inner vector is shorter, leave an empty space
            }
        }

        // Newline at the end of the row
        outputFile << "\n";
    }

    outputFile.close();      
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
    if (par.velocities)
      par.output_sizes = true;

    if (par.do_voronoi)
    {
      par.highT=false;
      CPM->Voronoi(par.sizex, par.sheet_depth, par.sheet_shift);
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

  vector<vector<double>> cooperativities(par.n_orgs);
  vector<vector<double>> dewetting_length(par.n_orgs);
  vector<int> depin_time(par.n_orgs);

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();
    dishes[i].CPM->start_network(network_list.at(i));
    dishes[i].CPM->Set_evoJ(par.J_stem_diff);
    dishes[i].CPM->SetAreas(par.cell_areas);

    bool stayed_together=true;

    // equilibriate cells with high T
    dishes[i].CPM->CopyProb(par.T);

    int t=0;

    bool depin=false;

    for (; t < par.mcs; t++)
    {  
      if (t==0 && (par.lambda_perimeter > 0 || par.lambda_perimeter_phase>0))
      {
        // cout << par.cell_addition_rate << '\t' << par.J_med << '\t' << par.lambda_perimeter << endl;
        par.H_perim = true;
        dishes[i].CPM->SetPerims(par.ptarget_perimeter);
        dishes[i].CPM->MeasureCellPerimeters();
      }

      if (t == par.init_wetting)
      {
        dishes[i].CPM->WetTopCells(par.dewet_length, par.dewet_cell_depth);
      }

      if (t > par.init_wetting && t % 100 == 0)
      {
        if (!depin)
        {
          bool check = dishes[i].CPM->WettingDepinned();
          if (check)
          {
            depin_time[i] = t;
            depin = true;
          }
        }
        double dl = dishes[i].CPM->WettingRatio();
        dewetting_length[i].push_back(dl);
      }




      if (true && t>= par.begin_network)
      {
        if (t == par.begin_network)
        {
          dishes[i].CPM->StartWettingNetwork();
        }

        if (t % par.update_freq == 0)
        {
          dishes[i].CPM->update_phase_network(t);
          dishes[i].AverageChemCell();  

          // speed up initial PDE diffusion
          for (int r=0;r<par.program_its;r++) 
          {
            dishes[i].PDEfield->Secrete(dishes[i].CPM);
            dishes[i].PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
          } 
        }
      }




      if (t > 200 && par.measure_time_order_params && t % 1 == 0)
      {
        dishes[i].CPM->PhaseShapeIndex(t);
        dishes[i].CPM->PhaseHexaticOrder(t);
      }

      if (par.velocities && t % 1 == 0)
      {
        dishes[i].CPM->RecordMasses(true);
        if (t > par.coop_start)
        {
          double coop = dishes[i].CPM->Cooperativity(1);
          cooperativities[i].push_back(coop);
        }
      }

      // dishes[i].CPM->DiscreteGrowthAndDivision(t);
      if (t % par.cell_addition_rate == 0 && t > 200 && par.add_cells)
      {
        int cnum = dishes[i].CPM->FindHighestCell();
        int mnum = dishes[i].CPM->TopStalk();
        int check_n_points = dishes[i].CPM->CheckAddPoints();
        int counter = 0;
        bool set=false;
        if (check_n_points < 60)
        {
          t = par.mcs;            
        }
        else
        {
          while (!set)
          {
            pair<int,int> val = dishes[i].CPM->ChooseAddPoint(mnum);
            set = dishes[i].CPM->SpawnCell(val.first, val.second, cnum, t);
            ++counter;
            if (counter > 4)
            {
              set = true;
              cerr << "Error in spawn cell with phase count " << dishes[i].CPM->CountPhaseOnCells() << endl;
              t=par.mcs;
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

      if (par.pics_for_opt && t % 500 == 0)
      {
        string dirn = par.pic_dir;
        if (mkdir(dirn.c_str(), 0777) != -1)
        {
          cout << "Directory created." << endl;
        }

        for (int org=0; org < par.n_orgs; ++org)
        {
          dishes[i].CPM->ColourCells(par.phase_evolution);
          fft new_org(par.sizex,par.sizey);
          new_org.ImportCPM(dishes[org].get_cpm());
          string f2 = "org-";
          string n2 = to_string(org);
          string ftype = ".png";
          string foutput = dirn + "/" + f2 + n2 + "-" + to_string(t) + ftype;
          new_org.cpmOutput(foutput);
        }
      }


    }


  }

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;

  ostringstream stream;
  stream << fixed << setprecision(2) << par.J_stem; // Setting precision to 2 decimal points
  string formatted_value = stream.str();

  int t_shape_count{};
  int t_hex_count{};

  if (par.measure_time_order_params)
  {
    vector<map<int, vector<pair<int,double>>>> hexdata;
    vector<map<int, vector<pair<int,double>>>> shapedata;



    for (int i = 0; i < par.n_orgs;++i)
    {
      hexdata.push_back(dishes[i].CPM->Get_time_hexatic_order());
      shapedata.push_back(dishes[i].CPM->Get_time_shape_index());

    }


    string oname = par.data_file + "/hex_time-" + formatted_value + ".dat";
    t_hex_count = WriteData(hexdata, oname);

    oname = par.data_file + "/shape_time-" + formatted_value + ".dat";
    t_shape_count = WriteData(shapedata, oname);
  }

  string coopname = par.data_file + "/coop-" + formatted_value + ".dat";
  OutputCooperativities(cooperativities, coopname);

  double avg_phase_remained = 0;
  for (int i=0; i < par.n_orgs;++i)
  {

    int n_phase = dishes[i].CPM->CountPhaseOnCells();
    avg_phase_remained += n_phase;

    // int empty_amount = dishes[i].CPM->EmptySpace();
    // empty_spaces[i] = empty_amount;
  }

  
  // double avg_empty_space = std::accumulate(empty_spaces.begin() + start, empty_spaces.begin() + half, 0.0) / half;
  avg_phase_remained = avg_phase_remained / par.n_orgs;

  double avg_depin = std::accumulate(depin_time.begin(), depin_time.end(), 0.0);
  avg_depin /= double(par.n_orgs);

  string fname = par.data_file + "/dewetting-" + formatted_value + ".dat";
  OutputColumnData(dewetting_length, fname);

  ofstream outfile;
  string infoname = par.data_file + "/info.txt";
  outfile.open(infoname, ios::app);  // Append mode
  outfile << par.J_stem << '\t' << avg_depin << '\t' << t_hex_count << '\t' << t_shape_count << '\t' << double(n_times_apart) / double(par.n_orgs) << endl;
  outfile.close();

  delete[] dishes;

}




int main(int argc, char *argv[])  
{
  par.pics_for_opt = true;

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
  par.gene_record=false;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  par.velocities=true;
  par.output_sizes = false;
  par.measure_time_order_params=true;
  Parameter();
  
  par.phase_evolution = true;
  par.min_phase_cells=4;
  par.mcs = 4000;
  par.sheet_hex=false;
  par.n_orgs = 1;
  par.do_voronoi = true;
  par.add_cells = false;

  par.coop_wtime=3000;
  par.coop_stime=0;
  par.coop_start=1000;


  par.sizex=300;
  par.sizey=200;
  par.begin_network = par.mcs;

  bool perimeter_model = false;

  vector<vector<vector<int>>> networks{};
  for (int i = 0; i < par.n_orgs; ++i)
  {
    networks.push_back(par.start_matrix);
  }
  par.J_stem = 7;
  while (par.J_stem < 10)
  {
    
    par.J_diff = par.J_stem + 8.;
    par.J_med = par.J_diff / 2 + 0.25;
    par.J_med2 = par.J_med;
    par.J_stem_diff = par.J_diff;
    process_population(networks);
    par.J_stem+=0.25;
  }

  
  // finished
  par.CleanUp();

  return 0;
}
