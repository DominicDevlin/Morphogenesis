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


void OutputIntColumnData(vector<vector<int>> &odata, string fname)
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


    par.highT=false;
    int xtoshift = par.sizex/2 - par.dewet_length/2;
    int ytoshift = par.sizey/2 - par.L2/2;
    CPM->Voronoi(par.dewet_length,round(par.L2+5), ytoshift, xtoshift);
    
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

  vector<vector<double>> cooperativities(par.n_orgs);
  // should be normalised to zero.
  vector<vector<double>> dewetting_ratio(par.n_orgs);
  vector<vector<int>> dewetting_length(par.n_orgs);

  vector<vector<double>> shape_proportions(par.n_orgs);
  // vector<vector<double>> nbh_exchange_rates(par.n_orgs);

  ostringstream makefll;
  makefll << fixed << setprecision(2) << par.J_L; // Setting precision to 2 decimal points
  string fnamer = par.data_file + "/" + makefll.str();
  if (par.record_transitions)
  {
    if (mkdir(fnamer.c_str(), 0777) == -1)
      cerr << "Error : " << strerror(errno) << endl;
    else
      cout << "Directory created." << endl;
  }



  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();
    dishes[i].CPM->start_network(network_list.at(i));
    dishes[i].CPM->Set_evoJ(par.J_SL);
    dishes[i].CPM->SetAreas(par.cell_areas);

    dishes[i].CPM->WetAllCells();
    // equilibriate cells with high T
    dishes[i].CPM->CopyProb(par.T);



    int t=0;

    for (; t < par.mcs; t++)
    {  
      dishes[i].CPM->SetCellCenters();

      if (t==10 && par.MakeEpithelia)
      {
        dishes[i].CPM->AddEpithelialLayer();
      }

      if (t == 10)
      {
        par.init_wet_length = dishes[i].CPM->WettingLength();
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

      if (t > 0 && t % 100 == 0)
      {
        double dl = dishes[i].CPM->WettingRatio();
        dewetting_ratio[i].push_back(dl);
        dewetting_length[i].push_back(dishes[i].CPM->WettingLength());
      }



      if (t > par.struct_avg_interval-1 && par.measure_time_order_params && t % par.measure_interval == 0)
      {
        dishes[i].CPM->MeasureHexaticOrder();
        dishes[i].CPM->MeasureShapeIndex();    
      }
      if (par.measure_time_order_params && t % par.struct_avg_interval == 0 && t > par.struct_avg_interval)
      {
        dishes[i].CPM->AverageShapeIndex();
        dishes[i].CPM->AverageHexaticOrder();
      } 


      if (par.record_transitions && t>10 && t % 10 == 0)
      {
        vector<vector<double>> shared_centres = dishes[i].CPM->find_shared_centres();
        if (shared_centres.size() > 0)
        {
          ostringstream stream;
          stream << fixed << setprecision(2) << par.J_L; // Setting precision to 2 decimal points
          string formatted_value = stream.str();
          string fnamee = fnamer + "/transitions-" + formatted_value + "-" + to_string(i+1) +".dat";
          ofstream outfile;
          outfile.open(fnamee, ios::app);  // Append mode
          outfile << fixed << setprecision(3);
          for (auto &vv : shared_centres)
          {
            outfile << t << '\t' << vv[4] << '\t' << vv[5] << '\t' << vv[0] << '\t' << vv[1] <<'\t' << vv[2] <<'\t' << vv[3] << vv[6] << endl;
          }
          outfile.close();
        }
      }


      // if (t >= par.init_wetting && t % 1000 == 0)
      // {
      //   double nbh_exchange = dishes[i].CPM->NeighbourExchangeRate();
      //   if (nbh_exchange >= 0)
      //     nbh_exchange_rates[i].push_back(nbh_exchange);
      // }
      // if (t > 100 + par.measure_interval && par.measure_time_order_params && t % par.measure_interval == 0)
      // {
      //   double shape_pr = dishes[i].CPM->ReturnShapeProportion();

      //   shape_proportions[i].push_back(shape_pr);
      // }
     

      dishes[i].CPM->AmoebaeMove(t);

      // if (t % 1000 == 0 && t > 0)
      // {
      //   dishes[i].CPM->RemoveUnconnectedCells();
      // }

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
  

  if (par.measure_time_order_params)
  {
    int container_size = par.mcs / par.struct_avg_interval - 2;
    vector<double> shape_index_output(container_size, 0.);
    vector<int> shape_counts(container_size, 0);
    vector<int> hex_counts(container_size, 0);
    vector<double> hex_order_output(container_size, 0.);
    for (int i = 0; i < par.n_orgs; ++i)
    {
      vector<double>& org_shapes = dishes[i].CPM->ReturnShapeIndex();
      vector<double>& org_hexes = dishes[i].CPM->ReturnHexaticOrder();
      for (int j = 0; j < container_size; ++j)
      {
        if (org_shapes[j] > 0.01)
        {
          shape_index_output[j] += org_shapes[j];
          shape_counts[j] += 1;
        }
        

        if (org_hexes[j] > 0.01)
        {
          hex_order_output[j] += org_hexes[j];
          hex_counts[j] += 1;
        }
      }
    }
    for (int i = 0; i < container_size; ++i)
    {
      if (shape_counts[i]==0)
      {
        shape_index_output[i] = 0.;
      }
      else
      {
        shape_index_output[i] /= shape_counts[i];
      }

      if (hex_counts[i]==0)
      {
        hex_order_output[i] = 0.;
      }
      else
      {
        hex_order_output[i] /= hex_counts[i];
      }
    }
    
    ostringstream stream;
    stream << fixed << setprecision(2) << par.J_L; // Setting precision to 2 decimal points
    string formatted_value = stream.str();
    string oname = par.data_file + "/hex_time-" + formatted_value + ".dat";
    ofstream outfile;
    outfile.open(oname, ios::app);  // Append mode
    outfile << fixed << setprecision(3);
    for (int i = 0; i < container_size; ++i)
    {
      outfile << i*par.struct_avg_interval + 2*par.struct_avg_interval << '\t' << hex_order_output[i] << endl;
    }
    outfile.close();

    oname = par.data_file + "/shape_time-" + formatted_value + ".dat";
    outfile.open(oname, ios::app);  // Append mode
    outfile << fixed << setprecision(3);
    for (int i = 0; i < container_size; ++i)
    {
      outfile << i*par.struct_avg_interval + 2*par.struct_avg_interval << '\t' << shape_index_output[i] << endl;
    }
    outfile.close();

  }





  ostringstream stream;
  stream << fixed << setprecision(2) << par.J_L; // Setting precision to 2 decimal points
  string formatted_value = stream.str();

  int t_shape_count{};
  int t_hex_count{};

  if (par.velocities)
  {
    string coopname = par.data_file + "/coop-" + formatted_value + ".dat";
    OutputCooperativities(cooperativities, coopname);
  }

  string fname = par.data_file + "/dewetting.ratio-" + formatted_value + ".dat";
  OutputColumnData(dewetting_ratio, fname);

  // fname = par.data_file + "/neighbour.exchange-" + formatted_value + ".dat";
  // OutputColumnData(nbh_exchange_rates, fname);

  fname = par.data_file + "/dewetting.length-" + formatted_value + ".dat";
  OutputIntColumnData(dewetting_length, fname);

  ofstream outfile;
  string infoname = par.data_file + "/info.txt";
  outfile.open(infoname, ios::app);  // Append mode
  outfile << par.J_L << '\t' << t_hex_count << '\t' << t_shape_count << endl;
  outfile.close();

  delete[] dishes;

}




int main(int argc, char *argv[])  
{
  par.pics_for_opt = false;

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

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;
  Parameter();
  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record=false;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  par.velocities=false;
  par.output_sizes = false;
  par.measure_time_order_params=true;
  par.record_transitions=false;
  par.MakeEpithelia=true;
  
  
  par.phase_evolution = true;
  par.mcs = 400000;
  par.sheet_hex=false;
  par.n_orgs = 120;
  par.do_voronoi = true;

  par.coop_wtime=3000;
  par.coop_stime=0;
  par.coop_start=1000;

  par.sizex=300;
  par.sizey=250;

  
  par.dewet_cell_depth=2;
  par.conserved_dewet_distance = 150;
  // double tmp_length = (sizex - 100 - 2 * sqrt((1240 * dewet_cell_depth ) / M_PI)) / 2.;

  par.L2 = sqrt((sqrt(3.)/2) * par.cell_areas ) * par.dewet_cell_depth;
  double tmp_length = par.L2 + (sqrt(pow(par.L2,2) + pow(M_PI * par.conserved_dewet_distance,2) )) / M_PI;
  par.dewet_length=round(tmp_length);
  par.theoretical_diameter = 2 * sqrt((par.dewet_length * par.L2)/M_PI);
  
ddd

  vector<vector<vector<int>>> networks{};
  for (int i = 0; i < par.n_orgs; ++i)
  {
    networks.push_back(par.start_matrix);
  }
  par.J_L = 1;
  while (par.J_L < 9.1)
  {
    // DO NOT CHANGE THESE NUMBERS! Keeps surface tension at 4.25
    par.J_S = par.J_L + 8.;
    par.J_med = par.J_S / 2 + 0.25;
    par.J_med2 = par.J_med;
    par.J_SL = par.J_S;

    if (par.MakeEpithelia)
    {
      par.epiJ=1;
      par.epiM=1;
      par.gamma_circle=4.25;
      par.epiJelse = par.gamma_circle + par.J_L / 2 + par.epiJ / 2;
    }


    process_population(networks);
    par.J_L+=0.25;
  }

  
  // finished
  par.CleanUp();

  return 0;
}
