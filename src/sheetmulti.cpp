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
#include <sys/stat.h>

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

#include <sys/stat.h>
#include <cstring>

using namespace std;


void WriteData(const vector<pair<int, double>>& shapedata, const string& oname)
{
  ofstream outfile;
  outfile.open(oname, ios::app);  // Append mode

  int phase_counts{};

  // First, find the maximum number of rows required
  int max_rows = 0;
  vector<int> rows{};

  for (const auto& [index, value] : shapedata) 
  {
    if (index + 1 > max_rows) 
    {
        max_rows = index + 1;
        rows.push_back(index);
    }     
  }

  // Write the header
  outfile << fixed << setprecision(6);
  
  for (int &row : rows) 
  {
    
    outfile << row;

    // Calculate the average for this row if there are matching pairs
    vector<double> values;
    for (const auto& [index, value] : shapedata) 
    {
      if (index == row) 
      {
        values.push_back(value);
      }
    }

    if (!values.empty()) 
    {
      double sum = accumulate(values.begin(), values.end(), 0.0);
      double mean = sum / values.size();
      outfile << "\t" << mean;  // Output the mean in the second column
    } 
    else 
    {
      cout << "Error in time output" << endl;
      outfile << "\t";  // No data for this row, leave empty
    }

    outfile << endl;  // Newline after each row
  }

  outfile.close(); 
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

    CPM->FillGrid();
    CPM->ConstructInitCells(*this);
    
    // for (int i=0;i<par.divisions;i++) {
    //   CPM->DivideCells();
    // }

    if (par.do_voronoi)
    {
      par.highT=false;
      CPM->Voronoi();
      par.start_sheet_measure=0;
    }
    else
    {
      CPM->Voronoi();
    }
    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();

    
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

void process_population()
{

  vector<vector<double>> hex_order(par.n_orgs);
  vector<vector<double>> shape_index(par.n_orgs);  
  vector<vector<double>> cell_displacements;

  Dish *dishes = new Dish[par.n_orgs];

  omp_set_num_threads(par.n_orgs);
  #pragma omp parallel for
  for (int i = 0; i < par.n_orgs; ++i)
  {
    dishes[i].CPM->set_num(i + 1);
    // does init block above.
    dishes[i].Init();



    if (par.sheetmix)
    {
      dishes[i].CPM->StartSheetTypes();

    }
    if (par.highT)
    {
      // double diff = par.highT_temp - par.T;
      // double toT = par.highT_temp - diff * (double(t) / double(par.highT_time));
      double toT = par.highT_temp;
      dishes[i].CPM->CopyProb(toT);
      dishes[i].CPM->Set_J(0.5);
      dishes[i].CPM->set_mixJ(1);      
    }
    else
    {
      dishes[i].CPM->CopyProb(par.T);
      dishes[i].CPM->Set_J(par.sheet_J);
      dishes[i].CPM->set_mixJ(par.sheetmixJ);
    }


    int t;

    for (t = 0; t < par.mcs; t++)
    {              
      // if (par.highT && t < par.highT_time)
      // {
      //   // double diff = par.highT_temp - par.T;
      //   // double toT = par.highT_temp - diff * (double(t) / double(par.highT_time));
      //   double toT = par.highT_temp;
      //   dishes[i].CPM->CopyProb(toT);
      // }
      if (par.highT && t==par.highT_time)
      {
        dishes[i].CPM->CopyProb(par.T);
        dishes[i].CPM->Set_J(par.sheet_J);
        dishes[i].CPM->set_mixJ(par.sheetmixJ);
      }

      if (par.velocities)
      {
        dishes[i].CPM->RecordMasses();
      }

      if (par.sheetmix)
      {
        dishes[i].CPM->RandomSheetType();
      }

      // if (t % 10 == 0 && t >= par.start_sheet_measure && t<= par.end_sheet_measure && par.sheet_hex)
      // {
      //   dishes[i].CPM->initVolume();
      //   dishes[i].CPM->adjustPerimeters();
      //   // vector<double> tperims = dishes[i].CPM->PerimitersRadiusN(sqrt(13));
      //   vector<double> tperims = dishes[i].CPM->TruePerimeters();
      //   vector<double> volumes = dishes[i].CPM->GetVolumes();

      //   vector<double> tadhesion = dishes[i].CPM->TrueAdhesion();

      //   double avg{};
      //   for (int j = 0; j < tperims.size(); ++j)
      //   {
      //     double sindex = tperims[j] / sqrt(volumes[j]);
      //     // cout << i << '\t';
      //     avg+=sindex;
      //     shape_index[i].push_back(sindex);
      //   }
      //   avg/=tperims.size();

      //   // now do hexatic order
      //   vector<double> hexes = dishes[i].CPM->GetHexes();
      //   for (auto &j : hexes)
      //   {
      //     hex_order[i].push_back(j);
      //   }
      // }

      if (t % 1 == 0 && t >= par.start_sheet_measure && t<= par.end_sheet_measure && par.sheet_hex)
      {
        dishes[i].CPM->ShapeOrder(t);
        dishes[i].CPM->HexaticOrder(t);
      }



      dishes[i].CPM->AmoebaeMove(t);
    }
  }

  if (mkdir(par.data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;




  if (par.output_sizes)
  {
    for (int i = 0; i < par.n_orgs; ++i)
    {
      vector<vector<double>> displc = dishes[i].CPM->ReturnMSD();
      for (auto i : displc)
        cell_displacements.push_back(i);
    }
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << par.T << "-" << par.sheet_J;
    string s = stream.str();

    string var_name = par.data_file + "/msd" + s + ".dat"; 
    ofstream outfile;
    outfile.open(var_name, ios::app);  

    int timer = 1;
    for (auto j=0;j<par.mcs-par.equilibriate-1;++j)
    {
      double msd=0;
      int n = 0;

      for (auto i=0;i<cell_displacements.size();++i)
      {
        msd += cell_displacements[i][j];
        ++n;
      }
      // outfile << timer << "\t" << msd/((double)n) << endl;
      outfile << (msd)/((double)n) << endl;

      ++timer;
    }

    outfile.close();
  }

  // if (par.sheet_hex)
  // {

  //   std::stringstream stream;
  //   stream << std::fixed << std::setprecision(1) << par.sheet_J;
  //   string s = stream.str();
  //   string var_name = par.data_file + "/shapes-" + s + ".dat"; 
  //   ofstream outfile;
  //   outfile.open(var_name, ios::app);  

  //   // Print the vector of vectors as columns
  //   for (size_t i = 0; i < shape_index.size(); i++) 
  //   {
  //     for (size_t j = 0; j < shape_index[i].size(); j++) 
  //     {

  //       outfile << shape_index[i][j] << endl;
  //     }
  //   }
  //   outfile.close();  


  //   // now do same for hexatic order
  //   var_name = par.data_file + "/hex-order-" + s + ".dat";
  //   outfile.open(var_name, ios::app);  
  //   // Print the vector of vectors as columns
  //   for (size_t i = 0; i < hex_order.size(); i++) 
  //   {
  //     for (size_t j = 0; j < hex_order[i].size(); j++) 
  //     {

  //       outfile << hex_order[i][j] << endl;
  //     }
  //   }
  //   outfile.close();  
  // }

  if (par.sheet_hex)
  {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << par.sheet_J;
    string s = stream.str();
    vector<pair<int,double>> hexdata = dishes[0].CPM->Get_sheet_hexatic_order();
    vector<pair<int,double>> shapedata = dishes[0].CPM->Get_sheet_shape_index();

    for (int i = 1; i < par.n_orgs;++i)
    {
      vector<pair<int,double>> next = dishes[i].CPM->Get_sheet_hexatic_order();
      for (auto&kv : next)
      {
        hexdata.push_back(kv);
      }
      vector<pair<int,double>> shape_next = dishes[i].CPM->Get_sheet_shape_index();
      for (auto&kv : shape_next)
      {
        shapedata.push_back(kv);
      }
    }
    string oname = par.data_file + "/hex_time-" + s + ".dat";
    WriteData(hexdata, oname);

    oname = par.data_file + "/shape_time-" + s + ".dat";
    WriteData(shapedata, oname);       
  }




  delete[] dishes;

}




int main(int argc, char *argv[]) {

  
  par.graphics=false;
  par.contours=false;
  par.print_fitness=true;
  par.randomise=false;
  par.gene_output=false;
  par.gene_record=false;
  // par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);
  par.mcs=100000 + par.equilibriate;
  par.sizex=200;
  par.sizey=200;
  par.end_program=0;
  par.sheet=true;
  par.periodic_boundaries=true;

  par.n_orgs = 120;

  par.velocities=false;
  par.sheet_hex=true;
  par.measure_interval=10;
  if (par.velocities)
    par.output_sizes = true;
  else
    par.output_sizes = false;
  Parameter();


  if (par.sheet_hex)
  {
    par.mcs = par.end_sheet_measure + par.equilibriate;
  }

  par.periodic_boundaries = true;
  par.flush_cells = true;


  bool single=false;
  if (single)
  {
    par.sheet_J = 4;
    process_population();
  }
  else
  {
    vector<double> Jlist;
    int n_trials = ceil((par.sheet_maxJ-par.sheet_minJ)/par.J_width)+1;
    for (int i = 0; i < n_trials; ++i)
    {
      double J = par.sheet_minJ + par.J_width*i;
      Jlist.push_back(J);
    }

    vector<vector<double>> hex(par.n_orgs);
    vector<vector<double>> shape(par.n_orgs);

    vector<vector<vector<double>>> hex_order{};
    vector<vector<vector<double>>> shape_index{};

    for (int i = 0; i < n_trials; ++i)
    {
      hex_order.push_back(hex);
      shape_index.push_back(shape);
      par.sheet_J = Jlist[i];
      if (par.sheetmix)
        par.sheetmixJ = 2*par.sheet_J;
      cout << par.sheet_J << endl;
      process_population();
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
