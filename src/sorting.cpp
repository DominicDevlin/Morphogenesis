/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
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
#include "storage.h"
#include "connections.h"
#include "fft.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

#include <sys/stat.h>
#include <cstring>

using namespace std;




INIT 
{
  try 
  {
    CPM->set_seed();
    CPM->set_datafile(par.data_file);
    // Define initial distribution of cells
    if (par.make_sheet)
    {
      CPM->ConstructSheet(par.sheetx,par.sheety);
      par.divisions = 6;
    }
      
    else
      CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
    CPM->ConstructInitCells(*this);
    if (par.velocities)
      par.output_sizes = true;
    
    // If we have only one big cell and divide it a few times
    // we start with a nice initial clump of cells. 
    // 
    // The behavior can be changed in the parameter file using 
    // parameters n_init_cells, size_init_cells and divisions
    for (int i=0;i<par.divisions;i++) {
      CPM->DivideCells();
    }
    
    // Assign a random type to each of the cells
    CPM->SetRandomTypes();

    CPM->start_network(par.start_matrix);


    par.print_fitness = true;
    par.node_threshold = 0;// int(floor((par.mcs - par.adult_begins) / 40) * 2 * 10);

    if (par.set_colours)
    {
      CPM->SetColours();
    }


    if (par.store)
    {
      if (mkdir("data_film", 0777) == -1)
        cerr << "Data film 2 " << strerror(errno) << endl;
      else
        cout << "data_film created." << endl;       
    }

    
  } 
  catch(const char* error) 
  {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);

  }

}


void addEdge(map<int, set<int>>& graph, int u, int v) {
    graph[u].insert(v);
    // graph[v].insert(u); // Remove this line if the graph is directed
}

bool dfs(map<int, set<int>>& graph, int current, int target, set<int>& visited) {
    if (current == target) return true; // target node found
    visited.insert(current); // mark the current node as visited
    for (int neighbor : graph[current]) {
        if (visited.find(neighbor) == visited.end()) { // if neighbor hasn't been visited
            if (dfs(graph, neighbor, target, visited)) return true;
        }
    }
    return false;
}

bool isPathExists(vector<pair<int, int>>& edges, vector<int>& nodes, int start, int end) {
    // Build the graph
    map<int, set<int>> graph;
    for (auto edge : edges) {
        addEdge(graph, edge.first, edge.second);
    }
    
    // Set of visited nodes
    set<int> visited;
    
    // Start DFS from the start node
    return dfs(graph, start, end, visited);
}





TIMESTEP { 
 
  try {

    static int t=0;
 
    static Dish *dish=new Dish();
    
    if (t < 1)
    { 
      cout << "calling init" << endl;
      dish->Init();
      dish->CPM->CopyProb(par.T);
      
    }
    static Info *info=new Info(*dish, *this);
    
    // record initial expression state. This occurs before any time step updates. 
    if (t == 100)
    {
      if (par.flush_cells)
      {
        dish->CPM->SetAllStates();
        dish->PDEfield->FlushGrid();
      }
      if (par.gene_output)
        dish->CPM->record_GRN();

      if (par.output_init_concs)
        dish->CPM->OutputInitConcs();
    }
      

    // programmed cell division section
    if (t < par.end_program)
    {

      if (t % par.div_freq == 0 && t <= par.div_end && !par.make_sheet)
      {
        dish->CPM->Programmed_Division(); // need to get the number of divisions right. 
      }


     
      if (t >= par.begin_network && t % par.update_freq == 0)
      {

        dish->CPM->update_network(t);
        dish->AverageChemCell(); 
        
        if (par.gene_output)
          dish->CPM->record_GRN();   

        // nop point recording velocities here?
        // if (par.velocities)
        // {
        //   dish->CPM->RecordMasses();
        // }
        if (par.output_gamma)
          dish->CPM->RecordGamma(); 

        // speed up initial PDE diffusion
        for (int r=0;r<par.program_its;r++) 
        {
          dish->PDEfield->Secrete(dish->CPM);
          dish->PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        } 
   
      }
    }
    else
    {
      if (t % par.update_freq == 0)
      {
        dish->CPM->update_network(t);
        // if (par.noise && t > par.noise_start)
          // dish->CPM->add_noise();
          
        dish->AverageChemCell(); 
        if (par.gene_output)
        {
          dish->CPM->record_GRN();
          dish->CPM->CountTypesTime();
        }

        if (par.output_gamma)
        {
          dish->CPM->RecordGamma();
        }

      }
      // if (par.insitu_shapes && t % 500 == 0)
      // {
      //   dish->CPM->ShapeIndexByState();
      //   dish->CPM->AdhesionByState();
      // }
      
      if (par.velocities)
      {
        dish->CPM->RecordMasses();
      }

      if (par.output_sizes)
      {
        dish->CPM->RecordSizes();
      }

      // if (par.gene_record)
      // {
      //   dish->CPM->RecordTypes();
      // }

      // if (par.gene_record && t > par.adult_begins)
      // {
      //   dish->CPM->RecordAdultTypes();
      // }
      for (int r=0;r<par.pde_its;r++) 
      {
        if (!par.hold_morph_constant)
        {
          dish->PDEfield->Secrete(dish->CPM);
          dish->PDEfield->Diffuse(1); // might need to do more diffussion steps ? 
        }

      }

      // print individual chemical concentrations. 
      // if (t % 100 == 0 && t > par.end_program)
      // {
      //   dish->PDEfield->print_concentrations(dish->CPM);
      // }

      dish->CPM->CellGrowthAndDivision(t);
    }
    dish->CPM->AmoebaeMove(t);
    if (t == par.mcs-1 && par.gene_output)
    {

      if (mkdir(par.data_file.c_str(), 0777) == -1)
        cerr << "Error : " << strerror(errno) << endl;
      else
        cout << "Directory created." << endl;  
      dish->CPM->print_cell_GRN();
    }


    // fitness stuff
    if (t > par.mcs * par.fitness_begin && t % par.fitness_typerate == 0)
    {
      dish->CPM->update_fitness();
    }
    


    if (t == par.mcs - 1)
    {

      cout << "n stable types: " << dish->CPM->CountStableTypes() << endl;

      if (par.output_gamma)
        dish->CPM->OutputGamma();

      if (par.output_sizes)
      {
        dish->CPM->OutputSizes();
        dish->CPM->Vectorfield();
      }
        
      // if (par.umap)
      dish->CPM->ColourIndex();

      dish->CPM->get_fitness();

      map<int, int> phens = dish->CPM->get_phenotype_time();
      map<int, int> types = dish->CPM->get_AdultTypes();  

      map<pair<int,int>,int> edge_tally{};
      
      // check if there are super long cycles. Need to account for this tiny edge case where there is a >3000 mcs cycle (very annoying)
      bool cycling = false;//dish->CPM->CycleCheck();
      // if (cycling && par.cycle_check)
      // {
      //   cout << "There is cycling!!" << endl;
      //   dish->CPM->set_long_switches(edge_tally);
      // }
      // else
      // {
      dish->CPM->set_switches(edge_tally);


      vector<vector<int>> scc;
      if (par.insitu_shapes)
      {
        par.node_threshold = 0;
        par.prune_edges = true;
      
        auto edge_copy = edge_tally;
        for (auto it = edge_copy.begin(); it != edge_copy.end();)
        {
          if (it->second < par.prune_amount)
          {
            it = edge_copy.erase(it);
          }
          else
          {
            ++it;
          }
        }

        map<int,int>subcomps{};
        Graph ungraph(types.size());

        subcomps = ungraph.CreateUnGraph(phens, types, edge_copy);
        scc = ungraph.GetComps(types, 2000);
        for (auto i : scc)
        {
            cout << "component: ";
            for (int j : i)
            {
                cout << j << "  ";
            }
            cout << std::endl;
        }
        vector<int> remaining_nodes;
        for (vector<int> &i : scc)
        {
          if (i.size() < 1)
            cerr << "error in scc size!\n";
          remaining_nodes.push_back(i[0]);
        }
        vector<vector<int>> dendrogram(remaining_nodes.size());
          
        vector<pair<int,int>> edges{};
        vector<int> nodes{};
        for (auto i : edge_copy)
        {
          if (i.second > 10)
          {
            edges.push_back(i.first);
            if (std::find(nodes.begin(), nodes.end(), i.first.first) == nodes.end())
            {
              nodes.push_back(i.first.first);
            }
            if (std::find(nodes.begin(), nodes.end(), i.first.second) == nodes.end())
            {
              nodes.push_back(i.first.second);
            }
          }

        }

        for (int i = 0; i < remaining_nodes.size(); i++)
        {
          for (int j = 0; j < remaining_nodes.size(); j++)
          {
            if (j != i)
            {
              // cout << i << '\t' << j << endl;
              int start = remaining_nodes[i];
              int end = remaining_nodes[j];
              // cout << start << '\t' << end << endl;
              // Check if there is a path from start to end
              if (isPathExists(edges, nodes, start, end)) 
              {
                dendrogram[i].push_back(j);
                // std::cout << "There is a path between " << start << " and " << end << std::endl;
              } 
              else 
              {
                // std::cout << "No path exists between " << start << " and " << end << std::endl;
              }
            }
          }
        }
        vector<int> dead_ends{};
        for (int i = 0; i < dendrogram.size(); i++)
        {
          if (dendrogram[i].size() == 0)
          {
            cout << "Dead end at " << remaining_nodes[i] << endl;
            dead_ends.push_back(i);
          }
        }
        set<int> dead_ends_set(dead_ends.begin(), dead_ends.end());

        // Filter dendrogram
        auto new_end = std::remove_if(dendrogram.begin(), dendrogram.end(), [&dead_ends_set](const std::vector<int>& pair) 
        {
          // Check if any element of the pair is not in dead_ends_set
          for (int elem : pair) 
          {
            if (dead_ends_set.find(elem) == dead_ends_set.end()) 
            {
              return true; // remove the element in dendrogram
            }
          }
          return false; // keep the element from dendrogram
        });

        // Erase the removed elements
        dendrogram.erase(new_end, dendrogram.end());

        int max_diffs=0;

        for (int i = 0; i < dendrogram.size(); ++i)
        {
          cout << "SCC differentiates into: ";
          for (int j = 0; j < dendrogram[i].size(); ++j)
          {
            cout << remaining_nodes[dendrogram[i][j]] << '\t';
          }
          cout << endl;
          if (dendrogram[i].size() > max_diffs)
            max_diffs = dendrogram[i].size();
        }
      }



      if (par.potency_edges)
      {
        // entire program is run from ungraph now
        map<int,int>subcomps{};
        if (cycling && par.cycle_check)
        {
          Graph ungraph(phens.size());
          subcomps = ungraph.CreateUnGraph(phens, phens, edge_tally, 1, true);          
        }
        else
        {
          Graph ungraph(types.size());
          par.node_threshold = int(floor((par.mcs - par.adult_begins) / 40) * 20);
          par.prune_edges = false;
          subcomps = ungraph.CreateUnGraph(phens, types, edge_tally);
          vector<vector<int>> result = ungraph.GetComps(types);

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
        // dish->CPM->CellVelocities();
        // dish->CPM->scc_momenta(scc);
        // dish->CPM->momenta();
        // dish->CPM->diff_anisotropy(scc);
        if (par.division_anisotropy)
          dish->CPM->division_anisotropy(scc);
      }
        
      // dish->CPM->SpecialVelocity();
      if (par.record_directions)
      {
        dish->CPM->Directionality();
        // dish->CPM->SingleCellDirection();
      }
   
    }

     
    //printing every 1000 steps. Do other debugging things here as well. 
    if (t % 1000 == 0)
    {
      dish->CPM->print_random_cell();
      cout << "Number of cell types: " << dish->CPM->get_ntypes() << endl;
      cout << t << " TIME STEPS HAVE PASSED." << endl;
      dish->CPM->PrintPhenotypes();
      // dish->CPM->WhiteSpace();
      // dish->CPM->DeviationFromCircle();
      cout << "AVG BINDING: " << dish->CPM->AverageBinding() << endl;
      cout << "NUMBER OF MEDIUM PROTEINS ON AVG: " << dish->CPM->AvgMedsOn() << endl;
      cout << "ORG MASS IS: " << dish->CPM->Mass() << endl;
      double center[] = {0.0,0.0};
      dish->CPM->get_center(center);
      cout << "x center: " << center[0] << "   y center: " << center[1] << endl;
      dish->CPM->PrintColours();
      // dish->CPM->prop_success();
    }

    // used to create morphogen stuff
    if (t==3000)
    {
      // dish->PDEfield->PrintAxisConcentrations(true, 125);
      // dish->CPM->OutputProteinNorms();
    }


    if (t >= 6000 && t < 8000 && t % 40 == 0 && par.scramble)
      dish->CPM->swap_cells();


    // for spawning a lot of morphogen at a specific point, or changing cell types + morphogen

    if (par.convert_cells && par.convert_time == t)
    {
      // dish->CPM->ConvertToStem(par.convert_x,par.convert_y,par.convert_size,par.convert_to_type, dish->PDEfield, true, par.clear_radius);
      // dish->CPM->ConvertToStem(100,230,par.convert_size,par.convert_to_type, dish->PDEfield, true, par.clear_radius);  
    }


    if (t == 6998)
    {
      // dish->CPM->ConvertToStem(125,95,40,11907, dish->PDEfield, true, 45);
      // dish->IntroduceMorphogen(1, 120, 90);
    }

    // for removing cells. 
    if (t == 12000)
    {
      // dish->CPM->DestroyCellsByRadius(34.);
    }

    if (t == 12000)
    {
      // fft test;
      // test.AllocateGrid(par.sizex, par.sizey);
      // test.ImportGrid(dish->CPM->ReturnGrid());//, dish->CPM);
      // test.PolarTransform();
      // test.PolarToOutput();

      // test.ShiftGrid(test.GetPolar(), 8);

      // test.PolarToOutput("polar-shift.png");

      // test.ReflectGrid(test.GetPolar());
      // test.PolarToOutput("polar3.png");

    }





    static vector<bool> con_plot(4, false);
    static vector<int> con_colours{5, 7, 14, 37};

    //cerr << "Done\n";
    if (par.graphics && t%5==0)// !(t%par.screen_freq)) 
    {
      
      BeginScene();
      ClearImage();

      // Plot the dish. 
      dish->Plot(this);
      
      // static vector<array<int,2>> perim;
      // static vector<int> pcells;
      // if (t > 1000 && t % 500 == 0)
      // {
      //   // perim = dish->CPM->PerimeterCAC();
      //   // pcells = dish->CPM->CellsFromCAC(perim);
      //   // pcells = dish->CPM->LinkPerimeter();
      //   // dish->CPM->AngleCurvature();
      // }
      
      // if (t > 1500)
      // {
      //   // dish->CPM->DrawListofCAC(this, perim);
      //   // dish->CPM->DrawPerimeter(this, pcells);
      // }


      // static bool c4 = false;

      if (t>0 && t % par.begin_network == 0)
      {
        for (int c=0;c<par.n_diffusers;++c)
          con_plot[c] = dish->PDEfield->CheckSecreting(c);
      }

      if (t>par.end_program && par.contours)
      {
        
        for (int c=0;c<par.n_diffusers;++c)
        {
          if (con_plot[c])
            dish->PDEfield->ContourPlot(this,c,con_colours[c]);
        }

        // if (c1)
        //   dish->PDEfield->ContourPlot(this,0,5);
        // if (c2)
        //   dish->PDEfield->ContourPlot(this,1,7);
        // if (par.n_diffusers > 2)
        // {
        //   if (c3)
        //     dish->PDEfield->ContourPlot(this,2,14);
          // if (c4)
          //   dish->PDEfield->ContourPlot(this,3,37);
        // }


        // this function plots a shade for the PDE field and is very computationally costly. Isn't needed. 
        // dish->PDEfield->Plot(this,0);
        // dish->PDEfield->Plot(this,1);
        // dish->PDEfield->Plot(this,2);
      }


    
      char title[400];
      snprintf(title,399,"CellularPotts: %.2f hr",dish->PDEfield->TheTime()/3600);      



      //char title[400];
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      EndScene();
      info->Menu();
     
    }
  
    // storage function. 
    if (par.store && !(t%par.storage_stride))//  || t == 3041) 
    {
      char fname[200];
      sprintf(fname,"%s/extend%07d.png",par.datadir,t);
    
      BeginScene();
      ClearImage();    
      dish->Plot(this);

      if (t>par.end_program && par.contours)
      {

        if (t>0 && t % par.begin_network == 0)
        {
          for (int c=0;c<par.n_diffusers;++c)
            con_plot[c] = dish->PDEfield->CheckSecreting(c);
        }

        if (t>par.end_program && par.contours)
        {
          
          for (int c=0;c<par.n_diffusers;++c)
          {
            if (con_plot[c])
              dish->PDEfield->ContourPlot(this,c,con_colours[c]);
          }
            
        }

        // c1 = dish->PDEfield->CheckSecreting(0);
        // c2 = dish->PDEfield->CheckSecreting(1);
        // if (par.n_diffusers > 2)
        // {
        //   c3 = dish->PDEfield->CheckSecreting(2);
        //   // c4 = dish->PDEfield->CheckSecreting(3);
        // }
        // if (c1)
        //   dish->PDEfield->ContourPlot(this,0,293);
        // if (c2)
        //   dish->PDEfield->ContourPlot(this,1,291);
        // if (c3)
        //   dish->PDEfield->ContourPlot(this,2,292);
      }

      
      EndScene();
    
      Write(fname);
        
    }

    t++;
  } 
  catch(const char* error) 
  {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }


}



int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}




int main(int argc, char *argv[]) {
  
	try 
  {

#ifdef QTGRAPHICS
    QApplication a(argc, argv);
#endif
        
    Parameter();
    // Read parameters
    bool read = false;
    if (read)
      par.Read(argv[1]);
    // Seed(par.rseed);
    
    //QMainWindow mainwindow w;
#ifdef QTGRAPHICS
    QtGraphics g(par.sizex*2,par.sizey*2);
    //a.setMainWidget( &g );
    a.connect(&g, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );

    if (par.graphics)
      g.show();
    
    a.exec();
#else
    X11Graphics g(par.sizex*2,par.sizey*2);
    int t;

    for (t=0;t<par.mcs;t++) {

      g.TimeStep();
    
    }
#endif
    
  } catch(const char* error) {
    std::cerr << error << "\n";
    exit(1);
  }
  catch(...) {
    std::cerr << "An unknown exception was caught\n";
  }
  return 0;
}
