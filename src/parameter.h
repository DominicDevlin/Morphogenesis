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
#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <iostream>
#include <vector>
using namespace std;
class Parameter {
  
 public: 
  Parameter();
  ~Parameter();
  // void inject();
  void CleanUp(void);
  void Read(const char *filename);
  void Write(ostream &os) const;
  double T;
  int target_area;
  int target_length;
  double lambda;
  double lambda2;
  char * Jtable;
  int conn_diss;
  bool vecadherinknockout;
  bool extensiononly;
  int chemotaxis;
  int border_energy;
  int neighbours;
  bool periodic_boundaries;
  int n_chem;
  double * diff_coeff; //DIFFUSER CORRECTION NEEDED
  double * decay_rate;
  double * secr_rate;
  double saturation;
  double dt;
  double dx;
  int pde_its;
  int n_init_cells;
  int size_init_cells;
  int sizex;
  int sizey;
  int divisions;
  int mcs;
  int evs;
  int rseed;
  double subfield;
  int relaxation;
  int storage_stride;
  bool graphics;
  bool contours;
  bool screen_freq;
  bool store;
  char * datadir;
  double tlength1;
  double tlength2;


  int n_orgs;
  int n_replicates;

  int n_pred;
  int n_mutations;
  double mut_rate;
  double polm_rate;


  //programmed division parameters
  int end_program;
  int begin_network;
  int div_freq;
  int div_end;
  int update_freq;
  double eT;
  double lT;

  bool morphogen;
  double morphdecay;

  // when to collect fitness
  double fitness_begin;
  int fitness_typerate;
  int fitness_somrate;
  bool print_fitness;

  bool all_divide;

  bool potency_edges;
  bool prune_edges;
  double node_percent;
  int node_threshold;
  int adult_begins;

  bool insert_randoms;

  int n_diffusers;
  int n_locks;
  int med_table[5] = {5,4,3,2,1};

  int div_threshold;
  int gthresh;
  int stem_gthresh;

  int shrink;
  int s_shrink;
  

  bool scramble;

  int max_on; //max number of proteins on before a cell can't divide.
  int max_cells; //max cells as a cap (slows down sims, dont expect to reach this number)

  vector<vector<int>> start_matrix; 
  vector<bool> start_polarity;
  bool polarity_on;
  
  int n_genes; // number of genes. not all genes interact with gene network (EXCLUDES LOCKS AND KEYS)
  int n_activators; //number of genes that can activate network (<= n_genes)
  int n_functional;
  int gene_vector_size;
  int n_lockandkey; // number of lock and keys (==), stored in separate vector for ease
  int n_length_genes;
  int n_mediums;
  int n_TF; 
  int n_stem;
  int n_MF;
  int cycle_size; // number of gene sets allowed in the vector "cycle" to check for gene cycling. 
  int cycle_threshold; // number of different cycles allowed before it can count for fitness. 

  int tloc1;
  int tloc2;

  int n_in;
  int n_gr;

  int gloc;
  int shrink_loc;
  bool shrink_on;

  bool enzymes;
  int e1_loc;
  int e2_loc;
  double reaction_rate;

  bool randomise;

  double theta;
  double d_rate;

  vector<double> new_g; 

  bool gene_output;
  bool gene_record;
  int start_record;
  
  int min_contig;

  double morphogen_dose;


  int div1;
  int div2;
  int begin_movement;
  int program_its;

  bool asymmetry_selection;
  bool asym_only;
  double swap_selection;
  bool growth_selection;

  bool starter;
  vector<vector<int>> start_n;

  bool set_colours;

  uint64_t pickseed;

 private:
};

ostream &operator<<(ostream &os, Parameter &p);
const char *sbool(const bool &p);


#endif
