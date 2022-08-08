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


#include "parameter.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include "output.h"
#include "parse.h"
#include <cmath>

Parameter::Parameter()
{
  // init params
  target_area = 4900;
  size_init_cells = 70; // this is equal to the radius of the circle (done by eden growth). 
  n_init_cells = 1;
  divisions = 0;
  
  //most important parameters.
  // T is defunct in current simulation, using eT and lT instead. 
  eT = 2; // temperature during programmed divisions 
  lT = 2; // temperature during development
  T = 2;
  evs = 10000;
  mcs = 12100;
  target_length = 0;
  lambda = 0.5;
  lambda2 = 0.1; // I think having this too high increases the volume too quickly.
  gthresh = 4; // tau used by Paulien. Want growth to be by squeezing and not temperature fluctuations. 
  stem_gthresh = 3; // going to have a lower growth threshold for stem cells. 

  graphics = true;
  contours = false;

  // Make sure this is off for evolution!!
  randomise = false;
  // output the entire network states over the course of a single simulation (off for evolution!)
  gene_output = false;

  // let all cells be capable of cell division. MUST BE FALSE FOR stem cell evolution testing. False promotes stem cell evolution?
  all_divide = false;

  // KEEP THIS TO FALSE FOR EVOLUTION
  print_fitness = true;

  // morphogen wave at the end of programmed division. 
  // The morphogen occupies a hidden spot at the back of the cell "genes" vector. This allows it to decay but not be increased. 
  morphogen = false;


  // If starting with a specific matrix for evolution, use the function in main(). 
  start_matrix = { { 0, -1, 0, 0, 0, 0, 0, 0 }, { -1, 0, 0, 0, -1, 0, -1, 0 }, { 0, 0, 0, 0, 0, 1, -1, 0 }, { -1, 0, 0, 0, -1, 1, 1, 0 }, { 1, -1, 0, 0, 1, 0, 1, 0 }, { 1, -1, 1, 1, -1, 0, 0, 0 }, { -1, 0, 0, 0, -1, 1, 1, 1 }, { 0, 1, 0, 1, 0, -1, -1, 0 }, { 0, -1, 0, 2, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 1, 0, -1 }, { -1, -1, 1, -1, 0, 0, 1, -2 }, { 0, -1, 0, 0, 1, 0, -1, 0 }, { -1, 0, 0, 0, 0, 1, 0, 1 }, { 0, 0, 0, 0, 0, 2, -1, 0 }, { 0, 0, 1, 1, 0, 0, 1, -1 }, { 0, 1, 0, -1, -1, 0, -1, -1 }, { 0, 0, 0, 1, 0, 1, 1, -1 }, { 0, 0, 0, 0, 1, 0, 0, 2 }, { 0, 0, 0, 0, -1, 0, -1, 0 }, { 0, 1, 0, 0, 0, 2, -1, -1 }, };
  





  // polarity matrix. Used for transcription factors being at one side of the cell upon reproduction. 
  start_polarity = { 0, 0, 1, 0 };
  polarity_on = true;


  // params for evolution
  n_orgs = 16; // MUST BE MULTIPLE OF 4 !!
  n_mutations = 1;
  // mut rate for gene network
  mut_rate = 0.5;
  // mutation rate for polarities
  polm_rate = 0.2;

  //programmed division parameters
  end_program = 800;
  begin_network = 400;
  div_freq = 50;
  div_end = 300;
  // gene network update frequency
  update_freq = 40;


  // when to collect fitness
  fitness_begin = 0.8;
  // frequency of time steps fitness is collected
  fitness_typerate = 200;
  // frequency of time steps to calculate rate of somatic cell production
  fitness_somrate = 2000;




  // target length with 1 gene or 2 genes on. These are multipliers (area / tlength = true target length)
  tlength1 = 3;
  tlength2 = 2;


  //growing and dividng
  max_cells = 500;
  div_threshold = 100;
  // thresholds which cell has to be GREATER THAN before its target volume shifts to its actual volume. 
  
  shrink = -10;


  //basic simulation parameters. 
  sizex = 400;
  sizey = 400;
  rseed = -1;



  // gene network parameters
  n_diffusers = 2;
  n_lockandkey = 10; // number of lock and keys (==), stored in separate vector for ease
  n_locks = 5; // must be half lockandkey. locks = keys
  // number of transcription factors. The first two transcription factors must be ON for stem cell identity. 
  n_TF = 4; 
  n_stem = 2;
  n_length_genes = 2;
  n_MF = 2;

  n_genes = n_diffusers + n_lockandkey + n_TF + n_length_genes +n_MF; // number of genes. All gene types must sum to this value (except if using morphogenwave, then activators is +1).

  n_activators = n_diffusers + n_TF+n_MF + morphogen; //number of genes that can activate network (<= n_genes)

  // location of target length genes in genome. 
  tloc1 = n_diffusers + n_MF + n_TF;
  tloc2 = tloc1+1;

    // maximum number of somatic genes on before a cell cannot divide. HAVE CHANGED THIS TO TF GENES
  // max_on = 3;
  
  cycle_size = 8; // number of past states held by a cell for determining cell types. Selecting against cycling cell types. 
  cycle_threshold = 2; // Dom thinks it is okay to increase this from 2 to 3.  Change back if needed.   


  // start concentrations are (in order): diffusers(external), maternal factors, 
  // transcription factors, target length genes  -- the other networks hold lock then keys --

  // dont need this if starting all at 0. However, TF need to start at 1 to have some input into network for MF 0,0, otherwise will always be single cells at beginning.
  // new_g = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // 

  theta = -0.3;
  d_rate = exp(-0.29);
  morphdecay = exp(-0.2);





  /// Dom has deprecated. 
  Jtable = strdup("J.dat");

  // keep this at 2-3
  neighbours = 2;
  // not sure what this does
  conn_diss = 2000;

  // DONT TOUCH!
  vecadherinknockout = false;
  extensiononly = false;
  chemotaxis = 0;
  border_energy = 100;
  periodic_boundaries = false;

  // PDE field parameters. I have defined separate rates for each diffuser, however they are the same for now.
  // This may change in the future if I make one an evolvable parameter. 
  n_chem = 0; // Dom not currently using, instead using n_diffusers


  // Do not need to allocate memory from these because not calling from file (memory always allocated regardless, no need for pointer)
  // diff_coeff = new double[1];
  // decay_rate = new double[1];
  // secr_rate = new double[1];


  diff_coeff[0] = 8e-7; // Keeping it at this for now. Maybe this could be evolvable. 
  diff_coeff[1] = 8e-7;
  
  decay_rate[0] = 2e-3;
  decay_rate[1] = 2e-3;
  
  secr_rate[0] = 2.4e-3;
  secr_rate[1] = 2.4e-3;
  saturation = 0;
  dt = 1.0;
  dx = 1/((double)sizex);
  pde_its = 1;

  // n_chem = 1;
  // diff_coeff = new double[1];
  // diff_coeff[0] = 0.5e-13;
  // decay_rate = new double[1];
  // decay_rate[0] = 1.2e-4;
  // secr_rate = new double[1];
  // secr_rate[0] = 0.5e-4;
  // saturation = 0;
  // dt = 2.0;
  // dx = 2.0e-6;
  // pde_its = 15;


  
  subfield = 1.0;
  relaxation = 0;

  // storing images. 
  storage_stride = 20;
  screen_freq = 100;
  store = false;
  datadir = strdup("data_film2");
}

Parameter::~Parameter() {
  // destruct parameter object
  // free string parameter
  CleanUp();
}

void Parameter::CleanUp(void) {
  if (Jtable) 
     free(Jtable);
  // if (diff_coeff) 
  //    free(diff_coeff);
  // if (decay_rate) 
  //    free(decay_rate);
  // if (secr_rate) 
  //    free(secr_rate);
  if (datadir) 
     free(datadir);

}

void Parameter::Read(const char *filename) {
  
  static bool ReadP=false;

  if (ReadP) {

    //throw "Run Time Error in parameter.cpp: Please Read parameter file only once!!";
    CleanUp();
    
  } else
    ReadP=true;

  FILE *fp=OpenReadFile(filename);


  T = fgetpar(fp, "T", 50., true);
  target_area = igetpar(fp, "target_area", 100, true);
  target_length = igetpar(fp, "target_length", 60, true);
  lambda = fgetpar(fp, "lambda", 50, true);
  lambda2 = fgetpar(fp, "lambda2", 5.0, true);
  Jtable = sgetpar(fp, "Jtable", "J.dat", true);
  conn_diss = igetpar(fp, "conn_diss", 2000, true);
  vecadherinknockout = bgetpar(fp, "vecadherinknockout", false, true);
  extensiononly = bgetpar(fp, "extensiononly", false, true);
  chemotaxis = igetpar(fp, "chemotaxis", 1000, true);
  border_energy = igetpar(fp, "border_energy", 100, true);
  neighbours = igetpar(fp, "neighbours", 2, true);
  periodic_boundaries = bgetpar(fp, "periodic_boundaries", false, true);
  // n_chem = igetpar(fp, "n_chem", 1, true);
  // diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, true);
  // decay_rate = dgetparlist(fp, "decay_rate", n_chem, true);
  // secr_rate = dgetparlist(fp, "secr_rate", n_chem, true);
  // saturation = fgetpar(fp, "saturation", 0, true);
  // dt = fgetpar(fp, "dt", 2.0, true);
  // dx = fgetpar(fp, "dx", 2.0e-6, true);
  // pde_its = igetpar(fp, "pde_its", 15, true);
  n_init_cells = igetpar(fp, "n_init_cells", 100, true);
  size_init_cells = igetpar(fp, "size_init_cells", 10, true);
  sizex = igetpar(fp, "sizex", 200, true);
  sizey = igetpar(fp, "sizey", 200, true);
  divisions = igetpar(fp, "divisions", 0, true);
  mcs = igetpar(fp, "mcs", 10000, true);
  rseed = igetpar(fp, "rseed", -1, true);
  subfield = fgetpar(fp, "subfield", 1.0, true);
  relaxation = igetpar(fp, "relaxation", 0, true);
  storage_stride = igetpar(fp, "storage_stride", 10, true);
  graphics = bgetpar(fp, "graphics", true, true);
  store = bgetpar(fp, "store", false, true);
  datadir = sgetpar(fp, "datadir", "data_film", true);

}

const char *sbool(const bool &p) {

  const char *true_str="true";
  const char *false_str="false";
  if (p)
    return true_str;
  else
    return false_str;
}

void Parameter::Write(ostream &os) const {
  setlocale(LC_NUMERIC, "C");

  os << " T = " << T << endl;
  os << " target_area = " << target_area << endl;
  os << " target_length = " << target_length << endl;
  os << " lambda = " << lambda << endl;
  os << " lambda2 = " << lambda2 << endl;

  if (Jtable) 
    os << " Jtable = " << Jtable << endl;
  os << " conn_diss = " << conn_diss << endl;
  os << " vecadherinknockout = " << sbool(vecadherinknockout) << endl;
  os << " extensiononly = " << sbool(extensiononly) << endl;
  os << " chemotaxis = " << chemotaxis << endl;
  os << " border_energy = " << border_energy << endl;
  os << " neighbours = " << neighbours << endl;
  os << " periodic_boundaries = " << sbool(periodic_boundaries) << endl;
  os << " n_chem = " << n_chem << endl;
  os << " diff_coeff = "<< diff_coeff[0] << endl;
  os << " decay_rate = "<< decay_rate[0] << endl;
  os << " secr_rate = "<< secr_rate[0] << endl;
  os << " saturation = " << saturation << endl;
  os << " dt = " << dt << endl;
  os << " dx = " << dx << endl;
  os << " pde_its = " << pde_its << endl;
  os << " n_init_cells = " << n_init_cells << endl;
  os << " size_init_cells = " << size_init_cells << endl;
  os << " sizex = " << sizex << endl;
  os << " sizey = " << sizey << endl;
  os << " divisions = " << divisions << endl;
  os << " mcs = " << mcs << endl;
  os << " rseed = " << rseed << endl;
  os << " subfield = " << subfield << endl;
  os << " relaxation = " << relaxation << endl;
  os << " storage_stride = " << storage_stride << endl;
  os << " graphics = " << sbool(graphics) << endl;
  os << " store = " << sbool(store) << endl;

  if (datadir) 
    os << " datadir = " << datadir << endl;
}


ostream &operator<<(ostream &os, Parameter &p) {
  p.Write(os);
  return os;
}

Parameter par;
