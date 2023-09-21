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
    //basic grid parameters. 
    sizex = 250;
    sizey = 250;
    mcs = 12100;

    // show on screen
    graphics = false;
    // show morphogen gradients
    contours = false;

    // Generate a random genome
    randomise = false;

    // ANALYSIS PARAMS: note that there is slow down when these are turned on. 
    // output data for analysis (connectivity, gene expression, state transitions)
    gene_output = true;
    // gene record needs to always be on to test network connectivity. 
    gene_record = true;
    //for umap
    umap = true;

    // record velocities for all cells
    velocities = true;
    record_directions = false;

    // record cell sizes
    output_sizes = true;

    // record gamma parameter
    output_gamma = false;

    // read genomes from file
    file_genomes = false;

    //for storing images
    store = true;
    //set specific colours (SHOULD ALWAYS BE FALSE UNLESS NEEDED)
    set_colours = true;
    // Start from specific seed. USE 0 for random seed. (SHOULD ALWAYS be 0 unless need specific seed.)
    pickseed=13041561340170254908;
    rseed = -1;

    // KEEP THIS TO FALSE FOR EVOLUTION
    print_fitness = true; 


    // This start matrix is for sorting, overlap and transitions. For evolution start matrix, see start_n below 
    start_matrix = { { 1, 1, -1, 1, 0, -1, 0, 1, 1 }, { 1, 1, 0, 0, 0, 0, -1, 0, 1 }, { 0, 2, -1, 0, 2, 0, 0, 0, 0 }, { 0, 0, 0, 2, 0, 0, 0, 0, 0 }, { 0, 0, 1, -1, -1, 0, 1, 0, -1 }, { 0, 1, 1, 0, -1, 0, 0, 0, 0 }, { -1, -1, 1, 1, 1, 1, 1, 2, 0 }, { 1, -1, 0, 0, 0, 0, 0, 0, 2 }, { 1, 1, -1, 0, -1, 0, -1, 0, 0 }, { 2, 1, 1, 1, 0, 0, 0, 0, 0 }, { 0, 0, 1, 0, 0, -1, 1, 0, 0 }, { 0, 1, 0, -1, 1, -2, 0, 0, 0 }, { -2, 0, 0, 0, 0, 0, 1, 1, 1 }, { -1, 0, 0, 0, 0, 0, 1, -2, -1 }, { -1, -1, 0, 2, -1, -1, 2, 0, 0 }, { 1, -1, 0, 0, -1, 1, 0, 0, -1 }, { 0, 0, 0, 0, -1, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, -2, -2, -2 }, { 0, 1, 0, 0, 2, 1, 0, 0, 0 }, { 2, 0, 0, 0, 0, -1, -2, 1, 0 }, { 1, 0, 1, 0, 1, 2, 0, -1, 1 }, { 1, 1, -1, 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 0, -1, 0, 1, 0 }, { 0, 0, 0, 0, 0, -1, 0, 1, 0 }, { 0, -1, 0, -1, 0, 0, 0, 1, 0 }, { -1, -2, 1, 0, 0, 0, -1, 2, 2 }, { 0, 0, 0, 0, 0, 1, 0, 1, 2 }, };
    



    n_orgs = 60; // MUST BE MULTIPLE OF 4 !! OR 10 FOR POTENCY TESTn_orgs = 16; // MUST BE MULTIPLE OF 4 !! OR 10 FOR POTENCY TEST
    n_replicates = 2;
    // edges and nodes only at end of simulation.
    potency_edges = true;
    // what mcs to start measuring adult types & differentiation. set to 8000 for all results
    adult_begins = 8000;
    // prune tiny edges (<1 per org) from graph ( we keep this to false )
    prune_edges = false;
    // scramble cells
    scramble = false;
    // flat threshold for nodes
    node_threshold = int(floor((mcs - adult_begins) / 40) * 2 * 10 * n_orgs);
    // DEPRACATED - prune nodes below this percent. Should probably set this as a minimum value (i.e. 10 cells equivalent)
    node_percent = 0.03;

/*Conditions for evolution */
    // select for movement of cells towards one side of the boundary.
    asymmetry_selection = false; 
    asym_only = false;
    swap_selection = 240.; // the average fitness of population needed to switch from asym_only to asymmetry selection. 
    // start from a certain network
    growth_selection=false;
    elongation_selection = false;
    starter = false;
    
    start_n = { { 2, 1, -2, 0, -1, 0, 0, 2, 0 }, { 2, 0, 0, 0, 0, 0, 0, 0, -1 }, { 0, 0, -1, 0, -1, 1, -2, 0, 0 }, { 2, -1, 1, 1, 0, 0, 1, -2, -2 }, { -1, 0, 0, 1, 1, 1, 0, -1, 0 }, { 1, 2, 0, -1, 1, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 1, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 1, 0, -2 }, { 0, 1, 0, 0, -1, 0, 1, -1, 1 }, { 1, 1, 0, 0, 0, 1, 0, 0, 0 }, { 1, 0, 0, 0, 1, 0, 0, -1, 0 }, { 0, 1, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 1, -2, -1, 0, 0, 0 }, { 1, 1, 1, -1, 0, 0, 0, 0, -2 }, { 1, 0, 1, 2, 0, 0, 0, 0, 1 }, { 0, 0, -2, 1, 0, 0, 0, 0, 0 }, { 0, -1, 1, 0, 2, 0, 0, 0, 1 }, { 0, 0, 2, 0, 0, 2, -1, 0, 1 }, { 0, 1, 2, -1, 1, 0, 0, 1, 1 }, { 1, -2, 0, 0, -1, 0, -2, 0, 0 }, { 0, -1, 0, 1, 0, 1, 0, -2, 2 }, { 0, -1, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, -1, 0, 0, 1, 1, 2, 0 }, { 0, 2, -2, 0, -1, 1, 1, 0, 0 }, { -1, -1, 0, 0, 1, 0, 1, 0, 0 }, { 0, 1, -1, 0, 0, 1, 0, 1, 0 }, { 1, 0, 1, -2, 0, -2, -1, 1, -1 }, };


  

    evs = 10000;
    insert_randoms = true;
    n_mutations = 1;
    // mut rate for gene network
    mut_rate = 0.5;
    // mutation rate for polarities
    polm_rate = 0.2;
    n_pred = n_orgs / 2;


    // init params for organisms
    target_area = 4900;
    size_init_cells = 70; // this is equal to the radius(diameter?) of the circle (done by eden growth). 
    n_init_cells = 1;
    divisions = 0;
    
    // gene network update frequency
    update_freq = 40;

    //programmed division parameters
    end_program = 1600;
    begin_network = 75;
    div_freq = 50;
    begin_movement=1200;
    program_its = 8; // we are doing more PDE iterations during the program. 
    div_end = 250;

    // add noise to regulatory network 
    noise = false;
    // noise amount
    noise_dose=0.012;
    noise_start = 6000;

    colour_index = { {83351, 76}, {90351, 75}, {91375, 74}, {66943, 71}, {115983, 70}, {66847, 69}, {66831, 68}, {67047, 66}, {1527, 63}, {91623, 62}, {124167, 61}, {91407, 59}, {67063, 65}, {2039, 27}, {9623, 23}, {9631, 22}, {25991, 18}, {116103, 16}, {83431, 3}, {17799, 14}, {67583, 10}, {83331, 2}, {83447, 4}, {83211, 49}, {66823, 67}, {50567, 15}, {116099, 1}, {75159, 72}, {10135, 39}, {34183, 24}, {115971, 0}, {99715, 43}, {83967, 9}, {9607, 12}, {58759, 56}, {83311, 7}, {75263, 73}, {10239, 20}, {1415, 13}, {83203, 8}, {115975, 17}, {99719, 25}, {83327, 6}, {67071, 19}, {2047, 11}, {66975, 33}, {83231, 41}, {124295, 57}, {1431, 29}, {83455, 5}, {34179, 30}, {50563, 31}, {99591, 26}, {83335, 32}, {9727, 36}, {10143, 21}, {75167, 37}, {67575, 64}, {1535, 46}, {1943, 28}, {66967, 38}, {91399, 60}, {10119, 40}, {115979, 42}, {92143, 54}, {10231, 44}, {91527, 58}, {26503, 45}, {66959, 34}, {83215, 47}, {66951, 35}, {83207, 48}, {26095, 51}, {25999, 52}, {91535, 53}, {91631, 50}, {26607, 55}, };




    // show output of all comparisons for overlap. Only use when comparing a small number of organisms. 
    // NOTE - MUST RUN with tag: "-platform offscreen" when using cluster (there is no display).
    overlap_images = false;
    overlap_orgs = 60;
    // true = compare different genomes, false = compare same genomes
    between_org_overlap = false;

    // Basic Cellular Potts parameters
    eT = 3; // temperature during programmed divisions 
    lT = 3; // temperature during development
    T = 3;
    target_length = 0;
    lambda = 0.5;
    lambda2 = 0.1;


    // DEPRACATED: let all cells be capable of cell division. MUST BE FALSE FOR stem cell evolution testing. False promotes stem cell evolution?
    all_divide = true;


    //name of data file
    data_file = "org-data";

    //count stem and differentiated cells
    stem_counts = false;

    count_bud_cells = false;

    // print single cell proteins
    single_cell = false;
    single_type = 108034;
    // turn all cells into this state
    single_states = { 0.993083, 1.17567e-24, 0.99268, 0.993102, 0.374408, 1.15058e-11, 0.00707728, 1.29675e-08, 0.993101, 0.993102, 0.00164919, 1.35233e-07, 0.993102, 0.993102, 0.00190157, 0.993102, 2.67398e-06, 1.83131e-09, 0.658413, 0.98249, 1.73228e-11, 0.0127765, 7.57861e-08, 5.68932e-07, 5.82723e-12, 5.82716e-12, 3.21722e-09, 0.922665, 1.59416e-09, 0.763114,  };

    flush_cells = false;


    // when to collect fitness
    fitness_begin = 0.9;
    // frequency of time steps fitness is collected
    fitness_typerate = 100;
    // frequency of time steps to calculate rate of somatic cell production
    fitness_somrate = 2000;
    // select for increasing cell types
    


    // target length with 1 gene or 2 genes on. These are multipliers (area / tlength = true target length)
    tlength1 = 3;
    tlength2 = 2;


    //growing and dividng
    max_cells = 500;
    div_threshold = 100;
    // thresholds which cell has to be GREATER THAN before its target volume shifts to its actual volume. 
    
    // shrink gene is neutral for simulations because it has no effect. Good for comparison to neutral rate of evolution
    gthresh = 2; // tau used by Paulien. Want growth to be by squeezing and not temperature fluctuations. 
    shrink = -16;
    s_shrink = -16;
    shrink_on = true;


    // Number of morphogens, Does changes depending on sim
    n_diffusers = 3;

    // enzymes that can break down the morphogen
    enzymes = false;

    // gene network parameters
    
    n_lockandkey = 10; // number of lock and keys (==), stored in separate vector for ease
    n_locks = 5; // must be half lockandkey. locks = keys
    
    //adding new medium genes to release constraint on keys
    n_mediums = 5;


    // number of transcription factors. The first two transcription factors must be ON for stem cell identity. 
    n_TF = 4; 
    
    n_length_genes = 2;
    n_MF = 2;

    // number of genes. All gene types must sum to this value (except if using morphogenwave, then activators is +1).
    n_genes = n_diffusers + n_lockandkey + n_mediums + n_TF + n_length_genes +n_MF + shrink_on + (enzymes * n_diffusers); 

    n_activators = n_diffusers + n_TF+n_MF; //number of genes that can activate network (<= n_genes)

    
    n_functional = n_lockandkey + n_length_genes + n_mediums;

    gene_vector_size = n_diffusers + n_TF + n_length_genes +n_MF + n_gr + n_in + shrink_on + (enzymes * n_diffusers);

    // location of target length genes in genome. 
    tloc1 = n_diffusers + n_MF + n_TF;
    tloc2 = tloc1+1;

    shrink_loc = tloc2+1;

    e1_loc=shrink_loc+1;
    e2_loc=e1_loc+1;


    // n_gr = 1; // increment target area of cell (growth gene)
    // n_in = 1; // decrement target area of cell 
    // gloc = tloc2 + 1;
    // + n_gr + n_in // add to n_genes
    // stem_gthresh = 2; // going to have a lower growth threshold for stem cells. 
    n_stem = 2; // not in use


    // maximum number of somatic genes on before a cell cannot divide. HAVE CHANGED THIS TO TF GENES
    // max_on = 3;
    
    cycle_size = 8; // number of past states held by a cell for determining cell types. Selecting against cycling cell types. 
    cycle_threshold = 1; // Dom thinks it is okay to increase this from 2 to 3.  Change back if needed.   


    // start concentrations are (in order): diffusers(external), maternal factors, 
    // transcription factors, target length genes  -- the other networks hold lock then keys --

    // dont need this if starting all at 0. However, TF need to start at 1 to have some input into network for MF 0,0, otherwise will always be single cells at beginning.
    // new_g = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // 

    theta = -0.3;
    
    delta_t = 0.25;
    // decay rate should be 1 - delta_t, but I accidentally used this.
    // It has no effect, as the difference between exp(-0.29) and 1-delta_t is very small.
    d_rate = exp(-0.29);
    morphdecay = exp(-0.2);

    morphogen_dose = 40.;


    min_contig = 25;


    /// Dom has deprecated. 
    Jtable = strdup("J.dat");

    // keep this at 2-3
    neighbours = 2;
    // high value ensures cells are never broken apart by copy attempts.
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


    // Do not need to allocate memory from these because not calling from file (memory always allocated regardless, no need for pointer); 
    
    diff_coeff = new double[n_diffusers];
    decay_rate = new double[n_diffusers];
    secr_rate = new double[n_diffusers];


    saturation = 0;
    dt = 1.0;
    dx = 1/((double)sizex);
    pde_its = 1;

    diff_coeff[0] = 8e-7; // Keeping it at this for now. Maybe this could be evolvable. 
    diff_coeff[1] = 8e-7;

    decay_rate[0] = 2e-3;
    decay_rate[1] = 2e-3;
    
    secr_rate[0] = 2.4e-3;
    secr_rate[1] = 2.4e-3;

    reaction_rate = 5e-3; // small rate = 5e-3; // large rate = 1e-2

    if (n_diffusers > 2)
    {
      
      diff_coeff[2] = 8e-7; 
      decay_rate[2] = 2e-3;
      secr_rate[2] = 2.4e-3;
      
      // Morphogens with shorter range 

      // diff_coeff[2] = 8e-7;
      // decay_rate[2] = 5e-3;
      // secr_rate[2] = 5.5e-3;


      // Morphogens with longer range

      // diff_coeff[2] = 4e-6;
      // decay_rate[2] = 1e-3;
      // secr_rate[2] = 1.5e-3;


    }


    // morphogen wave at the end of programmed division. 
    // The morphogen occupies a hidden spot at the back of the cell "genes" vector. This allows it to decay but not be increased.
    // Depracated
    morphogen = false;


    // polarity matrix. Used for transcription factors being at one side of the cell upon reproduction. Depracated 
    start_polarity = { 0, 0, 0, 0 };
    polarity_on = false;

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
    storage_stride = 500;
    // for some reason this isn't working.
    screen_freq = 200;
    // using a new one until I figure out why
    n_screen_freq = 100;


    datadir = strdup("data_film");
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




    // T = fgetpar(fp, "T", 50., true);
    // target_area = igetpar(fp, "target_area", 100, true);
    // target_length = igetpar(fp, "target_length", 60, true);
    // lambda = fgetpar(fp, "lambda", 50, true);
    // lambda2 = fgetpar(fp, "lambda2", 5.0, true);
    // Jtable = sgetpar(fp, "Jtable", "J.dat", true);
    // conn_diss = igetpar(fp, "conn_diss", 2000, true);
    // vecadherinknockout = bgetpar(fp, "vecadherinknockout", false, true);
    // extensiononly = bgetpar(fp, "extensiononly", false, true);
    // chemotaxis = igetpar(fp, "chemotaxis", 1000, true);
    // border_energy = igetpar(fp, "border_energy", 100, true);
    // neighbours = igetpar(fp, "neighbours", 2, true);
    // periodic_boundaries = bgetpar(fp, "periodic_boundaries", false, true);
    // n_chem = igetpar(fp, "n_chem", 1, true);
    // diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, true);
    // decay_rate = dgetparlist(fp, "decay_rate", n_chem, true);
    // secr_rate = dgetparlist(fp, "secr_rate", n_chem, true);
    // saturation = fgetpar(fp, "saturation", 0, true);
    // dt = fgetpar(fp, "dt", 2.0, true);
    // dx = fgetpar(fp, "dx", 2.0e-6, true);
    // pde_its = igetpar(fp, "pde_its", 15, true);
    // n_init_cells = igetpar(fp, "n_init_cells", 100, true);
    // size_init_cells = igetpar(fp, "size_init_cells", 10, true);
    // sizex = igetpar(fp, "sizex", 200, true);
    // sizey = igetpar(fp, "sizey", 200, true);
    // divisions = igetpar(fp, "divisions", 0, true);
    // mcs = igetpar(fp, "mcs", 10000, true);
    // rseed = igetpar(fp, "rseed", -1, true);
    // subfield = fgetpar(fp, "subfield", 1.0, true);
    // relaxation = igetpar(fp, "relaxation", 0, true);
    // storage_stride = igetpar(fp, "storage_stride", 10, true);
    // graphics = bgetpar(fp, "graphics", true, true);
    // store = bgetpar(fp, "store", false, true);
    // datadir = sgetpar(fp, "datadir", "data_film", true);

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
