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
    // show on screen
    graphics = true;
    // show morphogen gradients
    contours = true;
    // draw cell displacement paths
    draw_paths = false;

    // Generate a random genome
    randomise = false;

    // ANALYSIS PARAMS: note that there is slow down when these are turned on. 
    // output data for analysis (connectivity, gene expression, state transitions)
    gene_output = true;
    // gene record needs to always be on to test network connectivity. 
    gene_record = true;
    // include regulatory proteins in the state space 
    max_statespace = true;

    //for umap
    umap = false;

    // record momenta for all cells etc
    velocities = false;
    record_directions = false;

    // record cell sizes
    output_sizes = false;

    // record gamma parameter
    output_gamma = false;

    // read genomes from file
    file_genomes = true;

    //name of data file
    data_file = "org-data";

    //record shape index of cells
    record_shape = false;

    //for storing images
    store = true;

    // Start from specific seed. USE 0 for random seed. (Should be 0 unless need specific seed.)
    pickseed=0;//4626157915171642161;
    rseed = -1;

    // KEEP THIS TO FALSE FOR EVOLUTION
    print_fitness = true; 

    // This start matrix is for sorting, overlap and transitions. For evolution start matrix, see start_n below 
    start_matrix = { { -1, 0, 1.5, 1, 2, -1, 0, 0, 0, 0 }, { -1.5, -0.5, -1.5, -1, -2, -2, 0.5, -0.5, -2, 0 }, { 2, 2, -1, 0, -0.5, 1, -1.5, -0.5, 1, -1 }, { 0.5, -1.5, -0.5, -1.5, -1, -1, -2, -0.5, 1, -2 }, { 2, -1, -1, 1, 0.5, 0, 0, 0, 0, 0 }, { 0, -1, 2, 2, -1.5, 1.5, 1, 2, 2, 1.5 }, { 0.5, 2, -1, -2, -1, 0.5, 0.5, 0.5, 2, 0 }, { 0, 0.5, 0.5, 1, 0, 2, 1, -0.5, 1.5, -1 }, { -0.5, -2, 1, 0.5, 0.5, 0.5, 0, 0, 1, -1 }, { 0, 1.5, -1.5, 0, 1.5, -0.5, 0, 1, 2, 2 }, { -2, -0.5, -0.5, 1, -2, 1, -0.5, 0, 0, 0 }, };



/* Cellular Potts parameters */
    sizex = 150;
    sizey = 250;
    mcs = 18000;
    T = 3;
    target_length = 0;
    lambda = 0.5;
    lambda2 = 0; // WARNING - do not move from 0 (deltaH function has changed)
    div_threshold = 100;
    // thresholds which cell has to be GREATER THAN before its target volume shifts to its actual volume. 
  
    // shrink gene is neutral for simulations because it has no effect. Good for comparison to neutral rate of evolution
    
    periodic_boundaries = false;
    // keep this at 2= moore neighbourhood. 2 used in simulations. 
    neighbours = 2;
    // high value ensures cells are never broken apart by copy attempts.
    conn_diss = 2000;
    shrink_on=false;

/* adhesion params */

/*phase params*/ 
    phase_evolution=true;
    J_stem=3.5;
    J_diff=12;
    J_med=6.25;//0.5*J_diff+0.5;//0.25 + 0.5*J_diff;//0.5+0.5*J_diff;
    J_stem_diff=12;//J_diff + 0.5;//(J_diff - J_stem);
    // J_med=8;
    J_med2=J_med;//0.5*J_diff+0.5;
    Vs_max = 0.5; // 1;
    Vd_max = 0; // 1; 

    n_diffusers=4;
    // morphogen parameters
    secr_rate = new double[n_diffusers];
    diff_coeff = new double[n_diffusers];
    decay_rate = new double[n_diffusers];

start_matrix = { { 0, 0, 1, -0.5, 2, -1, 0.5, -0.5, 0, 0.5 }, { 0, 1.5, 0, -1.5, -1, -1.5, -0.5, 0, 0.5, 0 }, { 0.5, 1, -1, -0.5, 0, -1.5, -0.5, 0.5, 1.5, -2 }, { -2, -0.5, -1, 2, 1.5, -1, 1, 0.5, 1, 0 }, { 2, 2, 0, 0.5, 0, 0, 0, -1, 0, 0.5 }, { 1, -1.5, 1.5, -1.5, -2, 2, 0, 0, -0.5, 1.5 }, { 2, -0.5, 0.5, 0, 0, 0, -1, 1.5, 0, 1.5 }, { 2, 0, 1, 1, -1, -2, -1, 0.5, -2, -1.5 }, { 0.5, -1, -1, 0, 0, 0.5, -1, 2, 0, -2 }, { -2, -1, -1.5, 0, 1, -1.5, 0, 1, 2, -1.5 }, { -2, 0, 2, 0.5, -1.5, 1, -1, -2, -0.5, -1.5 }, };
secr_rate[0] = 0.0045879; decay_rate[0] = 0.003; diff_coeff[0] = 2.61153e-07; secr_rate[1] = 0.0005; decay_rate[1] = 0.003; diff_coeff[1] = 4.1761e-08; secr_rate[2] = 0.0005; decay_rate[2] = 0.003; diff_coeff[2] = 7.03154e-09; secr_rate[3] = 0.000757631; decay_rate[3] = 0.003; diff_coeff[3] = 4.618e-08;

    // GRN params
    n_TF = 4; 
    n_length_genes = 0;
    n_MF = 2;
    gthresh = 0;



    melting_adhesion = false;
    tip_max = 50;
    tip_min = 0;
    melt = -30;
    slope = 4;
    addition_rate=1;
    if (melting_adhesion)
    {
      gthresh = 50;
    v_melt = -30;
    v_slope = -4;
    penalty=250;
    optimization_replicates = 6;
    pics_for_opt = false;
    pics_for_opt_interval = 100;
    max_div_time = 20000;

    
/*iterators */
    shrink = -16;
    // this is a neutrally evolving gene if needed
    shrink_on = false;
    // difference between maximum and minimum cell J
    interval1 = maxJ-minJ;
    // addition of J for each lock and key pair
    interval2 = interval1 / (double)n_lockandkey;

    n_genes = n_diffusers + n_TF + n_MF + !phase_evolution * (n_lockandkey + n_mediums + shrink_on + n_length_genes + (enzymes * n_diffusers)) + phase_evolution*(1);
    // number of genes. All gene types must sum to this value (except if using morphogenwave, then activators is +1).
    // n_genes = n_diffusers + n_lockandkey + n_mediums + n_TF + n_length_genes + n_MF + shrink_on + (enzymes * n_diffusers); 
    n_activators = n_diffusers + n_TF+n_MF; //number of genes that can activate network (<= n_genes)
    n_functional = phase_evolution * (1) + !phase_evolution*(n_lockandkey + n_length_genes + n_mediums);
    gene_vector_size = n_diffusers + n_TF + n_MF +  !phase_evolution * (n_length_genes + n_MF + shrink_on + (enzymes * n_diffusers));
    

    //location of maternal factors in genome
    mfloc1 = n_diffusers;
    mfloc2 = mfloc1 + 1;

    // location of target length genes in genome. 
    tloc1 = n_diffusers + n_MF + n_TF;
    tloc2 = tloc1+1;

    shrink_loc = tloc2+1;
    // DEPRACATED
    e1_loc=shrink_loc+1;
    e2_loc=e1_loc+1;
    e1_loc=shrink_loc+1;
    e2_loc=e1_loc+1;


/* sheet related parameters */
    sheet=false;
    sheet_hex = false;
    sheet_J = 6;
    sheet_minJ=0.5;
    sheet_maxJ=12.5;
    J_width=0.5;

    sheetmix=true;
    sheetmixJ=2*sheet_J;
    sheetcol1=4;
    sheetcol2=6;
    mix_swaprate=0.;


    // sheet anisotropy adhesion;
    lambda3=1;


    highT=true;
    highT_time = 10000;
    highT_temp = 5.;

    start_sheet_measure = highT_time + 100;
    end_sheet_measure = start_sheet_measure + 100;

    // diffusion parameters
    waiting_time = 2000;
    // equilibrate MUST be bigger than highT_time (pref 2000, 1000)
    equilibriate = highT_time + 1000;

    // if temperature or adhesion energies etc. are not integers, we need to set this to false, so that
    // dH is calculated on the fly. 
    IntegerHamiltonian = false;

    insitu_shapes=true;

/* differentiation parameters */

    // edges and nodes only at end of simulation (always true).
    potency_edges = true;
    // what mcs to start measuring adult types & differentiation. set to 6000 for all results
    adult_begins = 1000;

    // prune tiny edges (<1 per org) from graph, but true for separating stem and differentiated or nearly stem where necessary)
    prune_edges = true;
    prune_amount = 5;

    cycle_check = false;

    // flat threshold for nodes
    node_threshold = 0; // int(floor((mcs - adult_begins) / 40) * 2 * 10 * n_orgs);

    // DEPRACATED - prune nodes below this percent. Should probably set this as a minimum value (i.e. 10 cells equivalent)
    // using node threshold above
    node_percent = 0.03;

    cycle_size = 8; // number of past states held by a cell for determining cell types. Selecting against cycling cell types. 
    cycle_threshold = 1; // Dom thinks it is okay to increase this from 2 to 3.  Change back if needed.   



/*Conditions for evolution */
    // select for movement of cells towards one side of the boundary.
    asymmetry_selection = true; 
    asym_only = true;
    swap_selection = 240.; // the average fitness of population needed to switch from asym_only to asymmetry selection. 
    type_selection=false;
    swap_selection2 = 110;
    growth_selection=false;
    elongation_selection = false;
    starter = false;
    n_orgs = 60; // should be multiple of 4, 60 used for evolution
    // start from a certain network
    start_n = { { 0, 2, -1 }, { 2, 0, 0 }, { 0, -2, 2 }, { -1, -1, 1 } };
    if (starter)
    {
      start_n = { { 0, -0.5, 0, 2, -1, -0.5, -0.5, 0, 0 }, { 1, 0.5, -1, 0, 0, -0.5, 1.5, -1.5, 1.5 }, { 0, -0.5, -1.5, 0, 0, 0, 2, 1, -1 }, { 2, 0.5, 0.5, 0, 0, 0, -0.5, 0, 0 }, { 1, 0, 0, -2, 2, 0, -0.5, -0.5, 0 }, { 1.5, -0.5, 1.5, 0, 0, -1, 1, 0, 0 }, { 2, -0.5, -0.5, -0.5, -2, 0, 0, -1, 0 }, { 0, 0.5, 1, 0, 0.5, -1, -1, 0, -1 }, { -0.5, -1, -2, 1, -1, 0.5, 0, 2, 0 }, { -1.5, 0.5, -0.5, -1, 1, -0.5, -1.5, -0.5, -0.5 }, };
      secr_rate[0] = 0.00375233; decay_rate[0] = 3e-3; diff_coeff[0] = 4.23754e-07; secr_rate[1] = 0.00121324; decay_rate[1] = 3e-3; diff_coeff[1] = 9.25668e-09; secr_rate[2] = 0.00119461; decay_rate[2] = 3e-3; diff_coeff[2] = 1.15174e-08;      
    }

    select_switch = false;
    // fluctuating selection interval
    fluctuate_interval1 = 100;
    fluctuate_interval2 = 20;

    evo_pics = true;
    pic_gen_interval = 100; 
    pic_dir = "images";
  
    evs = 10000;
    insert_randoms = true;
    n_mutations = 1;
    // mut rate for gene network
    mut_rate = 0.5;
    // mutate for J stem diff (no longer in use)
    J_mutate_probability=0.1;

    // morphogen mutations
    morph_mean_mut=0.;
    morph_stdev_mut=0.1;


    // mutation rate for polarities
    polm_rate = 0.2;
    n_pred = n_orgs / 2;
    min_contig = 25;

    // when to collect fitness
    fitness_begin = 0.9;
    // frequency of time steps fitness is collected
    fitness_typerate = 100;



/* init conditions and so forth */
    // init params for organisms
    target_area = 4900;
    size_init_cells = 70; // this is equal to the radius(diameter?) of the circle (done by eden growth). 
    n_init_cells = 1;
    divisions = 0;
    offset = 75;//75

    //programmed division parameters
    end_program = 100;
    begin_network = 50;
    div_freq = 10;
    // begin_movement=1200;
    program_its = 1; // we are doing more PDE iterations during the program. 
    div_end = 50;

    if (melting_adhesion)
    {
      end_program=0;
    }

/* GRN */
    update_freq = 40;
    theta = -0.3;
    delta_t = 0.25;
    d_rate = 1;

    // add noise to regulatory network 
    noise = false;
    // noise amount
    noise_dose=0.02;
    noise_start = 3000;



    // show output of all comparisons for overlap. Only use when comparing a small number of organisms. 
    // NOTE - MUST RUN with tag: "-platform offscreen" when using cluster (there is no display).
    overlap_images = false;
    overlap_orgs = 10;
    // true = compare different genomes, false = compare same genomes
    between_org_overlap = false;
    // do translations
    translation = true;
    t_interval = 10;
    nt_intervals = 5;
    n_replicates = 2;

    /* colours */
    set_colours = true;
    use_colour_index = true;
    colour_index = { {0, 10}, {1, 11},};

    //record location of cell divisions
    division_anisotropy = false;



    //count stem and differentiated cells
    stem_counts = false;

    count_bud_cells = false;

    // print single cell proteins
    single_cell = false;
    // the phenotype number to return
    single_type = 350;
    // start all cells from "single state" initial condition for ex vivo
    flush_cells = false;
   // turn all cells into this state at beginning of development
    flush_states = { 4.64542e-09, 1, 0.692659, 1, 1, 5.10909e-12, 6.54489e-28, 3.35757e-22, 1.34902e-36, 0.999999, 1.72289e-16, 1, 6.20555e-07, 1.16941, 0.829774,  };



    // convert cells at certain time point to square with radius as shown (radius is half length of square).
    convert_cells = false;
    convert_time = 3000;
    choose_alive_cell = false;
    convert_x = 45;
    convert_y = 125;
    convert_size = 25;
    clear_radius = 40;
    convert_to_type = 97539;
    convert_states = { 0.00289074, 2.95394e-32, 4.73742e-06, 0.993102, 0.993102, 4.00988e-09, 0.993101, 0.00231804, 0.678508, 0.992729, 0.993101, 0.542056, 0.993102, 0.993102, 0.992443, 0.993102, 6.89652e-12, 0.00205639, 0.00223481, 6.97606e-12, 1.37808e-20, 0.993101, 0.512162, 0.538056, 3.55936e-09, 3.39506e-09, 0.00218893, 0.0127123, 6.26531e-11, 0.00402398,  };
   
    // print list of concentrations, one for each cell type (used for future input)
    print_type_concentrations = true;

    // print initial cell protein concentrations
    output_init_concs = false;

    // used to make a sheet of cells. Use in combination with flush_cells and convert cells
    make_sheet=false;
    sheetx=250;
    sheety=100;
    triangle_x=175;
    triangle_y=75;


    // limit morphogen to amount (prevents differentiation in some ex vivo organisms)
    limit_morph = false;
    limit_amount = 0.8;

    morphogen_dose = 40.;

    hold_morph_constant = false;
    // scramble cells
    scramble = false;

    // record copies
    // recordcopies = true;
    // mintype = 1300;
    // maxtype = 2000;

    // DONT TOUCH!
    vecadherinknockout = false;
    extensiononly = false;
    chemotaxis = 0;
    border_energy = 1;//1000


    // n_gr = 1; // increment target area of cell (growth gene)
    // n_in = 1; // decrement target area of cell 
    // gloc = tloc2 + 1;
    // + n_gr + n_in // add to n_genes

    // maximum number of somatic genes on before a cell cannot divide. HAVE CHANGED THIS TO TF GENES
    // max_on = 3;
    
    // start concentrations are (in order): diffusers(external), maternal factors, 
    // transcription factors, target length genes  -- the other networks hold lock then keys --

    // dont need this if starting all at 0. However, TF need to start at 1 to have some input into network for MF 0,0, otherwise will always be single cells at beginning.
    // new_g = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // 

    /* morphogens */
    // I have defined separate rates for each diffuser, however they are the same for now.
    // This may change in the future if I make one an evolvable parameter. 
    n_chem = 0; // Dom not currently using, instead using n_diffusers


    subfield = 1.0;
    relaxation = 0;

    pde_divisor=2;
    pde_sx=sizex/pde_divisor;
    pde_sy=sizey/pde_divisor;

    saturation = 0;
    dt = 1;
    dx = double(1)/double(125);// 1/((double)sizex);
    pde_its = 1;


    
    // depracated - enzymes that can break down the morphogen
    enzymes = false;
    reaction_rate = 5e-3; 

    // DEPRACATED: let all cells be capable of cell division.
    all_divide = true;
    // polarity matrix. Used for transcription factors being at one side of the cell upon reproduction. Depracated 
    start_polarity = { 0, 0, 0, 0 };
    polarity_on = false;

    // storing images.
    storage_stride = 500;
    // for some reason this isn't working. Hard code in sorting if necessary. 
    screen_freq = 200;


    datadir = strdup("data_film");

      /// Dom has deprecated. 
      Jtable = strdup("J.dat");
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

/*new params*/
    // n_lockandkey = 10; // number of lock and keys (==), stored in separate vector for ease
    // n_locks = n_lockandkey / 2;
    // n_TF = 4; 
    // n_length_genes = 2;
    // minJ=4;
    // maxJ=24;
    // n_mediums=5;
    // med_table = new int[n_mediums];
    // med_table[0] = 8;//8;
    // med_table[1] = 5;//5;
    // med_table[2] = 3;
    // med_table[3] = 1;//1;
    // med_table[4] = 1;
    // n_diffusers=4;
    // n_MF=2;
    // minM=6;
    // gthresh = 2; 

    // tlength1 = 2; // target length with 1 gene or 2 genes on. These are multipliers. 2 is approximately circle.
    // tlength2 = 6;

    // // morphogen parameters
    // secr_rate = new double[n_diffusers];
    // diff_coeff = new double[n_diffusers];
    // decay_rate = new double[n_diffusers];
  
    // secr_rate[0] = 5e-3;
    // decay_rate[0] = 3e-3;
    // diff_coeff[0] = 4e-7; 

    // secr_rate[1] = 5e-3;
    // decay_rate[1] = 3e-3;
    // diff_coeff[1] = 4e-7;

    // secr_rate[2] = 5e-3;
    // decay_rate[2] = 3e-3;
    // diff_coeff[2] = 4e-7;

    // secr_rate[3] = 5e-3;
    // decay_rate[3] = 3e-3;
    // diff_coeff[3] = 4e-7;

/*stem-cell system project params*/
    // phase_evolution=false;
    // n_lockandkey = 10; // number of lock and keys (==), stored in separate vector for ease
    // n_locks = n_lockandkey / 2;
    // n_TF = 4; 
    // n_length_genes = 2;
    // minJ=4;
    // maxJ=24;
    // n_mediums=5;
    // n_diffusers=3;
    // med_table = new int[n_mediums];
    // med_table[0] = 5;//8;
    // med_table[1] = 4;//5;
    // med_table[2] = 3;
    // med_table[3] = 2;//1;
    // med_table[4] = 1;
    // n_MF=2;
    // minM=6;
    // gthresh = 2; 
    // shrink_on=true;
        
    // tlength1 = 2; // target length with 1 gene or 2 genes on. These are multipliers. 2 is approximately circle.
    // tlength2 = 6;

    // secr_rate = new double[n_diffusers];
    // diff_coeff = new double[n_diffusers];
    // decay_rate = new double[n_diffusers];
  
    // secr_rate[0] = 2.4e-3;
    // decay_rate[0] = 2e-3;
    // diff_coeff[0] = 8e-7; 

    // secr_rate[1] = 2.4e-3;
    // decay_rate[1] = 2e-3;
    // diff_coeff[1] = 8e-7; 

    // secr_rate[2] = 2.4e-3;
    // decay_rate[2] = 2e-3;
    // diff_coeff[2] = 8e-7;  