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


  #include "../model/parameter.h"
  #include <cstdio>
  #include <cstring>
  #include <cstdlib>
  #include <cerrno>
  #include <iostream>
  #include "../model/output.h"
  #include "../model/parse.h"
  #include <cmath>


  Parameter::Parameter()
  {
    // show on screen
    graphics = true;
    // show morphogen gradients
    contours = false;
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
    max_statespace = false;

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
    pickseed=0;//15120834811895428147//4626157915171642161;//4766666018663198866used seed for tagaki
    rseed = -1;

    // KEEP THIS TO FALSE FOR EVOLUTION
    print_fitness = true; 


/* Cellular Potts parameters */
    sizex = 300;// was using 300 x 200 for wetting, 200x300 for elongation. Testing 512x200 with dewet length of 36
    sizey = 300;
    mcs = 2000001;
    // NOTE - TEMPERATURE CURRENTLY DEFUNCT SINCE IT IS SET TO 1!
    T = 1;
    // NOTE: lambda must be divided by A_0 to maintain constant force

    periodic_boundaries = false;
    // copy neighbourhood 2 used in old simulations.
    // NOTE - FOR DETAILED BALANCE WE NEED COPY NEIGHBOURHOOD = 1 (see Durand 2016)
    // NOTE - ADHESION AND PERIM NEIGHBOURHOOD MUST BE EQUAL (unless one energy is non-existent)
    adhesion_neighbourhood=4;
    perimeter_neighbourhood=adhesion_neighbourhood;
    copy_neighbourhood=1;
    neigh_multipliers={1, 3, 5, 11, 15, 18, 26};
    neigh_multiplier=double(neigh_multipliers[adhesion_neighbourhood-1]);

    bulk_modulus = 13;
    cell_target_area = 100;
    lambda = bulk_modulus / cell_target_area;// 130;
    div_threshold = 150;
    synthetic_max_area=150;
    synthetic_min_area=75;

    H_perim = true;
    elastic_modulus = 1;
    ptarget_perimeter = 42;
    ptarget_perimeter = ptarget_perimeter * (neigh_multipliers[perimeter_neighbourhood-1]);
    // Note - value must be divided by P_0 to maintain constant force if P_0 is to change.
    lambda_perimeter = elastic_modulus / ptarget_perimeter;// 8;
    Ecad_elastic_change=5;
      

    // active motion should depend on E/P/N cadherin
    // P/N cadherin binding shoudlnt change active motion. 
    // E cadherin should decrease with E cadherin binding
    active_motion = true;
    motility_strength = 0.25;
    Ecadherin_bound_motility_loss=0.15;
    Ncadherin_bound_motility_loss=0.15;
    Pcadherin_bound_motility_loss=0.15;
    persistence_time = 40.;

    

    // note - currently active motion must be on for gravity.
    add_gravity=true;
    lambda_gravity=0.4/sizex;

    // high value ensures cells are never broken apart by copy attempts.
    // This value is only used in the slightly faster CPM implementation where 
    // detailed balance is not ensured.
    conn_diss = 2000;


/* synthetic params */
    make_synthetic = true;
    synthetic_update_step=100;
    check_cell_bindings_step=25;

    production_rate_synNotch=0.02;
    decay_synNotch_bound=0.02;

    binding_rate_CD19_synNotch = 0.5;
    decay_synNotch_unbound=0.02;
    decay_synNotch_intra=0.04;

    GFP_production_rate=0.01;
    decay_GFP=0.002;
    lo_cadherin_production_rate=0.01;


    E_cadherin_production_rate=0.04;
    E_cadherin_saturation_constant=0.25;
    hill_coefficient=3.0;
    decay_E_cadherin_unbound=0.01;
    decay_E_cadherin_bound=0.005;
    c_max = 1.5;

    // this concentration is too high, needs to come down i think (and get scaling right)
    random_binding_protein_production=0.003;
    decay_random_binding_protein_bound=0.001;
    decay_random_binding_protein_unbound=0.01;

    synthetic_dt=0.3;

    synthetic_Jm=1.3;

    // not using atm
    Jmed_scaling=0;

    synthetic_Jcell_baseline = 1;
    JEcadherin_scaling=1.5;//5.2;
    JPcadherin_scaling=1.5;//2.6;
    JNcadherin_scaling=1.5;//2.6;

    Jrandom_scaling_E=0.5;
    Jrandom_scaling_N=0.5;
    Jrandom_scaling_P=0.5;
    
    init_random_binding=1.;
    // IMPORTANT NOTE!!
    // CORTICAL TENSION SHOULD GO DOWN FOR ALL CELLS!!! AS THEY BIND MORE TO OTHER CELLS
    // = CELLS AT PERIPHERY WILL BE CIRCULAR, CELLS INSIDE WILL BE FLOPPY
    // NOTE - EFFECT WILL BE ENHANCED WITH CADHERINS BUT NOT LIMITED To
    // Note - i changed this to simply change lambdaP

    proportion_celltype2 =0.43; // i used 0.7 for 3 layer and 0.47? for asymmetric
    start_radius=150;
    start_density=0.2;
    

    // Here we decide the genes of c1 and c2.
    // first = E_cadherin high, second = E_cadherin low, third = P_cadherin, fourth = N_cadherin, 5 = CD19, 6=GFP, 7=mCherry
    spheroid_const={0,0,1,0,0,0,0};
    spheroid_GFP_induced={0,0,0,0,0,0,0};
    spheroid_mCherry_induced={0,0,0,0,0,0,0};
    spheroid_CD19_induced={0,0,0,0,0,0,0};

    make_spheroid=false;
    make_sparse_cells=true;

    // three layered structure
    // c1_const={0,0,0,0,0,0,0};
    // c2_const={0,0,0,0,1,0,0};
    // c1_GFP_induced={0,0,0,0,0,0,0};
    // c2_GFP_induced={0,1,0,0,0,0,1};
    // c1_mCherry_induced={0,0,0,0,0,0,0};
    // c2_mCherry_induced={0,0,0,0,0,0,0};
    // c1_CD19_induced={1,0,0,0,0,1,0};
    // c2_CD19_induced={0,0,0,0,0,0,0};


    // two layered structure CD19 + Ecad
    // c1_const={0,0,0,0,0,0,0};
    // c2_const={0,0,0,0,1,0,0};
    // c1_GFP_induced={0,0,0,0,0,0,0};
    // c2_GFP_induced={0,0,0,0,0,0,0};
    // c1_mCherry_induced={0,0,0,0,0,0,0};
    // c2_mCherry_induced={0,0,0,0,0,0,0};
    // c1_CD19_induced={1,0,0,0,0,1,0};
    // c2_CD19_induced={0,0,0,0,0,0,0};

    // asymmetric
    c1_const={0,0,0,0,0,0,0};
    c2_const={0,0,0,0,1,0,0};
    c1_GFP_induced={0,0,0,0,0,0,0};
    c2_GFP_induced={0,0,1,0,0,0,1};
    c1_mCherry_induced={0,0,0,0,0,0,0};
    c2_mCherry_induced={0,0,0,0,0,0,0};
    c1_CD19_induced={0,0,0,1,0,1,0};
    c2_CD19_induced={0,0,0,0,0,0,0};

    // for spheroid stuff
    // c1_const={0,0,0,0,0,0,0};
    // c2_const={0,0,0,0,0,0,0};
    // c1_GFP_induced={0,0,0,0,0,0,1};
    // c2_GFP_induced={0,0,0,0,0,0,1};
    // c1_mCherry_induced={0,0,0,0,0,0,0};
    // c2_mCherry_induced={0,0,0,0,0,0,0};
    // c1_CD19_induced={0,0,0,0,0,0,0};
    // c2_CD19_induced={0,0,0,0,0,0,0};

    // we have GFP, mcherry and cd19 (in that order). This vector decides whether
    // they are morphogens or not. 1=morph, 2=surface.
    if (make_spheroid)
      morph_or_surface={1,0,0};
    else
      morph_or_surface={0,0,0};



    // morphogen stuff.
    n_diffusers = 3; // morphogens (cant be less than one)
    secr_rate = new double[n_diffusers];
    diff_coeff = new double[n_diffusers];
    diff_coeff_cell = new double [n_diffusers];
    decay_rate = new double[n_diffusers];
    decay_rate_cell = new double[n_diffusers];
    subfield = 1.0;
    relaxation = 0;
    saturation = 0;
    dt = 1.0;
    dx = double(1)/double(250);// 1/((double)sizex);
    pde_its = 4;
    
    // GFP
    diff_coeff[0] = 2e-6;
    diff_coeff_cell[0]=1e-7;
    decay_rate[0] = 0.03e-3;
    decay_rate_cell[0]=0.09e-3;
    secr_rate[0] = 1e-3;

    //mCHERRY
    diff_coeff[1] = 2e-6;
    diff_coeff_cell[1]=1e-6;
    decay_rate[1] = 0.03e-3;
    decay_rate_cell[1]=0.4e-3;
    secr_rate[1] = 1e-3;

    //CD19 (probably not needed)
    diff_coeff[2] = 4e-6;
    diff_coeff_cell[2]=2e-6;
    decay_rate[2] = 0.03e-3;
    decay_rate_cell[2]=0.4e-3;
    secr_rate[2] = 1e-3;

    // for debugging
    thetime=0;



/*Conditions for evolution */
    // select for movement of cells towards one side of the boundary.
    asymmetry_selection = true; 
    asym_only = false;
    swap_selection = 240.; // the average fitness of population needed to switch from asym_only to asymmetry selection. 
    // start from a certain network
    growth_selection=false;
    elongation_selection = false;
    starter = false;
    n_orgs = 60; // should be multiple of 4, 60 used for evolution

    start_n = { { 0, 2, -1 }, { 1, 0, 0 }, { 0, -2, 2 }, { -1, -1, 1 } };
    evo_pics = false;
    pic_gen_interval = 1;
    pic_dir = "images";
  
    evs = 10000;
    insert_randoms = false;
    n_mutations = 1;
    // mut rate for gene network
    mut_rate = 0.25;
    // mutate for J stem diff (no longer in use)
    J_mutate_probability=0.1;
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
    init_area = 10240;
    size_init_cells = 80; // this is equal to the radius(diameter?) of the circle (done by eden growth). 
    eden_growth=false;
    n_init_cells = 1;
    divisions = 0;

    //programmed division parameters
    end_program = 7;
    begin_network = 1000;
    div_freq = 1;
    // begin_movement=1200;
    program_its = 1; // we are doing more PDE iterations during the program. 
    div_end = 6;

    /* colours */
    set_colours = true;
    use_colour_index = false;
    colour_index = { {0, 10}, {1, 11},};

    //record location of cell divisions
    division_anisotropy = false;


    // DONT TOUCH!
    vecadherinknockout = false;
    extensiononly = false;
    chemotaxis = 0;
    border_energy = 1;//1000

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

  void Parameter::CleanUp(void) 
  {
    if (Jtable) 
      free(Jtable);
    if (diff_coeff) delete[] diff_coeff;
    if (decay_rate) delete[] decay_rate;
    if (secr_rate) delete[] secr_rate;
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
    // init_area = igetpar(fp, "init_area", 100, true);
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
    os << " init_area = " << init_area << endl;
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
    os << " neighbours = " << copy_neighbourhood << endl;
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
