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

    // record momenta for all cells etc
    velocities = false;
    record_directions = false;

    // read genomes from file
    file_genomes = true;

    //name of data file
    data_file = "org-data";

    //for storing images
    store = true;

    // Start from specific seed. USE 0 for random seed. (Should be 0 unless need specific seed.)
    pickseed=0;//16045985250248971749;
    rseed = -1;

    // KEEP THIS TO FALSE FOR EVOLUTION
    print_fitness = true; 

/* Cellular Potts parameters */
    sizex = 300;// was using 300 x 200 for wetting, 200x300 for elongation. Testing 512x200 with dewet length of 36
    sizey = 300;
    mcs = 30001;
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

    bulk_modulus = 5;
    cell_target_area = 400;
    lambda = bulk_modulus / cell_target_area;// 130;
    div_threshold = 150;
    synthetic_max_area=cell_target_area+2;
    synthetic_min_area=cell_target_area-2;

    H_perim = true;
    elastic_modulus = 1;
    ptarget_perimeter = 114;
    perim_offset = 16;
    ptarget_perimeter = ptarget_perimeter * (neigh_multipliers[perimeter_neighbourhood-1]);
    perim_offset = perim_offset * (neigh_multipliers[perimeter_neighbourhood-1]);
    // Note - value must be divided by P_0 to maintain constant force if P_0 is to change.
    lambda_perimeter = elastic_modulus / ptarget_perimeter;// 8;

    // active motion should depend on E/P/N cadherin
    // P/N cadherin binding shoudlnt change active motion. 
    // E cadherin should decrease with E cadherin binding
    active_motion = true;
    motility_strength = 0.25; // not that this term depends on the cell size (1/sqrt(area))
    motility_zero = motility_strength * sqrt(cell_target_area);
    persistence_time = 40.;

    initialise_sox_time=800;
    time_till_full_expression=12000;

    // smaller this is the smoother the curve between losers and winners (this is important)
    switch_like=2000.;

    set_loser_colours=false;

    apop_signal_noise=1.; // was 1.5
    apop_noise_tau=0.1;  // was 0.1
    apop_dt=0.05;
    death_decay_rate=1.; // was 0.15

    // the important params
    loser_sox2_adhesion=0.7; //-0.1;
    loser_loser_adhesion=0.7;// -0.7;
    loser_sox17_adhesion=0.6;//-0.1;
    apop_threshold=20;

    // high value ensures cells are never broken apart by copy attempts.
    // This value is only used in the slightly faster CPM implementation where 
    // detailed balance is not ensured.
    conn_diss = 2000;


    starting_fraction_losers=0.33;
    target_sox2_prob=0.6;
    loser_perim_increase=0.;
/* adhesion params */

    // baseline J value for adhesion between cells and blastocoel
    Jblasto=0.5; //0.5
    // modulation of sox17 expressing cell to medium
    sox17_blasto_adhesion=0.;
    // modulation of sox2 expressing cell to medium
    sox2_blasto_adhesion=-0.;
    loser_blasto_adhesion=-0.;

    // baseline J value between cells
    J_cell_baseline=1.2; //1.2
    // binding of sox2 to sox2
    sox2binding=0.7;
    // binding of sox17 to sox17 =
    sox17binding=0.5;
    // binding between sox2 and sox17
    sox2vs17binding=0.6;


    // J cell zona is the same for all zona. Sticky part has different form non sticky just for specific adhesions.
    J_cell_zona = 1.2; //1.2
    Jzona_sox2 = 0.0;
    Jzona_sox17 = 0.2;
    Jzona_loser=0;
    // added zona adhesion for sox2 sox17 for sticky part
    J_cell_zona_sticky=2.0; //2.0
    Jzona_sticky_sox2extra=1.4;
    Jzona_sticky_sox17extra=0.;

    init_blasto=1.5;
    init_zona=3.0;
    init_zona_sticky=1.0;
    init_cellcell=1.0;

    // end of adhesion params
    adhesion_multiplier=1.5;


    J_cell_baseline=J_cell_baseline*adhesion_multiplier;
    sox2binding=sox2binding*adhesion_multiplier;
    sox17binding=sox17binding*adhesion_multiplier;
    sox2vs17binding=sox2vs17binding*adhesion_multiplier;
    loser_loser_adhesion=loser_loser_adhesion*adhesion_multiplier;
    loser_sox2_adhesion=loser_sox2_adhesion*adhesion_multiplier;
    loser_sox17_adhesion=loser_sox17_adhesion*adhesion_multiplier;
    J_cell_zona=J_cell_zona*adhesion_multiplier;
    Jzona_sox2=Jzona_sox2*adhesion_multiplier;
    Jzona_sox17=Jzona_sox17*adhesion_multiplier;
    Jzona_loser=Jzona_loser*adhesion_multiplier;
    J_cell_zona_sticky=J_cell_zona_sticky*adhesion_multiplier;
    Jzona_sticky_sox2extra=Jzona_sticky_sox2extra*adhesion_multiplier;
    Jzona_sticky_sox17extra=Jzona_sticky_sox17extra*adhesion_multiplier;
    Jblasto=Jblasto*adhesion_multiplier;
    sox17_blasto_adhesion=sox17_blasto_adhesion*adhesion_multiplier;
    sox2_blasto_adhesion=sox2_blasto_adhesion*adhesion_multiplier;
    loser_blasto_adhesion=loser_blasto_adhesion*adhesion_multiplier;



    sox_threshold=0.2;

    // make an oval zona pellucida that does not move
    make_zona_pellucida = true;


    E_cadherin_production_rate=0.04;
    E_cadherin_saturation_constant=0.25;
    hill_coefficient=3.0;
    decay_E_cadherin_unbound=0.01;
    decay_E_cadherin_bound=0.005;
    c_max = 1.5;

    JEcadherin_scaling=1.5;//5.2;
    
    // IMPORTANT NOTE!!
    // CORTICAL TENSION SHOULD GO DOWN FOR ALL CELLS!!! AS THEY BIND MORE TO OTHER CELLS
    // = CELLS AT PERIPHERY WILL BE CIRCULAR, CELLS INSIDE WILL BE FLOPPY
    // NOTE - EFFECT WILL BE ENHANCED WITH CADHERINS BUT NOT LIMITED To
    // Note - i changed this to simply change lambdaP
    make_sparse_cells=true;
    start_radius=80;
    start_density=1.0;
    
    


/* init conditions and so forth */
    // init params for organisms
    init_area = 10240;
    size_init_cells = 80; // this is equal to the radius(diameter?) of the circle (done by eden growth). 
    eden_growth=false;
    n_init_cells = 4;
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
    {
      free(Jtable);
      Jtable = nullptr;
    }
    if (datadir)
    {
      free(datadir);
      datadir = nullptr;
    }
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
