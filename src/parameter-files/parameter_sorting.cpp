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
    sizex = 250;// was using 300 x 200 for wetting, 200x300 for elongation. Testing 512x200 with dewet length of 36
    sizey = 250;
    mcs = 2000001;
    // NOTE - TEMPERATURE CURRENTLY DEFUNCT SINCE IT IS SET TO 1!
    T = 1;
    // NOTE: lambda must be divided by A_0 to maintain constant force
    // copy neighbourhood 2 used in old simulations.
    // NOTE - FOR DETAILED BALANCE WE NEED COPY NEIGHBOURHOOD = 1 (see Durand 2016)
    // NOTE - ADHESION AND PERIM NEIGHBOURHOOD MUST BE EQUAL (unless one energy is non-existent)
    adhesion_neighbourhood=5;
    perimeter_neighbourhood=adhesion_neighbourhood;
    copy_neighbourhood=1;
    neigh_multipliers={1, 3, 5, 11, 15, 18, 26};
    neigh_multiplier=double(neigh_multipliers[adhesion_neighbourhood-1]);

    bulk_modulus = 13;
    cell_target_area = 100;
    lambda = bulk_modulus / cell_target_area;// 130;
    div_threshold = 100;

    H_perim = true;
    elastic_modulus = 2;
    ptarget_perimeter = 40;
    ptarget_perimeter = ptarget_perimeter * (neigh_multipliers[perimeter_neighbourhood-1]);
    // Note - value must be divided by P_0 to maintain constant force if P_0 is to change.
    lambda_perimeter =( elastic_modulus) / (ptarget_perimeter );// 8;
      

    // high value ensures cells are never broken apart by copy attempts.
    // This value is only used in the slightly faster CPM implementation where 
    // detailed balance is not ensured.
    conn_diss = 2000;
    

    // sorting parameters
    startingAproportion=0.5; // A=0(false), B=1(true)
    init_J=-0.;
    dynJmed=0.;
    Jdyndiff=-0.1 / neigh_multiplier;
    AstaticJ=-1.5 / neigh_multiplier;
    BstaticJ=-1.5 / neigh_multiplier;
    
    // timescaler=0.000001;
    // //note this needs to be half (there are two meetings each recording)
    // timeadd_ifmet=5;


    make_sparse_cells=false;
    do_voronoi=true;
    periodic_boundaries = true;

    active_motion = false;
    motility_strength = 0.4;
    persistence_time = 500.;



    // for debugging
    thetime=0;


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

  /*leftover shit do not touch*/
    n_lockandkey = 4; // Locks+keys. number of lock = keys, stored in separate vectors. 
    n_locks = n_lockandkey / 2; // must be half lockandkey. 
    n_mediums = 2;
    med_table = new int[n_mediums]; // J values for cell with medium
    // morphogen stuff.
    n_diffusers = 3; // morphogens (cant be less than one)
    secr_rate = new double[n_diffusers];
    diff_coeff = new double[n_diffusers];
    diff_coeff_cell = new double [n_diffusers];
    decay_rate = new double[n_diffusers];
    decay_rate_cell = new double[n_diffusers];

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
    if (med_table) delete[] med_table;
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
