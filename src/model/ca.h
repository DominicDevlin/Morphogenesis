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

// mainpage.h contains no C++ code, it is for the main page of the
// documentation
#include "mainpage.h"

/*! Implementation of the Glazier & Graner cellular Potts model **/
#ifndef _CA_HH_
#define _CA_HH_
#include <vector>
#include <stdio.h>
#include "graph.h"
#include "pde.h"
//#include "dish.h"
#include "cell.h"
#include <array>
#include <map>
#include <set>

class Dish;

class Dir {
  
  /* To store a celldirection matrix */
  friend class CellularPotts;
public:

  Dir() {
    aa1=0.; aa2=0.;
    bb1=0.; bb2=0.;
    lb1=0.; lb2=0.;
  }

  double  aa1,aa2;
  double  bb1,bb2;
  double  lb1,lb2;
};

//! Counts of living cells broken down by differentiation state
//! (see CellularPotts::CountCellTypes).
struct CellTypeCounts {
  int zona_pellucida{};
  int sox2_high{};
  int sox17_high{};
  int loser{};
  int undifferentiated{};

  int total() const {
    return zona_pellucida + sox2_high + sox17_high + loser + undifferentiated;
  }
};

class CellularPotts {

  friend class Info;
  friend class Morphometry;
    
public:
  //! \brief Constructs a CA field. This should be done in "Dish".
  CellularPotts(std::vector<Cell> *cells, const int sizex=200, 
		const int sizey=200 );
  // empty constructor
  // (necessary for derivation)
  CellularPotts(void);


  /*! \brief Stretched induced cell growth and division.
    
  See Hogeweg (2000), Journal of Theoretical Biology.

  Find stretched cells, and increase their target area.
  Find enlarged cells, and divide them.*/
  void CellGrowthAndDivision(int time);

  int CountCells(void) const;

  //! \brief Classifies living cells into zona pellucida, Sox2-high (epiblast-like),
  //! Sox17-high (primitive-endoderm-like), loser (both markers above threshold,
  //! mutually inhibited) and undifferentiated (neither marker above threshold),
  //! and returns the count for each category.
  CellTypeCounts CountCellTypes(void) const;

  void set_num(int in);

  void set_seed(void);


  void set_datafile(string file);


  void update_cell_velocities_MCS();

  void FillGrid();

  void FractureSheet();

  void FractureSheet(int n_cells);

  double SumEnergy();

  void SetMediumArea();

  bool IsLocallyConnected(int* nbs, int check_val);

  bool SpawnCell(int x, int y, int cp_sigma, int time);

  int GetNewPerimeterIfXYWereAdded(int sxyp, int x, int y, const int* neighbor_spins);

  int GetNewPerimeterIfXYWereRemoved(int sxy, int x, int y, const int* neighbor_spins);

  void MeasureCellPerimeters();

  void MeasureSinglePerimeter(int targetsigma);


  void SetPerims(int tperim=0);

  void SetSoxColours(double tfrac);

  void ApoptoseDeadCells();

  void PopulateSparseCells(double density, double R, int shiftx, int shifty);

  void PopulateDenseCells(double density, double R, int shiftx, int shifty);


  /* SYNTHETIC MULTICELLULAR STRUCTURE METHODS */
  void SyntheticNetwork();

  void StartSyntheticNetwork(int start_point=0);

  void StartSyntheticNetwork(Cell &newcell);

  // this is going to be quite costly and iterate over the whole grid so we should be 
  // efficien and try and double things up.
  void SurfaceBindings();

  void OutputSyntheticNetwork(int thetime);

  void SyntheticGrowth(int t=0);

  void UpdateSyntheticCellConstraints();

  void UpdateActiveMotion();

  void MakeSpheroid(int centerx, int centery, int radius);

  void DivideCellsNoGrid(vector<bool> which_cells);

  void ClearGrid();

  void PopulateDenseCellsInZonaRadius(double density, double R, int shiftx, int shifty, double h, double k, double a, double b, double n);


  inline double GetDynamicAdhesion(int i, int j) const
  {
    return DynamicAdhesions[i][j];
    // int row = std::max(i, j);
    // int col = std::min(i, j);
    // return DynamicAdhesions[(row * (row + 1) >> 1) + col];
  }

  inline double DynamicAdhesionDiff(int i, int j) const
  {
    // if (i==j)
    //   return 0;
    // else if (i == 0)
    //   return par.dynJmed;
    // else if (j == 0)
    //   return par.dynJmed;
    if (i != j)
    {
      return par.Jdyndiff;
    }
    else if (i==0)
    {
      return par.AstaticJ;
      // cout << "return 2: " << GetDynamicAdhesion(i, j) << endl;
      // return GetDynamicAdhesion(i, j);
    }
    else
    {
      return par.BstaticJ;
    }
  }

  void SetAreas(int tarea);

  void MakeZonaPellucida(double h, double k, double a, double b, double n);

  void DifferentiateZonaPellucida();

  void SetMotilityStrengths();

  void Vectorfield();


  // personal random numbers for xoshiro RNG (each grid has its own state)
  uint64_t s_val[4]{1,1,1,1};

  // Keyword virtual means, that derived classed (cppvmCellularPotts) can override
  // this function and carry out the memory allocation in their preferred way
  // Every time AllocateSigma is called in the base class methods
  // the function belonging the actual type will be called
  virtual void AllocateSigma(int sx, int sy);
  
  // destructor must also be virtual
  virtual ~CellularPotts();

  /*! \brief Plots the dish to the screen or to a movie and searches the
   neighbours. 

   These distinct tasks have been lumped together in the
   same method because both for drawing the black lines between the
   cells and for searching the neighbours the cell borders have to be
   determined. */
  int **SearchNandPlot(Graphics *g=0, bool get_neighbours=true);


  vector<vector<int>> SearchNforVertices();

  vector<pair<int,int>> SearchNforEdges();
  

  //! Plot the dish to Graphics window g
  inline void Plot(Graphics *g) {
    SearchNandPlot(g, false);
  }
  
  //! Searches the cells' neighbors without plotting
  inline int **SearchNeighbours(void) {
    return SearchNandPlot(0, true);
  }

  //! Return the total area occupied by the cells
  inline int Mass(void) {
    int mass=0;
    for (int i=0;i<sizex*sizey;i++) {
      if (sigma[0][i]>0) mass++;
    }
    return mass;
  }

  inline void Set_J(double J)
  {
    internal_J = J;
  }

  inline void set_mixJ(double mixJ)
  {
    internal_mixJ=mixJ;
  }

  inline void Set_evoJ(double J)
  {
    evo_J = J;
  }


  void InitialiseRandomSoxValues();


  void ToxictoLonelyCells();

  void NeighbourBasedActiveMotion(double tfrac);

  void InnerCellMassDivisions(int t);

  void DrawDivisionTimes();

  void CheckIfDivisionHit(int t);

  void NeighbourBasedApoptosis();




  /*! Plot the cells according to their cell identity, not their type.
    
  The black lines are omitted.
  */
  void PlotSigma(Graphics *g, int mag=2);

  
  //! Divide all cells.
    void DivideCells(void) {
	  std::vector<bool> tmp;
    DivideCells(tmp);
  }
  
    /*! Divide all cells marked "true" in which_cells.
      
    \param which_cells is a vector<bool> with the same number of
    elements as the number of cells. It is a mask indicating which
    cells should be divided; each cell marked true will be divided.
      
     If which_cells is empty, this method divides all cells.
    */
    void DivideCells(std::vector<bool> which_cells, int t=0);
    
    /*! Implements the core CPM algorithm. Carries out one MCS.
      \return Total energy change during MCS.
    */
    int AmoebaeMove(long tsteps, PDE *PDEfield=0);

    int AmoebaeMoveLegacy(long tsteps, PDE *PDEfield=0);
  
    /*! \brief Read initial cell shape from XPM file.
      Reads the initial cell shape from an 
      include xpm picture called "ZYGXPM(ZYGOTE)",
      and it allocates enough cells for it to the Dish */
    // void ReadZygotePicture(void);

    
    // int BlobCounting(void); (not implemented?)
    
    // Used internally to assign a Cell object to each "blob"
    // "blobs" should already have different indices.
    //
    // (i.e. to start from a binary image you'd probably first want
    // to apply a blob counting algorithm
    void ConstructInitCells(Dish &beast);
    
    //! Returns the number of completed Monte Carlo steps.
    inline int Time() const {
      return thetime;
    }
  

    // not currently used? In Critter implementation (see Hogeweg
    // 2000) this was used to have cells divide at double their original area.
  inline int ZygoteArea() const {
    return zygote_area;
  }

  //! \brief Return the horizontal size of the CA plane.
  inline int SizeX() const {
    return sizex;
  }
  
  //! \brief Return the vertical size of the CA plane.
  inline int SizeY() const {
    return sizey;
  }
  
  /*! \brief Return the value of lattice site (x,y).

  i.e. This will return the index of the cell which occupies site (x,y). */
  inline int Sigma(const int x, const int y) const {
    return sigma[x][y];
  }
  
  // Was used to make it possible to enlarge the Graphics window in
  // X11 and replace the contents interactively. Not currently supported.
  void Replace(Graphics *g);

  /*! In this method the principal axes of the cells are computed using
   the method described in "Biometry", box 15.5 
   \return a pointer to a "new[]"ed array containing the directions.
   The memory has to be freed afterwards using the delete[] operator
  */
  Dir *FindCellDirections(void) const;

  /*! \brief Initialize the CA plane with n circular cells fitting in
    a cellsize^2 square.*/
  int ThrowInCells(int n, int cellsize);

  /*! \brief Initialize the CA plane with n cells using an Eden growth algorithm.

  \param n: Number of cells.
  \param cellsize: Number of Eden growth iterations.
  \param subfield: Defines a centered frame of size (size/subfield)^2 in which all cell will be positioned. 
  \return Index of last cell inserted.
  */
  int GrowInCells(int n_cells, int cellsize, double subfield=1.);
  int GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y);
  
  //! \brief Adds a new Cell and returns a reference to it.
  inline Cell &AddCell(Dish &beast) {
    cell->push_back(Cell(beast));
    return cell->back();
  }
  /*! \brief Display the division planes returned by FindCellDirections.
    
  \param g: Graphics window
  \param celldir: cell axes as returned by FindCellDirections.
  */
  void ShowDirections(Graphics &g, const Dir *celldir) const;
  
  //! \brief Returns the mean area of the cells. 
  double MeanCellArea(void) const;

  double MeanCellPerimeter(void) const;
  
  /*! \brief Returns the cell density.

  Cell density is defined as the area occupied by cells divided by the size of the field.
  */
  double CellDensity(void) const; 
  
  //! \brief Set target lengths of all cells to the value given in parameter file.
  void ResetTargetLengths(void);
 
  int spins_converted;
  
  /*! \brief Give each cell a random cell type.
    
  The number of cell types is defined by the J parameter file. (See
  Jtable in parameter file).
  */
  void SetRandomTypes(void);

  void SetColours(void);

  /*! Cells grow until twice their original target_length, then
    divide, with rate "growth_rate"
  */
  void GrowAndDivideCells(int growth_rate);
  
  inline Cell &getCell(int c) {
    return (*cell)[c];
  }

  /*! Draw convex hull around all cells.
    \return The area of the convex hull in lattice sites.
  */
  double DrawConvexHull(Graphics *g, int color=1);
  
  /*! Calculate compactness (summed_area/hull_area) of all cells.
    This is a good measure for the density.
    \return Compactness.
  */
  double Compactness(double *res_compactness = 0, 
		     double *res_area = 0, 
		     double *res_cell_area = 0);
  
  void CopyProb(double T);

private:
  void IndexShuffle(void);
  // double DeltaH(int x,int y, int sxyp, int tsteps, PDE *PDEfield=0);
  double DeltaH(int x,int y, int sxyp, int tsteps, const int* neighbor_spins, PDE *PDEfield=0);
  bool Probability(int DH);
  void ConvertSpin(int x,int y,int kp);
  void ConvertSpinPerim(int x, int y, int kp, const int* neighbor_spins);
  void GetNeighborsSafe(int x, int y, int* nbs); 
  void SprayMedium(void);
  int CopyvProb(double DH,  double stiff);
  void FreezeAmoebae(void);
  void MeasureCellSizes(void);
  void MeasureCellSize(Cell &c);
  
  bool ConnectivityPreservedP(int x, int y);
  




  // little debugging function to print the site and its neighbourhood
  inline void PrintSite(int x,int y) {
	  std::cerr << "--------\n";
	  std::cerr << "[" << sigma[x-1][y-1] << " " << sigma[x][y-1] << " " << sigma[x+1][y-1] << "]\n";
	  std::cerr << "[" << sigma[x-1][y] << " " << sigma[x][y] << " " << sigma[x+1][y] << "]\n";
	  std::cerr << "[" << sigma[x-1][y+1] << " " << sigma[x][y+1] << " " << sigma[x+1][y+1] << "]\n";
  }

protected:
	void BaseInitialisation(std::vector<Cell> *cell);
  
protected:
  int **sigma;
  int **inside_elipse;
  int sizex;
  int sizey;


  int **old_nbhs;
  int old_cell_count{};

  vector<int> old_med_nbhs;
  bool exchange_encounter=false;

  vector<double> prev_x_med;
  vector<double> prev_y_med;
  vector<int> prev_sigmas;

  int medp_success=0;
  int medp_count=0;

  std::map<int, std::set< std::pair<int, int>>> cellVolumeList;
  map<int, int> vlist;
  std::map<int, std::set< std::pair<int, int>>> cellPerimeterList;


  

private:
  bool frozen;
  static const int nx[37], ny[37];
  static const int nbh_level[8];
  static int shuffleindex[9];
  std::vector<Cell> *cell;
  int zygote_area;
  int thetime;
  int n_nb;
  int n_nb_adh;
  int n_nb_perim;

  int stack[8]; // stack to count number of different surrounding cells, CHANGE TO MEMBER FUNCTION

  int init_colours=-20;


  int zona_sigma;
  int zona_sigma_sticky;

  double internal_T;

  double internal_J;

  double internal_mixJ;

  double evo_J;

  map<int,vector<double>> state_shape_index;

  map<int,vector<double>> state_hexatic_order;

  map<int, vector<pair<int,double>>> time_hexatic_order;

  map<int, vector<pair<int,double>>> time_shape_index;

  vector<pair<int,double>> sheet_hexatic_order;

  vector<pair<int,double>> sheet_shape_order;

  map<int,vector<double>> state_adhesion;

  double tmp_hex_order;
  double transition_point;
  double tmp_avg_shape=0;

  double proportion_over_transition;
  int proportion_counter;

  int start_width;
  int opt_starty;

  double hexatic_tally{};
  double hexatic_counter{};
  vector<double> hex_vec{};

  double shape_tally{};
  double shape_counter{};
  vector<double> shape_vec{};

  int rows;
  int cols;

  double n_grads[5] = {tan(-(M_PI)/12.), tan(-(2*M_PI)/12.), tan(-(3*M_PI)/12.), tan(-(4*M_PI)/12.), tan(-(5*M_PI)/12.)};
  double p_grads[5] = {tan((M_PI)/12.), tan((2*M_PI)/12.), tan((3*M_PI)/12.), tan((4*M_PI)/12.), tan((5*M_PI)/12.)};

  // used for calculating the shane in shape over time in the lattice. 
  // int ***Shape;
  // int scount=0;
  bool ShapeMaintained=true;

  //for calculating organism fitness.
  vector<int> som_cell_list;
  vector<double> type_fitness_list;
  vector<double> shape_fitness_list;

  map<int, bool> touching_medium;
  // long flip_true{};
  // long flip_false{};
  // double dH_tally{};
  // double dH_neg{};

  double leftover_mass_stem{};
  double leftover_mass_diff{};

  int** prev_nbs;
  vector<vector<double>> DynamicAdhesions;
  vector<int> DynamicMeeting;


  vector<map<int,int>> TypeCounts;


  // list of non-cycling cell types
  map<int,int> type_list;

  vector<vector<int>> matrix;
  vector<bool> polarity;
  int org_num=1;
  map<int,int> transition_cooldown_list;


  //count grid hits for awkward gradients
  long griditcount;

  int n_segments;

  long rad_count;

  bool early_contig;

  double previous_angle=0;
  double oldx=0;
  double oldy=0;
  double oldxcen=0;
  double oldycen=0;

  string data_file;


};


#endif
