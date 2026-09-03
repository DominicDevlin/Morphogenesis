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
#ifndef _CELL_HH_
#define _CELL_HH_

#include "parameter.h"
//#define EMPTY -1
#include <math.h>
#include <unordered_map>
#include <utility>
#include <map>
#include <deque>
#include <numeric>
#include <algorithm>


extern Parameter par;
class Dish;

class Cell
{ 

  friend class Dish;
  friend class CellularPotts;
  friend class Info;

public:
  /*! \brief Constructor to insert a cell into Dish "who"
    
  Used to add a new Cell to the dish: new Cell(dish,
  celtype).
  */
  Cell(Dish &who, int
  settau=1) {
    owner=&who;
    ConstructorBody(settau);
  }

  Cell(void) {
    chem = new double[par.n_chem];
    diffs = new double[par.n_diffusers];
  };
  
  ~Cell(void);
  
  
  //! Default copy constructor.
  Cell(const Cell &src) {
        
    // make an exact copy (for internal use)
    sigma=src.sigma;
    area=src.area;
    target_area=src.target_area;
    length=src.length;
    target_length=src.target_length;
    mother=src.mother;
    daughter=src.daughter;
    times_divided=src.times_divided;
    date_of_birth=src.date_of_birth;
    colour_of_birth=src.colour_of_birth;
    tau=src.tau;
    alive=src.alive;
    v[0]=src.v[0];
    v[1]=src.v[1];
    n_copies=src.n_copies;
    sum_x=src.sum_x;
    sum_y=src.sum_y;
    sum_xx=src.sum_xx;
    sum_yy=src.sum_yy;
    sum_xy=src.sum_xy;
    owner=src.owner;

    velocity_histories_x = src.velocity_histories_x;
    velocity_histories_y = src.velocity_histories_y;
    prev_com_x = src.prev_com_x;
    prev_com_y = src.prev_com_y;
    avg_vx= src.avg_vx;
    avg_vy= src.avg_vy;

    lambda_2 = src.lambda_2;
    lambda = src.lambda;
    c_type=src.c_type;

    xcen = src.xcen;
    ycen = src.ycen;

    xcens = src.xcens;
    ycens = src.ycens;

    time_created = src.time_created;


    perimeter = src.perimeter;
    target_perimeter = src.target_perimeter;

    apop_noise_state=src.apop_noise_state;

    diffs = new double[par.n_diffusers];

    for (int i=0;i<par.n_diffusers;i++)
    {
      diffs[i]=src.diffs[i];
    }

    
    chem = new double[par.n_chem];
    for (int ch=0;ch<par.n_chem;ch++)
      chem[ch]=src.chem[ch];


    Sox2_concentration = src.Sox2_concentration;
    sox2_internal_adhesion = src.sox2_internal_adhesion;
    Sox17_concentration = src.Sox17_concentration;
    sox17_internal_adhesion = src.sox17_internal_adhesion;
    lonely_cell = src.lonely_cell;

    div_times = src.div_times;

    cell_perim_constraint = src.cell_perim_constraint;

    epithelial = src.epithelial;
    synNotch_bound = src.synNotch_bound;
    synNotch_unbound = src.synNotch_unbound;
    synNotch_intra = src.synNotch_intra;
    E_cadherin = src.E_cadherin;
    CD19 = src.CD19;
    random_binding_proteins = src.random_binding_proteins;
    touching_med = src.touching_med;
    mCherry=src.mCherry;
    GFP=src.GFP;
    opposing_GFP=src.opposing_GFP;
    P_cadherin=src.P_cadherin;
    N_cadherin=src.N_cadherin;
    opposing_P_cadherin=src.opposing_P_cadherin;
    opposing_N_cadherin=src.opposing_N_cadherin;
    spheroid_cell=src.spheroid_cell;
    leftover_area=src.leftover_area;

    opposing_CD19 = src.opposing_CD19;
    opposing_E_cadherin = src.opposing_E_cadherin;
    opposing_mCherry=src.opposing_mCherry;
    opposing_GFP=src.opposing_GFP;
    opposing_N_cadherin=src.opposing_N_cadherin;
    opposing_P_cadherin=src.opposing_P_cadherin;

    f_opposing_GFP = src.f_opposing_GFP;
    f_opposing_CD19 = src.f_opposing_CD19;
    f_opposing_mCherry=src.f_opposing_mCherry;
    f_opposing_E_cad = src.f_opposing_E_cad;
    f_opposing_N_cad = src.f_opposing_N_cad;
    f_opposing_P_cad = src.f_opposing_P_cad;

    constitutives=src.constitutives;
    GFP_induced=src.GFP_induced;
    mCherry_induced=src.mCherry_induced;
    CD19_induced=src.CD19_induced;

    centerx = src.centerx;
    centery = src.centery;

    cell_elastic_mod=src.cell_elastic_mod;
    motility_strength=src.motility_strength;

    accumulated_death_signals = src.accumulated_death_signals;

    medium_touch_count=src.medium_touch_count;
    
  }
  
  /*! \brief Add a new cell to the dish.
   
     Call it as: new Cell(parent, true); mother will be modified for
     ancestry administration!  

     \param settau 
     Cell type of daughter cell.
  */
  void CellBirth(Cell &mother);
 
  /*! \brief Assignment operator.

  Called if one cell is assigned to another. Remember to change both
  assignment operator and copy constructor when adding new attributes
  to Cell.
  */
  inline Cell &operator=(const Cell &src) {
    colour=src.colour;
    alive=src.alive;
    sigma=src.sigma;
    area=src.area;
    tau=src.tau;
    target_area=src.target_area;
    v[0]=src.v[0];
    v[1]=src.v[1];
    n_copies=src.n_copies;
    
    sum_x=src.sum_x;
    sum_y=src.sum_y;
    sum_xx=src.sum_xx;
    sum_yy=src.sum_yy;
    sum_xy=src.sum_xy;

    velocity_histories_x = src.velocity_histories_x;
    velocity_histories_y = src.velocity_histories_y;
    prev_com_x = src.prev_com_x;
    prev_com_y = src.prev_com_y;
    avg_vx= src.avg_vx;
    avg_vy= src.avg_vy;
    
    length=src.length;
    target_length=src.target_length;
    owner=src.owner;

    lambda_2 = src.lambda_2;
    lambda = src.lambda;
    c_type=src.c_type;

    time_created = src.time_created;

    xcen = src.xcen;
    ycen = src.ycen;
    xcens = src.xcens;
    ycens = src.ycens;

    div_times = src.div_times;

    perimeter = src.perimeter;
    target_perimeter = src.target_perimeter;

    cell_perim_constraint = src.cell_perim_constraint;
    apop_noise_state=src.apop_noise_state;


    Sox2_concentration = src.Sox2_concentration;
    sox2_internal_adhesion = src.sox2_internal_adhesion;
    Sox17_concentration = src.Sox17_concentration;
    sox17_internal_adhesion = src.sox17_internal_adhesion;
    lonely_cell = src.lonely_cell;

    epithelial = src.epithelial;
    synNotch_bound = src.synNotch_bound;
    synNotch_unbound = src.synNotch_unbound;
    synNotch_intra = src.synNotch_intra;
    E_cadherin = src.E_cadherin;
    CD19 = src.CD19;
    opposing_CD19 = src.opposing_CD19;
    random_binding_proteins = src.random_binding_proteins;
    touching_med = src.touching_med;
    mCherry=src.mCherry;
    GFP=src.GFP;
    opposing_GFP=src.opposing_GFP;
    P_cadherin=src.P_cadherin;
    N_cadherin=src.N_cadherin;

    opposing_CD19 = src.opposing_CD19;
    opposing_E_cadherin = src.opposing_E_cadherin;
    opposing_mCherry=src.opposing_mCherry;
    opposing_GFP=src.opposing_GFP;
    opposing_N_cadherin=src.opposing_N_cadherin;
    opposing_P_cadherin=src.opposing_P_cadherin;

    spheroid_cell=src.spheroid_cell;
    constitutives=src.constitutives;
    GFP_induced=src.GFP_induced;
    mCherry_induced=src.mCherry_induced;
    CD19_induced=src.CD19_induced;

    cell_elastic_mod=src.cell_elastic_mod;
    motility_strength=src.motility_strength;

    leftover_area=src.leftover_area;

    f_opposing_GFP = src.f_opposing_GFP;
    f_opposing_CD19 = src.f_opposing_CD19;
    f_opposing_mCherry=src.f_opposing_mCherry;
    f_opposing_E_cad = src.f_opposing_E_cad;
    f_opposing_N_cad = src.f_opposing_N_cad;
    f_opposing_P_cad = src.f_opposing_P_cad;
  
    accumulated_death_signals = src.accumulated_death_signals;
    medium_touch_count=src.medium_touch_count;


    diffs = new double[par.n_diffusers];

    for (int i=0;i<par.n_diffusers;i++)
    {
      diffs[i]=src.diffs[i];
    }


    chem = new double[par.n_chem];
    for (int ch=0;ch<par.n_chem;ch++)
      chem[ch]=src.chem[ch];
    
    return *this;

    centerx = src.centerx;
    centery = src.centery;

  }

  /*! \brief Returns false if Cell has apoptosed (vanished). */
  inline bool AliveP(void) const {
    return alive;
  }
  
  inline void makeAlive(void) {
    alive = true;
  }

  //! Returns the cell colour.
  inline int Colour(void) const {
   
    /* if (par.dynamicJ) 
      return colour;
      else */
      // return tau+1;
      // return colour;
      return c_type;
    
  };

  //! Set cell type of this Cell.
  inline void setTau(int settau) {
    tau=settau;
  }
  
  //! Get cell type of this Cell.
  inline int getTau(void) {
    return tau;
  }



  
  //! Set color of this cell to new_colour, irrespective of type.
  inline int SetColour(const int new_colour) {
    return colour=new_colour;
  }


  inline void set_ctype(const int col)
  {
    // cout << sigma << '\t' << c_type << endl;
    c_type = col;
    // if (sigma == 97)
    // {
    //   c_type = 2;
    //   cout << xcen << '\t' << ycen << endl;
    // }

  }


  /* \brief Returns the energy between this cell and cell2. 

  Called from CellularPotts::DeltaH.
  **/
  // int EnergyDifference(Cell &cell2) const;
  double EnergyDifference(Cell &cell2); // DOM CHANGED TO NON-CONST

  double EnergyDifference(Cell &cell2, bool phase, double Jstemdiff=par.J_SL);

  double Melt(Cell &cell2, int x);

  //! Return Cell's actual area.
  inline int Area() const {
    return area;
  }
  
  //! Return Cell's target area.
  inline int TargetArea() const {
    return target_area;
  }
  

  //! Debugging function used to print the cell's current inertia tensor (as used for calculations of the length )
  inline void PrintInertia(void) {
    
    double ixx=(double)sum_xx-(double)sum_x*sum_x/(double)area;
    double iyy=(double)sum_yy-(double)sum_y*sum_y/(double)area;
    double ixy=(double)sum_xy-(double)sum_x*sum_y/(double)area;

    cerr << "ixx = " << ixx << "\n";
    cerr << "iyy = " << iyy << "\n";
    cerr << "ixy = " << ixy << "\n";
    
  }

  // return the current length
  inline double Length(void) {
    return length;
  }

/*! \brief Clears the table of J's.  

This is only important for a
feature called "DynamicJ's", where J-values depend on internal states
of the cells (such as a genetic network; see e.g. Hogeweg et
al. 2000). The current version of TST does not include such functionality.
*/
  // static void ClearJ(void);
  double polarvec[9]; 
  void RenormPolarVec(void);
  
  /*! \brief Returns the maximum cell identity number in the Dish.
    This would normally be the number of cells in the Dish, although
   the number includes apoptosed cells.
  */
  // static inline int MaxSigma() {
  //   return maxsigma;
  // }

  //! Returns the cell's cell identity number.
  inline int Sigma() const {
    return sigma;
  }
  
  //! Sets the target area of the cell.
  inline int SetTargetArea(const int new_area, double lpi, double lsi)
  {
    target_area=new_area;
    UpdatePerimeterConstraint(lpi, lsi);
    return target_area;
  }
  
  //! Sends the current cell into apoptosis
  inline void Apoptose() {
    alive=false;
  }
  
  //! Decrement the cell's target area by one unit.
  inline int IncrementTargetArea() {
    return ++target_area;
  }
  
  //! Increment the cell's target area by one unit.
  inline int DecrementTargetArea() {
    return --target_area;
  }

  //! Cell lineage tracking: get the cell's parent
  inline int Mother(void) const { return mother; }
  
  //! Cell lineage tracking: get the cell's daughter
  inline int Daughter(void) const { return daughter; }
  
  //! Returns a counter keeping track of the number of divisions
  inline int TimesDivided(void) const { return times_divided; }
  
  //! Returns Monte Carlo Step (MCS) when this cell originated.
  inline int DateOfBirth(void) const { return date_of_birth; }
  
  //! Returns the cell type at the time of birth. 
  inline int ColourOfBirth(void) const { return colour_of_birth; }

  //! Returns the bond energy J between this cell and cell c2.
  // inline int GetJ(const Cell &c2) const {
  //   return J[sigma][c2.sigma];
  // }

  double EmbryoEnergy(Cell &cell2, int zona_sigma, int zona_sigma_sticky, double t, bool debug=false);
  //! Sets bond energy J between cell type t1 and t2 to val
  // inline static int SetJ(int t1,int t2, int val) {
  //   return J[t2][t1]=J[t1][t2]=val;
  // }
  double EquilibrateEnergy(Cell &cell2, int zona_sigma, int zona_sigma_sticky, double t);




  // Deal with gradient measurements:

  //! Set the current gradient of the cell to g. Currently not in use.
  inline double* SetGrad(double *g) {
    grad[0]=g[0];
    grad[1]=g[1];
    return grad;
  }
  
  //! Returns the cell's measured gradient. Currently not in use.
  inline const double* GetGrad(void) const {
    return grad;
  } 
  
  //! Returns the cell's measured gradient. Currently not in use.
  inline double GradX() const {
    return grad[0];
  }

  //! Returns the cell's measured gradient. Currently not in use.
  inline double GradY() const {
    return grad[1];
  }

  //! Currently not in use (remove?)
  inline double* AddToGrad(double *g) {
    grad[0]+=g[0];
    grad[1]+=g[1];
    return grad;
  }
   
  //! Currently not in use (remove?)
  inline void ClearGrad(void) {
    grad[0]=0.;
    grad[1]=0.;
  }  
  
  /*! After introducing a new Cell (e.g. with GrowInCell)
    call this function to set the moments and areas right.
  */
  void MeasureCellSize(Cell &c);

private:
  /*! \brief Read a table of static Js.
    First line: number of types (including medium)
    Next lines: diagonal matrix, starting with 1 element (0 0)
    ending with n elements */
  static void ReadStaticJTable(const char *fname);

  // used internally by dish in "CellGrowthAndDivision"
  // inline int GrowthThreshold(void) const { return growth_threshold; }
  
  // used internally by class CellularPotts
  inline void CleanMoments(void) {
    sum_x = sum_y = sum_xx = sum_xy = sum_yy = area = target_area  = 0;
  }
  // used internally by class CellularPotts
  inline double AddSiteToMoments(int x,int y, double new_l=-1. ) {
    
    // Add a site to the raw moments, then update and return the
    // length of the cell
    
    // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
    // Eventually this function may be used to carry
    // out all necessary adminstration at once

    if (par.periodic_boundaries)
    {

      double Lx = par.sizex - 2;
      double Ly = par.sizey - 2;

      double com_x = sum_x / (double)area;
      double com_y = sum_y / (double)area;

      double dx = (double)x - com_x;
      double dy = (double)y - com_y;

      // Find the "unfolded" coordinate relative to current COM
      x -= (int)(Lx * round(dx / Lx));
      y -= (int)(Ly * round(dy / Ly));
    }
    sum_x+=x;
    sum_y+=y;
    sum_xx+=x*x;
    sum_yy+=y*y;
    sum_xy+=x*y;
    
    // update length (see appendix. A, Zajac.jtb03), if length is not given
    // NB. 24 NOV 2004. Found mistake in Zajac's paper. See remarks in
    // method "Length(..)". 
    if (new_l<0.) {
      length=Length(sum_x,sum_y,sum_xx,sum_yy,sum_xy,area);
    } else {
      length=new_l;
    }
    return length;
  }

  // used internally by class CellularPotts
  inline double RemoveSiteFromMoments(int x,int y, double new_l=-1.) {
    
    // Remove a site from the raw moments, then update and return the
    // length of the cell
    
    // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
    // Eventually this function may be used to carry
    // out all necessary adminstration at once

    if (par.periodic_boundaries)
    {

      double Lx = par.sizex - 2;
      double Ly = par.sizey - 2;

      double com_x = sum_x / (double)area;
      double com_y = sum_y / (double)area;

      double dx = (double)x - com_x;
      double dy = (double)y - com_y;

      // Find the "unfolded" coordinate relative to current COM
      x -= (int)(Lx * round(dx / Lx));
      y -= (int)(Ly * round(dy / Ly));
    }

    sum_x-=x;
    sum_y-=y;
    sum_xx-=x*x;
    sum_yy-=y*y;
    sum_xy-=x*y;
    
    // update length (see app. A, Zajac.jtb03), if length is not given
    if (new_l<0.) {
      length=Length(sum_x,sum_y,sum_xx,sum_yy,sum_xy,area);
    } else {
      length=new_l;
    }
    return length;
  }

  //! \brief Calculates the length based on the given inertia tensor
  //components (used internally)
  inline double Length(long int s_x,long int s_y,long int s_xx,
				long int s_yy,long int s_xy,long int n) {
    
    // inertia tensor (constructed from the raw momenta, see notebook)
    double iyy=(double)s_xx-(double)s_x*s_x/(double)n;
    double ixx=(double)s_yy-(double)s_y*s_y/(double)n;
    double ixy=-(double)s_xy+(double)s_x*s_y/(double)n;
        
    double rhs1=(ixx+iyy)/2., rhs2=sqrt( (ixx-iyy)*(ixx-iyy)+4*ixy*ixy )/2.;

    double lambda_b=rhs1+rhs2;
    //double lambda_a=rhs1-rhs2;
    
    // according to Zajac et al. 2003:
    //return 2*sqrt(lambda_b);
    // Grumble, this is not right!!!
    // Must divide by mass!!!!!!

    // see: http://scienceworld.wolfram.com/physics/MomentofInertiaEllipse.html
    //    cerr << "n = " << n << "\n";
    return 4*sqrt(lambda_b/n);

    // 2*sqrt(lambda_b/n) give semimajor axis. We want the length.

  }
  
  // return the new length that the cell would have
  // if site (x,y) were added.
  // used internally by CellularPotts
  inline double GetNewLengthIfXYWereAdded(int x, int y) {
    
    return Length(sum_x+x,sum_y+y,
			   sum_xx+x*x,sum_yy+y*y,sum_xy+x*y,area+1);
    
  }

  // return the new length that the cell would have
  // if site (x,y) were removed
  // used internally by CellularPotts
  inline double GetNewLengthIfXYWereRemoved(int x, int y) {
    if (area == 1)
    {
      return 0;
    }

    return Length(sum_x-x,sum_y-y,
          sum_xx-x*x,sum_yy-y*y,sum_xy-x*y,area-1);

  }

  // ! Return Cell's perimeter
  inline int Perimeter() 
  { 
    return perimeter; 
  }

  // ! Return Cell's target perimeter
  inline int TargetPerimeter() { return target_perimeter; }
  // ! Set Cell's target perimeter
  inline int SetTargetPerimeter(const int new_perimeter) {
    return target_perimeter = new_perimeter;
  }
  // ! Set Cell's perimeter
  inline int SetPerimeter(const int new_perimeter) {
    return perimeter = new_perimeter;
  }

  //! Increments the cell's actual perimeter by 1 unit.
  inline int IncrementPerimeter() { return ++perimeter; }

  //! Decrements the cell's actual perimeter by 1 unit.
  inline int DecrementPerimeter() { return ++perimeter; }

  inline int IncrementTargetPerimeter() { return ++target_perimeter; }

  inline int DecrementTargetPerimeter() { return --target_perimeter; }

  inline void average_chem_synthetic()
  {
    if (par.morph_or_surface[0]==true)
      opposing_GFP = opposing_GFP / (double)(area);
    if (par.morph_or_surface[1]==true)
      opposing_mCherry = opposing_mCherry / (double)(area);
    if (par.morph_or_surface[2]==true)
      opposing_CD19 = opposing_CD19 / (double)(area);
  }


  inline double chem_conc(int d)
  {
    return diffs[d];
  }



  double CalculateJfromKeyLock(vector<bool>& key2, vector<bool>& lock2 );


  // static vector<bool> spare;

  double PhaseJ(bool &phase, double &Jstemdiff, bool &epith);

  double PhaseJwithMed();

  double phaseJfromMed();

  double J_equation(int x); 

  inline void set_lambda_2(double l)
  {
    lambda_2 = l;
  }

  inline double get_lambda_2()
  {
    return lambda_2;
  }

  inline void set_lambda(double l)
  {
    lambda = l;
  }

  void SetLambdaByBulk()
  {
    lambda = par.bulk_modulus / target_area;
  }

  inline double get_lambda()
  {
    return lambda;
  }

  inline void set_xcen()
  {
    xcen = double(sum_x) / area;
  }

  inline void set_ycen()
  {
    ycen = double(sum_y) / area;
  }


  inline void set_xcen(double x)
  {
    xcen = x;
  }

  inline void set_ycen(double y)
  {
    ycen = y;
  }

  inline double get_xcen()
  {
    return xcen;
  }

  inline double get_ycen()
  {
    return ycen;
  }

  inline void RecordMass()
  {
    xcens.push_back(xcen);
    ycens.push_back(ycen);
    cout << xcen << '\t' << ycen << '\t' << double(sum_x) / area << '\t' << double(sum_y) / area << endl;

  }


  inline vector<double>& get_xcens()
  {
    return xcens;
  }

  inline vector<double>& get_ycens()
  {
    return ycens;
  }


  inline void cellmed()
  {
    ++medium_contact;
  }

  inline void cellcell()
  {
    ++cell_contact;
  }
  inline void cellperim()
  {
    ++cell_perim;
  }


  inline int GetPhenotype()
  {
    return phenotype;
  }


  inline int get_time_created()
  {
    return time_created;
  }


  inline void SetTimeCreated(int &time)
  {
    time_created = time;
  }


  inline void MassToList()
  {
    mass_list.push_back(area);
  }

  inline vector<double>& GetMassList()
  {
    return mass_list;
  }




  inline void ResetActiveMotion()
  {
    velocity_initialised=false;
    prev_com_x=0;
    prev_com_y=0;
    velocity_histories_x.clear();
    velocity_histories_y.clear();
    velocity_histories_x.assign(par.persistence_time, 0.);
    velocity_histories_y.assign(par.persistence_time, 0.);
  }

  /* active matter methods */

  // called every time step to determine the direction of cell motion over last N steps.
  inline void update_velocity()
  {
    // calculate velocity here

    double com_x = double(sum_x) / double(area);
    double com_y = double(sum_y) / double(area);

    // if (prev_com_x < 0.0001)
    // {
    //   prev_com_x = com_x;
    //   prev_com_y = com_y;
    // }
    if (!velocity_initialised)
    {
      prev_com_x = com_x;
      prev_com_y = com_y;
      velocity_initialised = true;
      // Initialize histories with 0 so the cell starts neutral
      return; 
    }


    double v_x = com_x - prev_com_x;
    double v_y = com_y - prev_com_y;

    // double Lx = par.sizex -2;
    // double Ly = par.sizey -2;
    // if (par.periodic_boundaries)
    // {
    //     if (v_x >  Lx / 2.0) v_x -= Lx;
    //     if (v_x < -Lx / 2.0) v_x += Lx;
    //     if (v_y >  Ly / 2.0) v_y -= Ly;
    //     if (v_y < -Ly / 2.0) v_y += Ly;
    // }
    // need function 
    avg_vx -= velocity_histories_x.back() / par.persistence_time;
    avg_vy -= velocity_histories_y.back() / par.persistence_time;
    avg_vx += v_x / par.persistence_time;
    avg_vy += v_y / par.persistence_time;

    // cout << "velocity debugging: " << avg_vx << '\t' << avg_vy << endl;
    
    velocity_histories_x.push_front(v_x);
    velocity_histories_x.pop_back();

    velocity_histories_y.push_front(v_y);
    velocity_histories_y.pop_back();

    prev_com_x = com_x;
    prev_com_y = com_y;
  }

  inline double cell_velx()
  {
    return avg_vx;
  }
  inline double cell_vely()
  {
    return avg_vy;
  }

  inline double ActiveDotProduct_added(int x, int y)
  {
    double Lx = par.sizex - 2;
    double Ly = par.sizey - 2;

    com_x = (double)sum_x / area;
    com_y = (double)sum_y / area;

    double dx = (double)x - com_x;
    double dy = (double)y - com_y;


    
    if (par.periodic_boundaries) 
    {
        dx -= Lx * round(dx / Lx);
        dy -= Ly * round(dy / Ly);

        // if (dx >  Lx / 2.0) dx -= Lx;
        // if (dx < -Lx / 2.0) dx += Lx;
        // if (dy >  Ly / 2.0) dy -= Ly;
        // if (dy < -Ly / 2.0) dy += Ly;
    }

    // Displacement of COM: dCOM = (x - COM_old) / (Area + 1)
    com_shiftx = dx / (double)(area + 1);
    com_shifty = dy / (double)(area + 1);

    // Energy contribution: area * (dCOM . Velocity)
    return (double)area * (com_shiftx * avg_vx + com_shifty * avg_vy);


    // double dirx = double(sum_x+x)/double(area+1) - double(sum_x)/double(area);
    // double diry = double(sum_y+y)/double(area+1) - double(sum_y)/double(area);

    // double dot_product = area * (dirx * avg_vx + diry * avg_vy);
    // return dot_product; 
  }

  inline double ActiveDotProduct_removed(int x, int y)
  {
    double Lx = par.sizex - 2;
    double Ly = par.sizey - 2;

    com_x = (double)sum_x / area;
    com_y = (double)sum_y / area;

    // Vector from COM to the pixel being removed
    double dx = (double)x - com_x;
    double dy = (double)y - com_y;

    if (par.periodic_boundaries) 
    {
        dx -= Lx * round(dx / Lx);
        dy -= Ly * round(dy / Ly);
        // if (dx >  Lx / 2.0) dx -= Lx;
        // if (dx < -Lx / 2.0) dx += Lx;
        // if (dy >  Ly / 2.0) dy -= Ly;
        // if (dy < -Ly / 2.0) dy += Ly;
    }

    // cout << dx << endl;
    // if (abs(dx) > 20)
    //   cout << x << '\t' << y << '\t' << com_x << '\t' << com_y << '\t' << dx << endl;

    // Displacement of COM: dCOM = (COM_old - x) / (Area - 1)
    // Note: Removing a pixel moves the COM in the opposite direction
    com_shiftx = -dx / (double)(area - 1);
    com_shifty = -dy / (double)(area - 1);
    // if (abs(toreturn) > 1)
    // {
    //   cout << shift_x << '\t' << avg_vx << '\t' << shift_y << '\t' << avg_vy << '\t' << toreturn << endl;
    // }
      // cout << toreturn << endl;
    // cout << (double)area * (shift_x * avg_vx + shift_y * avg_vy) << endl;
    return (double)area * (com_shiftx * avg_vx + com_shifty * avg_vy);
}  

  double Gravity()
  {
    // x^2 gravity, we typically use coeffiecnt of..
    double newcom_x = com_x + com_shiftx;
    double newcom_y = com_y + com_shifty;
    double delta_x = com_x - centerx;
    double delta_y = com_y - centery;
    double delta_xnew = newcom_x - centerx;
    double delta_ynew = newcom_y - centery;
    // cout << centerx << '\t' << newcom_x << '\t' << com_x << endl;
    double old_energy = par.lambda_gravity * (delta_x * delta_x + delta_y * delta_y);
    double new_energy = par.lambda_gravity * (delta_xnew * delta_xnew + delta_ynew * delta_ynew);
    return new_energy - old_energy;

    // x^4 gravity, coefficient of approx 0.00000008;
      // 1. Calculate the squared distance from the center for the current position
    // double delta_x = com_x - centerx;
    // double delta_y = com_y - centery;
    // double dist_sq_old = delta_x * delta_x + delta_y * delta_y;

    // // 2. Calculate the squared distance for the proposed new position
    // double newcom_x = com_x + com_shiftx;
    // double newcom_y = com_y + com_shifty;
    // double delta_xnew = newcom_x - centerx;
    // double delta_ynew = newcom_y - centery;
    // double dist_sq_new = delta_xnew * delta_xnew + delta_ynew * delta_ynew;
    // // cout << centerx << '\t' << newcom_x << '\t' << com_x << endl;

    // // 3. Compute quartic energy: E = lambda * (r^2)^2 = lambda * r^4
    // // This creates a flat bottom and steep walls.
    // double old_energy = par.lambda_gravity * (dist_sq_old * dist_sq_old);
    // double new_energy = par.lambda_gravity * (dist_sq_new * dist_sq_new);

    // return new_energy - old_energy;
  }


/* synthetic structure methods */


void setsynNotch_bound(double new_value)
{
  synNotch_bound = new_value;
}
double& getsynNotch_bound()
{
  return synNotch_bound;
}

void setsynNotch_unbound(double new_value)
{
  synNotch_unbound = new_value;
}
double& getsynNotch_unbound()
{
  return synNotch_unbound;
}

void setsynNotch_intra(double new_value)
{
  synNotch_intra = new_value;
}
double& getsynNotch_intra()
{
  return synNotch_intra;
}

void setE_cadherin(double new_value)
{
  E_cadherin = new_value;
}
double& getE_cadherin()
{
  return E_cadherin;
}

void setRandomBindingProteins(double new_value)
{
  random_binding_proteins = new_value;
}
double& getRandomBindingProteins()
{
  return random_binding_proteins;
}

void setmCherry(double new_value)
{
  mCherry = new_value;
}
double& getmCherry()
{
  return mCherry;
}
void setGFP(double new_value)
{
  GFP = new_value;
}
double& getGFP()
{
  return GFP;
}

double& getOppositeGFP()
{
  return f_opposing_GFP;
}

double& getN_cadherin()
{
  return N_cadherin;
}

void setN_cadherin(double new_value)
{
  N_cadherin = new_value;
}
double& getP_cadherin()
{
  return P_cadherin;
}

void setP_cadherin(double new_value)
{
  P_cadherin = new_value;
}

double& getOpposingN_cadherin()
{
  return f_opposing_N_cad;
}

void setOpposingN_cadherin(double new_value)
{
  opposing_N_cadherin = new_value;
}
double& getOpposingP_cadherin()
{
  return f_opposing_P_cad;
}

void setOpposingP_cadherin(double new_value)
{
  opposing_P_cadherin = new_value;
}



void setCD19(double new_value)
{
  CD19 = new_value;
}
double& getCD19()
{
  return CD19;
}

double& getOpposingCD19()
{
  return f_opposing_CD19;
}

double& getOpposing_E_cadherin()
{
  return f_opposing_E_cad;
}

void ResetTempSurfaceBindings()
{
  if (!par.morph_or_surface[0])
    opposing_GFP=0;
  if (!par.morph_or_surface[1])
    opposing_mCherry=0;    
  if (!par.morph_or_surface[2])
    opposing_CD19=0;

  opposing_E_cadherin=0;
  opposing_N_cadherin=0;
  opposing_P_cadherin=0;
  
}


void ResetFinalSurfaceBindings()
{
  f_opposing_GFP=0;
  f_opposing_CD19=0;
  f_opposing_E_cad=0;
  f_opposing_N_cad=0;
  f_opposing_P_cad=0;
  f_opposing_mCherry=0;
  
}

void AddtoSurfaces(bool bcd19, double bE_cad, double bGFP, double bP_cad, double bN_cad)
{
  if (!par.morph_or_surface[0])
  {
    if (bGFP > 1)
      opposing_GFP = opposing_GFP + 1.;
    else
      opposing_GFP = opposing_GFP + bGFP;
  }
  if (!par.morph_or_surface[2])
  {
    opposing_CD19=opposing_CD19 + bcd19;
  }

  if (bE_cad > 1)
    opposing_E_cadherin = opposing_E_cadherin + 1.;
  else
    opposing_E_cadherin = opposing_E_cadherin + bE_cad;

  if (bP_cad > 1)
    opposing_P_cadherin = opposing_P_cadherin + 1.;
  else
    opposing_P_cadherin = opposing_P_cadherin + bP_cad;
  
  if (bN_cad)
    opposing_N_cadherin = opposing_N_cadherin + 1.;
  else
    opposing_N_cadherin = opposing_N_cadherin + bN_cad;
  // cout << "ADDING:" << opposing_CD19 << '\t' << opposing_E_cadherin << endl;
}

void AverageSurfaceBindings()
{
  // cout << opposing_CD19 << '\t' << perimeter << endl;
  if (!par.morph_or_surface[0])
    opposing_GFP = opposing_GFP / double(perimeter);
  if (!par.morph_or_surface[2])
    opposing_CD19 = opposing_CD19 / double(perimeter);

  opposing_E_cadherin = opposing_E_cadherin / double(perimeter);
  opposing_N_cadherin = opposing_N_cadherin / double(perimeter);
  opposing_P_cadherin = opposing_P_cadherin / double(perimeter);
  
  if (opposing_GFP > f_opposing_GFP) f_opposing_GFP=opposing_GFP;
  if (opposing_CD19 > f_opposing_CD19) f_opposing_CD19=opposing_CD19;
  if (opposing_E_cadherin > f_opposing_E_cad) f_opposing_E_cad=opposing_E_cadherin;
  if (opposing_N_cadherin > f_opposing_N_cad) f_opposing_N_cad=opposing_N_cadherin;
  if (opposing_P_cadherin > f_opposing_P_cad) f_opposing_P_cad=opposing_P_cadherin;

  ResetTempSurfaceBindings();

}


bool& isSpheroid()
{
  return spheroid_cell;
}

void setSpheroid(bool S)
{
  spheroid_cell = S;
}

vector<bool> GetConstitutives()
{
  return constitutives;
}
vector<bool>GetGFP_induced()
{
  return GFP_induced;
}
vector<bool>GetMcherry_induced()
{
  return mCherry_induced;
}
vector<bool>GetCD19_induced()
{
  return CD19_induced;
}

void SetConstitutives(vector<bool> incoming)
{
  constitutives=incoming;
}
void SetGFP_induced(vector<bool> incoming)
{
  GFP_induced=incoming;
}
void SetMcherry_induced(vector<bool> incoming)
{
  mCherry_induced=incoming;
}
void SetCD19_induced(vector<bool> incoming)
{
  CD19_induced=incoming;
}

void SetMotilityStrength(double mots)
{
  motility_strength=mots;
}

double& GetMotilityStrength()
{
  return motility_strength;
}

void SetElasticMod(double emod)
{
  cell_elastic_mod=emod;
}

double& GetElasticMod()
{
  return cell_elastic_mod;
}

void LeftoverTargetArea(double fta)
{
  leftover_area+=fta;
  while (leftover_area > 1)
  {
    ++target_area;
    leftover_area-=1;
  }
}

double& GetFauxTargetArea()
{
  return leftover_area;
}


inline void UpdateAdhesions() 
{
  double t2  = 1.0 / (1.0 + std::exp(-par.switch_like * (Sox2_concentration  - par.sox_threshold)));
  double t17 = 1.0 / (1.0 + std::exp(-par.switch_like * (Sox17_concentration - par.sox_threshold)));

  // 2. Apply mutual inhibition so that if BOTH are ~1, both adhesions become ~0
  sox2_internal_adhesion  = t2  * (1.0 - t17);
  sox17_internal_adhesion = t17 * (1.0 - t2);
}

inline double getSox2adhesion()
{
  return sox2_internal_adhesion;
}

inline double getSox17adhesion()
{
  return sox17_internal_adhesion;
}

inline double getSox2()
{
  return Sox2_concentration;
}

inline void setSox2(double newsox, double perim_increase, double sox17perim_increase)
{
  Sox2_concentration=newsox;
  UpdatePerimeterConstraint(perim_increase, sox17perim_increase);

  UpdateAdhesions();

}

inline double getSox17()
{
  return Sox17_concentration;
}

inline void setSox17(double newsox, double perim_increase, double sox17perim_increase)
{
  Sox17_concentration=newsox;
  UpdatePerimeterConstraint(perim_increase, sox17perim_increase);

  UpdateAdhesions();

  // cout << sox2_internal_adhesion << '\t' << sox17_internal_adhesion << endl;
  
}

inline void OutputPerim()
{
  cout << "here: " << sigma << '\t' << perimeter << '\t' << target_perimeter << '\t' << cell_perim_constraint << '\t' << area << '\t' << target_area << '\t' << lambda << endl;
}

inline void SetSoxColour(double t)
{
    // 1. Fixed Sox2 / Sox17 baseline colour (Index 2 to 202)
    double weight = 0.5 * (sox2_internal_adhesion - sox17_internal_adhesion + 1.0);
    double target_offset = (weight - 0.5) * 200.0;
    
    int index = 102 + static_cast<int>(std::round(target_offset));
    index = std::clamp(index, 2, 202);

    set_ctype(index);
    

    double is_loser = std::max(sox2_internal_adhesion * sox17_internal_adhesion, 
                                (1.0 - sox2_internal_adhesion) * (1.0 - sox17_internal_adhesion));
    if (is_loser > 0.9)
    {
        double t_clamped = std::clamp(t, 0.0, 1.0);
        int start=203;
        if (par.set_loser_colours)
          start=280;        
        // t = 0 -> index 203 (Blue)
        // t = 0.5 -> index 253 (Yellow)
        // t = 1 -> index 303 (Red)
        int loser_index = start + static_cast<int>(std::round(t_clamped * 100.0));
        loser_index = std::clamp(loser_index, 203, 303);

        set_ctype(loser_index);
    }
}

inline void MakeLonely(bool lonely)
{
  lonely_cell=true;
}

inline bool IsLonely()
{
  return lonely_cell;
}

inline double& GetDeathSignals()
{
  return accumulated_death_signals;
}

inline void SetShapeIndex(double sindex)
{
  shape_index=sindex;
}

inline double GetShapeIndex()
{
  return shape_index;
}

// Which mechanism started killing this cell - lonely/blastocoel extrusion
// (ToxictoLonelyCells) or neighbour-competition signalling
// (NeighbourBasedApoptosis). Both shrink the cell identically, so this has
// to be recorded at the point the shrink is triggered; by the time the CPM
// dynamics actually bring the cell's area to 0 (ConvertSpin), the cause is
// no longer distinguishable from the mechanics alone.
static const int DEATH_CAUSE_NONE = 0;
static const int DEATH_CAUSE_LONELY = 1;
static const int DEATH_CAUSE_SIGNAL = 2;

inline int GetDeathCause() const
{
  return death_cause;
}

//! Records why this cell started dying, unless a cause is already recorded
//! (first trigger wins - a cell can become lonely and cross the signal
//! threshold on different steps).
inline void MarkDeathCause(int cause)
{
  if (death_cause == DEATH_CAUSE_NONE)
    death_cause = cause;
}

inline void ClearDeathCause()
{
  death_cause = DEATH_CAUSE_NONE;
}



inline void UpdatePerimeterConstraint(double loser_perim_inc, double sox17_perim_inc)
{
  double is_looser = max(sox2_internal_adhesion * sox17_internal_adhesion, (1. - sox2_internal_adhesion) * (1. - sox17_internal_adhesion));
  double added_perim1 = loser_perim_inc * static_cast<int>(round((par.ptarget_perimeter) * sqrt(double(target_area) / double(par.cell_target_area)))) * is_looser;

  // hypoblast cells must have higher perimeter (maybe make this dynamic eventually?)
  double added_perim2 = sox17_internal_adhesion * sox17_perim_inc * static_cast<int>(round((par.ptarget_perimeter) * sqrt(double(target_area) / double(par.cell_target_area))));


  target_perimeter = static_cast<int>(round((par.ptarget_perimeter) *
      sqrt(double(target_area) / double(par.cell_target_area)))) + added_perim1 + added_perim2 - par.perim_offset;

  // if (sox17_perim_inc > 0.1 && sox17_internal_adhesion > 0.5)
  //   cout << "sox17: " << added_perim2 << endl;
  
  // if (loser_perim_inc > 0.1 && is_looser > 0.2)
  //   cout << "is looser: " << is_looser << '\t' << added_perim1 << endl;

  // cout << sqrt(double(target_area) / double(par.cell_target_area)) << endl;
  // out << target_perimeter << '\t' << target_area << endl;
  if (target_perimeter < 20)
    target_perimeter=1;

  cell_elastic_mod = par.elastic_modulus;
  cell_perim_constraint = cell_elastic_mod / double(target_perimeter);


}

void setDivisionTimes(vector<int> tt)
{
  div_times=tt;
}

vector<int> getDivisionTimes()
{
  return div_times;
}



void setTouchingMed(bool is)
{
  touching_med = is;
}

bool& getTouchingMed()
{
  return touching_med;
}

void setPerimConstraint(double is)
{
  cell_perim_constraint = is;
}

double& getPerimConstraint()
{
  return cell_perim_constraint;
}

double& GetApopNoiseState()
{
  return apop_noise_state;
}

bool CheckLooser()
{
  return looser_cell;
}

inline void SetLooser()
{
  double is_looser = max(sox2_internal_adhesion * sox17_internal_adhesion, (1. - sox2_internal_adhesion) * (1. - sox17_internal_adhesion));
  if (is_looser>0.5)
  {
    looser_cell=true;
  }
  else
  {
    looser_cell=false;
  }
}

inline void IncrementMediumTouchCount() 
{
    ++medium_touch_count;
}

inline int GetMediumTouchCount() const 
{
    return medium_touch_count;
}

inline double CheckLooserValue()
{
  return max(sox2_internal_adhesion * sox17_internal_adhesion, (1. - sox2_internal_adhesion) * (1. - sox17_internal_adhesion));
 
}

inline void AddMedCount()
{
  ++medcount;
}

inline void AddNonMedCount()
{
  ++notmedcount;
}

inline void ResetMedCounts()
{
  medcount=0;
  notmedcount=0;
}

inline double GetProportionMed()
{
  if ((notmedcount+medcount)==0)
  {
    return 0;
  }
  return static_cast<double>(medcount) / static_cast<double>(medcount + notmedcount);
}

private:
//! Increments the cell's actual area by 1 unit.
  inline int IncrementArea() {
    return ++area;
  }

  //! Decrements the cell's actual area by 1 unit.
  inline int DecrementArea() {
    return --area;
  }

  
  /*! \brief Sets target area to actual area, to remove "pressure".

  This is useful when reading an initial condition from an image.
  */
  inline int SetAreaToTarget(void) {
    return area=target_area;
  }

  //! Called whenever a cell is constructed, from constructor
  void ConstructorBody(int settau=1);
  
  // returns the maximum cell type index
  // (depends on Jtable)
  // static int MaxTau(void) {
  //   return maxtau;
  // }

protected:
  int colour;
  bool alive;
  int sigma; // cell identity, 0 if medium
  int tau; // Cell type, when dynamicJ's are not used

  bool epithelial=false;

  // Two dimensional (square) array of ints, containing the J's.
  double length; // length of the cell;
  double target_length;

  // Dynamically increased when cells are added to the system
  // unless a static Jtable is used (currently this is the default situation)
  // Amount: the number of Cell instantations, INCLUDING copies
  // For internal use only.
  // Reading amount is NOT the way to get the number of cells!!
  // static int  **J;
  // static int maxtau;
  // static int amount; // Dom removed static
  // static int capacity;
  // static int maxsigma; // the last cell identity number given out


  double Sox2_concentration;
  double sox2_internal_adhesion;
  double Sox17_concentration;
  double sox17_internal_adhesion;
  bool lonely_cell;

  int medium_touch_count{};

  //current state of the cell
  int phenotype;

  vector<double> mass_list;

  double apop_noise_state{};
  
  int medcount{};
  int notmedcount{};

  double shape_index{};

  int perimeter;        // amount of cell's membrane
  int target_perimeter; // cell's target membrane length  

  bool exposed{true};

  int c_type{2};

  /* parameters for synthetic structures */
  // we are going to change this so that concentrations depend on cell size
  double synNotch_bound{};
  double synNotch_unbound{};
  double synNotch_intra{};
  double E_cadherin{};
  double N_cadherin{};
  double P_cadherin{};  
  double CD19{};

  double opposing_CD19{};
  double opposing_E_cadherin{};
  double opposing_mCherry{};
  double opposing_GFP{};
  double opposing_N_cadherin{};
  double opposing_P_cadherin{};

  double random_binding_proteins{};
  double GFP{};
  double mCherry{};

  bool spheroid_cell{};

  double f_opposing_GFP;
  double f_opposing_CD19;
  double f_opposing_E_cad;
  double f_opposing_N_cad;
  double f_opposing_P_cad;
  double f_opposing_mCherry;

  vector<bool>constitutives;
  vector<bool>GFP_induced;
  vector<bool>mCherry_induced;
  vector<bool>CD19_induced;
  
  bool touching_med{};

  double cell_perim_constraint;

  double motility_strength;
  double leftover_area;

  vector<int> div_times;

  bool looser_cell;

  vector<double> xcens;
  vector<double> ycens;
  int time_created=0;



  // determine energy change to reach target length
  double lambda_2{};

  double lambda{};

  // mass centers. NOTE: THESE ARE NOT INHERITED!
  double xcen{};
  double ycen{};

  int medium_contact{};
  int cell_contact{};
  int cell_perim{};

  double accumulated_death_signals{};
  int death_cause{};

  // indices of mother and daughter
  // (Note: no pointers, cells may be relocated)
  int mother;
  int daughter;
  int times_divided;
  int date_of_birth;
  int colour_of_birth;
  
  int area;
  int target_area;

  /* for energy calculations shift in c.o.m*/
  double com_x;
  double com_y;
  double com_shiftx;
  double com_shifty;
  double centerx;
  double centery;

  double cell_elastic_mod;
  

  double v[2];
  int n_copies; // number of expansions of this cell
  // gradient of a chemical (to be extended to the total number chemicals)
  double grad[2];
  
  double *diffs; // concentration of diffusers based on PDE field. 

  double *chem;
  // Raw moments of the cells
  // Are used to calculate minor and major axes
  // and center of mass
  // are locally adjusted, so axes are easily
  // and quickly calculated!


  // active motion terms
  deque<double> velocity_histories_x;
  deque<double> velocity_histories_y;
  double prev_com_x = 0;
  double prev_com_y = 0;
  double avg_vx=0;
  double avg_vy=0;
  bool velocity_initialised=false;
  
  // N.B: N is area!
  
  long int sum_x;
  long int sum_y;
  long int sum_xx;
  long int sum_yy;
  long int sum_xy;
  
  Dish *owner; // pointer to owner of cell. Dom changed from const so that amount can be counted. 

};

#endif
