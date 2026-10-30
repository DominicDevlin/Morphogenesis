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


    shrinker = src.shrinker;


    phenotype = src.phenotype;

    perimeter = src.perimeter;
    target_perimeter = src.target_perimeter;

    phase_protein_conc = src.phase_protein_conc;
    phase_state = src.phase_state;
    medium_protein_conc = src.medium_protein_conc;
    medium_state = src.medium_state;


  

    epithelial = src.epithelial;

    centerx = src.centerx;
    centery = src.centery;
    
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

    phenotype = src.phenotype;


    shrinker = src.shrinker;

    xcen = src.xcen;
    ycen = src.ycen;


    phase_protein_conc = src.phase_protein_conc;
    phase_state = src.phase_state;
    medium_protein_conc = src.medium_protein_conc;
    medium_state = src.medium_state;


    perimeter = src.perimeter;
    target_perimeter = src.target_perimeter;

    epithelial = src.epithelial;
    
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
  
  /*! \brief Return Cell's target length.
    
  Length constraint is documented in Merks et al. 2006, Dev. Biol. 
  */
  inline double TargetLength() const 
  {
    if (phase_state)
    {
      return 2* sqrt(area/M_PI);
      // return length;
    }
    else
      return target_length;// sqrt(area)*par.target_length;
    // return target_length;
  }

  //! Set the Cell's target length
  inline double SetTargetLength(double l) {
    return target_length=l;
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
  inline int SetTargetArea(const int new_area) {
    return target_area=new_area;
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

  
  //! Sets bond energy J between cell type t1 and t2 to val
  // inline static int SetJ(int t1,int t2, int val) {
  //   return J[t2][t1]=J[t1][t2]=val;
  // }

  
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


  inline void TransformPhase(bool ptype)
  {
    if (ptype)
    {
      phase_protein_conc=1.;
      phase_state = true;
    }
    else
    {
      phase_protein_conc=0.;
      phase_state = false;
    }
  }

  inline double& get_phase_J(void)
  {
    return phase_protein_conc;
  }

  inline double& get_phase_M(void)
  {
    return medium_protein_conc;
  }

  inline bool IsEpithelia(void)
  {
    return epithelial;
  }

  inline bool& getmJ(void)
  {
    return medium_state;
  }  




  double CalculateJfromKeyLock(vector<bool>& key2, vector<bool>& lock2 );

  double CalculateJwithMed(void);

  // static vector<bool> spare;

  double CalculateJfromMed(vector<bool>& medp2);

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

  inline double get_lambda()
  {
    return lambda;
  }



  inline void set_shrink(bool fate)
  {
    shrinker = fate;
  }

  inline bool get_shrink(void)
  {
    return shrinker;
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

  inline void resetCAC()
  {
    medium_contact=0;
    cell_contact=0;
    cell_perim=0;
  }

  void Phenotype();

  int RegPhenotype();


  inline int GetPhenotype()
  {
    return phenotype;
  }

  void RecordSwitch(vector<bool> &v1, uint64_t rndm);

  void RecordLongSwitch(vector<bool> &v1, uint64_t rndm);




  inline void SetShapeIndex(double ind)
  {
    shape_index = ind;
  }

  inline double& GetShapeIndex()
  {
    return shape_index;
  }


  inline bool& GetSortingType()
  {
    return sorting_type;
  }

  inline void SetSortingType(bool t)
  {
    sorting_type = t;
    if (sorting_type)
      c_type = 4;
    else
      c_type = 90;
  }



  int LocksKeysScore(vector<bool>& locks, vector<bool>& keys);

  int CheckMedsOn();

  
  inline bool& GetPhase()
  {
    return phase_state;
  }

  void SetSheetType(bool tp)
  {
    sheet_type = tp;
    if (sheet_type)
      c_type = par.sheetcol1;
    else
      c_type = par.sheetcol2;
  }

  void SwapSheetType()
  {
    if (sheet_type)
    {
      sheet_type = false;
      c_type = par.sheetcol1;
    }
    else
    {
      sheet_type = true;
      c_type = par.sheetcol2;
    }
  }  

  bool GetSheetType()
  {
    return sheet_type;
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






inline double Gravity()
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

  //current state of the cell
  int phenotype;




  double phase_protein_conc;
  bool phase_state=false;
  double medium_protein_conc;
  bool medium_state;


  int perimeter;        // amount of cell's membrane
  int target_perimeter; // cell's target membrane length  

  bool exposed{true};

  int c_type{2};

  double EnDif(Cell &cell2);

  double SheetDif(Cell &cell2, double &sJ=par.sheet_J, double &sheetmix=par.sheetmixJ);

  double SyntheticEnergy(Cell &cell2);

  // if cell shrinks more easily
  bool shrinker{false};

  bool sorting_type{};
  // determine energy change to reach target length
  double lambda_2{};

  double lambda{};

  // mass centers. NOTE: THESE ARE NOT INHERITED!
  double xcen{};
  double ycen{};

  int medium_contact{};
  int cell_contact{};
  int cell_perim{};


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

  double lambda_perim;

  

  double v[2];
  int n_copies; // number of expansions of this cell
  // gradient of a chemical (to be extended to the total number chemicals)
  
  bool sheet_type=false;


  double *chem;
  // Raw moments of the cells
  // Are used to calculate minor and major axes
  // and center of mass
  // are locally adjusted, so axes are easily
  // and quickly calculated!

  double theta{};

  // active motion terms
  deque<double> velocity_histories_x;
  deque<double> velocity_histories_y;
  double prev_com_x = 0;
  double prev_com_y = 0;
  double avg_vx=0;
  double avg_vy=0;
  bool velocity_initialised=false;
  
  // N.B: N is area!
  double shape_index;
  
  long int sum_x;
  long int sum_y;
  long int sum_xx;
  long int sum_yy;
  long int sum_xy;
  
  Dish *owner; // pointer to owner of cell. Dom changed from const so that amount can be counted. 

};

#endif
