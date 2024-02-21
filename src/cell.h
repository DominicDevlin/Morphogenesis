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

    fitness=src.fitness;
    genes=src.genes;
    diff_genes=src.diff_genes;
    lambda_2 = src.lambda_2;
    lambda = src.lambda;
    c_type=src.c_type;

    // more network things
    locks=src.locks;
    locks_bool=src.locks_bool;
    keys=src.keys;
    keys_bool=src.keys_bool;
    medp=src.medp;
    medp_bool=src.medp_bool;
    full_set=src.full_set;
    cycles=src.cycles;
    gene_recordings=src.gene_recordings;

    stemness = src.stemness;
    shrinker = src.shrinker;

    div_time = src.div_time;
    div_phen = src.div_phen;
    phenotype_history = src.phenotype_history;
    phentime = src.phentime;
    adulttime = src.adulttime;
    phenotype = src.phenotype;
    switches = src.switches;
    long_switches = src.switches;

    xcen = src.xcen;
    ycen = src.ycen;

    xcens = src.xcens;
    ycens = src.ycens;
    vel_phens = src.vel_phens;

    time_created = src.time_created;


    gamma_list = src.gamma_list;
    mass_list = src.mass_list;


    diffs = new double[par.n_diffusers];

    for (int i=0;i<par.n_diffusers;i++)
    {
      diffs[i]=src.diffs[i];
    }

    
    chem = new double[par.n_chem];
    for (int ch=0;ch<par.n_chem;ch++)
      chem[ch]=src.chem[ch];
    
    
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
    
    length=src.length;
    target_length=src.target_length;
    owner=src.owner;

    fitness=src.fitness;
    genes=src.genes;
    diff_genes=src.diff_genes;
    lambda_2 = src.lambda_2;
    lambda = src.lambda;
    c_type=src.c_type;


    locks=src.locks;
    locks_bool=src.locks_bool;
    keys=src.keys;
    keys_bool=src.keys_bool;
    medp=src.medp;
    medp_bool=src.medp_bool;
    full_set=src.full_set;
    cycles=src.cycles;
    gene_recordings=src.gene_recordings;

    div_time = src.div_time;
    div_phen = src.div_phen;
    phenotype_history = src.phenotype_history;
    phentime = src.phentime;
    adulttime = src.adulttime;
    phenotype = src.phenotype;
    switches = src.switches;
    long_switches = src.switches;

    time_created = src.time_created;

    stemness = src.stemness;
    shrinker = src.shrinker;

    xcen = src.xcen;
    ycen = src.ycen;
    xcens = src.xcens;
    ycens = src.ycens;
    vel_phens = src.vel_phens;

    gamma_list = src.gamma_list;
    mass_list = src.mass_list;



    diffs = new double[par.n_diffusers];

    for (int i=0;i<par.n_diffusers;i++)
    {
      diffs[i]=src.diffs[i];
    }


    chem = new double[par.n_chem];
    for (int ch=0;ch<par.n_chem;ch++)
      chem[ch]=src.chem[ch];
    
    return *this;

  }

  /*! \brief Returns false if Cell has apoptosed (vanished). */
  inline bool AliveP(void) const {
    return alive;
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

  inline int GetColour(void)
  {
    return colour;
  }



  inline void set_ctype(const int col)
  {
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
  inline double TargetLength() const {
    return sqrt(area);
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

  inline double get_fitness(void)
  {
    return fitness;
  }

  inline void set_fitness(double fit)
  {
    fitness = fit;
  }

  inline void set_genes(vector<double> new_g)
  {
    genes = new_g;
  }

  inline vector<double>& get_genes(void)
  {
    return genes;
  }

  inline void print_genes(void)
  {
    cout << "Concentrations: ";
    for (int i = 0; i < genes.size(); ++i)
    {
      if (i < par.n_diffusers)
        cout << " - Gene " << i+1 << " is: " << diff_genes.at(i);
      else
        cout << " - Gene " << i+1 << " is: " << genes.at(i);
    }
    cout << endl << "Diffuser proteins: ";
    for (int i = 0; i < par.n_diffusers; ++i)
      cout << " - Diffuser " << i+1 << " is: " << genes.at(i);
    cout << endl;
    cout << "lock bools: ";
    for (int i = 0; i < par.n_locks; ++i) 
      cout << locks_bool.at(i) << "  ";
    cout << "key bools: ";
    for (int i = 0; i < par.n_locks; ++i) 
      cout << keys_bool.at(i) << "  ";

    cout << "medp bools: ";
    for (int i = 0; i < par.n_mediums; ++i) 
      cout << medp_bool.at(i) << "  ";
    cout << endl; 


    cout << "Target length is: " << target_length << endl;

  }


  inline vector<double>& get_diffusers(void)
  {
    return diff_genes;
  }

  inline void set_diffusers(vector<double> new_diff)
  {
    diff_genes = new_diff;
  }


  inline void set_locks(vector<double>& nlock)
  {
    locks = nlock;
  }
  inline void set_keys(vector<double>& nkeys)
  {
    keys = nkeys;
  }

  inline vector<double>& get_locks(void)
  {
    return locks;
  }
  inline vector<double>& get_keys(void)
  {
    return keys;
  }
  inline vector<bool>& get_locks_bool(void)
  {
    return locks_bool;
  }
  inline vector<bool>& get_keys_bool(void)
  {
    return keys_bool;
  }

  inline void set_meds(vector<double>& nmeds)
  {
    medp = nmeds;
  }
  inline vector<double>& get_medp(void)
  {
    return medp;
  }
  inline vector<bool>& get_medp_bool(void)
  {
    return medp_bool;
  }


  inline vector<bool>& get_set(void)
  {
    return full_set;
  }

  inline void set_lists(void)
  {
    locks.resize(par.n_locks);
    keys.resize(par.n_locks);
    locks_bool.resize(par.n_locks);
    keys_bool.resize(par.n_locks);
    medp.resize(par.n_mediums);
    medp_bool.resize(par.n_mediums);
    full_set.resize(par.n_lockandkey + par.n_length_genes + par.n_mediums);
    cycles.resize(par.cycle_size);
    gene_recordings.resize(par.n_genes+par.n_diffusers);
  }

  

  inline void average_chem()
  {
    for (int i=0; i<par.n_diffusers;++i)
    {
      diffs[i] = diffs[i] / (double)(area);
      if (par.limit_morph)
      {
        if (diffs[i] > par.limit_amount)
          diffs[i] = par.limit_amount;
      }
      genes[i] = diffs[i];
    }
  }

  inline double chem_conc(int d)
  {
    return diffs[d];
  }



  double CalculateJfromKeyLock(vector<bool>& key2, vector<bool>& lock2 );

  double CalculateJwithMed(void);

  // static vector<bool> spare;

  double CalculateJfromMed(vector<bool>& medp2);


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

  // used for checking cycles in network. Holds the booleanised network from the last 4 update steps. Depracated.
  inline void add_to_cycle()
  {
    if (cycles.size() > par.cycle_size - 1)
      cycles.erase(cycles.begin());
    
    cycles.push_back(full_set);
  }

  bool checkforcycles(int max);


  void add_to_vectors();


  inline vector<vector<double>>& get_history(void)
  {
    return gene_recordings;
  }


  inline void set_stem(bool fate)
  {
    stemness = fate;
  }

  inline bool get_stem(void)
  {
    return stemness;
  }

  inline void set_shrink(bool fate)
  {
    shrinker = fate;
  }

  inline bool get_shrink(void)
  {
    return shrinker;
  }

  inline void set_xcen(double x)
  {
    xcen = x;
  }

  inline void set_ycen(double y)
  {
    ycen = y;
  }

  inline void RecordMass()
  {
    xcens.push_back(xcen);
    ycens.push_back(ycen);
  }

  inline vector<double>& get_xcens()
  {
    return xcens;
  }

  inline vector<double>& get_ycens()
  {
    return ycens;
  }

  inline vector<int>& get_velphens()
  {
    return vel_phens;
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



  inline vector<tuple<int,int,uint64_t>>& get_switches()
  {
    return switches; 
  }

  inline vector<tuple<int,int,uint64_t>>& get_long_switches()
  {
    return long_switches; 
  }



  inline void AddPhenotype() // how long each cell spends in each phentoype
  {
    phentime[phenotype] += 1;
  }

  inline unordered_map<int,int>& PhenTime()
  {
    return phentime;
  }

  inline void reset_recordings()
  {
    phentime.clear();
    adulttime.clear();
  }

  inline void AddType()
  {
    adulttime[phenotype] += 1;
  }

  inline unordered_map<int,int>& AdultTime()
  {
    return adulttime;
  }

  bool limit_cycle();

  inline int get_time_created()
  {
    return time_created;
  }


  inline vector<double>& get_div_x()
  {
    return div_x_cen;
  }

  inline vector<double>& get_div_y()
  {
    return div_y_cen;
  }

  inline vector<int>& get_mass_div_time()
  {
    return mass_div_time;
  }


  inline void RecordDivision(int time)
  {
    if (time > par.end_program)
    {
      pair<int, int> newp = {time, sigma};
      pair<int, int> newt = {phenotype, time};
      div_time.push_back(newp);
      div_phen.push_back(newt);
      // cout << "division recorded with phenotype: " << phenotype << endl;
    }
    time_created = time;
    
    if (par.division_anisotropy && time > par.end_program)
    {
      div_x_cen.push_back(xcens.back());
      div_y_cen.push_back(ycens.back());
      mass_div_time.push_back(xcens.size());
    }

  }


  inline vector<pair<int, int>>& get_divisions()
  {
    return div_time;
  }

  inline vector<int>& TypeHistory()
  {
    return phenotype_history;
  }

  inline vector<pair<int,int>>& DivisionPhenotype()
  {
    return div_phen;
  }

  inline void set_death_tag(bool n)
  {
    death_tag = n;
  }

  inline bool get_death_tag()
  {
    return death_tag;
  }


  inline double return_enzyme(int loc)
  {
    return genes[loc];
  }

  int LocksKeysScore(vector<bool>& locks, vector<bool>& keys);

  int CheckMedsOn();

  inline void AddGamma(double g)
  {
    gamma_list.push_back(g);
  }

  inline vector<double>& GetGamma()
  {
    return gamma_list;
  }

  inline void MassToList()
  {
    mass_list.push_back(area);
  }

  inline vector<double>& GetMassList()
  {
    return mass_list;
  }

  inline void MaxSet()
  {
    full_set.resize(par.n_functional+par.n_activators, 0);
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

  vector<pair<int,int>> div_time;
  vector<pair<int,int>> div_phen;

  vector<int> phenotype_history;

  unordered_map<int, int> phentime;
  unordered_map<int, int> adulttime;

  vector<double> gamma_list; 
  vector<double> mass_list;


  bool exposed{true};

  int c_type{2};

  double EnDif(Cell &cell2);


  // static int maxsigma; // the last cell identity number given out, Dom removed
  double fitness{};

  // I will make the gene network I vector for now.
  vector<double> genes; // should initialise to appropriate value

  vector<double> diff_genes;


  vector<double> locks{};
  vector<bool> locks_bool{};

  vector<double> keys{};
  vector<bool> keys_bool{};

  vector<double> medp{};
  vector<bool> medp_bool{};

  vector<bool> full_set{};

  vector<vector<bool>> cycles; 

  vector<vector<double>> gene_recordings;


  bool death_tag=false;

  vector<double> xcens;
  vector<double> ycens;
  vector<int> vel_phens;
  int time_created=0;

  vector<double> div_x_cen;
  vector<double> div_y_cen;
  vector<int> mass_div_time;


  vector<tuple<int,int, uint64_t>> switches;
  vector<tuple<int,int, uint64_t>> long_switches;

  // if cell is stem cell or not
  bool stemness{false};

  // if cell shrinks more easily
  bool shrinker{false};


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
  
  // N.B: N is area!
  
  long int sum_x;
  long int sum_y;
  long int sum_xx;
  long int sum_yy;
  long int sum_xy;
  
  Dish *owner; // pointer to owner of cell. Dom changed from const so that amount can be counted. 

};

#endif
