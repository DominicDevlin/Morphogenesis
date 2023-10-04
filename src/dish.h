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

/*! \class Dish
  \brief The virtual Petri dish.
  
  Hosts the cells with states and the CA-plane.
*/

#ifndef CRITTER_H_
#define CRITTER_H_
#include <vector>
#include "graph.h"
#include "random.h"
#include "pde.h"
#include "cell.h"
#include "ca.h"

namespace ColourMode {
  enum { State,CellType,Sigma,Auxilliary };
}

class Dish {

  friend class Info;

public:

  // Dummy constructor
  Dish();
  
  /*! \brief Init defines the initial state of the virtual
    cell culture.
    
    Define Init() in your main file describing the simulation set up,
    within the block INIT { }. See for examples vessel.cpp and
    sorting.cpp.
  */
  void Init(void);

  void ConstructorBody(void);
  
  virtual ~Dish();
 
  /*! \brief Plot the Dish to graphics window g.

  Simply calls CPM->Plot.
  */
  void Plot(Graphics *g);



  //DOM FUNCTIONS

  void start_network(vector<vector<int>> start_matrix, vector<bool> start_pol);

  void randomise_network(void);

  void AverageChemCell(void);

  void IntroduceMorphogen(int num, int xloc, int yloc);




  
  int ZygoteArea(void) const;
  
  //! Returns the number of completed Monte Carlo Steps.
  int Time(void) const;
  
  //! Returns the number of cells in the dish, excluding apoptosed cells.
  int CountCells(void) const;
   
  //! \brief. Returns the summed area of all cells in the dish
  int Area(void) const;

  //! \brief Returns the summed of all cells target area in the dish
  int TargetArea(void) const;

  //! \brief Returns the horizontal size of the dish.
  int SizeX(void);
  
  //! \brief Returns the horizontal size of the dish.
  int SizeY(void);

  //! \brief Returns a reference to cell number "c"
  inline Cell &getCell(int c) {
    return cell[c];
  }

  inline CellularPotts* get_cpm()
  {
    return CPM;
  }

  // Dom has changed these static variables to functions 
  void add_maxsigma(void);
  int get_maxsigma(void);
  void set_maxsigma(int);
  
  PDE *PDEfield;
  CellularPotts *CPM;

  // Was used for gradient measurements, not functional now.
  void ClearGrads(void);

  void MeasureChemConcentrations(void);
protected:
  //! Assign a the cell to the current Dish
  void SetCellOwner(Cell &which_cell);
  // int amount;
  int maxsigma;

private:
  bool CellLonelyP(const Cell &c, int **neighbours) const;

protected:
  //! The cells in the Petri dish; accessible to derived classes
  std::vector<Cell> cell;
};

#define INIT void Dish::Init(void)

#endif
