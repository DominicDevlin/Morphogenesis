/* 

Copyright 1995-2006 Roeland Merks, Nick Savill

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


/* CA.cpp: implementation of Glazier & Graner's Cellular Potts Model */

// This code derives from a Cellular Potts implementation written around 1995
// by Nick Savill

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include "sticky.h"
#include "random.h"
#include "ca.h"
#include "parameter.h"
#include "dish.h"
#include "sqr.h"
#include "crash.h"
#include "hull.h"
#include "connections.h"
#include <set>
#include <bits/stdc++.h>
#include <array>
#include <unordered_set>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include "fft.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <map>

using namespace std;


// #define FILESYSTEM
#ifdef FILESYSTEM
#include <filesystem>
#endif



#define ZYGFILE(Z) <Z.xpm>
#define XPM(Z) Z ## _xpm
#define ZYGXPM(Z) XPM(Z)

/* define default zygote */
/* NOTE: ZYGOTE is normally defined in Makefile!!!!!! */
#ifndef ZYGOTE
#define ZYGOTE init
#include "init.xpm"
// #else
// #include ZYGFILE(ZYGOTE)
#endif

/* STATIC DATA MEMBER INITIALISATION */
double copyprob[BOLTZMANN]; 


// X coordinates for neighbors (0 to 36)
const int CellularPotts::nx[37] = {
    0,                      // Self (r^2 = 0)
    0, 1, 0,-1,             // Level 1: r^2 = 1  (4 cells)
    1, 1,-1,-1,             // Level 2: r^2 = 2  (4 cells)
    0, 2, 0,-2,             // Level 3: r^2 = 4  (4 cells)
    1, 2, 2, 1,-1,-2,-2,-1, // Level 4: r^2 = 5  (8 cells)
    2, 2,-2,-2,             // Level 5: r^2 = 8  (4 cells)
    0, 3, 0,-3,             // Level 6: r^2 = 9  (4 cells)
    1, 3, 3, 1,-1,-3,-3,-1  // Level 7: r^2 = 10 (8 cells)
};

// Y coordinates for neighbors (0 to 36)
const int CellularPotts::ny[37] = {
    0,                      // Self
   -1, 0, 1, 0,             // Level 1
   -1, 1, 1,-1,             // Level 2
   -2, 0, 2, 0,             // Level 3
   -2,-1, 1, 2, 2, 1,-1,-2, // Level 4
   -2, 2, 2,-2,             // Level 5 
   -3, 0, 3, 0,             // Level 6
   -3,-1, 1, 3, 3, 1,-1,-3  // Level 7
};

const int CellularPotts::nbh_level[8] = { 0, 4, 8, 12, 20, 24, 28, 36 };
int CellularPotts::shuffleindex[9]={0,1,2,3,4,5,6,7,8};

extern Parameter par;


/** PRIVATE **/

using namespace std;
void CellularPotts::BaseInitialisation(vector<Cell> *cells) {
  CopyProb(par.T);
  cell=cells;
  if (par.copy_neighbourhood>=1 && par.copy_neighbourhood<=4)
    n_nb=nbh_level[par.copy_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4]).";

  if (par.adhesion_neighbourhood>=1 && par.adhesion_neighbourhood<=7)
    n_nb_adh=nbh_level[par.adhesion_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-7])";

  if (par.perimeter_neighbourhood>=1 && par.perimeter_neighbourhood<=7)
    n_nb_perim=nbh_level[par.perimeter_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-7])";
  
}

CellularPotts::CellularPotts(vector<Cell> *cells,
			     const int sx, const int sy ) {
  
  sigma=0;
  outside=0;
  frozen=false;
  thetime=0;
  zygote_area=0;

  
  BaseInitialisation(cells);
  sizex=sx;
  sizey=sy;
  rows = sizex-1;
  cols = sizey-1;

  old_nbhs = nullptr;

  AllocateSigma(sx,sy);
  
  
  // fill borders with special border state
  for (int x=0;x<sizex;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }
  
  if (par.copy_neighbourhood>=1 && par.copy_neighbourhood<=4)
    n_nb=nbh_level[par.copy_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";

  if (par.adhesion_neighbourhood>=1 && par.adhesion_neighbourhood<=7)
    n_nb_adh=nbh_level[par.adhesion_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-7])";

  if (par.perimeter_neighbourhood>=1 && par.perimeter_neighbourhood<=7)
    n_nb_perim=nbh_level[par.perimeter_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-7])";
}

CellularPotts::CellularPotts(void) {

  sigma=0;
  outside = 0;
  sizex=0; sizey=0;
  frozen=false;
  thetime=0;
  zygote_area=0;

  CopyProb(par.T);

  // fill borders with special border state
  for (int x=0;x<sizex;x++) 
  {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey;y++) 
  {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }
  if (par.copy_neighbourhood>=1 && par.copy_neighbourhood<=4)
    n_nb=nbh_level[par.copy_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";

  if (par.adhesion_neighbourhood>=1 && par.adhesion_neighbourhood<=4)
    n_nb_adh=nbh_level[par.adhesion_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";

  if (par.perimeter_neighbourhood>=1 && par.perimeter_neighbourhood<=4)
    n_nb_perim=nbh_level[par.perimeter_neighbourhood];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";
}

// destructor (virtual)
CellularPotts::~CellularPotts(void) {
  if (sigma) {
    free(sigma[0]);
    free(sigma);
    sigma=0;
  }

  if (outside)
  {
    free(outside[0]);
    free(outside);
    outside=0;
  }

  if (old_nbhs)
  {
    free(old_nbhs[0]);
    free(old_nbhs);
    old_nbhs=0;
  }

}

void CellularPotts::AllocateSigma(int sx, int sy) {
  
  sizex=sx; sizey=sy;
  
  sigma=(int **)malloc(sizex*sizeof(int *));
  if (sigma==NULL)
    MemoryWarning();
  
  sigma[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (sigma[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<sizex;i++) 
    sigma[i]=sigma[i-1]+sizey;}
  
  /* Clear CA plane */
   {for (int i=0;i<sizex*sizey;i++) 
     sigma[0][i]=0; }


  // do the same for outside plane
  outside=(int **)malloc(sizex*sizeof(int *));
  if (outside==NULL)
    MemoryWarning();
  
  outside[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (outside[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<sizex;i++) 
    outside[i]=outside[i-1]+sizey;}
  
  /* Clear CA plane */
   {for (int i=0;i<sizex*sizey;i++) 
     outside[0][i]=0; }

  {for (int i=1;i<sizex;i++) 
    outside[i]=outside[i-1]+sizey;}

  
}

void CellularPotts::IndexShuffle() {

  int i;
  int temp;
  int index1,index2;
  
  for (i=0;i<9;i++) {
    
    index1=RandomNumber(8, s_val);
    index2=RandomNumber(8, s_val);

    temp=shuffleindex[index1];
    shuffleindex[index1]=shuffleindex[index2];
    shuffleindex[index2]=temp;

  }
}

// Utility functions for grid iterations. 
int minu(int a, int b)
{
  return (a < b) ? a : b;
}

int nmin(int a, int b, int c)
{
    return minu(minu(a, b), c);
}

int nmax(int a, int b)
{
    return (a > b) ? a : b;
}


double sat(double x) {
  
  return x/(par.saturation*x+1.);
  //return x;

}

void CellularPotts::SetMediumArea()
{
  int totmed=0;
  for (int x = 1; x < sizex-1; ++x)
  {
    for (int y = 1; y < sizey-1; ++y)
    {
      if (!sigma[x][y])
      {
        ++totmed;
      }
    }
  }
  cell->at(0).SetTargetArea(totmed);
  cell->at(0).SetAreaToTarget();
}
  




// double CellularPotts::DeltaH(int x,int y, int sxyp, const int tsteps, PDE *PDEfield)       
// {
//   double DH = 0;
//   int i, sxy;
//   int neighsite;

//   /* Compute energydifference *IF* the copying were to occur */
//   sxy = sigma[x][y];

//   Cell& cell_sxy = (*cell)[sxy];
//   Cell& cell_sxyp = (*cell)[sxyp];
//   double Jen=0;


//   if (par.dynamic_sorting)
//   {
//     int type_sxy  = (sxy == 0)  ? 0 : cell_sxy.GetSortingType();
//     int type_sxyp = (sxyp == 0) ? 0 : cell_sxyp.GetSortingType();
//     for (i=1;i<=n_nb_adh;i++) 
//     {
//       int xp2,yp2;
//       xp2=x+nx[i]; yp2=y+ny[i];
//       if (par.periodic_boundaries) 
//       {
        
//         // since we are asynchronic, we cannot just copy the borders once 
//         // every MCS
        
//         if (xp2<=0)
//           xp2=sizex-2+xp2;
//         if (yp2<=0)
//           yp2=sizey-2+yp2;
//         if (xp2>=sizex-1)
//           xp2=xp2-sizex+2;
//         if (yp2>=sizey-1)
//           yp2=yp2-sizey+2;
      
//         neighsite=sigma[xp2][yp2];
    
//       } 
//       else 
//       {
//         if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1)
//           neighsite=-1;
//         else
//           neighsite=sigma[xp2][yp2];
//       }
      
//       if (neighsite==-1) 
//       { // border 
//         Jen += (sxyp==0?0:par.border_energy)-(sxy==0?0:par.border_energy);
//       } 
//       else 
//       {
//         int type_neigh = (neighsite == 0) ? 0 : (*cell)[neighsite].GetSortingType();
//         Jen += DynamicAdhesionDiff(type_sxyp, type_neigh) - DynamicAdhesionDiff(type_sxy, neighsite);
//       }
//     }
//   }
//   else if (par.make_synthetic)
//   {
//     for (i=1;i<=n_nb_adh;i++) 
//     {
//       int xp2,yp2;
//       xp2=x+nx[i]; yp2=y+ny[i];
//       if (par.periodic_boundaries) 
//       {
        
//         // since we are asynchronic, we cannot just copy the borders once 
//         // every MCS
        
//         if (xp2<=0)
//           xp2=sizex-2+xp2;
//         if (yp2<=0)
//           yp2=sizey-2+yp2;
//         if (xp2>=sizex-1)
//           xp2=xp2-sizex+2;
//         if (yp2>=sizey-1)
//           yp2=yp2-sizey+2;
      
//         neighsite=sigma[xp2][yp2];
    
//       } 
//       else 
//       {
//         if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1)
//           neighsite=-1;
//         else
//           neighsite=sigma[xp2][yp2];
//       }
      
//       if (neighsite==-1) 
//       { // border 
//         Jen += (sxyp==0?0:par.border_energy)-(sxy==0?0:par.border_energy);
//       } 
//       else 
//       {
//         Jen += (*cell)[sxyp].SyntheticEnergy((*cell)[neighsite]) - (*cell)[sxy].SyntheticEnergy((*cell)[neighsite]);
//       }
//     }
//   }
//   DH += Jen;// / (par.neigh_multiplier);

  
//   // lambda is determined by chemical 0
//   double lambda = (*cell)[sxy].get_lambda();
    
//   //cerr << "[" << lambda << "]";
//   if (par.medium_area_constraint)
//   {
//       DH += (int)((lambda * (2.+  2.  * (double) 
//               (  (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()
//               - (*cell)[sxy].Area() + (*cell)[sxy].TargetArea() )) ));
      
//       // cout << (*cell)[sxy].Area() << '\t' << (*cell)[sxy].TargetArea() << endl;
//   }
//   else
//   {
//     if ( sxyp == MEDIUM ) {
//       DH += (int)(lambda *  (1. - 2. *   
//               (double) ( (*cell)[sxy].Area() - (*cell)[sxy].TargetArea()) ));
//     }
//     else if ( sxy == MEDIUM ) {
//       DH += (int)((lambda * (1. + 2. *  
//               (double) ( (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()) )));
//     }
//     else
//       DH += (int)((lambda * (2.+  2.  * (double) 
//               (  (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()
//               - (*cell)[sxy].Area() + (*cell)[sxy].TargetArea() )) ));
//   }


//   /* Active motion term */
//   if (par.active_motion)
//   {
//     if ( sxyp == MEDIUM)
//     {
//       DH -= par.motility_strength * (*cell)[sxy].ActiveDotProduct_removed(x,y);
//       // cout << "active: " << par.motility_strength * (*cell)[sxy].ActiveDotProduct_removed(x,y) << endl;
//     }
//     else if (sxy == MEDIUM)
//     {
//       DH -= par.motility_strength * (*cell)[sxyp].ActiveDotProduct_added(x,y);

//     }
//     else
//     {
//       // cout << "dot product with cell: " << (*cell)[sxy].ActiveDotProduct_removed(x,y);
//       DH -= par.motility_strength * (*cell)[sxyp].ActiveDotProduct_added(x,y);
//       DH -= par.motility_strength * (*cell)[sxy].ActiveDotProduct_removed(x,y);
//     }
//   }

//   // gravity term for synthetic structures. I think this might be wrong doing medium this way?
//   if (par.add_gravity)
//   {
//     if ( sxyp == MEDIUM)
//     {
//       DH += (*cell)[sxy].Gravity();
//     }
//     else if (sxy == MEDIUM)
//     {
//       DH += (*cell)[sxyp].Gravity();
//     }
//     else
//     {
//       DH += (*cell)[sxy].Gravity();
//       DH += (*cell)[sxyp].Gravity(); 
//     }
//   }


  
//   if (par.H_perim) 
//   {
//     double DH_perimeter = 0;
//     if (sxyp == MEDIUM) 
//     {

//       DH_perimeter -=
//           (*cell)[sxy].getPerimConstraint() *
//           (DSQR((*cell)[sxy].Perimeter() - (*cell)[sxy].TargetPerimeter()) -
//           DSQR(GetNewPerimeterIfXYWereRemoved(sxy, x, y) -
//                 (*cell)[sxy].TargetPerimeter()));      
//     } 
//     else if (sxy == MEDIUM) 
//     {

//       DH_perimeter -=
//           (*cell)[sxyp].getPerimConstraint() *
//           (DSQR((*cell)[sxyp].Perimeter() - (*cell)[sxyp].TargetPerimeter()) -
//           DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y) -
//                 (*cell)[sxyp].TargetPerimeter()));      
//     }
//     // they're both cells
//     else 
//     {
//       DH_perimeter -=
//           (*cell)[sxyp].getPerimConstraint() *
//           ((DSQR((*cell)[sxyp].Perimeter() - (*cell)[sxyp].TargetPerimeter()) -
//             DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y) -
//                 (*cell)[sxyp].TargetPerimeter())));
//       DH_perimeter -=
//           (*cell)[sxy].getPerimConstraint() *
//           (DSQR((*cell)[sxy].Perimeter() - (*cell)[sxy].TargetPerimeter()) -
//             DSQR(GetNewPerimeterIfXYWereRemoved(sxy, x, y) -
//                 (*cell)[sxy].TargetPerimeter()));      
//     }
//     DH += DH_perimeter;
//   }
  
//   return DH;

// }



double CellularPotts::DeltaH(int x, int y, int sxyp, const int tsteps, const int* neighbor_spins)       
{
  double DH = 0;
  int i, sxy;
  int neighsite;

  /* Compute energy difference *IF* the copying were to occur */
  sxy = sigma[x][y];

  // Pre-fetch references to save standard vector lookups
  Cell& cell_sxy = (*cell)[sxy];
  Cell& cell_sxyp = (*cell)[sxyp];
  double Jen = 0;

  // ==========================================
  // ADHESION ENERGY CALCULATION
  // ==========================================
  if (par.dynamic_sorting)
  {
    int type_sxy  = (sxy == 0)  ? 0 : cell_sxy.GetSortingType();
    int type_sxyp = (sxyp == 0) ? 0 : cell_sxyp.GetSortingType();

    for (i = 1; i <= n_nb_adh; i++) 
    {
      neighsite = neighbor_spins[i];
      
      if (neighsite == -1) 
      { // out-of-bounds border 
        Jen += (sxyp == 0 ? 0 : par.border_energy) - (sxy == 0 ? 0 : par.border_energy);
      } 
      else 
      {
        int type_neigh = (neighsite == 0) ? 0 : (*cell)[neighsite].GetSortingType();
        
        // Only calculate adhesion if the neighbor belongs to a DIFFERENT cell.
        // If neighsite == sxyp (proposed state), the boundary is internal to the cell, so energy is 0.
        double E_final = (sxyp == neighsite) ? 0.0 : DynamicAdhesionDiff(type_sxyp, type_neigh);
        
        // If neighsite == sxy (initial state), the boundary was internal to the cell, so energy was 0.
        double E_initial = (sxy == neighsite) ? 0.0 : DynamicAdhesionDiff(type_sxy, type_neigh);

        Jen += (E_final - E_initial);
      }
    }
  }
  
  // NOTE: If you have other modes (par.sheet, par.melting_adhesion) re-add them 
  // here following the EXACT SAME pattern (just use neighsite = neighbor_spins[i])!

  DH += Jen; // / (par.neigh_multiplier);


  // ==========================================
  // AREA CONSTRAINT
  // ==========================================
    
  if (par.medium_area_constraint)
  {
      DH += (int)((par.lambda * (2. + 2. * (double) 
              (cell_sxyp.Area() - cell_sxyp.TargetArea()
             - cell_sxy.Area()  + cell_sxy.TargetArea()))));
  }
  else
  {
    if (sxyp == MEDIUM) {
      DH += (int)(par.lambda * (1. - 2. * (double) (cell_sxy.Area() - cell_sxy.TargetArea())));
    }
    else if (sxy == MEDIUM) {
      DH += (int)((par.lambda * (1. + 2. * (double) (cell_sxyp.Area() - cell_sxyp.TargetArea()))));
    }
    else {
      DH += (int)((par.lambda * (2. + 2. * (double) 
              (cell_sxyp.Area() - cell_sxyp.TargetArea()
             - cell_sxy.Area()  + cell_sxy.TargetArea()))));
    }
  }


  // ==========================================
  // ACTIVE MOTION TERM
  // ==========================================
  if (par.active_motion)
  {
    if (sxyp == MEDIUM)
    {
      DH -= par.motility_strength * cell_sxy.ActiveDotProduct_removed(x, y);
    }
    else if (sxy == MEDIUM)
    {
      DH -= par.motility_strength * cell_sxyp.ActiveDotProduct_added(x, y);
    }
    else
    {
      DH -= par.motility_strength * cell_sxyp.ActiveDotProduct_added(x, y);
      DH -= par.motility_strength * cell_sxy.ActiveDotProduct_removed(x, y);
    }
  }



  // ==========================================
  // PERIMETER CONSTRAINT
  // ==========================================
  if (par.H_perim) 
  {
    double DH_perimeter = 0;
    if (sxyp == MEDIUM) 
    {
      DH_perimeter -= par.lambda_perimeter *
          (DSQR(cell_sxy.Perimeter() - cell_sxy.TargetPerimeter()) -
           DSQR(GetNewPerimeterIfXYWereRemoved(sxy, x, y, neighbor_spins) - cell_sxy.TargetPerimeter()));      
    } 
    else if (sxy == MEDIUM) 
    {
      DH_perimeter -= par.lambda_perimeter *
          (DSQR(cell_sxyp.Perimeter() - cell_sxyp.TargetPerimeter()) -
           DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y, neighbor_spins) - cell_sxyp.TargetPerimeter()));      
    }
    else 
    {
      // they're both cells
      DH_perimeter -= par.lambda_perimeter *
          ((DSQR(cell_sxyp.Perimeter() - cell_sxyp.TargetPerimeter()) -
            DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y, neighbor_spins) - cell_sxyp.TargetPerimeter())));
            
      DH_perimeter -= par.lambda_perimeter *
          (DSQR(cell_sxy.Perimeter() - cell_sxy.TargetPerimeter()) -
           DSQR(GetNewPerimeterIfXYWereRemoved(sxy, x, y, neighbor_spins) - cell_sxy.TargetPerimeter()));      
    }
    DH += DH_perimeter;
  }
  
  return DH;
}







bool CellularPotts::Probability(int DH)
{
  if ( DH > BOLTZMANN-1 )
    return false;
  else if ( DH < 0 || RANDOM(s_val) < copyprob[DH] )
    return true; 
   return false; 
}



void CellularPotts::update_cell_velocities_MCS()
{

  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      c->update_velocity();
    }
  }
}

void CellularPotts::MeasureSinglePerimeter(int targetsigma)
{
  cell->at(targetsigma).SetPerimeter(0);
  
  for (int x = 1; x < sizex - 1; x++) 
  {
    for (int y = 1; y < sizey - 1; y++) 
    {
      if (sigma[x][y] == targetsigma) {
        for (int i = 1; i <= n_nb_perim; i++) {
          int xp2, yp2;
          xp2 = x + nx[i];
          yp2 = y + ny[i];
          if (par.periodic_boundaries) {
            if (xp2 <= 0)
              xp2 = sizex - 2 + xp2;
            if (yp2 <= 0)
              yp2 = sizey - 2 + yp2;
            if (xp2 >= sizex - 1)
              xp2 = xp2 - sizex + 2;
            if (yp2 >= sizey - 1)
              yp2 = yp2 - sizey + 2;
          }
          // did we find a border?
          if (sigma[xp2][yp2] != sigma[x][y]) {
            // add to the perimeter of the cell
            // (*cell)[sigma[x][y]].IncrementTargetPerimeter();
            (*cell)[sigma[x][y]].IncrementPerimeter();
          }
        }
      }
    }
  }
}

void CellularPotts::MeasureCellPerimeters() 
{


  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      c->SetPerimeter(0);
    }
  }

  for (int x = 1; x < sizex - 1; x++) 
  {
    for (int y = 1; y < sizey - 1; y++) {
      if (sigma[x][y] > 0) {
        for (int i = 1; i <= n_nb_perim; i++) {
          int xp2, yp2;
          xp2 = x + nx[i];
          yp2 = y + ny[i];
          if (par.periodic_boundaries) {
            if (xp2 <= 0)
              xp2 = sizex - 2 + xp2;
            if (yp2 <= 0)
              yp2 = sizey - 2 + yp2;
            if (xp2 >= sizex - 1)
              xp2 = xp2 - sizex + 2;
            if (yp2 >= sizey - 1)
              yp2 = yp2 - sizey + 2;
          }
          else if (xp2 >= sizex-1 || xp2 < 1 || yp2 >= sizey-1 || yp2 < 1)
          {
            continue;
          }
          // did we find a border?
          if (sigma[xp2][yp2] != sigma[x][y]) {
            // add to the perimeter of the cell
            // (*cell)[sigma[x][y]].IncrementTargetPerimeter();
            (*cell)[sigma[x][y]].IncrementPerimeter();
          }
        }
      }
    }
  }

  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      c->SetTargetPerimeter(par.ptarget_perimeter);
    }
  }  

  // for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) 
  // {
  //   if (c->AliveP())
  //   {
  //     cout << c->Perimeter() << endl;
  //   }
  // }


}

void CellularPotts::ConvertSpin(int x,int y,int kp)
{
  int tmpcell;
  if ( (tmpcell=sigma[x][y]) ) 
  { // if tmpcell is not MEDIUM
    (*cell)[tmpcell].DecrementArea();
    (*cell)[tmpcell].RemoveSiteFromMoments(x,y);
        
    if (!(*cell)[tmpcell].Area()) 
    {
      (*cell)[tmpcell].Apoptose();
      // cerr << "Cell " << tmpcell << " apoptosed\n";
    }
  }
  
  if ( (tmpcell = kp) ) 
  {// if tmpcell is not MEDIUM
    (*cell)[tmpcell].IncrementArea();
    (*cell)[tmpcell].AddSiteToMoments(x,y);
  }
  sigma[x][y] = kp;

}


void CellularPotts::ConvertSpinPerim(int x, int y, int kp, const int* neighbor_spins) 
{
  int tmpcell;

  if ((tmpcell = sigma[x][y])) { // if tmpcell is not MEDIUM
    (*cell)[tmpcell].DecrementArea();
    (*cell)[tmpcell].RemoveSiteFromMoments(x, y);
    
    // Pass neighbor_spins array here!
    (*cell)[tmpcell].SetPerimeter(
        GetNewPerimeterIfXYWereRemoved(tmpcell, x, y, neighbor_spins));
        
    if (!(*cell)[tmpcell].Area()) {
      (*cell)[tmpcell].Apoptose();
    }
  }

  if ((tmpcell = kp)) { // if tmpcell is not MEDIUM
    (*cell)[tmpcell].IncrementArea();
    (*cell)[tmpcell].AddSiteToMoments(x, y);
    
    // Pass neighbor_spins array here!
    (*cell)[tmpcell].SetPerimeter(
        GetNewPerimeterIfXYWereAdded(tmpcell, x, y, neighbor_spins));
  }
  sigma[x][y] = kp;
}




/** PUBLIC **/
int CellularPotts::CopyvProb(double DH,  double stiff) {

  double dd; 
  int s;
  s=(int)stiff;
  if (DH<=-s) 
    return 2;


  // if (par.IntegerHamiltonian)
  // {
  //   // if DH becomes extremely large, calculate probability on-the-fly
  //   if (DH+s > BOLTZMANN-1)
  //     dd=exp( -( (double)(DH+s)/par.T ));
  //   else
  //     dd=copyprob[DH+s]; 
  // }
  // else
  // {
  // we are slowing down sim by dong on the fly probabilities so that we can do doubles. 
  dd=exp( -( (double)(DH+s)/internal_T ));
  // }




  if (RANDOM(s_val)<dd) 
    return 1; 
  else 
    return 0;
} 

void CellularPotts::CopyProb(double T) {
  int i;
  internal_T = T;
  for ( i = 0; i < BOLTZMANN; i++ )
    copyprob[i] = exp( -( (double)(i)/internal_T ) );
}

void CellularPotts::FreezeAmoebae(void) 
{
  if (frozen) 
    frozen=FALSE;
  else
    frozen=TRUE;
}

// testing perimeter constraint here.



int CellularPotts::GetNewPerimeterIfXYWereAdded(int sxyp, int x, int y, const int* neighbor_spins) 
{
  int perim = (*cell)[sxyp].Perimeter();

  /* the cell with sigma sxyp wants to extend by adding lattice site (x, y). */
  for (int i = 1; i <= n_nb_perim; i++) {
    
    int neighsite = neighbor_spins[i];
    
    // Skip out-of-bounds just like your original code did
    if (neighsite == -1) {
      continue;
    }

    if (neighsite == sxyp) {
      perim--;
    } else {
      perim++;
    }
  }
  return perim;
}

int CellularPotts::GetNewPerimeterIfXYWereRemoved(int sxy, int x, int y, const int* neighbor_spins) {
  int perim = (*cell)[sxy].Perimeter();
  
  /* the cell with sigma sxy loses xy */
  for (int i = 1; i <= n_nb_perim; i++) {
    
    int neighsite = neighbor_spins[i];
    
    // -1 means it was an out-of-bounds border. Your original code skipped these 
    // with 'continue', so we do the exact same thing here.
    if (neighsite == -1) {
      continue;
    }

    if (neighsite == sxy) {
      perim++;
    } else {
      perim--;
    }
  }
  return perim;
}


void CellularPotts::GetNeighborsSafe(int x, int y, int* nbs) 
{
  // Clockwise offsets starting from top-left (-1, -1)
  const int cyc_nnx[8] = {-1,  0,  1,  1,  1,  0, -1, -1};
  const int cyc_nny[8] = {-1, -1, -1,  0,  1,  1,  1,  0};
  
  for (int i = 0; i < 8; i++) 
  {
    int tx = x + cyc_nnx[i];
    int ty = y + cyc_nny[i];

    if (par.periodic_boundaries) 
    {
      // Mathematically identical to your old periodic wrap
      if (tx <= 0) tx += sizex - 2;
      else if (tx >= sizex - 1) tx -= sizex - 2;
      
      if (ty <= 0) ty += sizey - 2;
      else if (ty >= sizey - 1) ty -= sizey - 2;
    }
    
    // If not periodic, it naturally reads the 0 / sizex-1 ghost cells exactly like the old code did
    nbs[i] = sigma[tx][ty];
  }
}



int CellularPotts::AmoebaeMove(long tsteps)
{ 
  int loop, p;
  thetime++;
  int SumDH = 0;
  
  if (frozen) 
    return 0;

  const int sx_inner = sizex - 2;
  const int sy_inner = sizey - 2;

  loop = sx_inner * sy_inner;

  // Use a fixed size (e.g., 64) to avoid Variable Length Array compiler errors. 
  // It is more than enough for up to 5th+ order neighborhoods.
  int present_states[64]; 
  int distinct_count = 0;

  // We need enough neighbors to satisfy BOTH Adhesion and Perimeter checks
  int max_nb = n_nb_adh;
  int neighbor_spins[64]; // Fixed size covers max_nb easily, allows 1-based indexing safely
 
  for (int i = 0; i < loop; i++) 
  {  
    // FAST random site selection (replaces expensive % and / modulo operations)
    int x = 1 + (int)(RANDOM(s_val) * sx_inner);
    int y = 1 + (int)(RANDOM(s_val) * sy_inner);
    
    int k = sigma[x][y];

    // =============================================================
    // 1. PRE-FETCH ALL REQUIRED NEIGHBORS EXACTLY ONCE
    // =============================================================
    for (int j = 1; j <= max_nb; j++) 
    {
      int tx = nx[j] + x;
      int ty = ny[j] + y;

      if (par.periodic_boundaries) {
          if (tx <= 0) tx += sx_inner;
          else if (tx >= sizex - 1) tx -= sx_inner;
          if (ty <= 0) ty += sy_inner;
          else if (ty >= sizey - 1) ty -= sy_inner;
          
          neighbor_spins[j] = sigma[tx][ty];
      } 
      else 
      {
          if (tx <= 0 || ty <= 0 || tx >= sizex - 1 || ty >= sizey - 1) {
              neighbor_spins[j] = -1; // -1 means out-of-bounds border
          } else {
              neighbor_spins[j] = sigma[tx][ty];
          }
      }
    }

    // =============================================================
    // 2. IDENTIFY TARGET STATES (using only copy neighborhood: n_nb)
    // =============================================================
    distinct_count = 0;
    for (int j = 1; j <= n_nb; j++)
    {
      int neighbor_val = neighbor_spins[j];
      
      if (neighbor_val == -1) 
          continue; // Skip out-of-bounds border states

      // Check if this state is already in our list
      bool seen = false;
      for (int u = 0; u < distinct_count; u++) {
          if (present_states[u] == neighbor_val) {
              seen = true;
              break;
          }
      }

      // If new unique state, add it
      if (!seen) 
      {
          present_states[distinct_count] = neighbor_val;
          distinct_count++;
      }
    }

    // 3. Randomly select a TARGET STATE from the unique list
    if (distinct_count == 0) 
      continue; // Should not happen unless isolated
    
    int rand_idx = (int)(distinct_count * RANDOM(s_val));
    int kp = present_states[rand_idx];

    // =============================================================
    // 4. CONNECTIVITY CHECK (Kept separate because it relies on 
    //    strict Clockwise 8-neighborhood ordering for Moore Check)
    // =============================================================
    int nbs[8];
    // FAST PATH: 95%+ of the time, the cell is entirely inside the grid
    if (!par.periodic_boundaries || (x > 1 && y > 1 && x < sizex - 2 && y < sizey - 2)) 
    {
        nbs[0] = sigma[x-1][y-1];
        nbs[1] = sigma[x  ][y-1];
        nbs[2] = sigma[x+1][y-1];
        nbs[3] = sigma[x+1][y  ];
        nbs[4] = sigma[x+1][y+1];
        nbs[5] = sigma[x  ][y+1];
        nbs[6] = sigma[x-1][y+1];
        nbs[7] = sigma[x-1][y  ];
    }
    else
    {
      GetNeighborsSafe(x, y, nbs);
    }

    bool check1 = (k == 0)  ? true : IsLocallyConnected(nbs, k);
    bool check2 = (kp == 0) ? true : IsLocallyConnected(nbs, kp);

    // =============================================================
    // 5. ATTEMPT COPY
    // =============================================================
    if (kp != -1 && check1 == true && check2 == true && k != kp) 
    {  
      int H_diss = 0;

      // PASS IN THE PRE-CALCULATED ARRAY HERE!
      double D_H = DeltaH(x, y, kp, tsteps, neighbor_spins);

      if ((p = CopyvProb(D_H, H_diss)) > 0) 
      {
        ++par.tmpcounter;
        if (par.H_perim) {
          // PASS IN THE PRE-CALCULATED ARRAY HERE TOO!
          ConvertSpinPerim(x, y, kp, neighbor_spins);
        } else {
          ConvertSpin(x, y, kp);
        }
      }
    } 
  }
  return SumDH;
}

// int CellularPotts::AmoebaeMoveLegacy(long tsteps, PDE *PDEfield)
// {
//   int loop,p;
//   //int updated=0;
//   thetime++;
//   int SumDH=0;
  
//   if (frozen) 
//     return 0;

//   loop=(sizex-2)*(sizey-2);
 
//   for (int i=0;i<loop;i++) 
//   {  
//     // take a random site
//     int xy = (int)(RANDOM(s_val)*(sizex-2)*(sizey-2));
//     int x = xy%(sizex-2)+1;
//     int y = xy/(sizex-2)+1; 
    
//     // take a random neighbour
//     int xyp=(int)(n_nb*RANDOM(s_val)+1);
//     int xp = nx[xyp]+x;
//     int yp = ny[xyp]+y;
    
//     int k=sigma[x][y];
    
//     int kp;
//     if (par.periodic_boundaries) 
//     {
//       // since we are asynchronic, we cannot just copy the borders once 
//       // every MCS
//       if (xp<=0)
// 	      xp=sizex-2+xp;
//       if (yp<=0)
// 	      yp=sizey-2+yp;
//       if (xp>=sizex-1)
// 	      xp=xp-sizex+2;
//       if (yp>=sizey-1)
// 	      yp=yp-sizey+2;
      
//       kp=sigma[xp][yp];
      
//     } 
//     else 
//     {
//       if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1)
// 	      kp=-1;
//       else
// 	      kp=sigma[xp][yp];
//     }
//     // int type1 = (*cell)[sigma[xp][yp]].GetPhenotype();    
//     // int type2 = (*cell)[sigma[xp][yp]].GetPhenotype();    

//     // test for border state (relevant only if we do not use 
//     // periodic boundaries)
//     if (kp!=-1) 
//     {  
//       // Don't even think of copying the special border state into you!
    
//       if ( k  != kp ) 
//       {
//         /* Try to copy if sites do not belong to the same cell */
//         // connectivity dissipation:
//         int H_diss=0;
//         if (!ConnectivityPreservedP(x,y)) 
//           H_diss=par.conn_diss;
        
//         double D_H=DeltaH(x,y,kp, tsteps, PDEfield);
        
//         // dH_tally += D_H;
//         // if ((type1 > par.mintype && type1 < par.maxtype) || (type2 > par.mintype && type2 < par.maxtype))
//         //   cout << D_H << endl;
//         // bool is_med_attempt = false;
//         // if (sigma[x][y] == 0 && (*cell)[sigma[xp][yp]].GetPhase() == true || sigma[xp][yp] == 0 && (*cell)[sigma[x][y]].GetPhase() == true)
//         // {
//         //   is_med_attempt = true;
//         //   ++medp_count;
//         // }
//         if ((p=CopyvProb(D_H,H_diss))>0) 
//         {
//           if (par.H_perim)
//             ConvertSpinPerim( x,y,kp );
//           else
//           {
//             ConvertSpin( x,y,kp );
//           }  
//         }
//         //   if (par.recordcopies)
//         //   {
//         //     if ((type1 > par.mintype && type1 < par.maxtype) || (type2 > par.mintype && type2 < par.maxtype))
//         //     {
//         //       ++flip_true;
//         //       SumDH+=D_H;
//         //       dH_neg+=D_H;
//         //     }
//         //   }
          
//         // }
//         // else
//         // {
//         //   if (par.recordcopies)
//         //     if ((type1 > par.mintype && type1 < par.maxtype) || (type2 > par.mintype && type2 < par.maxtype))
//         //     {
//         //       ++flip_false;
//         //     }
//         // }
//         // if (Probability(D_H)) 
//         // {
//         //   ConvertSpin( x,y,xp,yp );
//         //   SumDH+=D_H;
//         // }
//       }
//     } 
//   }
//   return SumDH;
  
// }




void CellularPotts::StartDynamicAdhesion()
{
  prev_nbs = GetNeighbourArray();

  // note that we should make it +1 bigger because of medium
  int ncells = (*cell).size();
  int total_elements = (ncells * (ncells + 1)) / 2;
  vector<double> onevec(ncells, par.init_J);
  for (int i = 0; i < ncells; ++i)
  {
    DynamicAdhesions.push_back(onevec);
  }
  // DynamicAdhesions.resize(total_elements, par.init_J);
  DynamicMeeting.resize(total_elements, 0);

}

void CellularPotts::AddtoMeeting(int i, int j)
{
  int row = std::max(i, j);
  int col = std::min(i, j);
  DynamicMeeting[(row * (row + 1)) / 2 + col] += par.timeadd_ifmet;
}

void CellularPotts::SnapMeeting(int i, int j)
{
  int row = std::max(i, j);
  int col = std::min(i, j);
  DynamicMeeting[(row * (row + 1)) / 2 + col] = 0;
}

int CellularPotts::GetMeeting(int i, int j)
{
  int row = std::max(i, j);
  int col = std::min(i, j);
  return DynamicMeeting[(row * (row + 1)) / 2 + col];
}


void CellularPotts::UpdateDynamicAdhesion()
{
  // for starters we will just do neighbour or not to increase/snap adhesion.
  // We can add in a cell joining threshold later

  // first, we find neighbours

  int **nbs = GetNeighbourArray();
  
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      int sig = c->Sigma();
      int j=0;
      while(nbs[sig][j] != EMPTY && nbs[sig][j] > 0)
      {
        int nbh = nbs[sig][j];
        if (nbh > 0)
        {
          if (c->GetSortingType() == (*cell)[nbh].GetSortingType())
          {
            AddtoMeeting(sig, nbh);
          }
        }
        // cout << nbh << '\t';
        ++j;
      }         
      // cout << endl;
      // SNAP condition
      int k = 0; 
      while (prev_nbs[sig][k] != EMPTY && prev_nbs[sig][k] > 0)
      {
        int prev_nbh = prev_nbs[sig][k];
        if (prev_nbh > 0)
        {
          bool found_in_current = false;
          int i = 0;
          // cout << prev_nbh << '\t';
          // Search for prev_nbh in the current nbs array
          while (nbs[sig][i] != EMPTY)
          {
            if (nbs[sig][i] == prev_nbh)
            {
                found_in_current = true;
                break; // Stop searching once we find it
            }
            ++i;
          }
          // Handle the result
          if (!found_in_current)
          {
            // cout << GetMeeting(sig, prev_nbh) << '\t' << sig << '\t' << prev_nbh << '\t' << (*cell)[sig].GetSortingType() << '\t' << (*cell)[prev_nbh].GetSortingType() << endl;
            SnapMeeting(sig, prev_nbh);
          }
        }
        ++k; // Move to the next previous neighbor
      }
      // cout << endl;
    }
  }
  int vecsize = DynamicAdhesions.size();
  for (int i = 1; i < vecsize; ++i)
  {
    bool ctype = (*cell)[i].GetSortingType();
    if (ctype == true)
    {
      if (par.Ahascortex)
      {
        for (int j = 1; j<vecsize; ++j)
        {
          int cortex_time = GetMeeting(i, j);
          DynamicAdhesions[i][j] = par.AdynJmin + ( par.AdynJmax -  par.AdynJmin ) * exp(-par.timescaler * pow(double(cortex_time), 2 ) ); 
          // if (cortex_time > 100)
            // cout << cortex_time << '\t' << DynamicAdhesions[i][j] << endl;
        }
      }
      else
      {
        for (int j = 1; j<vecsize; ++j)
        {
          DynamicAdhesions[i][j] = par.AstaticJ;
        }
      }
    }
    else
    {
      if (par.Bhascortex)
      {
        for (int j = 1; j<vecsize; ++j)
        {
          int cortex_time = GetMeeting(i, j);
          DynamicAdhesions[i][j] = par.BdynJmin + ( par.BdynJmax -  par.BdynJmin ) * exp(-par.timescaler * pow(cortex_time, 2 ) ); 
        }
      }
      else
      {
        for (int j = 1; j<vecsize; ++j)
        {
          DynamicAdhesions[i][j] = par.BstaticJ;
        }
      }
    }
  }


  free(prev_nbs[0]);
  free(prev_nbs);
  prev_nbs=0;

  // cout << "neighbours swapped and total: " << n_swapped_neighbours << '\t' << counter << endl;
  prev_nbs=(int **)malloc((cell->size()+1)*sizeof(int *));
  if (prev_nbs==NULL) 
    MemoryWarning();
  
  prev_nbs[0]=(int *)malloc((cell->size()+1)*(cell->size()+1)*sizeof(int));
  if (prev_nbs[0]==NULL)
    MemoryWarning();
  
  for (int i=1;i<(int)cell->size()+1;i++)
    prev_nbs[i]=prev_nbs[i-1]+(cell->size()+1);
  
  for (int i=0;i<((int)cell->size()+1)*((int)cell->size()+1);i++)
    prev_nbs[0][i]=nbs[0][i]; 

  free(nbs[0]);
  free(nbs);

}


void CellularPotts::SetSortingTypesRandomly()
{

  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {

      double rand = RANDOM(s_val);
      if (rand > par.startingAproportion)
      {
        c->SetSortingType(true);
      }
      else
      {
        c->SetSortingType(false);
      }
    }
  }

}

double CellularPotts::MeasureDomainSizeR(std::vector<double>* out_C_r)
{
    const int sx_inner = sizex - 2;
    const int sy_inner = sizey - 2;
    // we dont really need the correlation function for values larger than 400.
    int max_r = std::min(std::min(sx_inner, sy_inner) / 2  , 400);
    if (max_r < 1) return 0.0;
    
    // ------------------------------------------------------------------------
    // 1. PRECOMPUTE TAU GRID (Eliminates pointer chasing)
    // ------------------------------------------------------------------------
    std::vector<std::vector<double>> tau_grid(sizex, std::vector<double>(sizey, 0.0));
    std::vector<std::vector<bool>> valid_grid(sizex, std::vector<bool>(sizey, false));
    
    double sum_tau = 0.0;
    long long count_p = 0;
    
    for (int x = 1; x <= sx_inner; x++) {
        for (int y = 1; y <= sy_inner; y++) {
            int id = sigma[x][y];
            if (id != 0) {
                double t = (*cell)[id].GetSortingType() ? 1.0 : -1.0;
                tau_grid[x][y] = t;
                valid_grid[x][y] = true;
                sum_tau += t;
                count_p++;
            }
        }
    }
    
    if (count_p == 0) return 0.0;
    double mean_tau_sq = (sum_tau / count_p) * (sum_tau / count_p);

    // Setup output array if requested
    if (out_C_r) {
        out_C_r->assign(max_r + 1, 0.0);
        // C(0) is perfectly correlated with itself: (1.0 * 1.0) - mean_tau_sq
        (*out_C_r)[0] = 1.0 - mean_tau_sq; 
    }

    // ------------------------------------------------------------------------
    // 2. PRECOMPUTE RINGS
    // ------------------------------------------------------------------------
    std::vector<std::vector<std::pair<int, int>>> rings(max_r + 1);
    for (int dx = -max_r; dx <= max_r; dx++) {
        for (int dy = -max_r; dy <= max_r; dy++) {
            int r_sq = dx*dx + dy*dy;
            if (r_sq == 0 || r_sq > max_r*max_r) continue;
            int r = (int)std::round(std::sqrt(r_sq));
            if (r <= max_r) {
                rings[r].push_back({dx, dy});
            }
        }
    }

    // ------------------------------------------------------------------------
    // 3. CACHE-OPTIMIZED RING EVALUATION
    // ------------------------------------------------------------------------
    int stride = 4; 
    double prev_C = 0.0;
    bool has_prev = false;
    
    double R_t = -1.0;
    bool found_Rt = false;

    for (int r = 1; r <= max_r; r++) 
    {
        if (rings[r].empty()) continue;

        double current_sum = 0.0;
        long long current_count = 0;

        // INVERTED LOOPS: Offset -> X -> Y for maximum CPU Cache speed
        for (const auto& offset : rings[r]) 
        {
            int dx = offset.first;
            int dy = offset.second;
            
            for (int x = 1; x <= sx_inner; x += stride) 
            {
                int tx = x + dx;
                
                // X-boundary evaluated OUTSIDE the y-loop!
                if (par.periodic_boundaries) {
                    if (tx <= 0) tx += sx_inner;
                    else if (tx >= sizex - 1) tx -= sx_inner;
                } else {
                    if (tx <= 0 || tx >= sizex - 1) continue; 
                }
                
                // Y-loop accesses sequential memory in C++ vectors
                for (int y = 1; y <= sy_inner; y += stride) 
                {
                    if (!valid_grid[x][y]) continue;

                    int ty = y + dy;
                    if (par.periodic_boundaries) {
                        if (ty <= 0) ty += sy_inner;
                        else if (ty >= sizey - 1) ty -= sy_inner;
                    } else {
                        if (ty <= 0 || ty >= sizey - 1) continue; 
                    }
                    
                    if (valid_grid[tx][ty]) {
                        current_sum += tau_grid[x][y] * tau_grid[tx][ty];
                        current_count++;
                    }
                }
            }
        }

        if (current_count > 0) 
        {
            double C_r = (current_sum / current_count) - mean_tau_sq;
            
            if (out_C_r) {
                (*out_C_r)[r] = C_r;
            }

            // Check for zero-crossing
            if (!found_Rt && has_prev && prev_C >= 0.0 && C_r < 0.0) 
            {
                double slope = C_r - prev_C;
                R_t = (slope != 0.0) ? ((r - 1) - (prev_C / slope)) : (double)r;
                found_Rt = true;
                
                // EARLY EXIT: If they didn't ask for the full curve, stop calculating!
                if (!out_C_r) {
                    return R_t;
                }
            }
            prev_C = C_r;
            has_prev = true;
        }
    }
    
    return R_t;
}





//! Check if the set of neighbors with value 'check_val' forms a single connected component.
//! This ensures that removing a pixel (Candidate) doesn't split a cell, 
//! and adding a pixel (Target) doesn't create a handle/hole.

bool CellularPotts::IsLocallyConnected(int *nbs, int check_val) 
{
  
    int n_borders = 0;
    
    for (int i = 0; i < 8; i++) 
    {
        int s_nb = nbs[i];
        int s_next_nb = nbs[(i + 1) & 7]; // Bitwise AND wraps 7->0 instantly
        
        if ((s_nb == check_val || s_next_nb == check_val) && (s_nb != s_next_nb)) 
        {
            n_borders++;
            // EARLY EXIT: Don't waste CPU cycles checking the rest
            if (n_borders > 2) 
              return false; 
        }
    }
    return true; // if we survived the loop, it's connected
  
  // OLD CODE BELOW!!
  // Neighbor offsets for the 8-neighbor cycle (with overlap for i and i+1)
  // const int cyc_nnx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1 };
  // const int cyc_nny[10] = {0, -1,-1,-1, 0, 1, 1,  1,  0, -1 };

  // // int n_borders = 0;

  // for (int i = 1; i <= 8; i++) {
    
  //   // Calculate raw coordinates for the current neighbor (i) 
  //   // and the next neighbor in the cycle (i+1)
  //   int nx1 = x + cyc_nnx[i];
  //   int ny1 = y + cyc_nny[i];
    
  //   int nx2 = x + cyc_nnx[i+1];
  //   int ny2 = y + cyc_nny[i+1];

  //   // Apply Periodic Boundaries if enabled
  //   // Logic adapted strictly from GetNewPerimeterIfXYWereRemoved
  //   if (par.periodic_boundaries) {
      
  //     // Wrap current neighbor (nx1, ny1)
  //     if (nx1 <= 0) nx1 = sizex - 2 + nx1;
  //     if (ny1 <= 0) ny1 = sizey - 2 + ny1;
  //     if (nx1 >= sizex - 1) nx1 = nx1 - sizex + 2;
  //     if (ny1 >= sizey - 1) ny1 = ny1 - sizey + 2;

  //     // Wrap next neighbor (nx2, ny2)
  //     if (nx2 <= 0) nx2 = sizex - 2 + nx2;
  //     if (ny2 <= 0) ny2 = sizey - 2 + ny2;
  //     if (nx2 >= sizex - 1) nx2 = nx2 - sizex + 2;
  //     if (ny2 >= sizey - 1) ny2 = ny2 - sizey + 2;
  //   }

  //   int s_nb = sigma[nx1][ny1];
  //   int s_next_nb = sigma[nx2][ny2];
    
  //   if ((s_nb == check_val || s_next_nb == check_val) && (s_nb != s_next_nb)) 
  //   {
  //     // check whether s_nb is adjacent to non-identical site, count it
  //     n_borders++;
  //   }
  // }

  // // Connectivity check: In a locally connected grid on a square lattice, 
  // // there should be no more than 2 transitions entering/leaving the cell region.
  // if (n_borders > 2) 
  // {
  //   return false;
  // }
  // else 
  // {
  //   return true;
  // }
}


// Predicate returns true when connectivity is locally preserved
// if the value of the central site would be changed
bool CellularPotts::ConnectivityPreservedP(int x, int y) {
  
  // Use local nx and ny in a cyclic order (starts at upper left corner)
  // first site is repeated, for easier looping
  const int cyc_nx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1 };
  const int cyc_ny[10] = {0, -1,-1,-1, 0, 1, 1,  1,  0, -1 };
  
  int sxy=sigma[x][y]; // the central site
  if (sxy==0) return true;

  int n_borders=0; // to count the amount of sites in state sxy bordering a site !=sxy

  int stackp=-1;
  bool one_of_neighbors_medium=false;
  
  for (int i=1;i<=8;i++) {
    
    int s_nb=sigma[x+cyc_nx[i]][y+cyc_ny[i]];
    int s_next_nb=sigma[x+cyc_nx[i+1]][y+cyc_ny[i+1]];
    
    if ((s_nb==sxy || s_next_nb==sxy) && (s_nb!=s_next_nb)) {
      
      // check whether s_nb is adjacent to non-identical site,
      // count it
      n_borders++;
    }
    int j;
    bool on_stack_p=false;
    
    // we need the next heuristic to prevent stalling at
    // cell-cell borders
    // do not enforce constraint at two cell interface(no medium)
    if (s_nb) 
    {
      for (j=stackp;j>=0;j--) 
      {
        if (s_nb==stack[j]) 
        {
          on_stack_p=true;
          break;
        }
      }
      if (!on_stack_p) 
      {
	      if (stackp>6) 
        {
	        cerr << "Stack overflow, stackp=" << stackp << "\n";
	      }
	      stack[++stackp]=s_nb;
      }
    }
    else 
    {
      one_of_neighbors_medium=true;
    }
  }
  
  // number of different neighbours is stackp+1;
  if (n_borders>2 && ( (stackp+1)>2 || one_of_neighbors_medium) ) 
  {
    return false;
  }
  else 
    return true;

}



/** A simple method to plot all sigma's in window
    without the black lines */
void CellularPotts::PlotSigma(Graphics *g, int mag) {
  
  for (int x=0;x<sizex;x++) 
    for (int y=0;y<sizey;y++) {
      for (int xm=0;xm<mag;xm++)
	for (int ym=0;ym<mag;ym++)
      g->Point( sigma[x][y], mag*x+xm, mag*y+ym);
  }
  
}





int **CellularPotts::GetNeighbourArray()
{
  int i, j, q;
  int num_cells = (int)cell->size();
  
  /* Allocate neighbour matrix */
  int **neighbours = (int **)malloc((num_cells + 1) * sizeof(int *));
  if (neighbours == NULL) {
    MemoryWarning();
    return NULL;
  }
  
  neighbours[0] = (int *)malloc((num_cells + 1) * (num_cells + 1) * sizeof(int));
  if (neighbours[0] == NULL) {
    free(neighbours);
    MemoryWarning();
    return NULL;
  }
 
  for (i = 1; i < num_cells + 1; i++) {
    neighbours[i] = neighbours[i - 1] + (num_cells + 1);
  }
  
  /* Clear this matrix */
  for (i = 0; i < (num_cells + 1) * (num_cells + 1); i++) {
    neighbours[0][i] = EMPTY;  
  }

  /* Scan grid and find neighbors with periodic boundary conditions */
  for (i = 1; i < sizex-1; i++) {
    for (j = 1; j < sizey-1; j++) {
      
      int current = sigma[i][j];
      int next_i = (i + 1) % (sizex-2);
      int next_j = (j + 1) % (sizey-2);
      
      int right = sigma[next_i][j];
      int bottom = sigma[i][next_j];
      
      /* Compare with right neighbor */
      if (current != right) {
        if (current > 0) {
          for (q = 0; q < num_cells; q++) {
            if (neighbours[current][q] == EMPTY) { 
              neighbours[current][q] = right;  
              break;
            } else if (neighbours[current][q] == right) {
              break;
            }
          }
        }
        if (right > 0) {
          for (q = 0; q < num_cells; q++) {
            if (neighbours[right][q] == EMPTY) { 
              neighbours[right][q] = current; 
              break;
            } else if (neighbours[right][q] == current) {
              break;
            }
          }
        }
      }
      
      /* Compare with bottom neighbor */
      if (current != bottom) {
        if (current > 0) {
          for (q = 0; q < num_cells; q++) {
            if (neighbours[current][q] == EMPTY) { 
              neighbours[current][q] = bottom;  
              break; 
            } else if (neighbours[current][q] == bottom) {
              break;
            }
          }
        }
        if (bottom > 0) {
          for (q = 0; q < num_cells; q++) {
            if (neighbours[bottom][q] == EMPTY) { 
              neighbours[bottom][q] = current; 
              break;
            } else if (neighbours[bottom][q] == current) {
              break;
            }
          }
        }
      }
    }
  }
  
  return neighbours;
}



int **CellularPotts::SearchNandPlot(Graphics *g, bool get_neighbours)
{
  int i, j,q;
  int **neighbours=0;
  
  
  /* Allocate neighbour matrix */
  if (get_neighbours) {
    neighbours=(int **)malloc((cell->size()+1)*sizeof(int *));
    if (neighbours==NULL) 
      MemoryWarning();
    
    neighbours[0]=(int *)malloc((cell->size()+1)*(cell->size()+1)*sizeof(int));
    if (neighbours[0]==NULL)
      MemoryWarning();
   
    for (i=1;i<(int)cell->size()+1;i++)
      neighbours[i]=neighbours[i-1]+(cell->size()+1);
    
    /* Clear this matrix */
    for (i=0;i<((int)cell->size()+1)*((int)cell->size()+1);i++)
      neighbours[0][i]=EMPTY;  
  }

  for ( i = 0; i < sizex-1; i++ )
    for ( j = 0; j < sizey-1; j++ ) {
      

      int colour;
      if (sigma[i][j]<=0) {
	colour=0;
      } else {
	colour = (*cell)[sigma[i][j]].Colour();
      }
      
      if (g && sigma[i][j]>0)  /* if draw */ 
        g->Point( colour, 2*i, 2*j);
      
      if ( sigma[i][j] != sigma[i+1][j] )  /* if cellborder */ /* etc. etc. */
	{
	  if (g) 
	    g->Point( 1, 2*i+1, 2*j );
	  if (get_neighbours) {
	    if (sigma[i][j]>0) {
	      for (q=0;q<(int)cell->size();q++)
		if (neighbours[sigma[i][j]][q]==EMPTY) { 
		  neighbours[sigma[i][j]][q]=sigma[i+1][j];  
		  break;
		}
		else
		  if (neighbours[sigma[i][j]][q]==sigma[i+1][j]) 
		    break;
	    }
	    if (sigma[i+1][j]>0) {
	      for (q=0;q<(int)cell->size();q++)
		if (neighbours[sigma[i+1][j]][q]==EMPTY) { 
		  neighbours[sigma[i+1][j]][q]=sigma[i][j]; 
		  break;
		}
		else
		  if (neighbours[sigma[i+1][j]][q]==sigma[i][j]) 
		    break;
	    }
	  }
	} 
      else
        if (g && sigma[i][j]>0) 
          g->Point( colour, 2*i+1, 2*j );
      
      
      if ( sigma[i][j] != sigma[i][j+1] ) {
	
        if (g) 
	  g->Point( 1, 2*i, 2*j+1 );
	
	if (get_neighbours) {
	  if (sigma[i][j]>0) {
	    for (q=0;q<(int)cell->size();q++)
	      if (neighbours[sigma[i][j]][q]==EMPTY) { 
		neighbours[sigma[i][j]][q]=sigma[i][j+1];  
		break; 
	      }
	      else
		if (neighbours[sigma[i][j]][q]==sigma[i][j+1]) 
		  break;
	  }
	  
	  if (sigma[i][j+1]>0) {
	    
	    for (q=0;q<(int)cell->size();q++)
	      if (neighbours[sigma[i][j+1]][q]==EMPTY) { 
		neighbours[sigma[i][j+1]][q]=sigma[i][j]; 
		break;
	      }
	      else
		if (neighbours[sigma[i][j+1]][q]==sigma[i][j]) 
		  break;
	  }
	}
      } 
      else
        if (g && sigma[i][j]>0) 
          g->Point( colour, 2*i, 2*j+1 );
      
      /* Cells that touch eachother's corners are NO neighbours */ 
      
      if (sigma[i][j]!=sigma[i+1][j+1] 
	  || sigma[i+1][j]!=sigma[i][j+1] ) { 
        if (g) 
          g->Point( 1, 2*i+1, 2*j+1 ); 
      }
      else
        if (g && sigma[i][j]>0) 
          g->Point( colour, 2*i+1, 2*j+1 );
    }
  
  if (get_neighbours)
    return neighbours;
  else 
    return 0;

}


bool containsTargetVector(const vector<int>& target, const vector<std::vector<int>>& vec) 
{
  int tcount = target.size();
  for (const auto& subVec : vec) 
  {
    std::unordered_set<int> subVecSet(subVec.begin(), subVec.end());  // Convert subVec to a set for easy lookup
    bool containsAll = true;
    int counter=0;
    // Check if all elements in the target are in the current sub-vector
    for (int num : target) 
    {
      if (subVecSet.find(num) != subVecSet.end()) 
      {
        ++counter;
      }
    }
    if (counter == tcount)
    {
      return true;
    }
  }
  return false;  // No such vector found
}

double CellularPotts::CalculateABBoundaryLength()
{
    int total_ab_boundary = 0;
    int total_boundary = 0; // Track the total cell-cell boundary length

    // Iterate through the entire grid
    for (int i = 0; i < sizex; i++) {
        for (int j = 0; j < sizey; j++) {
            
            int current_id = sigma[i][j];
            
            // Skip if the current pixel is empty medium (<= 0)
            if (current_id <= 0) continue; 
            
            // Get the type of the current cell (false = A, true = B)
            bool current_type = (*cell)[current_id].GetSortingType();

            // 1. Check the RIGHT neighbor
            if (i + 1 < sizex) {
                int right_id = sigma[i + 1][j];
                
                // Check if right pixel is a valid cell and is a different cell
                if (right_id > 0 && current_id != right_id) {
                    total_boundary++; // It's a cell-cell boundary
                    
                    bool right_type = (*cell)[right_id].GetSortingType();
                    
                    // If the types are different (one is A, one is B)
                    if (current_type != right_type) {
                        total_ab_boundary++;
                    }
                }
            }

            // 2. Check the BOTTOM neighbor
            if (j + 1 < sizey) {
                int bottom_id = sigma[i][j + 1];
                
                // Check if bottom pixel is a valid cell and is a different cell
                if (bottom_id > 0 && current_id != bottom_id) {
                    total_boundary++; // It's a cell-cell boundary
                    
                    bool bottom_type = (*cell)[bottom_id].GetSortingType();
                    
                    // If the types are different (one is A, one is B)
                    if (current_type != bottom_type) {
                        total_ab_boundary++;
                    }
                }
            }
            
        }
    }

    // Prevent division by zero if there are no cell boundaries at all
    if (total_boundary == 0) {
        return 0.0;
    }

    // Return the relative boundary length as a double between 0.0 and 1.0
    return static_cast<double>(total_ab_boundary) / total_boundary;
}



vector<vector<int>> CellularPotts::SearchNforVertices()
{
  int x, y,q;
  int neighsite;
  vector<vector<int>> neighbours;
  int bound = 1;
  if (!par.periodic_boundaries)
    bound = 3;

  for ( x = bound; x < sizex-bound; x++ )
    for ( y = bound; y < sizey-bound; y++ ) 
    {
      int curcell=sigma[x][y];
      vector<int> tempn{};
      tempn.push_back(curcell);
      for (int i=1;i<=n_nb;i++) 
      {
        int xp2,yp2;
        xp2=x+nx[i]; yp2=y+ny[i];
        if (par.periodic_boundaries)
        {

          if (xp2<=0)
            xp2=sizex-2+xp2;
          if (yp2<=0)
            yp2=sizey-2+yp2;
          if (xp2>=sizex-1)
            xp2=xp2-sizex+2;
          if (yp2>=sizey-1)
            yp2=yp2-sizey+2;
        
          neighsite=sigma[xp2][yp2];
          
          // cout << "WHAT THE.." << neighsite << endl;
      
        } 
        else 
        {
          if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1)
            neighsite=-1;
          else
            neighsite=sigma[xp2][yp2];
        }
        if (curcell != neighsite)
        {
          if (find(tempn.begin(), tempn.end(), neighsite) == tempn.end())
          {
            tempn.push_back(neighsite);
          }
        }
      } 
      // must be at least three cells for it to be a vertex
      if (tempn.size() > 2)
      {
        bool is_in = containsTargetVector(tempn, neighbours);
        if (!is_in)
        {
          neighbours.push_back(tempn);
        }
      }
    }

  return neighbours;
}

#include <tuple>

// // Function to calculate the cross product of vectors ab and ac
// double crossProduct(double ax, double ay, double bx, double by, double cx, double cy) {
//     return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
// }

// // Function to reorder points in counterclockwise order
// void reorderCounterClockwise(double &ax, double &ay, double &bx, double &by, double &cx, double &cy) {
//     // Calculate the cross product of vectors ab and ac
//     double cross = crossProduct(ax, ay, bx, by, cx, cy);

//     // If the cross product is negative, the points are in clockwise order, so we swap b and c
//     if (cross < 0) {
//         std::swap(bx, cx);
//         std::swap(by, cy);
//     }
//     // If cross product is zero, the points are collinear, and no reordering is needed.
// }

struct cellPoint {
    double x, y;
};

// Function to calculate the cross product of two vectors (p0p1) and (p0p2)
double crossProduct(const cellPoint &p0, const cellPoint &p1, const cellPoint &p2) {
    return (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
}

// Function to find the centroid of a set of points
cellPoint calculateCentroid(const std::vector<cellPoint> &points) {
    double cx = 0, cy = 0;
    for (const auto &point : points) {
        cx += point.x;
        cy += point.y;
    }
    return { cx / points.size(), cy / points.size() };
}

// Function to reorder points in counterclockwise order
void reorderCounterClockwise(std::vector<cellPoint> &points) {
    if (points.size() <= 1) return; // No need to reorder if there's 1 or fewer points

    // Calculate the centroid of the points as the reference point
    cellPoint centroid = calculateCentroid(points);

    // Sort the points based on the angle relative to the centroid
    std::sort(points.begin(), points.end(), [&centroid](const cellPoint &a, const cellPoint &b) {
        double cross = crossProduct(centroid, a, b);
        return cross > 0;  // Counterclockwise order
    });
}



vector<vector<double>> InverseEquiTriangle(vector<vector<double>> &matrix)
{
  // Result matrix to store the multiplication result
  vector<vector<double>> result = {{0.,0.},{0.,0.}};

  // Define the inverse of the equilateral triangle shape tensor matrix
  double inverse_triangle_matrix[2][2] = {{1, -1 / std::sqrt(3)},
                                          {0, 2 / std::sqrt(3)}};

  // Perform matrix multiplication: result = shape_matrix * inverse_triangle_matrix
  for (int i = 0; i < 2; ++i) 
  {
    for (int j = 0; j < 2; ++j) 
    {
      for (int k = 0; k < 2; ++k) 
      {
          result[i][j] += matrix[i][k] * inverse_triangle_matrix[k][j];
      }
    }
  }
  return result;
}

// Function to multiply two 2x2 matrices
std::vector<std::vector<double>> multiplyMatrices(const vector<std::vector<double>>& a, const vector<std::vector<double>>& b) 
{
  std::vector<std::vector<double>> result(2, std::vector<double>(2, 0.0));
  result[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
  result[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
  result[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
  result[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];
  return result;
}



vector<vector<double>> CalculateQTensor(vector<cellPoint>& points)
{
  vector<vector<double>> n_matrix{{points[1].x - points[0].x, points[2].x-points[0].x}, {points[1].y - points[0].y, points[2].y-points[0].y}};
  vector<vector<double>> smatrix = InverseEquiTriangle(n_matrix);

  double smatrix_det = (smatrix[0][0]*smatrix[1][1]) - (smatrix[0][1] * smatrix[1][0]);

  // we now have the shape matrix, and need to compute 
  // the trace, symmetric, traceless part  and antisymmetric part

  // first compute trace:
  double trace = smatrix[0][0] + smatrix[1][1];
  double trace_part = trace/2;

  // traceless symmetric part:
  vector<vector<double>> tless_sym = {{(smatrix[0][0]-smatrix[1][1])/2, (smatrix[0][1]+smatrix[1][0])/2},
                                      {(smatrix[1][0]+smatrix[0][1])/2, (smatrix[1][1]-smatrix[0][0])/2}};

  double tless_sym_mag = sqrt( ( pow(tless_sym[0][0], 2) + pow(tless_sym[0][1], 2)));
  

  vector<vector<double>> antisym = {{0., (smatrix[0][1]-smatrix[1][0])/2},
                                    {(smatrix[1][0]-smatrix[0][1])/2, 0.}};

  // theta = arctan2 of (s^a_yx, txx)
  double theta = atan2(antisym[1][0], trace_part);
  vector<vector<double>> rotation_matrix = {{std::cos(-theta), std::sin(-theta)}, {-std::sin(-theta), std::cos(-theta)}};

  vector<vector<double>> q = multiplyMatrices(tless_sym, rotation_matrix);

  double prefactor = (1.0 / tless_sym_mag) * asinh(tless_sym_mag / sqrt(smatrix_det));

  q[0][0] *= prefactor;
  q[0][1] *= prefactor;
  q[1][0] *= prefactor;
  q[1][1] *= prefactor;
  return q;
}




vector<pair<int,int>> CellularPotts::SearchNforEdges()
{
  int x, y,q;
  int neighsite;
  vector<pair<int,int>> neighbours;
  int bound = 1;
  if (!par.periodic_boundaries)
  {
    bound = 3;
  }

  for ( x = bound; x < sizex-bound; x++ )
    for ( y = bound; y < sizey-bound; y++ ) 
    {
      int curcell=sigma[x][y];
      vector<int> tempn{};
      tempn.push_back(curcell);
      for (int i=1;i<=n_nb;i++) 
      {
        int xp2,yp2;
        xp2=x+nx[i]; yp2=y+ny[i];
        if (par.periodic_boundaries)
        {

          if (xp2<=0)
            xp2=sizex-2+xp2;
          if (yp2<=0)
            yp2=sizey-2+yp2;
          if (xp2>=sizex-1)
            xp2=xp2-sizex+2;
          if (yp2>=sizey-1)
            yp2=yp2-sizey+2;
        
          neighsite=sigma[xp2][yp2];
          
          // cout << "WHAT THE.." << neighsite << endl;
      
        } 
        else 
        {
          if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1)
            neighsite=-1;
          else
            neighsite=sigma[xp2][yp2];
        }
        if (curcell != neighsite && neighsite > 0)
        {
          if (find(tempn.begin(), tempn.end(), neighsite) == tempn.end())
          {
            tempn.push_back(neighsite);
          }
          }
        // must be at least three cells for it to be a vertex
        if (tempn.size() == 2)
        {
          pair<int,int> edge = {tempn[0], tempn[1]};
          pair<int,int> edge2 = {tempn[1], tempn[0]};
          if (find(neighbours.begin(), neighbours.end(), edge) == neighbours.end() && find(neighbours.begin(), neighbours.end(), edge2) == neighbours.end())
          {
            neighbours.push_back(edge);
          }
        }
      } 
    }

  return neighbours;
}

pair<double,double> CellularPotts::PhaseZValues()
{
  int phase_on_edges{};
  int phase_off_edges{};
  int phase_on_vertices{};
  int phase_off_vertices{};
  vector<pair<int,int>> neighbours = SearchNforEdges();
  for (auto &i : neighbours)
  {
    if (i.first != 0 && i.second != 0)
    {
      bool phase1 = cell->at(i.first).GetPhase();
      bool phase2 = cell->at(i.second).GetPhase();
      if (phase1==true && phase2==true)
      {
        ++phase_on_edges;
      }
      else if (phase1==false && phase2 == false)
      {
        ++phase_off_edges;
      }
    }
    else if (i.first == 0)
    {
      bool phase = cell->at(i.second).GetPhase();
      if (phase==true)
      {
        ++phase_on_edges;
      }
      else
      {
        ++phase_off_edges;
      }      
    }
    else if (i.second == 0)
    {
      bool phase = cell->at(i.first).GetPhase();
      if (phase==true)
      {
        ++phase_on_edges;
      }
      else
      {
        ++phase_off_edges;
      }      
    }
  }
  vector<vector<int>> vertices = SearchNforVertices();
  for (auto &i : vertices)
  {
    bool phase_on = true;
    bool phase_off = false;
    for (int &j : i)
    {
      if (j == 0)
      {
        continue;
      }
      bool phase = cell->at(j).GetPhase();
      if (phase != phase_on)
        phase_on = false;
      if (phase != phase_off)
        phase_off = true;
    }
    if (phase_on)
    {
      ++phase_on_vertices;
    }
    if (!phase_off)
      ++phase_off_vertices;
  }
  // cout << phase_on_edges << '\t' << phase_on_vertices << endl;
  double Z_on = 2 * double(phase_on_edges) / double(phase_on_vertices);
  double Z_off = 2 * double(phase_off_edges) / double(phase_off_vertices);
  pair<double,double> toreturn = {Z_off, Z_on};
  return toreturn;
}



// void CellularPotts::ReadZygotePicture(void) {
 
//   int pix,cells,i,j,c,p,checkx,checky;
//   char **pixelmap;
//   char pixel[3];

//   sscanf(ZYGXPM(ZYGOTE)[0],"%d %d %d %d",&checkx,&checky,&cells,&pix);

//   if ((checkx>sizex)||(checky>sizey)) { 
//     std::cerr <<  "ReadZygote: The included xpm picture is smaller than the grid!\n";
//     std::cerr << "\n Please adjust either the grid size or the picture size.\n";
//     std::cerr << sizex << "," << sizey << "," << checkx << "," << checky << "\n";
//     exit(1);
//   } 
  
//   pixelmap=(char **)malloc(cells*sizeof(char *));
//   if (pixelmap==NULL) MemoryWarning();

//   pixelmap[0]=(char *)malloc(cells*3*sizeof(char));
//   if (pixelmap[0]==NULL) MemoryWarning();
  
//   for(i=1;i<cells;i++) 
//     pixelmap[i]=pixelmap[i-1]+3;

//   for (i=0;i<cells;i++) {
//     for (j=0;j<pix;j++)
//       pixelmap[i][j]=ZYGXPM(ZYGOTE)[i+1][j];
//     pixelmap[i][pix]='\0';
//   }

//   for (i=0;i<sizex*sizey;i++) sigma[0][i]=0;
//   fprintf(stderr,"[%d %d]\n",checkx,checky);
  
//   int offs_x, offs_y;
//   offs_x=(sizex-checkx)/2;
//   offs_y=(sizey-checky)/2;
  
//   for (i=0;i<checkx;i++)
//     for (j=0;j<checky;j++) {
//       for (p=0;p<pix;p++)
//         pixel[p]=ZYGXPM(ZYGOTE)[cells+1+j][i*pix+p];
      
//       pixel[pix]='\0';

//       for (c=0;c<cells;c++) {
// 	if (!(strcmp(pixelmap[c],pixel))) {
// 	  if ( (sigma[offs_x+i][offs_y+j]=c) ) {
	
// 	    // if c is _NOT_ medium (then c=0)
// 	    // assign pixel values from "sigmamax"
// 	    sigma[offs_x+i][offs_y+j]+=(Cell::MaxSigma()-1);
// 	  }
// 	}
	
//       }
//     }

//   free(pixelmap[0]);
//   free(pixelmap);
// }


void CellularPotts::ConstructInitCells (Dish &beast) {
  
  // Get the maximum cell ID (mostly equal to the cell number)
  int loop=sizex*sizey;
  int cells=0;
  for (int i=0;i<loop;i++) {
    if (cells<sigma[0][i]) cells=sigma[0][i];
  }

  // cerr << "[ cells = " << cells << "]\n";

  // construct enough cells for the zygote.  "cells", contains the
  // number of colours (excluding background).
  { 
    for (int i=0; i<cells; i++) {
      cell->push_back(Cell(beast));
    }
  }
  
  // Set the area and target area of the cell
  // makes use of the pointer to the Cell pointer of Dish
  // which is a member of CellularPotts 
  MeasureCellSizes();
  
  // set zygote_area to mean cell area.
  int mean_area=0;
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    mean_area+=c->Area();
  }
  if (cells!=0) 
    mean_area/=cells;
  
  zygote_area=mean_area;

  // cout << "mean_area = " << mean_area << "\n";
  // set all cell areas to the mean area
  {
    for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
      if (par.init_area) {
	c->SetTargetArea(par.init_area);
      } else	 {
	c->SetTargetArea(mean_area);
      }
    }
  }
}

double CellularPotts::SumEnergy()
{
  double sumH=0;
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) 
  {
    if (c->AliveP()) 
    {
      int ca = c->Area();
      int ta = c->TargetArea();
      int perim = c->Perimeter();
      int tperim = c->TargetPerimeter();

      // cout << ca << '\t' << ta << '\t' << perim << '\t' << tperim << endl;

      int area_energy = par.lambda * pow((ca - ta),2);
      int perim_energy = par.lambda_perimeter * pow((perim - tperim),2);

      sumH += area_energy;
      sumH += perim_energy;
    } 
  }
  return sumH;
}







void CellularPotts::MeasureCellSizes(void) {
  
  // Clean areas of all cells, including medium

  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->CleanMoments();
  }
  
  // calculate the area of the cells
  for (int x=1;x<sizex-1;x++) {
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	(*cell)[sigma[x][y]].IncrementTargetArea();
	(*cell)[sigma[x][y]].IncrementArea();
	(*cell)[sigma[x][y]].AddSiteToMoments(x,y);

      }
    }
  }
  
}

void CellularPotts::MeasureCellSize(Cell &c) {
  
  c.CleanMoments();
  
  // calculate the area of the cell
  for (int x=1;x<sizex-1;x++) {
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y] == c.sigma) {
	(*cell)[sigma[x][y]].IncrementTargetArea();
	(*cell)[sigma[x][y]].IncrementArea();
	(*cell)[sigma[x][y]].AddSiteToMoments(x,y);

      }
    }
  }
  
//   // set the actual area to the target area
//   {
//   for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
//     c->SetAreaToTarget();

//   }

}

Dir *CellularPotts::FindCellDirections(void) const
{ 
  
  double *sumx=0,*sumy=0;
  double *sumxx=0,*sumxy=0,*sumyy=0;
  double *n=0;  

  double xmean=0,ymean=0,sxx=0,sxy=0,syy=0;
  double D,lb1=0,lb2=0;

  Dir *celldir;

  /* Allocation of sufficient memory space */
  if( (sumx= (double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning(); 
  else
    if( (sumy= (double *)malloc(cell->size()*sizeof(double)))==NULL) 
      MemoryWarning();
    else
      if ((sumxx=(double *)malloc(cell->size()*sizeof(double)))==NULL) 
	MemoryWarning();
      else
	if((sumxy=(double *)malloc(cell->size()*sizeof(double)))==NULL) 
	  MemoryWarning();
	else
	  if((sumyy=(double *)malloc(cell->size()*sizeof(double)))==NULL)
	    MemoryWarning();
	  else
	    if((n=(double *)malloc(cell->size()*sizeof(double)))==NULL) 
	      MemoryWarning();
  
  
  if ( !(celldir=new Dir[cell->size()]) )
    MemoryWarning();

  	
  /* Initialization of the variables */
   
  for (int i=0;i<(int)cell->size();i++) {
    
    sumx[i]=0.;
    sumy[i]=0.;
    sumxx[i]=0.;
    sumxy[i]=0.;
    sumyy[i]=0.;
    n[i]=0L;

  }


  /* Find sumx, sumy, sumxx and sumxy for all cells */

  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++) 
      if (sigma[x][y]>0) {
	sumx[0]+=(double)x;
	sumy[0]+=(double)y;
	sumxx[0]+=(double)x*x;
	sumxy[0]+=(double)x*y;
	sumyy[0]+=(double)y*y;
	
	n[0]++;
	
	sumx[sigma[x][y]]+=(double)x;
	sumy[sigma[x][y]]+=(double)y;
	
	sumxx[sigma[x][y]]+=(double)x*x;
	sumxy[sigma[x][y]]+=(double)x*y;
	sumyy[sigma[x][y]]+=(double)y*y;
	
	n[sigma[x][y]]++;
	
      }
  
  /* Compute the principal axes for all cells */
  
  {
    for (int i=0;i<(int)cell->size();i++) {
    
      if (n[i]>10) {
      
	xmean=((double)sumx[i])/((double)n[i]);
	ymean=((double)sumy[i])/((double)n[i]);

	sxx=(double)(sumxx[i])-((double)(sumx[i]*sumx[i]))/(double)n[i];
	sxx=sxx/(double)(n[i]-1);

	sxy=(double)(sumxy[i])-((double)(sumx[i]*sumy[i]))/(double)n[i];
	sxy=sxy/(double)(n[i]-1);

	syy=(double)(sumyy[i])-((double)(sumy[i]*sumy[i]))/(double)n[i];
	syy=syy/(double)(n[i]-1);

	D=sqrt( (sxx+syy)*(sxx+syy)-4.*(sxx*syy-sxy*sxy) );
	lb1=(sxx+syy+D)/2.;lb2=(sxx+syy-D)/2.;
	celldir[i].lb1=lb1; celldir[i].lb2=lb2; 
      }
      if (sxy==0.0)
	celldir[i].bb1=1.; 
      else
	celldir[i].bb1=sxy/(lb1-syy);
    
      if (fabs(celldir[i].bb1)<.00001) {
	if (celldir[i].bb1>0.) 
	  celldir[i].bb1=.00001;
	else 
	  celldir[i].bb1=-.00001;
      }
 
      celldir[i].aa1=ymean-xmean*celldir[i].bb1;
      celldir[i].bb2= (-1.)/celldir[i].bb1;
    
      celldir[i].aa2=ymean-celldir[i].bb2*xmean;     
    }
		  
  }

  /* bevrijd gealloceerd geheugen */
  free(sumx);
  free(sumy);
  free(sumxx);
  free(sumxy);
  free(sumyy);
  free(n);

  return celldir;
 
}

void CellularPotts::ShowDirections(Graphics &g, const Dir *celldir) const
{
  int i;
  
  if (cell->size()>1) 
    for (i=1;i<(int)cell->size();i++)
      g.Line(0,(int)(2*celldir[i].aa1),sizex*2,(int)((celldir[i].aa1+celldir[i].bb1*sizey)*2),2);
  
}



void CellularPotts::DivideCellsNoGrid(vector<bool> which_cells)
{

  int n_cells = which_cells.size();
  for (int i=0; i<n_cells;++i)
  {
    if (which_cells[i]==true)
    {
      Cell *motherp=&((*cell)[i]);
      Cell *daughterp;

      // add daughter cell, copying states of mother
      daughterp=new Cell(*(motherp->owner));
      daughterp->CellBirth(*motherp);
      cell->push_back(*daughterp);
      // renew pointer to mother
      motherp=&((*cell)[i]);
      delete daughterp;
      // array may be relocated after "push_back"
      // renew daughter pointers
      daughterp=&(cell->back());
    }
  }
}


// Multithreading method
void CellularPotts::set_seed()
{
  if (par.pickseed)
    s_val[0] = par.pickseed;
  else
  {
    s_val[0] = Seed(org_num);
    if (par.gene_output && par.gene_record)
      par.pickseed = s_val[0];
  }
  if (par.print_fitness)
  {
    cout << "Seed is: " << s_val[0] << endl;
  }

  
}

void CellularPotts::set_datafile(string file)
{
  data_file = file;
}

void CellularPotts::SetAreas(int tarea)
{
  vector<Cell>::iterator i;
  for ( (i=cell->begin(),i++); i!=cell->end(); i++) 
  {
    i->SetTargetArea(tarea);
  }
}

void CellularPotts::SetPerims(int tperim)
{
  vector<Cell>::iterator i;
  for ( (i=cell->begin(),i++); i!=cell->end(); i++) 
  {
    i->SetTargetPerimeter(par.ptarget_perimeter);
  }
}


void CellularPotts::DivideCells(vector<bool> which_cells, int t)
{
  // for the cell directions
  Dir *celldir=0;
  
  /* Allocate space for divisionflags */
  int *divflags=(int *)malloc((cell->size()*2+5)*sizeof(int));
  
  /* Clear divisionflags */
  for (int i=0;i<(int)(cell->size()*2+5);i++) 
    divflags[i]=0;
  
  
  if ( !(which_cells.size()==0 || which_cells.size()>=cell->size()) ) {
    throw "In CellularPotts::DivideCells, Too few elements in vector<int> which_cells.";
  }
  
  /* division */
  {
  for (int i=0;i<sizex;i++)
    for (int j=0;j<sizey;j++) 
      if (sigma[i][j]>0) // i.e. not medium and not border state (-1)
      { 

        // Pointer to mother. Warning: Renew pointer after a new
        // cell is added (push_back). Then, the array *cell is relocated and
        // the pointer will be lost...
        
        Cell *motherp=&((*cell)[sigma[i][j]]);
        Cell *daughterp;
        
        /* Divide if NOT medium and if DIV bit set or divide_always is set */
        // if which_cells is given, divide only if the cell
        // is marked in which_cells.
        if  ( !which_cells.size() || which_cells[motherp->sigma] )    
        {

          if (!(divflags[ motherp->Sigma() ]) ) 
          {
      
            // add daughter cell, copying states of mother
            daughterp=new Cell(*(motherp->owner));
            daughterp->CellBirth(*motherp);
            cell->push_back(*daughterp);
            // renew pointer to mother
            motherp=&((*cell)[sigma[i][j]]);

            divflags[ motherp->Sigma() ]=daughterp->Sigma();
            delete daughterp;
            // array may be relocated after "push_back"
            
            // renew daughter pointers
            daughterp=&(cell->back());


          /* administration on the onset of mitosis */
          
          /* Ancestry is taken care of in copy constructor of Cell 
            see cell.h: Cell(const Cell &src, bool newcellP=false) : Cytoplasm(src) {} */
          
          /* inherit  polarity of mother */
          // All that needs to be copied is copied in the copy constructor
          // of Cell and in the default copy constr. of its base class Cytoplasm
          // note: also the celltype is inherited
          } 
          else 
          {
            daughterp=&((*cell)[ divflags[motherp->Sigma()] ]);
          }
            
          /* Now the actual division takes place */
          /* If celldirections where not yet computed: do it now */
          if (!celldir) 
            celldir=FindCellDirections();
        
            /* if site is below the minor axis of the cell: sigma of new cell */
          if (j>((int)(celldir[motherp->sigma].aa2 + celldir[motherp->sigma].bb2*(double)i))) 
          { 
            motherp->DecrementArea();
            motherp->DecrementTargetArea();
            motherp->RemoveSiteFromMoments(i,j);
            sigma[i][j]=daughterp->Sigma();
            daughterp->AddSiteToMoments(i,j);
            daughterp->IncrementArea();
            daughterp->IncrementTargetArea();

          } 
        }
      }
  }  
  if (par.H_perim)
    MeasureCellPerimeters();

 
  if (celldir) 
    delete[] (celldir);
  
  if (divflags)
    free(divflags);
}   

int CellularPotts::CountCells(void) const 
{
  int amount=0;
  vector<Cell>::const_iterator i;
  for ( (i=cell->begin(),i++); i!=cell->end(); i++) 
  {
    if (i->AliveP()) 
    {
      amount++;
    }
  }
  return amount;
}


bool CellularPotts::SpawnCell(int x, int y, int cp_sigma, int time)
{
  if (x < 5 || x > sizex-5 || y < 5 || y > sizey-5)
  {
    cerr << "Spawned cell outside of grid. Simulation should now end.\n";
    return true;
  }
  else if (sigma[x][y] > 0)
  {
    cerr << "Spawned cell in bad spot.\n"; 
    cerr << "xy is: " << x << '\t' << y << '\n' << endl;
    return false;
  }
  else
  {
    Cell *new_cell;
    Cell *cp=&((*cell)[cp_sigma]);
    new_cell = new Cell(*(cp->owner));
    new_cell->CellBirth(*cp);
    cell->push_back(*new_cell);
    delete new_cell;
    // cp = &((*cell)[cp_sigma]);
    // new_cell=&(cell->back());
    

    queue<pair<int, int>> q; // Queue for BFS
    q.push({x, y});

    int cell_size = par.cell_target_area;
    int cell_sigma=cell->back().Sigma();
    // cout << "cp sigma is: " << cp_sigma << endl;
    // cout << "sigma is: " << cell_sigma << endl;
    // cout << x << '\t' << y << endl;
    
    int filledPixels = 0;
    while (!q.empty() && filledPixels < cell_size) 
    {
      auto [x, y] = q.front();
      q.pop();
      // Skip if this cell is already filled or out of bounds
      if (x < 1 || x >= sizex-1 || y < 1 || y >= sizey-1)
      {
        return true;
      }
      else if (sigma[x][y]>0 ) 
      {
        // cout << "FAILED TO ADD TO: " << x << '\t' << y << endl;
        continue;
      }
      // cout << "added site: " << x << '\t' << y << endl;
      // Fill the pixel
      sigma[x][y] = cell_sigma;
      cell->back().AddSiteToMoments(x,y);
      cell->back().IncrementArea();
      cell->back().IncrementTargetArea();
      filledPixels++;
      
      // Add neighboring cells (up, down, left, right) to the queue
      for (int i = 0; i < nbh_level[2]; i++) 
      {
        int newX = x + nx[i];
        int newY = y + ny[i];
        if (newX > 1 && newX < sizex-1 && newY >= 1 && newY < sizey-1 && sigma[newX][newY]==0) 
        {
          q.push({newX, newY});
        }
      }
    }
    if (par.H_perim)
      MeasureSinglePerimeter(cell->back().Sigma());
  }
  return true;
}





// Function to calculate the Euclidean distance between two points
double calculate_distance(const pair<int, int>& p1, const pair<int, int>& p2) {
    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}



void CellularPotts::ConstructSheet(int xm, int ym)
{    
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++)
    {
      // if (x < xm && y < ym)
      // {
      //   sigma[x][y] = 1;
      // }

      // this makes a triangle
      if (x < par.triangle_x && y > par.triangle_y + x)
      {
        sigma[x][y] = 1;
      }     
    }

}





/**! Fill the plane with initial cells 
 \return actual amount of cells (some are not draw due to overlap) */
int CellularPotts::ThrowInCells(int n,int cellsize) {
  
  //  int gapx=(sizex-nx*cellsize)/(nx+1);
  //int gapy=(sizey-ny*cellsize)/(ny+1);
  
  int cellnum=1;

  for (int i=0;i<n;i++) {
    
    // draw a circle at x0, y0
    int x0=RandomNumber(sizex, s_val);
    int y0=RandomNumber(sizey, s_val);
   
    bool overlap=false;
    
    // check overlap
    for (int x=0;x<cellsize;x++)
      for (int y=0;y<cellsize;y++)
	if ( ( 
	      ( (x-cellsize/2)*(x-cellsize/2)+(y-cellsize/2)*(y-cellsize/2) )<
	      ( (cellsize/2)*(cellsize/2))) &&
	     ( x0+x<sizex && y0+y<sizey ) )
	  if (sigma[x0+x][y0+y]) {
	    overlap=true;
	    break;
	  }
    
    if (!overlap) {
      for (int x=0;x<cellsize;x++)
	for (int y=0;y<cellsize;y++)
	  if ( ( 
		( (x-cellsize/2)*(x-cellsize/2)+(y-cellsize/2)*(y-cellsize/2) )<
		( (cellsize/2)*(cellsize/2))) &&
	       ( x0+x<sizex && y0+y<sizey ) )
	    sigma[x0+x][y0+y]=cellnum;
      
      cellnum++;
    }
  }
  cerr << "[ cellnum = " << cellnum << "]";

  // repair borders
  // fill borders with special border state
  for (int x=0;x<sizex-1;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey-1;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }

  {for (int x=1;x<sizex-2;x++) {
      sigma[x][1]=0;
      sigma[x][sizey-2]=0;
    }}
  {for (int y=1;y<sizey-2;y++) {
      sigma[1][y]=0;
      sigma[sizex-2][y]=0;
    }}
  return cellnum;
} 

  

// Function to fill grid with cell. 
void CellularPotts::FillGrid()
{
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) 
    {
      sigma[x][y]=1;
    }
}

double euclideanDistance(int x1, int y1, int x2, int y2, int sizex, int sizey) 
{
  // Calculate direct distances
  double dx = std::abs(x2 - x1);
  double dy = std::abs(y2 - y1);
  if (par.periodic_boundaries)
  {
    // Apply periodic boundary conditions
    if (dx > sizex / 2) {
        dx = sizex - dx - 2;  // Wrap around horizontally
    }
    if (dy > sizey / 2) {
        dy = sizey - dy - 2;  // Wrap around vertically
    }
  }

  
  // Return the Euclidean distance
  return std::sqrt(dx * dx + dy * dy);
}



void CellularPotts::ClearGrid()
{
  for (int x = 1; x < sizex - 1; ++x) 
  {
    for (int y = 1; y < sizey - 1; ++y) 
    {
      sigma[x][y] = 0;
    }
  }
}

// Added 'double R' to the function parameters
void CellularPotts::PopulateSparseCells(double density, double R, int shiftx, int shifty)
{

  int current_cells = CountCells();
  // 1. Define the usable active grid space
  int W = sizex - 2;
  int H = sizey - 2;
  
  double target_area = static_cast<double>(par.cell_target_area);
  
  // Guard against invalid parameters
  if (target_area <= 0 || density <= 0.0 || R <= 0.0) return;
  
  // 2. Calculate lattice spacing directly from density and target area.
  // This ensures the density inside the radius is exactly what was requested.
  double area_per_center = target_area / density;
  
  // Hexagonal area per center = 2 * sqrt(3) * r^2
  double r = std::sqrt(area_per_center / (2.0 * std::sqrt(3.0)));

  // 3. Determine lattice parameters to cover the grid
  int num_cols = static_cast<int>(std::round(W / (2.0 * r)));  
  int num_rows = static_cast<int>(std::round(H / (std::sqrt(3.0) * r)));
  if (num_cols <= 0) num_cols = 1;
  if (num_rows <= 0) num_rows = 1;

  // 4. Calculate actual lattice bounds to perfectly center the cells on the grid
  double max_x = 0, max_y = 0;
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      double cx = col * 2.0 * r;
      double cy = row * std::sqrt(3.0) * r;
      if (row % 2 == 1) cx += r; // Stagger odd rows
      
      if (cx > max_x) max_x = cx;
      if (cy > max_y) max_y = cy;
    }
  }

  // Offset ensures perfectly symmetric empty padding near the walls
  double offset_x = (W - max_x) / 2.0 + 1.0 + shiftx; 
  double offset_y = (H - max_y) / 2.0 + 1.0 + shifty;

  // Find the exact mathematical center of the grid
  double grid_center_x = sizex / 2.0 + shiftx;
  double grid_center_y = sizey / 2.0 + shifty;

  // 5. Generate the centers, filtering by distance to grid center <= R
  struct VPoint { double x, y; int id; };
  std::vector<VPoint> centers;

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      double cx = col * 2.0 * r;
      double cy = row * std::sqrt(3.0) * r;
      if (row % 2 == 1) cx += r;
      
      double final_cx = cx + offset_x;
      double final_cy = cy + offset_y;
      
      // Calculate distance from grid center. 
      // Re-using your euclideanDistance function ensures periodic boundaries are respected if applicable.
      double dist = euclideanDistance(final_cx, final_cy, grid_center_x, grid_center_y, sizex, sizey);
      
      if (dist <= R) {
          centers.push_back({final_cx, final_cy, -1});
      }
      if (final_cx < 5 || final_cx > sizex-5 || final_cy < 5 || final_cy > sizey-5)
        cerr << "warning: some centers are outside of domain\n";
    }
  }

  int final_ncells = centers.size();
  if (final_ncells == 0) {
      std::cout << "Density/Radius too low to generate any cells inside R." << std::endl;
      return;
  }

  // 6. Split sheet to prepare sufficient cell instances
  FractureSheet(final_ncells);

  // Map the new alive cell IDs to our spatial centers
  std::vector<int> valid_ids;
  for (auto c = std::next(cell->begin(), current_cells); c != cell->end(); ++c) {
    if (c == cell->begin()) continue; // Skip medium/background index 0
    if (c->AliveP()) {
      valid_ids.push_back(c->Sigma());
    }
  }

  for (size_t i = 0; i < centers.size(); ++i) {
      if (i < valid_ids.size()) {
          centers[i].id = valid_ids[i];
      }
  }


  // 8. Draw the Voronoi domains, strictly bounded to achieve `par.cell_target_area` size.
  // Formula for circle radius: R = sqrt(A / pi). 
  // We apply a slight 5% buffer to account for discrete pixelation artifacts cutting areas short.
  double radius_limit = std::sqrt(target_area / M_PI) * 1.05; 

  for (int x = 1; x < sizex - 1; ++x) {
      for (int y = 1; y < sizey - 1; ++y) {
        double minDistance = std::numeric_limits<double>::max();
        int closestCenter = -1;
        
        for (const auto& center : centers) {
          if (center.id == -1) continue;
          
          double dist = euclideanDistance(x, y, center.x, center.y, sizex, sizey);
          if (dist < minDistance) {
              minDistance = dist;
              closestCenter = center.id;
          }
        }
          
        if (minDistance < radius_limit) {
            sigma[x][y] = closestCenter;
        }
      }
  }

  // 9. Re-evaluate actual populated cell areas
  for (auto c = cell->begin(); c != cell->end(); ++c) {
    if (c == cell->begin()) 
      continue;
    if (c->AliveP()) 
      c->area = 0;
  }

  for (int x = 1; x < sizex - 1; ++x) {
    for (int y = 1; y < sizey - 1; ++y) {
      if (sigma[x][y] > 0) 
      {
        (*cell)[sigma[x][y]].area += 1;
        (*cell)[sigma[x][y]].makeAlive();
      }
    }   
  }
  
  // 10. Single unified cleanup and target initialization loop
  int deadcells = 0;
  for (auto c = std::next(cell->begin(), current_cells); c != cell->end(); ++c) 
  {
    if (c == cell->begin()) continue;
    if (c->AliveP()) 
    {
      if (c->area == 0) 
      {
        c->Apoptose();
        ++deadcells;
      } 
      else 
      {
        c->SetTargetArea(c->area);
      }
    }
  }
  MeasureCellSizes();

  std::cout << "Grid generated | Radius: " << R 
            << " | Density: " << density 
            << " | Cells populated: " << (final_ncells - deadcells)
            << " | Cells killed (0 area): " << deadcells << std::endl;
}



//split sheet into cells
void CellularPotts::FractureSheet()
{
  
  bool dividing = true;

  while (dividing)
  {
    vector<bool> which_cells(cell->size());
    dividing = false;
    vector<Cell>::iterator c;
    for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    {
      if (c->AliveP())
      {
        int area = c->Area();  
        if (area>par.div_threshold)
        {

          dividing = true;
          which_cells[c->Sigma()]=true;

        }
      }
    }
    if (dividing)
      DivideCells(which_cells);
  }
}

//split sheet into cells
void CellularPotts::FractureSheet(int n_cells)
{
  int counter=0;
  bool dividing = true;

  while (dividing)
  {
    vector<bool> which_cells(cell->size());
    vector<Cell>::iterator c;
    for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    {
      if (c->AliveP())
      {
        if (counter<=n_cells)
        {

          which_cells[c->Sigma()]=true;
          ++counter;
        }
      }
    }
    if (dividing)
      DivideCellsNoGrid(which_cells);
    if (counter >= n_cells)
    {
      dividing=false;
    }
  }
}



struct VPoint 
{
    double x, y;
    int id;
     // Corresponding Voronoi seed
};

vector<VPoint> HexaCenters(int m, int n, double r)
{
  vector<VPoint> centers;
  int num_cols = static_cast<int>(std::floor(m / (2 * r)));  
  int num_rows = static_cast<int>(std::floor(n / (sqrt(3) * r)));
  int center_count = 1;

    // Generate centers
    for (int row = 0; row < num_rows; ++row) 
    {
      for (int col = 0; col < num_cols; ++col) 
      {
          double x = col * 2 * r;
          double y = row * sqrt(3) * r;
          
          // Stagger odd rows
          if (row % 2 == 1) {
              x += r;  // Shift odd rows horizontally by r
          }
          
          // Ensure the center is within the grid bounds
          if (x < m+2 && y < n+2) {
              centers.push_back({x+2+(r/2), y+2+(r/2), center_count});
              center_count++;
          }
        }
    }
    // cout << "CENTRE COUNT IS: " << center_count << endl;
    return centers;
}


vector<VPoint> HexaCircleCenters(double circle_radius, double dist, int centerx, int centery, int starting_value)
{
    vector<VPoint> centers;
    
    // Calculate the max rows and columns needed to cover the circle's radius
    int max_row = static_cast<int>(std::floor(circle_radius / (sqrt(3) * dist)));
    int max_col = static_cast<int>(std::floor(circle_radius / (2 * dist))) + 1;
    
    int center_count = starting_value;
    double radius_squared = circle_radius * circle_radius;

    // Generate centers from -max to +max (centered at 0,0)
    for (int row = -max_row; row <= max_row; ++row) 
    {
        for (int col = -max_col; col <= max_col; ++col) 
        {
            double x = col * 2 * dist;
            double y = row * sqrt(3) * dist;
            
            // Stagger odd rows
            // (Using != 0 ensures negative odd rows are handled correctly in C++)
            if (row % 2 != 0) {
                x += dist;  // Shift odd rows horizontally by dist
            }
            
            // Ensure the center is within the circular bounds
            if ((x * x + y * y) <= radius_squared) {
                
                // Keep your original padding of + 2 + (dist/2)
                double final_x = x + 2 + (dist / 2);
                double final_y = y + 2 + (dist / 2);
                
                // uncomment the two lines below to shift the whole circle into the positive quadrant:
                final_x += centerx;
                final_y += centery;

                centers.push_back({final_x, final_y, center_count});
                center_count++;
            }
        }
    }
    
    // cout << "CENTRE COUNT IS: " << center_count - 1 << endl;
    return centers;
}


int HexaCounter(int m, int n, double r)
{
  // int num_rows = static_cast<int>(std::floor(m / (sqrt(3) * r)));
  // int num_cols = static_cast<int>(std::floor(n / (2 * r)));  

  int num_cols = static_cast<int>(std::floor(m / (2 * r)));  
  int num_rows = static_cast<int>(std::floor(n / (sqrt(3) * r)));


  int center_count = 0;
    // Generate centers
    for (int row = 0; row < num_rows; ++row) 
    {
      for (int col = 0; col < num_cols; ++col) 
      {
          double x = col * 2 * r;
          double y = row * sqrt(3) * r;
          
          // Stagger odd rows
          if (row % 2 == 1) {
              x += r;  // Shift odd rows horizontally by r
          }
          
          // Ensure the center is within the grid bounds
          if (x < m && y < n) {
              center_count++;
          }
        }
    }
    return center_count;
}

int CircleHexaCounter(double circle_radius, double dist)
{
  // Calculate how many rows/cols we need to check in the positive/negative directions.
  // We add 1 to the columns just to provide a safe buffer for the odd-row shift.
  int max_row = static_cast<int>(std::floor(circle_radius / (sqrt(3) * dist)));
  int max_col = static_cast<int>(std::floor(circle_radius / (2 * dist))) + 1;

  int center_count = 0;
  double radius_squared = circle_radius * circle_radius;

  // Generate centers from -max to +max to cover the whole circle centered at (0,0)
  for (int row = -max_row; row <= max_row; ++row) 
  {
    for (int col = -max_col; col <= max_col; ++col) 
    {
      double x = col * 2 * dist;
      double y = row * sqrt(3) * dist;
      
      // Stagger odd rows. 
      // NOTE: We use != 0 instead of == 1 because C++ modulo on negative numbers yields negative results.
      if (row % 2 != 0) 
      {
          x += dist;  // Shift odd rows horizontally by r
      }
      
      // Ensure the center is within the circle using x^2 + y^2 <= R^2
      // (Using squared values avoids a costly std::sqrt() call)
      if ((x * x + y * y) <= radius_squared) {
          center_count++;
      }
    }
  }
  return center_count;
}







void CellularPotts::Voronoi()
{
  // double total = sizex*sizey;
  // int ncells = round(total / 75.);
  // cout << ncells << endl;
  double A = double(par.cell_target_area);
  double distance = sqrt((A)/(2*sqrt(3)));
  double leftover = fmod(sizex-2, distance);
  int dividor = int(floor(double(sizex-2)/distance));
  // cout << "LEFTOVERS: " << leftover << '\t' << dividor << endl;
  distance += leftover/dividor;


  int ncells = HexaCounter(sizex-2,sizey-2,distance);
  FractureSheet(ncells);

  int periodic_length_x = sizex - 2;
  int periodic_length_y = sizey - 2;

  vector<VPoint> centers = HexaCenters(periodic_length_x, periodic_length_y, distance);
  // for (const auto& center : centers) 
  // {
  //     std::cout << "Center at (" << center.x << ", " << center.y << ")\n";
  // }


  for (int x = 1; x < sizex-1; ++x) {
      for (int y = 1; y < sizey-1; ++y) 
      {
        double minDistance = std::numeric_limits<double>::max();
        int closestCenter = -1;
        
        // Find the closest center to (i, j)
        for (const auto& center : centers) 
        {
          double dist = euclideanDistance(x, y, center.x, center.y, sizex, sizey);
          if (dist < minDistance) {
              minDistance = dist;
              closestCenter = center.id;
          }
        }
          
        // Assign the closest center id to the grid cell
        sigma[x][y] = closestCenter;
      }
  }

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->area = 0;
    }
  }

  for (int x=1; x<sizex; ++x)
    for (int y=1; y<sizey; ++y)
    {
      if (sigma[x][y] > 0)
      {
        (*cell)[sigma[x][y]].area +=1;
      }
    }   
  
  int deadcells{};
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      if (!c->area)
      {
        c->Apoptose();
        ++deadcells;
      }
      else
      {
        c->SetTargetArea(c->area);
        // cout << c->area << endl;
      }
    }
  }
}



void CellularPotts::VoronoiSeparated(int xlen, int ylen, int shift, int xshift, bool turnonphase)
{
  // 1. Calculate the base dimensions from the target area
  double A = double(par.cell_target_area);
  
  // This 'base_distance' represents the radius required for the target Area A 
  // in a packed configuration. We will use this to limit the drawn cell size.
  double base_distance = sqrt((A)/(2*sqrt(3)));

  // 2. Define a Spacing Factor to spread the points out.
  // A factor of 1.25 means the centers are 25% further apart than in the dense packing.
  // This creates the "separation" or gaps between cells.
  double spread_factor = 1.25; 
  double spacing_distance = base_distance * spread_factor;

  // 3. Adjust spacing to fit the grid perfectly (Leftover logic applied to spacing)
  // We apply the fitting logic to 'spacing_distance' because that determines the lattice.
  double leftover = fmod(xlen-2, spacing_distance);
  int dividor = int(floor(double(xlen-2)/spacing_distance));
  spacing_distance += leftover/dividor;

  // 4. Generate centers using the LARGER spacing_distance (Fewer cells, spread out)
  int ncells = HexaCounter(xlen-2,ylen-2, spacing_distance);
  FractureSheet(ncells);

  // Clear grid
  for (int x = 1; x < sizex-1; ++x) 
  {
    for (int y = 1; y < sizey-1; ++y) 
    {
      sigma[x][y]=0;
    }
  }

  int periodic_length_x = xlen - 2;
  int periodic_length_y = ylen - 2;

  // Get centers based on the spread out spacing
  vector<VPoint> centers = HexaCenters(periodic_length_x, periodic_length_y, spacing_distance);
  
  // Shift centers to the correct ROI
  for (auto& center : centers) 
  {
    center.x += sizex - xlen - xshift;
    center.y += sizey - ylen - shift;
  }

  // 5. Draw the cells, but LIMIT the radius to 'base_distance'
  // This ensures the cells are roughly the size of 'par.cell_target_area', 
  // but because the centers are spread out, there will be empty space (Medium 0) between them.
  
  // We use a slight multiplier (e.g., 1.1) on base_distance for the drawing limit 
  // to allow them to be slightly organic/hexagonal but definitely not touching 
  // neighbors calculated at 1.25x distance.
  double radius_limit = base_distance * 1.1; 

  for (int x = 1 + sizex - xlen - xshift; x < sizex - xshift - 1; ++x) {
      for (int y = 1 + sizey - ylen - shift; y < sizey - shift - 1; ++y) 
      {
        double minDistance = std::numeric_limits<double>::max();
        int closestCenter = -1;
        
        // Find the closest center
        for (const auto& center : centers) 
        {
          double dist = euclideanDistance(x, y, center.x, center.y, sizex, sizey);
          if (dist < minDistance) {
              minDistance = dist;
              closestCenter = center.id;
          }
        }
          
        // KEY CHANGE: Only assign if within the radius limit.
        // Since centers are spaced by `spacing_distance` (1.25x) and we cut off at 
        // `radius_limit` (~1.0x), gaps are guaranteed.
        if (minDistance < radius_limit) {
            sigma[x][y] = closestCenter;
        }
      }
  }

  // Standard cleanup and initialization logic from original function
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->area = 0;
    }
  }

  for (int x=1; x<sizex; ++x)
    for (int y=1; y<sizey; ++y)
    {
      if (sigma[x][y] > 0)
      {
        (*cell)[sigma[x][y]].area +=1;
      }
    }   
  
  int deadcells{};
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      if (!c->area)
      {
        c->Apoptose();
        ++deadcells;
      }
      else
      {
        c->SetTargetArea(c->area);
      }
    }
  }
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      if (!c->area)
      {
        c->Apoptose();
        ++deadcells;
      }
      else
      {
        c->SetTargetArea(c->area);
        c->makeAlive();
      }
    }
  }
  MeasureCellSizes();

  cout << "Total cells killed: " << deadcells << endl;
}


void CellularPotts::Voronoi(int xlen, int ylen, int shift, int xshift, bool turnonphase)
{


  // double total = sizex*sizey;
  // int ncells = round(total / 75.);
  // cout << ncells << endl;
  double A = double(par.cell_target_area);
  double distance = sqrt((A)/(2*sqrt(3)));
  double leftover = fmod(xlen-2, distance);
  int dividor = int(floor(double(xlen-2)/distance));
  // cout << "LEFTOVERS: " << leftover << '\t' << dividor << endl;
  distance += leftover/dividor;


  int ncells = HexaCounter(xlen-2,ylen-2,distance);
  // make lots of cells
  FractureSheet(ncells);
  // cout << ncells << endl;
  // first need to clear the grid:
  for (int x = 1; x < sizex-1; ++x) 
  {
    for (int y = 1; y < sizey-1; ++y) 
    {
      sigma[x][y]=0;
    }
  }


  int periodic_length_x = xlen - 2;
  int periodic_length_y = ylen - 2;

  vector<VPoint> centers = HexaCenters(periodic_length_x, periodic_length_y, distance);
  for (auto& center : centers) 
  {
    center.x += sizex - xlen - xshift;
    center.y += sizey - ylen - shift;
    // std::cout << "Center at (" << center.x << ", " << center.y << ")\n";
  }


  for (int x = 1 + sizex - xlen - xshift; x < sizex - xshift - 1; ++x) {
      for (int y = 1 + sizey - ylen - shift; y < sizey - shift - 1; ++y) 
      {
        double minDistance = std::numeric_limits<double>::max();
        int closestCenter = -1;
        
        // Find the closest center to (i, j)
        for (const auto& center : centers) 
        {
          double dist = euclideanDistance(x, y, center.x, center.y, sizex, sizey);
          if (dist < minDistance) {
              minDistance = dist;
              closestCenter = center.id;
          }
        }
          
        // Assign the closest center id to the grid cell
        sigma[x][y] = closestCenter;
      }
  }

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->area = 0;
    }
  }

  for (int x=1; x<sizex; ++x)
    for (int y=1; y<sizey; ++y)
    {
      if (sigma[x][y] > 0)
      {
        (*cell)[sigma[x][y]].area +=1;
      }
    }   
  
  int deadcells{};
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      if (!c->area)
      {
        c->Apoptose();
        ++deadcells;
      }
      else
      {
        c->SetTargetArea(c->area);
        c->makeAlive();
        // cout << c->area << endl;
      }
    }
  }
  MeasureCellSizes();

  cout << "Total cells killed: " << deadcells << endl;
}


void CellularPotts::ApoptoseDeadCells(void) {
  
  // Clean areas of all cells, including medium

  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->area = 0;
  }
  
  // calculate the area of the cells
  for (int x=1;x<sizex-1;x++) {
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	(*cell)[sigma[x][y]].IncrementArea();

      }
    }
  }
  
  // set the actual area to the target area
  {
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) 
  {
    if (c->Area() == 0)
    {
      cout << "Found dead cell/" << endl;
      c->Apoptose();
    }


  }
  }
}



int CellularPotts::GrowInCells(int n_cells, int cell_size, double subfield) {
  
  int sx = (int)((sizex-2)/subfield);
  int sy = (int)((sizey-2)/subfield);
  
  int offset_x = (sizex-2-sx)/2;
  int offset_y = (sizey-2-sy)/2;
  
  if (n_cells==1) 
  {
    return GrowInCells(1, cell_size, sizex/2, sizey/2, 0, 0);
  } else 
  {
    return GrowInCells(n_cells, cell_size, sx, sy, offset_x, offset_y);
  }
}

int CellularPotts::GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y) 
{
  // make initial cells using Eden Growth
  
  int **new_sigma=(int **)malloc(sizex*sizeof(int *));
  if (new_sigma==NULL)
    MemoryWarning();
  
  new_sigma[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (new_sigma[0]==NULL)  
    MemoryWarning();
  
  for (int i=1;i<sizex;i++) 
    new_sigma[i]=new_sigma[i-1]+sizey;
  
  /* Clear CA plane */
  { for (int i=0;i<sizex*sizey;i++) 
     new_sigma[0][i]=0; 
  }

  
  // scatter initial points, or place a cell in the middle 
  // if only one cell is desired
  int cellnum=cell->size()-1;

  if (n_cells>1) 
  {   
    { for (int i=0;i<n_cells;i++) 
    {
      sigma[RandomNumber(sx, s_val)+offset_x][RandomNumber(sy, s_val)+offset_y]=++cellnum;
    }}
  } 
  else 
  {
    sigma[sx+offset_x][sy+offset_y]=++cellnum;

  }


  if (par.eden_growth)
  // Do Eden growth for a number of time steps
  {for (int i=0;i<cell_size;i++) 
  {
    for (int x=1;x<sizex-1;x++)
      for (int y=1;y<sizey-1;y++) 
      {
	
        if (sigma[x][y]==0) {
          // take a random neighbour
          int xyp=(int)(8*RANDOM(s_val)+1);
          int xp = nx[xyp]+x;
          int yp = ny[xyp]+y;
          int kp;
          //  NB removing this border test yields interesting effects :-)
          // You get a ragged border, which you may like!
          if ((kp=sigma[xp][yp])!=-1)
            if (kp>(cellnum-n_cells))
              new_sigma[x][y]=kp;
            else
              new_sigma[x][y]=0;
          else
            new_sigma[x][y]=0;
          
        } else {
          new_sigma[x][y]=sigma[x][y];
        }
      }
    
    // copy sigma to new_sigma, but do not touch the border!
	{for (int x=1;x<sizex-1;x++) 
  {
    for (int y=1;y<sizey-1;y++) 
    {
	    sigma[x][y]=new_sigma[x][y];
    }
  }
  }}}
  else
  {
    double radius = sqrt(par.init_area / M_PI);

    // Iterate over the grid and fill the points within the circle
    for (int i = 0; i < sizex; ++i) {
        for (int j = 0; j < sizey; ++j) {
            // Calculate the distance from the center (x, y)
            double distance = sqrt(pow(i - sx-offset_x, 2) + pow(j - sy-offset_y, 2));

            // If the distance is less than or equal to the radius, mark the cell as part of the circle
            if (distance <= radius) {
                sigma[i][j] = 1;  // Mark cell inside the circle
            } else {
                sigma[i][j] = 0;  // Mark cell outside the circle
            }
        }
    }


  }



  free(new_sigma[0]);
  free(new_sigma);
  
  return cellnum;
}
  


double CellularPotts::CellDensity(void) const {
  
  // return the density of cells
  int sum=0;
  for (int i=0;i<sizex*sizey;i++) {
    if (sigma[0][i]) {
      sum++;
    }
  }
  return (double)sum/(double)(sizex*sizey);

}

double CellularPotts::MeanCellArea(void) const
{
  
  int sum_area=0, n=0;
  double sum_length=0.;
  vector<Cell>::iterator c=cell->begin(); ++c;
  
  for (; 
	c!=cell->end();
	c++) {
    
    sum_area+=c->Area();
    sum_length+=c->Length();
    n++;    
  }
  
  // cerr << "Mean cell length is " << sum_length/((double)n) << endl;
  return (double)sum_area/(double)n;
}

double CellularPotts::MeanCellPerimeter(void) const {
  
  int sum_perim=0, n=0;
  double sum_length=0.;
  vector<Cell>::iterator c=cell->begin(); ++c;
  
  for (; 
	c!=cell->end();
	c++) {
    
    sum_perim+=c->Perimeter();
    n++;    
  }
  return (double)sum_perim/(double)n;
}

void CellularPotts::ResetTargetLengths(void)  {
   vector<Cell>::iterator c=cell->begin(); ++c;

   for (;
        c!=cell->end();
        c++) {

     c->SetTargetLength(par.target_length);

} 

}

void CellularPotts::SetRandomTypes(void) {
  
  // each cell gets a random type 1..maxtau
  
  vector<Cell>::iterator c=cell->begin(); ++c;
  
  for (;c!=cell->end();c++) 
  {   
    c->setTau(1);
    c->set_ctype(2);
    c->SetTargetLength(0.0);
    
  } 
  
}

void CellularPotts::GrowAndDivideCells(int growth_rate) {

  vector<Cell>::iterator c=cell->begin(); ++c;
  vector<bool> which_cells(cell->size());

  for (;
       c!=cell->end();
       c++) {

    // only tumor cells grow and divide
    if (c->getTau()==2) {
     
      c->SetTargetArea(c->TargetArea()+growth_rate);
    
      if (c->Area()>par.init_area) {
	which_cells[c->Sigma()]=true;
      } else {
	which_cells[c->Sigma()]=false;
      }

      if (c->chem[1]<0.9) { //arbitrary oxygen threshold for the moment
	c->setTau(3);
      }
    } else {
      which_cells[c->Sigma()]=false;
    }

  }

  DivideCells(which_cells);

}

double CellularPotts::DrawConvexHull(Graphics *g, int color) {
  
  // Draw the convex hull of the cells
  // using Andrew's Monotone Chain Algorithm (see hull.cpp)

  // Step 1. Prepare data for 2D hull code
  
  // count number of points to determine size of array
  int np=0;
  for (int x=1;x<sizex-1;x++) 
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	np++;
      }
    }

  Point *p=new Point[np];
  
  int pc=0;
  for (int x=1;x<sizex-1;x++) 
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	p[pc++]=Point(x,y);
      }
    }

  // Step 2: call 2D Hull code
  Point *hull=new Point[np];
  int nph=chainHull_2D(p,np,hull);
  
  // Step 3: draw it
  for (int i=0;i<nph-1;i++) {
    g->Line(2*hull[i].x,2*hull[i].y,2*hull[i+1].x,2*hull[i+1].y, color);
  }

  
  // Step 4: calculate area of convex hull
  double hull_area=0.;
  for (int i=0;i<nph-1;i++) {
    hull_area+=hull[i].x*hull[i+1].y-hull[i+1].x*hull[i].y;
  }
  hull_area/=2.;

  //cerr << "Area = " << hull_area << "\n";
  
  delete[] p;
  delete[] hull;
  
  return hull_area;

}

double CellularPotts::Compactness(double *res_compactness, double *res_area, double *res_cell_area) {
  
  // Calculate compactness using the convex hull of the cells
  // We use Andrew's Monotone Chain Algorithm (see hull.cpp)

  // Step 1. Prepare data for 2D hull code
  
  // count number of points to determine size of array
  int np=0;
  for (int x=1;x<sizex-1;x++) 
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	np++;
      }
    }

  Point *p=new Point[np];
  
  int pc=0;
  for (int x=1;x<sizex-1;x++) 
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	p[pc++]=Point(x,y);
      }
    }

  // Step 2: call 2D Hull code
  Point *hull=new Point[np];
  int nph=chainHull_2D(p,np,hull);
  
  //// Step 3: draw it
  //for (int i=0;i<nph-1;i++) {
  //  g->Line(2*hull[i].x,2*hull[i].y,2*hull[i+1].x,2*hull[i+1].y, color);
  //}

  
  // Step 3: calculate area of convex hull
  double hull_area=0.;
  for (int i=0;i<nph-1;i++) {
    hull_area+=hull[i].x*hull[i+1].y-hull[i+1].x*hull[i].y;
  }
  hull_area/=2.;

  // Step 4: calculate total cell area
  double cell_area=0;

  vector<Cell>::const_iterator c;

  for ( (c=cell->begin(),c++);
       c!=cell->end();
       c++) {
    cell_area+=c->Area();
  }
  
  delete[] p;
  delete[] hull;


  // put intermediate results into optional pointers
  if (res_compactness) {
    *res_compactness = cell_area/hull_area;
  }
  if (res_area) {
    *res_area = hull_area;
  }
  if (res_cell_area) {
    *res_cell_area = cell_area;
  }

  // return compactness
  return cell_area/hull_area;

}




struct vec2d
{
  double x, y;
  vec2d(double x, double y) : x(x), y(y) {}

  // Magnitude of the vector
  double magnitude() const {
      return sqrt(x*x + y*y);
  }

  // Dot product of two vectors
  double dot(const vec2d& other) const {
      return x * other.x + y * other.y;
  }
};

// Comparator function to sort vectors clockwise
bool compareVec(const vec2d& v1, const vec2d& v2) {
    // Calculate angles
    double angle1 = atan2(v1.y, v1.x);
    double angle2 = atan2(v2.y, v2.x);

    // Return true if angle1 is less than angle2
    return angle1 < angle2; // Change to '<' for counter-clockwise
}


void CellularPotts::MeasureHexaticOrder()
{

  tmp_hex_order=0;
  int tmp_counter=0;
  int **ns = SearchNeighbours();
  int n_size = CountCells();
  for (int i = 1; i < n_size; ++i)
  {
    bool phaser = cell->at(i).GetPhase();
    bool isepi = cell->at(i).IsEpithelia();

    if (cell->at(i).AliveP() && phaser && !isepi)
    {
      double XCEN = cell->at(i).get_xcen();
      double YCEN = cell->at(i).get_ycen();
      vector<double> xcens{};
      vector<double> ycens{};
      int n_neighbours=0;
      bool med_check=false;
      int j = 0;
      while (ns[i][j] >= 0)
      {
        med_check = false;
        if (ns[i][j] == 0)
        {
          med_check=true;
          break;
        }
        else
        {
          double xc = cell->at(ns[i][j]).get_xcen();
          double yc = cell->at(ns[i][j]).get_ycen();
          xcens.push_back(xc);
          ycens.push_back(yc);
          ++n_neighbours;
        }
        ++j;
      }
      // cout << i << '\t' << n_neighbours << endl;

      if (med_check) 
        continue;

      vector<vec2d> com_vectors{};
      vec2d reference_axis(1.0,0.0);

      for (int n1 = 0; n1 < n_neighbours; ++n1)
      {
        double ABx = XCEN - xcens[n1];
        double ABy = YCEN - ycens[n1];
        vec2d newvec(ABx, ABy);
        com_vectors.push_back(newvec);
      }
      sort(com_vectors.begin(), com_vectors.end(), compareVec);
      // for (auto v : com_vectors)
      //   cout << v.x << '\t' << v.y << '\t';

      // cout << endl;
      vector<double> angles{};
      for (const auto& vec : com_vectors) 
      {
        double angle = atan2(vec.y, vec.x);
        angles.push_back(angle);
      }
      // Now use angles to calculate psi_6 for each particle
      complex<double> psi_sum(0,0);
      for (const auto& angle : angles) {
        psi_sum += std::exp(std::complex<double>(0, 6 * angle));
      }
      psi_sum /= static_cast<double>(angles.size());
      double psi_mag = std::abs(psi_sum);

      hexatic_tally += psi_mag;
      ++hexatic_counter;

    }
  }
  free(ns[0]);
  free(ns);
}

void CellularPotts::AverageHexaticOrder()
{
  if (hexatic_counter>0)
  {
    double to_add = hexatic_tally / hexatic_counter;
    hex_vec.push_back(to_add);
    hexatic_tally = 0;
    hexatic_counter = 0;
  }
  else
  {
    return hex_vec.push_back(-1);
  }
}

vector<double>& CellularPotts::ReturnHexaticOrder()
{
  return hex_vec;
}
