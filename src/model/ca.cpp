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
  



double CellularPotts::DeltaH(int x, int y, int sxyp, const int tsteps, const int* neighbor_spins, PDE *PDEfield)       
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

  for (i = 1; i <= n_nb_adh; i++) 
  {
    neighsite = neighbor_spins[i];
    if (neighsite == -1) 
    { // out-of-bounds border 
      Jen += (sxyp == 0 ? 0 : par.border_energy) - (sxy == 0 ? 0 : par.border_energy);
    } 
    else 
    {
      Jen += cell_sxyp.EmbryoEnergy((*cell)[neighsite], zona_sigma, zona_sigma_sticky) - cell_sxy.EmbryoEnergy((*cell)[neighsite], zona_sigma, zona_sigma_sticky);
    }
  }
  

  DH += Jen; // / (par.neigh_multiplier);


  // ==========================================
  // AREA CONSTRAINT
  // ==========================================
  double lambda = cell_sxy.get_lambda();
    
  if (par.medium_area_constraint)
  {
      DH += (int)((lambda * (2. + 2. * (double) 
              (cell_sxyp.Area() - cell_sxyp.TargetArea()
             - cell_sxy.Area()  + cell_sxy.TargetArea()))));
  }
  else
  {
    if (sxyp == MEDIUM) {
      DH += (int)(lambda * (1. - 2. * (double) (cell_sxy.Area() - cell_sxy.TargetArea())));
    }
    else if (sxy == MEDIUM) {
      DH += (int)((lambda * (1. + 2. * (double) (cell_sxyp.Area() - cell_sxyp.TargetArea()))));
    }
    else {
      DH += (int)((lambda * (2. + 2. * (double) 
              (cell_sxyp.Area() - cell_sxyp.TargetArea()
             - cell_sxy.Area()  + cell_sxy.TargetArea()))));
    }
  }


  // ==========================================
  // ACTIVE MOTION TERM
  // ==========================================
  if (par.active_motion && cell_sxy.Area() > 2)
  {
    double &mot_strength_sxy = cell_sxy.GetMotilityStrength();
    double &mot_strength_sxyp = cell_sxyp.GetMotilityStrength();
    if (sxyp == MEDIUM)
    {
      DH -= cell_sxy.GetMotilityStrength() * cell_sxy.ActiveDotProduct_removed(x, y);
    }
    else if (sxy == MEDIUM)
    {
      DH -= cell_sxyp.GetMotilityStrength() * cell_sxyp.ActiveDotProduct_added(x, y);
    }
    else
    {
      DH -= cell_sxyp.GetMotilityStrength() * cell_sxyp.ActiveDotProduct_added(x, y);
      DH -= cell_sxy.GetMotilityStrength() * cell_sxy.ActiveDotProduct_removed(x, y);
    }
    // if ((*cell)[sigma[x][y]].TargetArea()==0 && sigma[x][y] > 0)
    //   cout << cell_sxy.GetMotilityStrength() * cell_sxy.ActiveDotProduct_removed(x, y) << endl;
  }


  // ==========================================
  // PERIMETER CONSTRAINT
  // ==========================================
  if (par.H_perim && cell_sxy.Area() > 2) 
  {
    double DH_perimeter = 0;
    if (sxyp == MEDIUM) 
    {
      DH_perimeter -= cell_sxy.getPerimConstraint() *
          (DSQR(cell_sxy.Perimeter() - cell_sxy.TargetPerimeter()) -
           DSQR(GetNewPerimeterIfXYWereRemoved(sxy, x, y, neighbor_spins) - cell_sxy.TargetPerimeter()));      
    } 
    else if (sxy == MEDIUM) 
    {
      DH_perimeter -= cell_sxyp.getPerimConstraint() *
          (DSQR(cell_sxyp.Perimeter() - cell_sxyp.TargetPerimeter()) -
           DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y, neighbor_spins) - cell_sxyp.TargetPerimeter()));      
    }
    else 
    {
      // they're both cells
      DH_perimeter -= cell_sxyp.getPerimConstraint() *
          ((DSQR(cell_sxyp.Perimeter() - cell_sxyp.TargetPerimeter()) -
            DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y, neighbor_spins) - cell_sxyp.TargetPerimeter())));
            
      DH_perimeter -= cell_sxy.getPerimConstraint() *
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



int CellularPotts::AmoebaeMove(long tsteps, PDE *PDEfield)
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
          if (tx <= 0 || ty <= 0 || tx >= sizex - 1 || ty >= sizey - 1) 
          {
              neighbor_spins[j] = -1; // -1 means out-of-bounds border
          } 
          else 
          {
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

    if (distinct_count == 1 && present_states[0] == k)
      continue; 

    
    int rand_idx = (int)(distinct_count * RANDOM(s_val));
    int kp = present_states[rand_idx];

    if (k == kp)
      continue;

    if (par.make_zona_pellucida && (k==zona_sigma || kp==zona_sigma || k==zona_sigma_sticky || kp==zona_sigma_sticky))
    {
      continue;
    }
    // if (neighbor_spins[15] == zona_sigma)
    // {
    //   cout << "zona found: " << zona_sigma << endl;
    // }

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
      double D_H = DeltaH(x, y, kp, tsteps, neighbor_spins, PDEfield);

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

int CellularPotts::AmoebaeMoveLegacy(long tsteps, PDE *PDEfield)
{
  int loop,p;
  //int updated=0;
  thetime++;
  int SumDH=0;
  
  if (frozen) 
    return 0;

  const int sx_inner = sizex - 2;
  const int sy_inner = sizey - 2;

  int max_nb = n_nb_adh;
  int neighbor_spins[64];

  loop=(sizex-2)*(sizey-2);
 
  for (int i=0;i<loop;i++) 
  {  
    // take a random site
    int xy = (int)(RANDOM(s_val)*(sizex-2)*(sizey-2));
    int x = xy%(sizex-2)+1;
    int y = xy/(sizex-2)+1; 
    int k=sigma[x][y];
    

    // get all neighbours for adhesion/perim etc.
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

    // take a random neighbour
    int xyp=(int)(n_nb*RANDOM(s_val)+1);
    int xp = nx[xyp]+x;
    int yp = ny[xyp]+y;    
    
    int kp;
    if (par.periodic_boundaries) 
    {
      // since we are asynchronic, we cannot just copy the borders once 
      // every MCS
      if (xp<=0)
	      xp=sizex-2+xp;
      if (yp<=0)
	      yp=sizey-2+yp;
      if (xp>=sizex-1)
	      xp=xp-sizex+2;
      if (yp>=sizey-1)
	      yp=yp-sizey+2;
      
      kp=sigma[xp][yp];
      
    } 
    else 
    {
      if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1)
	      kp=-1;
      else
	      kp=sigma[xp][yp];
    }
    // int type1 = (*cell)[sigma[xp][yp]].GetPhenotype();    
    // int type2 = (*cell)[sigma[xp][yp]].GetPhenotype();    

    // test for border state (relevant only if we do not use 
    // periodic boundaries)
    if (kp!=-1) 
    {  
      // Don't even think of copying the special border state into you!
    
      if ( k  != kp ) 
      {
        /* Try to copy if sites do not belong to the same cell */
        // connectivity dissipation:
        int H_diss=0;
        if (!ConnectivityPreservedP(x,y)) 
          H_diss=par.conn_diss;
        
        double D_H=DeltaH(x,y,kp, tsteps, neighbor_spins);
        
        // dH_tally += D_H;
        // if ((type1 > par.mintype && type1 < par.maxtype) || (type2 > par.mintype && type2 < par.maxtype))
        //   cout << D_H << endl;
        // bool is_med_attempt = false;
        // if (sigma[x][y] == 0 && (*cell)[sigma[xp][yp]].GetPhase() == true || sigma[xp][yp] == 0 && (*cell)[sigma[x][y]].GetPhase() == true)
        // {
        //   is_med_attempt = true;
        //   ++medp_count;
        // }
        if ((p=CopyvProb(D_H,H_diss))>0) 
        {
          if (par.H_perim)
            ConvertSpinPerim( x,y,kp, neighbor_spins );
          else
          {
            ConvertSpin( x,y,kp );
          }  
        }
      }
    } 
  }
  return SumDH;
  
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


// double CellularPotts::CalculateABBoundaryLength()
// {
//     int total_ab_boundary = 0;
//     int total_boundary = 0; // Track the total cell-cell boundary length

//     // Iterate through the entire grid
//     for (int i = 0; i < sizex; i++) {
//         for (int j = 0; j < sizey; j++) {
            
//             int current_id = sigma[i][j];
            
//             // Skip if the current pixel is empty medium (<= 0)
//             if (current_id <= 0) continue; 
            
//             // Get the type of the current cell (false = A, true = B)
//             bool current_type = (*cell)[current_id].GetSortingType();

//             // 1. Check the RIGHT neighbor
//             if (i + 1 < sizex) {
//                 int right_id = sigma[i + 1][j];
                
//                 // Check if right pixel is a valid cell and is a different cell
//                 if (right_id > 0 && current_id != right_id) {
//                     total_boundary++; // It's a cell-cell boundary
                    
//                     bool right_type = (*cell)[right_id].GetSortingType();
                    
//                     // If the types are different (one is A, one is B)
//                     if (current_type != right_type) {
//                         total_ab_boundary++;
//                     }
//                 }
//             }

//             // 2. Check the BOTTOM neighbor
//             if (j + 1 < sizey) {
//                 int bottom_id = sigma[i][j + 1];
                
//                 // Check if bottom pixel is a valid cell and is a different cell
//                 if (bottom_id > 0 && current_id != bottom_id) {
//                     total_boundary++; // It's a cell-cell boundary
                    
//                     bool bottom_type = (*cell)[bottom_id].GetSortingType();
                    
//                     // If the types are different (one is A, one is B)
//                     if (current_type != bottom_type) {
//                         total_ab_boundary++;
//                     }
//                 }
//             }
            
//         }
//     }

//     // Prevent division by zero if there are no cell boundaries at all
//     if (total_boundary == 0) {
//         return 0.0;
//     }

//     // Return the relative boundary length as a double between 0.0 and 1.0
//     return static_cast<double>(total_ab_boundary) / total_boundary;
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

            daughterp->SetTimeCreated(t);
            motherp->SetTimeCreated(t);

            
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
    cell->back().SetTimeCreated(time);

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
}\


/* putting new methods here */

void CellularPotts::InitialiseRandomSoxValues()
{
  // --- NEW CONFIGURATION ---
  // Define your desired probability for a cell to be Sox2 dominant.
  // 0.8 means 80% Sox2, 20% Sox17. 0.5 is 50/50.
  // (You could easily make this a parameter like par.sox2_ratio)

  // Clamp the probability between 0.001 and 0.999 to prevent log(0) math errors
  if (par.target_sox2_prob < 0.001) par.target_sox2_prob = 0.001;
  if (par.target_sox2_prob > 0.999) par.target_sox2_prob = 0.999;
  
  double target_sox17_prob = 1.0 - par.target_sox2_prob;

  // The threshold upon which cell fate is decided. 
  // Based on your previous p=2.321928 math, this is 0.2.
  // (Ideally, set this to par.sox_threshold if it's accessible here)
  double threshold = 0.2; 

  // Mathematically calculate the exact exponents needed to shift the distributions
  // so that exactly 'target_sox2_prob' proportion of cells land above the threshold.
  double p2  = std::log(threshold) / std::log(1.0 - par.target_sox2_prob);
  double p17 = std::log(threshold) / std::log(1.0 - target_sox17_prob);
  // -------------------------

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP() && c->Sigma()!=zona_sigma)
    {
      // 1. Get two independent uniform random numbers between 0 and 1
      double u_a = RANDOM(s_val);
      double u_b = RANDOM(s_val);
      if (u_a <= 0.0) u_a = 0.0000001; 

      // 2. Box-Muller transform
      double radius = std::sqrt(-2.0 * std::log(u_a));
      double theta = 2.0 * M_PI * u_b;
      double z1 = radius * std::cos(theta);
      double z2 = radius * std::sin(theta);

      // 3. DEFINE YOUR 'x' HERE (Fraction of cases split across 0.2)
      // This maintains the mutual exclusion (if it's Sox2, it's NOT Sox17)
      double x = 1-par.starting_fraction_losers; 
      double rho = std::sin(M_PI * (0.5 - x)); 

      // Create negatively correlated variables
      double y1 = z1;
      double y2 = rho * z1 + std::sqrt(1.0 - rho * rho) * z2;

      // 4. Convert back to Uniform(0, 1)
      const double inv_sqrt2 = 0.7071067811865475;
      double u1 = 0.5 * std::erfc(-y1 * inv_sqrt2);
      double u2 = 0.5 * std::erfc(-y2 * inv_sqrt2);

      // 5. Apply the calculated skew exponents to achieve your desired ratio
      double sox2 = std::pow(u1, p2);
      double sox17 = std::pow(u2, p17);
      
      c->setSox2(sox2);
      c->setSox17(sox17);

      c->SetSoxColour();
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

void CellularPotts::PopulateDenseCellsInZonaRadius(double density, double R, int shiftx, int shifty, double h, double k, double a, double b, double n)
{
  int current_cells = CountCells();
  // 1. Define the usable active grid space
  int W = sizex - 2;
  int H = sizey - 2;
  
  double max_area = par.synthetic_max_area;
  double min_area = par.synthetic_min_area;
  double average_area = (max_area - min_area)/2. + min_area;
  
  // Guard against invalid parameters
  if (max_area <= 0 || density <= 0.0 || R <= 0.0 || a <= 0.0 || b <= 0.0) return;
  
  // 2. Calculate lattice spacing directly from density and target area.
  double area_per_center = average_area / density;
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

  // Offset ensures perfectly symmetric empty padding near the walls, shifted by user request
  double offset_x = (W - max_x) / 2.0 + 1.0 + shiftx; 
  double offset_y = (H - max_y) / 2.0 + 1.0 + shifty;

  // Find the exact mathematical center of the grid
  double grid_center_x = sizex / 2.0 + shiftx;
  double grid_center_y = sizey / 2.0 + shifty;

  // HELPER: Check if a point is strictly INSIDE the Zona Pellucida hollow cavity
  auto isInsideZonaCavity = [&](double px, double py) {
      double dx = px - h;
      double dy = py - k;
      
      // Ellipse equation f(x, y)
      double f = (dx * dx) / (a * a) + (dy * dy) / (b * b) - 1.0;
      
      // If f >= 0, we are on or outside the mathematical ellipse boundary
      if (f >= 0.0) return false; 
      
      // Calculate approximated distance to the ellipse curve
      double grad_x = (2.0 * dx) / (a * a);
      double grad_y = (2.0 * dy) / (b * b);
      double grad_mag = std::sqrt(grad_x * grad_x + grad_y * grad_y);
      
      double distance = (grad_mag == 0.0) ? std::min(a, b) : std::abs(f) / grad_mag;
      
      // Must be deeper inside the boundary than the zona thickness 'n'
      return distance > n; 
  };

  // 5. Generate the centers, filtering by BOTH grid center <= R AND ellipse cavity
  struct VPoint { 
    double x, y; 
    int id; 
    double target_area;
  };
  std::vector<VPoint> centers;

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      double cx = col * 2.0 * r;
      double cy = row * std::sqrt(3.0) * r;
      if (row % 2 == 1) cx += r;
      
      double final_cx = cx + offset_x;
      double final_cy = cy + offset_y;
      
      // Distance from the shifted grid center
      double dist = euclideanDistance(final_cx, final_cy, grid_center_x, grid_center_y, sizex, sizey);
      
      // Must be within radius R AND inside the Zona Pellucida cavity
      if (dist <= R && isInsideZonaCavity(final_cx, final_cy)) 
      {
        double rnd = RANDOM(s_val) * (max_area - min_area) + min_area;
        centers.push_back({final_cx, final_cy, -1, rnd});
      }
      if (final_cx < 5 || final_cx > sizex-5 || final_cy < 5 || final_cy > sizey-5)
        cerr << "warning: some centers are outside of domain\n";
    }
  }

  int final_ncells = centers.size();
  if (final_ncells == 0) {
      std::cout << "Conditions too restrictive to generate any cells (Check R and Ellipse bounds)." << std::endl;
      return;
  }

  // 6. Split sheet to prepare sufficient cell instances
  FractureSheet(final_ncells);

  // Map the new alive cell IDs to our spatial centers
  std::vector<int> valid_ids;
  for (auto c = std::next(cell->begin(), current_cells); c != cell->end(); ++c) {
    if (c == cell->begin()) continue; // Skip medium/background index 0
    if (c->AliveP()) 
    {
      valid_ids.push_back(c->Sigma());
    }
  }

  for (size_t i = 0; i < centers.size(); ++i) {
      if (i < valid_ids.size()) {
          centers[i].id = valid_ids[i];
      }
  }

  // 8. Draw the Voronoi domains
  for (int x = 1; x < sizex - 1; ++x) {
      for (int y = 1; y < sizey - 1; ++y) {
        
        // Safety check: Never overwrite the existing Zona Pellucida
        if (sigma[x][y] == zona_sigma) {
          continue;
        }

        // Distance from the shifted grid center
        double dist_to_grid_center = euclideanDistance(x, y, grid_center_x, grid_center_y, sizex, sizey);
        
        // Pixel must be within radius R AND strictly inside the hollow cavity
        if (dist_to_grid_center <= R && isInsideZonaCavity(x, y)) 
        {
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
            
          // Assign to closest center
          if (closestCenter != -1) {
            sigma[x][y] = closestCenter;
          }
        }
      }
  }

  // 9. Re-evaluate actual populated cell areas
  for (auto c = cell->begin(); c != cell->end(); ++c) {
    if (c == cell->begin()) 
      continue;
    if (c->AliveP() && c->Sigma() != zona_sigma) 
      c->area = 0;
  }

  for (int x = 1; x < sizex - 1; ++x) {
    for (int y = 1; y < sizey - 1; ++y) {
      if (sigma[x][y] > 0 && sigma[x][y] != zona_sigma) 
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
    if (c->Sigma() == zona_sigma) continue; // Skip Zona Pellucida cell

    if (c->AliveP()) 
    {
      if (c->area == 0) 
      {
        c->Apoptose();
        ++deadcells;
      } 
      else 
      {
        // Target area is now matched identically to the size the Voronoi algorithm generated.
        c->SetTargetArea(c->area); 
        c->setSpheroid(false);
      }
    }
  }
  MeasureCellSizes();

  std::cout << "Grid generated | Radius: " << R 
            << " | Density: " << density 
            << " | Cells populated: " << (final_ncells - deadcells)
            << " | Cells killed (0 area): " << deadcells << std::endl;
}


// Added 'double R' to the function parameters
void CellularPotts::PopulateSparseCells(double density, double R, int shiftx, int shifty)
{

  int current_cells = CountCells();
  // 1. Define the usable active grid space
  int W = sizex - 2;
  int H = sizey - 2;
  
  double max_area = par.synthetic_max_area;
  double min_area = par.synthetic_min_area;
  double average_area = (max_area - min_area)/2. + min_area;
  
  // Guard against invalid parameters
  if (max_area <= 0 || density <= 0.0 || R <= 0.0) return;
  
  // 2. Calculate lattice spacing directly from density and target area.
  // This ensures the density inside the radius is exactly what was requested.
  double area_per_center = average_area / density;
  
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
  struct VPoint 
  { double x, y; 
    int id; 
    double target_area;
    double radius_limit;
  };
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
      
      if (dist <= R) 
      {
        double rnd = RANDOM(s_val) * (max_area - min_area) + min_area;
        double r_limit = std::sqrt(rnd / M_PI) * 1.05; // 5% buffer applied individually
        centers.push_back({final_cx, final_cy, -1, rnd, r_limit});
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
    if (c->AliveP()) 
    {
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

  for (int x = 1; x < sizex - 1; ++x) {
      for (int y = 1; y < sizey - 1; ++y) {
        if (sigma[x][y] == zona_sigma) {
          continue;
        }

        double minDistance = std::numeric_limits<double>::max();
        int closestCenter = -1;
        double closestRadiusLimit=0.;
        
        for (const auto& center : centers) {
          if (center.id == -1) continue;
          
          double dist = euclideanDistance(x, y, center.x, center.y, sizex, sizey);
          if (dist < minDistance) {
              minDistance = dist;
              closestCenter = center.id;
              closestRadiusLimit=center.radius_limit;
          }
        }
          
        if (minDistance < closestRadiusLimit) 
        {
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
      if (sigma[x][y] > 0 && sigma[x][y] != zona_sigma) 
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
        c->setSpheroid(false);
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
      if (c->AliveP() && c->Sigma() != zona_sigma)
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
      if (c->AliveP() && c->Sigma() != zona_sigma)
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
	c++) 
  {
    // cout << c->Perimeter() << " and target is.. " << c->TargetPerimeter() << endl;
    sum_perim+=c->Perimeter();
    n++;    
  }
  return (double)sum_perim/(double)n;
}


void CellularPotts::SetRandomTypes(void) {
  
  // each cell gets a random type 1..maxtau
  
  vector<Cell>::iterator c=cell->begin(); ++c;
  
  for (;c!=cell->end();c++) 
  {   
    c->setTau(1);
    c->set_ctype(2);    
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


void CellularPotts::CellGrowthAndDivision(int time) 
{
  vector<bool> which_cells(cell->size());

  int cell_division=0;
  
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {

      int TA = c->TargetArea();
      int area = c->Area();
      int gthresh = par.gthresh;
      int sthresh=par.shrink;


      if ( (area-TA)>gthresh) // && area <= (double)(par.div_threshold) * 1.1) //  
      {
        int count= area-TA; //area-TA;
        while (count>0)
        {
          c->IncrementTargetArea();
          --count;
        }
      }
      else if ( (area-TA)<sthresh ) 
      {
        int count=TA-area;
        while (c->TargetArea() > 0 && count > 0)
        {
          c->DecrementTargetArea();
          --count;
        }
      }
      else if (area < 3)
      {
        c->SetTargetArea(0);
        c->set_lambda(100);
      }
      // else if (area > (double)(par.div_threshold) * 1.1)
      // {
      //   c->DecrementTargetArea();
      // }   
      
      if (area>par.div_threshold) // && c->checkforcycles(par.cycle_threshold) == false)
      {
        
        which_cells[c->Sigma()]=true;
        cell_division++;
      }
    }
  }
  // Divide scheduled cells
  if (cell_division) 
  {
    DivideCells(which_cells, time);
  }
  //Function that partitions the TF by polarity (mother, daughter, neither) is in divide cells  
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


// Function to plot a thick ellipse on a 2D grid
void CellularPotts::MakeZonaPellucida(double h, double k, double a, double b, double n) 
{                    
  // Basic bounds checking to ensure the grid exists

  //make room for a new cell!
  vector<bool> which_cells(cell->size());
  which_cells.back()=true;
  DivideCellsNoGrid(which_cells);
  zona_sigma = (*cell).back().Sigma();
  (*cell)[zona_sigma].set_ctype(1);
  
  double a2 = a * a;
  double b2 = b * b;
  int total_area=0;
  
  // Step 2: Iterate through the bounding box
  for (int x = 1; x <= sizex-1; ++x) 
  {
    for (int y = 1; y <= sizey-1; ++y) 
    {
      double dx = x - h;
      double dy = y - k;
      
      // The Ellipse Equation f(x, y)
      double f = (dx * dx) / a2 + (dy * dy) / b2 - 1.0;
      
      // The Gradient Magnitude |∇f(x,y)|
      double grad_x = (2.0 * dx) / a2;
      double grad_y = (2.0 * dy) / b2;
      double grad_mag = std::sqrt(grad_x * grad_x + grad_y * grad_y);
      
      // Step 3: Calculate approximated distance to the ellipse curve
      double distance;
      if (grad_mag == 0.0) {
          // Edge case: if we are exactly at the center (h, k), avoid dividing by zero.
          // The distance to the ellipse from the center is the minor axis length.
          distance = std::min(a, b); 
      } else {
          distance = std::abs(f) / grad_mag;
      }
      
      // Step 4: If the pixel is within distance 'n', modify its value in sigma
      if (distance <= n) 
      {
        sigma[x][y] = zona_sigma;
        ++total_area=0;
      }
    }
  }
  (*cell)[zona_sigma].SetTargetArea(total_area);
  (*cell)[zona_sigma].Apoptose();
}


void CellularPotts::DifferentiateZonaPellucida()
{
  //make room for a new cell!
  vector<bool> which_cells(cell->size());
  which_cells.back()=true;
  DivideCellsNoGrid(which_cells);
  zona_sigma_sticky = (*cell).back().Sigma();
  int total_area=0;
  (*cell)[zona_sigma_sticky].set_ctype(203);
  // Store modifications to prevent a chain-reaction in a single pass
  vector<std::pair<int, int>> to_change;
  int R = 4;
  // Step 2: Iterate through the grid
  for (int x = 1; x <= sizex-1; ++x) 
  {
    for (int y = 1; y <= sizey-1; ++y) 
    {
      if (sigma[x][y] == zona_sigma)
      {
        bool found = false;
        // Search neighboring pixels within the bounding box of radius R
        // std::max/std::min ensures we don't check outside the grid boundaries
        for (int nx = max(1, x - R); nx <= min(sizex - 1, x + R); ++nx) 
        {
          for (int ny = max(1, y - R); ny <= min(sizey - 1, y + R); ++ny) 
          {
            // Check if the neighbor is actually within the circular radius R
            if ((nx - x) * (nx - x) + (ny - y) * (ny - y) <= R * R) 
            {
              if (sigma[nx][ny] > 0 && sigma[nx][ny] != zona_sigma) 
              {
                found = true;
                break; // Stop searching if we found at least one
              }
            }
          }
          if (found) break; // Break out of the outer neighbor loop early
        }
        
        // If a valid pixel was found in radius R, mark this coordinate to be changed
        if (found) 
        {
          cout << x << '\t' << y << endl;
          pair<int,int>newchange={x,y};
          to_change.push_back(newchange);
          ++total_area;
        } 
      }
      
 
    }
  }
  
  for (auto pp : to_change)
  {
    sigma[pp.first][pp.second] = zona_sigma_sticky;
  }

  (*cell)[zona_sigma_sticky].SetTargetArea(total_area);
  (*cell)[zona_sigma_sticky].Apoptose();
}





void CellularPotts::SetMotilityStrengths()
{

  vector<Cell>::iterator c;

  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP()) 
    {
      c->SetMotilityStrength(par.motility_strength);
    }
  }
}




void CellularPotts::ToxictoLonelyCells()
{
  int **ns = SearchNeighbours();
  int n_size = CountCells();
  for (int i = 1; i < n_size; ++i)
  {
    if (true==true)//(cell->at(i).AliveP())
    {
      int nbh_count{};
      int j=0;
      int non_cell_count=0;
      while (ns[i][j] >= 0)
      {
        if (ns[i][j] != zona_sigma && ns[i][j] !=zona_sigma_sticky && ns[i][j] > 0)
          ++non_cell_count;
        ++j;
      }
      if (non_cell_count == 0)
      {
        cell->at(i).MakeLonely(true);
        if (cell->at(i).Area() < 50 && cell->at(i).TargetArea() > 0)
          cell->at(i).SetTargetArea(1);
        else if (cell->at(i).TargetArea() > 1)
          cell->at(i).DecrementTargetArea();

        double area_constraint = par.bulk_modulus / double(cell->at(i).TargetArea());
        cell->at(i).setAreaConstraint(area_constraint);
        int target_perim = round(double(par.ptarget_perimeter) * sqrt(double(cell->at(i).TargetArea())/double(par.cell_target_area)));
        if (target_perim < 2)
          target_perim=0;
        cell->at(i).SetTargetPerimeter(target_perim);
        
        double perim_constraint = (cell->at(i).GetElasticMod() / double(target_perim));
        cell->at(i).setPerimConstraint(perim_constraint);

        // cout << "cell number: " << cell->at(i).Sigma() << "  area: " << cell->at(i).Area() << '\t' << "  target area: " << cell->at(i).TargetArea() << "  perimeter: " << cell->at(i).Perimeter() << "   target perimeter: " << cell->at(i).TargetPerimeter() << endl;

      }
      else
      {
        cell->at(i).MakeLonely(false);
      }
    }
  }

  free(ns[0]);
  free(ns);

}

// still working on this
void CellularPotts::NeighbourBasedPerimeterConstraint()
{

  int **ns = SearchNeighbours();
  int n_size = CountCells();
  vector<double> neighbour_sox2_vals(n_size,0);
  vector<double> neighbour_sox17_vals(n_size,0);
  vector<int> cell_nbh_counts(n_size,0);
  for (int i = 1; i < n_size; ++i)
  {
    if (cell->at(i).AliveP())
    {
      int nbh_count{};
      int j=0;
      bool touching_med=false;
      // sox17 should maybe depend on blastocoel touching (do later)
      while (ns[i][j] >= 0)
      {
        if (ns[i][j] > 0)
        {
          ++cell_nbh_counts[i];
          neighbour_sox2_vals[i]+=(*cell)[ns[i][j]].getSox2adhesion();
          neighbour_sox17_vals[i]+=(*cell)[ns[i][j]].getSox17adhesion();
        }
        if (ns[i][j]==0)
          touching_med=true;
        ++j;
      }
      double stiffness_multiplier=1;
      if (cell_nbh_counts[i] > 0)
      {
        neighbour_sox2_vals[i] /= double(cell_nbh_counts[i]);
        neighbour_sox17_vals[i] /= double(cell_nbh_counts[i]);
        double p1 = neighbour_sox2_vals[i] * cell->at(i).getSox2adhesion() * 20;
        double p2 = neighbour_sox17_vals[i] * cell->at(i).getSox17adhesion() * 20;
        stiffness_multiplier += (p1+p2);
        cout << p1 << '\t' << p2 << endl;
      }
      if (touching_med==true)
      {
        double mot_strength = par.motility_strength - par.motility_strength * cell->at(i).getSox17adhesion();
        cell->at(i).SetMotilityStrength(mot_strength);
      }
      else
      {
        double mot_strength = par.motility_strength;
        cell->at(i).SetMotilityStrength(mot_strength);
      }
      int target_perim = cell->at(i).TargetPerimeter();
      double ideal_perim_constraint = (par.elastic_modulus / double(target_perim)) * stiffness_multiplier;

      double current_constraint = cell->at(i).getPerimConstraint();
      double smoothed_constraint = (current_constraint * 0.9) + (ideal_perim_constraint * 0.1);

      cell->at(i).setPerimConstraint(smoothed_constraint);
    }
  }

  free(ns[0]);
  free(ns);
}


















/* EVERYTHING BELOW IS SYNTHETIC NETWORK RELATED MAY BE USEFUL*/

double synNotch_bound_derivative(double c, double L, double V)
{
  return par.production_rate_synNotch - (par.decay_synNotch_bound*c) - par.binding_rate_CD19_synNotch * L * c - c/(3*V);
}

double synNotch_bound_rk4(double dt, double c, double L, double V)
{
  // the CD19 ligand is incorporated in L
  double k1 = synNotch_bound_derivative(c, L, V);
  double k2 = synNotch_bound_derivative(c + dt * k1/2.0, L, V);
  double k3 = synNotch_bound_derivative(c + dt*k2/2.0, L, V);
  double k4 = synNotch_bound_derivative(c + dt*k3, L, V);
  return c + (dt/6.0) * (k1  + 2.0*k2 + 2.0*k3 + k4);
}

double synNotch_unbound_derivative(double c, double cB, double L, double V)
{
  return par.binding_rate_CD19_synNotch * L * cB - par.decay_synNotch_unbound * c - c/(3*V);
}

double synNotch_unbound_rk4(double dt, double c, double cB, double L, double V)
{
  double k1 = synNotch_unbound_derivative(c, cB, L, V);
  double k2 = synNotch_unbound_derivative(c + dt * k1/2.0, cB, L, V);
  double k3 = synNotch_unbound_derivative(c + dt*k2/2.0, cB, L, V);
  double k4 = synNotch_unbound_derivative(c + dt*k3, cB, L, V);
  return c + (dt/6.0) * (k1  + 2.0*k2 + 2.0*k3 + k4);
}

double synNotch_intra_derivative(double c, double cB, double L, double V)
{
  return par.binding_rate_CD19_synNotch * L  * cB - par.decay_synNotch_intra * c - c/(3*V);
}


double synNotch_intra_rk4(double dt, double c, double cB, double L, double V)
{
  double k1 = synNotch_intra_derivative(c, cB, L, V);
  double k2 = synNotch_intra_derivative(c + dt * k1/2.0, cB, L, V);
  double k3 = synNotch_intra_derivative(c + dt*k2/2.0, cB, L, V);
  double k4 = synNotch_intra_derivative(c + dt*k3, cB, L, V);
  return c + (dt/6.0) * (k1  + 2.0*k2 + 2.0*k3 + k4);
}


double E_cadherin_derivative(double c, double I, double X, double prate, double V)
{
  // X is the proportion shared surface with cells also expressing E_cadherin (have to normalise to so that peak concentration = 1)
  return prate * (pow(I, par.hill_coefficient)/ (pow(par.E_cadherin_saturation_constant, par.hill_coefficient) + pow(I, par.hill_coefficient))) * (1. - c/par.c_max)
  - par.decay_E_cadherin_bound * c * X - par.decay_E_cadherin_unbound * c * (1-X) - c/(3*V);
}


double E_cadherin_rk4(double dt, double c, double I, double X, double prate, double V)
{
  double k1 = E_cadherin_derivative(c, I, X, prate, V);
  double k2 = E_cadherin_derivative(c + dt * k1/2.0, I, X, prate, V);
  double k3 = E_cadherin_derivative(c + dt*k2/2.0, I, X, prate, V);
  double k4 = E_cadherin_derivative(c + dt*k3, I, X, prate, V);
  return c + (dt/6.0) * (k1  + 2.0*k2 + 2.0*k3 + k4);
}

double random_binding_derivative(double c, double X, double V)
{
  return par.random_binding_protein_production
  - par.decay_random_binding_protein_bound * c * X - par.decay_random_binding_protein_unbound * c * (1-X) - c/(3*V);
}

double random_binding_rk4(double dt, double c, double X, double V)
{
  double k1 = random_binding_derivative(c, X, V);
  double k2 = random_binding_derivative(c + dt * k1/2.0, X, V);
  double k3 = random_binding_derivative(c + dt*k2/2.0, X, V);
  double k4 = random_binding_derivative(c + dt*k3, X, V);
  return c + (dt/6.0) * (k1  + 2.0*k2 + 2.0*k3 + k4);
}

double GFP_derivative(double c, double I, double V)
{
  return par.GFP_production_rate * (pow(I, par.hill_coefficient)/ (pow(par.E_cadherin_saturation_constant, par.hill_coefficient) + pow(I, par.hill_coefficient))) * (1. - c/par.c_max) - par.decay_GFP * c - c/(3*V);
}

double GFP_rk4(double dt, double c, double I, double V)
{
  double k1 = GFP_derivative(c, I, V);
  double k2 = GFP_derivative(c + dt * k1/2.0, I, V);
  double k3 = GFP_derivative(c + dt*k2/2.0, I, V);
  double k4 = GFP_derivative(c + dt*k3, I, V);
  return c + (dt/6.0) * (k1  + 2.0*k2 + 2.0*k3 + k4);
}


void CellularPotts::UpdateActiveMotion()
{

  int **ns = SearchNeighbours();
  int n_size = CountCells();
  for (int i = 1; i < n_size; ++i)
  {
    if (cell->at(i).AliveP())
    {
      double Ecad_conc = (*cell)[i].getE_cadherin();
      double Ncad_conc = (*cell)[i].getN_cadherin();
      double Pcad_conc = (*cell)[i].getP_cadherin();

      double E_avg_of_neighbours{};
      double N_avg_of_neighbours{};
      double P_avg_of_neighbours{};

      int nbh_count{};
      int j=0;
      while (ns[i][j] >= 0)
      {
        ++nbh_count;
        if (ns[i][j] == 0)
        {
          ++j;
          continue;
        }
        E_avg_of_neighbours += (*cell)[j].getE_cadherin();
        N_avg_of_neighbours += (*cell)[j].getN_cadherin();
        P_avg_of_neighbours += (*cell)[j].getP_cadherin();
        ++j;
      }
      if (E_avg_of_neighbours > 0)
        E_avg_of_neighbours/=double(nbh_count);
      if (N_avg_of_neighbours > 0)
        N_avg_of_neighbours/=double(nbh_count);
      if (P_avg_of_neighbours > 0)
        P_avg_of_neighbours/=double(nbh_count);

      // We assume that concentrations max out at 1.. hope this is okay...
      double Epart1 = Ecad_conc;
      double Epart2 = E_avg_of_neighbours;
      if (Epart1 > 1)
        Epart1=1;
      if (Epart2 > 1)
        Epart2=1;
      double Emot_strength = par.motility_strength - par.Ecadherin_bound_motility_loss * (Epart1 * Epart2);
      double Enew_elastic = par.elastic_modulus + par.Ecad_elastic_change * Epart1 * Epart2;


      double Npart1 = Ncad_conc;
      double Npart2 = N_avg_of_neighbours;
      if (Npart1 > 1)
        Npart1=1;
      if (Npart2 > 1)
        Npart2=1;
      double Nmot_strength = par.motility_strength - par.Ncadherin_bound_motility_loss * (Npart1 * Npart2);
      double Nnew_elastic = par.elastic_modulus + par.Ncad_elastic_change * Npart1 * Npart2;

      double Ppart1 = Pcad_conc;
      double Ppart2 = P_avg_of_neighbours;
      if (Ppart1 > 1)
        Ppart1=1;
      if (Ppart2 > 1)
        Ppart2=1;
      double Pmot_strength = par.motility_strength - par.Pcadherin_bound_motility_loss * (Ppart1 * Ppart2);
      double Pnew_elastic = par.elastic_modulus + par.Pcad_elastic_change * Ppart1 * Ppart2;



      double mot_strength = std::min({Nmot_strength, Emot_strength, Pmot_strength});
      double new_elastic = std::max({Enew_elastic, Nnew_elastic, Pnew_elastic});
      
      // if (mot_strength < 0.3)
      //   cout << Ecad_conc << '\t' << mot_strength << '\t' << new_elastic << endl;

      (*cell)[i].SetMotilityStrength(mot_strength);

      (*cell)[i].SetElasticMod(new_elastic);
      nbh_count=0;
    }
  }
  free(ns[0]);
  free(ns);
}


void CellularPotts::SurfaceBindings()
{
  // reset values
  // vector<Cell>::iterator c;
  // for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  // {
  //   if (c->AliveP())
  //   {
  //     c->ResetTempSurfaceBindings();
  //     c->setTouchingMed(false);
  //   }
  // }

  for (int x = 1; x < sizex - 1; x++) 
  {
    for (int y = 1; y < sizey - 1; y++) 
    {
      if (sigma[x][y] > 0) 
      {
        int current_cell = sigma[x][y];
        for (int i = 1; i <= n_nb_adh; i++) 
        {
          int xp2, yp2;
          xp2 = x + nx[i];
          yp2 = y + ny[i];
          if (xp2 > sizex-1 || xp2 < 1 || yp2 > sizey-1 || yp2 < 1)
            continue;
            
          if (par.periodic_boundaries) 
          {
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
          if (sigma[xp2][yp2] != sigma[x][y]) 
          {
            bool oppCD19 = (*cell)[sigma[xp2][yp2]].getCD19();
            double oppEcad = (*cell)[sigma[xp2][yp2]].getE_cadherin();
            double oppGFP = (*cell)[sigma[xp2][yp2]].getGFP();
            double oppNcad = (*cell)[sigma[xp2][yp2]].getN_cadherin();
            double oppPcad = (*cell)[sigma[xp2][yp2]].getP_cadherin();
            // cout << oppCD19 << '\t' << oppEcad << endl;
            (*cell)[current_cell].AddtoSurfaces(oppCD19, oppEcad, oppGFP, oppPcad, oppNcad);
            if (sigma[xp2][yp2]==0)
            {
              (*cell)[current_cell].setTouchingMed(true);
            }
          }
        }
      }
    }
  }
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      c->AverageSurfaceBindings();
    }
  }
}




void CellularPotts::StartSyntheticNetwork(int start_point)
{
  for (auto c=std::next(cell->begin(), start_point); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      c->SetMotilityStrength(par.motility_strength);
      c->SetElasticMod(par.elastic_modulus);
      double init_synNotch_bound = 1.0;
      double init_synNotch_unbound = 0.0;
      double init_synNotch_intra = 0.0;
      double init_E_cadherin = 0.0;
      double init_P_cadherin = 0.0;
      double init_N_cadherin = 0.0;
      double init_GFP = 0.0;
      double init_mCherry = 0.0;
      double init_random_binding_proteins = par.init_random_binding;

      c->setsynNotch_bound(init_synNotch_bound);
      c->setsynNotch_unbound(init_synNotch_unbound);
      c->setsynNotch_intra(init_synNotch_intra);
      c->setE_cadherin(init_E_cadherin);
      c->setP_cadherin(init_P_cadherin);
      c->setN_cadherin(init_N_cadherin);
      c->setGFP(init_GFP);
      c->setmCherry(init_mCherry);
      c->setRandomBindingProteins(init_random_binding_proteins);

      if (c->isSpheroid())
      {
        double init_P_cadherin=1.0;
        double init_random_binding_proteins = par.init_random_binding;
        c->setP_cadherin(init_P_cadherin);
        c->setRandomBindingProteins(init_random_binding_proteins);

        vector<double> init_diffusers{1.0};
        // c->set_diffusers(init_diffusers);
        cout << c->Sigma() << endl;

        c->SetConstitutives(par.spheroid_const);
        c->SetGFP_induced(par.spheroid_GFP_induced);
        c->SetMcherry_induced(par.spheroid_mCherry_induced);
        c->SetCD19_induced(par.spheroid_CD19_induced);


      }
      else
      {
        vector<double> init_diffusers{0.0};
        // c->set_diffusers(init_diffusers);
        // randomly make cell CD19 or not (move to different method eventually)
        double rand = RANDOM(s_val);
        if (rand < par.proportion_celltype2)
        {
          // assuming c2
          c->SetConstitutives(par.c2_const);
          c->SetGFP_induced(par.c2_GFP_induced);
          c->SetMcherry_induced(par.c2_mCherry_induced);
          c->SetCD19_induced(par.c2_CD19_induced);
        }
        else 
        {
          // assuming c1
          c->SetConstitutives(par.c1_const);
          c->SetGFP_induced(par.c1_GFP_induced);
          c->SetMcherry_induced(par.c1_mCherry_induced);
          c->SetCD19_induced(par.c1_CD19_induced);
        }
      }

    } 
  }
}

void CellularPotts::StartSyntheticNetwork(Cell &newcell)
{

  if (newcell.AliveP())
  {
    if (newcell.isSpheroid())
    {
      double init_P_cadherin=1.0;
      double init_random_binding_proteins = par.init_random_binding;
      newcell.setP_cadherin(init_P_cadherin);
      newcell.setRandomBindingProteins(init_random_binding_proteins);

      vector<double> init_diffusers{1.0};
      // newcell.set_diffusers(init_diffusers);
    }
    else
    {
      double init_synNotch_bound = 1.0;
      double init_synNotch_unbound = 0.0;
      double init_synNotch_intra = 0.0;
      double init_E_cadherin = 0.0;
      double init_P_cadherin = 0.0;
      double init_N_cadherin = 0.0;
      double init_GFP = 0.0;
      double init_mCherry = 0.0;
      double init_random_binding_proteins = par.init_random_binding;

      newcell.setsynNotch_bound(init_synNotch_bound);
      newcell.setsynNotch_unbound(init_synNotch_unbound);
      newcell.setsynNotch_intra(init_synNotch_intra);
      newcell.setE_cadherin(init_E_cadherin);
      newcell.setP_cadherin(init_P_cadherin);
      newcell.setN_cadherin(init_N_cadherin);
      newcell.setGFP(init_GFP);
      newcell.setmCherry(init_mCherry);
      newcell.setRandomBindingProteins(init_random_binding_proteins);

      vector<double> init_diffusers{0.0};
      // newcell.set_diffusers(init_diffusers);
    }

  } 
}

void CellularPotts::OutputSyntheticNetwork(int thetime)
{
  vector<Cell>::iterator c;

  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {  
      ofstream outfile;
      string out = data_file + "/cell-" + to_string(c->Sigma()) + ".dat";
      outfile.open(out, ios::app);
      outfile << c->getCD19() << '\t' << c->getsynNotch_bound() << '\t' << c->getsynNotch_unbound() << '\t' << c->getsynNotch_intra() << '\t' << c->getE_cadherin() << '\t' << '\t' << '\t' << c->getGFP() << '\t' << c->getmCherry() << endl;
      outfile.close();
    }
  }


}

void CellularPotts::UpdateSyntheticCellConstraints()
{
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      double area_constraint = par.bulk_modulus / double(c->TargetArea());
      c->setAreaConstraint(area_constraint);
      int target_perim = round(double(par.ptarget_perimeter) * sqrt(double(c->TargetArea())/double(par.cell_target_area)));
      c->SetTargetPerimeter(target_perim);
      
      double perim_constraint = (c->GetElasticMod() / double(target_perim));
      c->setPerimConstraint(perim_constraint);
      // cout << target_perim << '\t' << area_constraint << '\t' << perim_constraint << endl;

    }
  }
}



// NOTE: not called from embryo.cpp. This synNotch/GFP/mCherry/cadherin
// differentiation network (and the cellcolour it computes below) is legacy/
// unused in the current sox-driven sorting model, where cell colour comes
// from InitialiseRandomSoxValues() and adhesion from Cell::EmbryoEnergy().
// Kept here for reference in case this network is revisited later.
void CellularPotts::SyntheticNetwork()
{

  // WE MUST DO SURFACE BINDINGS FIRST!
  SurfaceBindings();

  vector<Cell>::iterator c;

  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      double dt = par.synthetic_dt;
      // averaging after surface bindings

      // cells either do or do not have the synethic network. So, we have to get
      // its network type and then decide how to update. For now, we say CD19.
      double& CD19 = c->getCD19();
      bool& spheroid_cell = c->isSpheroid();
      double& synNotch_bound = c->getsynNotch_bound();
      double& synNotch_unbound = c->getsynNotch_unbound();
      double& synNotch_intra = c->getsynNotch_intra();
      double& E_cadherin = c->getE_cadherin();
      double& GFP = c->getGFP();
      double& mCherry = c->getmCherry();
      double& N_cadherin = c->getN_cadherin();
      double& P_cadherin = c->getP_cadherin();

      double c_area = c->TargetArea();

      double opposite_Ecad = c->getOpposing_E_cadherin();
      double opposite_Pcad = c->getOpposingP_cadherin();
      double opposite_Ncad = c->getOpposingN_cadherin();

    
      double opposite_CD19=c->getOpposingCD19();
      double opposite_GFP=c->getOppositeGFP();

      /* all cells have random binding proteins*/
      // double& random_binding_proteins = c->getRandomBindingProteins();
      // random_binding_proteins = random_binding_rk4(dt, random_binding_proteins, opposite_Ecad);
      vector<bool> constitutives = c->GetConstitutives();
      vector<bool> GFP_induced = c->GetGFP_induced();
      vector<bool> mCherry_induced = c->GetMcherry_induced();
      vector<bool> CD19_induced = c->GetCD19_induced();

      if (constitutives[0]==true)
      {
        E_cadherin=1.;
      }
      if (constitutives[1]==true)
      {
        E_cadherin=0.5;
      }
      if (constitutives[2]==true)
      {
        P_cadherin==1.;
      }
      if (constitutives[3]==true)
      {
        N_cadherin=1.;
      }
      if (constitutives[4]==true)
      {
        CD19=1.;
      }
      if (constitutives[5]==true)
      {
        GFP=1.;
      }
      if (constitutives[6]==true)
      {
        mCherry=1.;
      }

      bool receptsGFP = std::find(GFP_induced.begin(), GFP_induced.end(), true) != GFP_induced.end();

      bool receptsCD19 = std::find(CD19_induced.begin(), CD19_induced.end(), true) != CD19_induced.end();
      
      if (receptsGFP)
      {
        synNotch_bound = synNotch_bound_rk4(dt, synNotch_bound, opposite_GFP, c_area);
        synNotch_unbound = synNotch_unbound_rk4(dt, synNotch_unbound, synNotch_bound, opposite_GFP, c_area);
        synNotch_intra = synNotch_intra_rk4(dt, synNotch_intra, synNotch_bound, opposite_GFP, c_area);
      }
      else if (receptsCD19)
      {
        synNotch_bound = synNotch_bound_rk4(dt, synNotch_bound, opposite_CD19, c_area);
        synNotch_unbound = synNotch_unbound_rk4(dt, synNotch_unbound, synNotch_bound, opposite_CD19, c_area);
        synNotch_intra = synNotch_intra_rk4(dt, synNotch_intra, synNotch_bound, opposite_CD19, c_area);
      }

      if (GFP_induced[0]==true || CD19_induced[0]==true)
      {
        E_cadherin = E_cadherin_rk4(dt, E_cadherin, synNotch_intra, opposite_Ecad, par.E_cadherin_production_rate, c_area);
      }
      if (GFP_induced[1]==true || CD19_induced[1]==true)
      {
        E_cadherin = E_cadherin_rk4(dt, E_cadherin, synNotch_intra, opposite_Ecad, par.lo_cadherin_production_rate, c_area);
      }
      if (GFP_induced[2]==true || CD19_induced[2]==true)
      {
        P_cadherin = E_cadherin_rk4(dt, P_cadherin, synNotch_intra, opposite_Pcad, par.E_cadherin_production_rate, c_area);
      }
      if (GFP_induced[3]==true || CD19_induced[3]==true)
      {
        N_cadherin = E_cadherin_rk4(dt, N_cadherin, synNotch_intra, opposite_Ncad, par.E_cadherin_production_rate, c_area);
      }
      if (GFP_induced[5]==true || CD19_induced[5]==true)
      {
        GFP = GFP_rk4(dt, GFP, synNotch_intra, c_area);
      }
      if (GFP_induced[6]==true || CD19_induced[6]==true)
      {
        mCherry = GFP_rk4(dt, mCherry, synNotch_intra, c_area);
      }
 
      // for colour output
      int rounded_GFP = round(GFP * 200);
      if (rounded_GFP > 100)
        rounded_GFP = 100;

      int rounded_cherry = round(mCherry * 200);
      if (rounded_cherry > 100)
        rounded_cherry = 100;

      int cellcolour{};
      if (rounded_cherry < 1)
      {
        cellcolour = rounded_GFP + 2;
      }
      else
      {
        cellcolour = rounded_cherry + 102;
      }
      c->set_ctype(cellcolour);
      
      if (c->isSpheroid())
      {
        c->set_ctype(203);
      }
      c->ResetFinalSurfaceBindings();
    }
  }
  UpdateSyntheticCellConstraints();
}

void CellularPotts::SyntheticGrowth(int t)
{

  vector<bool> which_cells(cell->size());
  int cell_division=0;
  // Note this function MUST come after synthetic network,
  // which sets the touching med variable

  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      double min_growth_rate=0.1;
      double max_growth_rate=0.2;
      
      double rand = RANDOM(s_val);

      int area=c->Area();

      double growth_rate = max_growth_rate-min_growth_rate * rand + min_growth_rate;
      c->LeftoverTargetArea(growth_rate);

      if (area>par.div_threshold) 
      {
        ofstream outfile;
        string out = data_file + "/division-times.dat";
        outfile.open(out, ios::app);
        int tc = t - c->get_time_created();
        outfile << tc << endl;
        outfile.close();
        which_cells[c->Sigma()]=true;
        cell_division++;
      }
    }
  }
  if (cell_division) 
  {
    DivideCells(which_cells, t);
  }

  // adjust constraints!
  UpdateSyntheticCellConstraints();
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



void CellularPotts::MakeSpheroid(int centerx, int centery, int radius)
{

  int current_cells = CountCells();

  if (radius > centerx/2 || radius > centery/2)
    cerr << "ERROR: SPHEROID RADIUS TOO LARGE!";
  // double total = sizex*sizey;
  // int ncells = round(total / 75.);
  // cout << ncells << endl;
  double A = double(par.cell_target_area);
  double distance = sqrt((A)/(2*sqrt(3)));
  double leftover = fmod(sizex-2, distance);
  int dividor = int(floor(double(sizex-2)/distance));
  // cout << "LEFTOVERS: " << leftover << '\t' << dividor << endl;
  distance += leftover/dividor;
  // cout << "distance between cells is: " << distance << endl;


  int ncells = CircleHexaCounter(radius, distance);
  FractureSheet(ncells);


  vector<VPoint> centers = HexaCircleCenters(radius, distance, centerx, centery, current_cells);


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
        int xdist = abs(x - centerx);
        int ydist = abs(y - centery);
        // Assign the closest center id to the grid cell if its inside the circle
        if ((xdist * xdist + ydist*ydist) < radius*radius)
        {
          sigma[x][y] = closestCenter;
        }
        else if (minDistance < distance / 2)
        {
          sigma[x][y] = closestCenter;

        }
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
  
  // 10. Single unified cleanup and target initialization loop
  int deadcells = 0;
  for (auto c = std::next(cell->begin(), current_cells); c != cell->end(); ++c) {
    if (c == cell->begin()) 
      continue;
    if (c->area == 0) 
    {
      c->Apoptose();
      ++deadcells;
    } 
    else 
    {
      c->SetTargetArea(c->area);
      c->makeAlive();
      c->setSpheroid(true);
    }
  }
  MeasureCellSizes();

  std::cout << "Grid generated | Radius: " << radius 
            << " | Cells populated: " << (ncells - deadcells)
            << " | Cells killed (0 area): " << deadcells << std::endl;

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
  // tperim is unused: the target perimeter is derived from par.ptarget_perimeter,
  // each cell's lineage and its *targeted* area (see Cell::UpdatePerimeterConstraint),
  // which is also kept in sync automatically whenever target area or Sox
  // concentrations change later on.
  (void)tperim;
  vector<Cell>::iterator i;
  for ( (i=cell->begin(),i++); i!=cell->end(); i++)
    i->UpdatePerimeterConstraint();
}


// Multithreading method
void CellularPotts::set_num(int in)
{
  org_num = in;
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


// may be useful
void CellularPotts::Vectorfield()
{
  if (par.velocities)
  {
    int i = 0;
    int interval = 1;
    vector<vector<double>> xdata{};
    vector<vector<double>> ydata{};

    for (; i < par.mcs-par.end_program;i+=interval)
    {
      vector<double> xpoint{};
      vector<double> ypoint{};

      vector<Cell>::iterator c;
      for ( (c=cell->begin(), c++);c!=cell->end();c++) 
      {
        if (c->AliveP())
        {
          vector<double>& xm = c->get_xcens();
          vector<double>& ym = c->get_ycens();

        
          // we want displacement from a while ago to account for back and forth motion
          double x2 = xm[i];
          double y2 = ym[i];

          xpoint.push_back(x2);
          ypoint.push_back(y2);

          
        }
      }
      xdata.push_back(xpoint);
      ydata.push_back(ypoint);
    }

    string var_name = data_file + "/xvector-data.dat";
    ofstream outfile;
    outfile.open(var_name, ios::app);

    for (vector<double>& i : xdata)
    {
      for (double &j : i)
      {
        outfile << j << '\t';
      }
      outfile << endl;
    }
    
    outfile.close();

    var_name = data_file + "/yvector-data.dat";
    outfile.open(var_name, ios::app);

    for (vector<double>& i : ydata)
    {
      for (double &j : i)
      {
        outfile << j << '\t';
      }
      outfile << endl;
    }

    outfile.close();
  }
  else 
  {
    cout << "called centroids without setting velocities to true" << endl;
  }

}