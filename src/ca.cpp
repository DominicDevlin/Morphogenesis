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


const int CellularPotts::nx[25] = {0, 0, 1, 0,-1, 1, 1,-1,-1, 0, 2, 0, -2, 1, 2, 2, 1,-1,-2,-2,-1, 0, 2, 0,-2 };
const int CellularPotts::ny[25] = {0,-1, 0, 1, 0,-1, 1, 1,-1,-2, 0, 2,  0,-2,-1, 1, 2, 2, 1,-1,-2,-2, 0, 2, 0 };

const int CellularPotts::nbh_level[5] = { 0, 4, 8, 20, 24 };
int CellularPotts::shuffleindex[9]={0,1,2,3,4,5,6,7,8};

extern Parameter par;


/** PRIVATE **/

using namespace std;
void CellularPotts::BaseInitialisation(vector<Cell> *cells) {
  CopyProb(par.T);
  cell=cells;
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4]).";
  
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
  
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";
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
  for (int x=0;x<sizex;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
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
  
int CellularPotts::DeltaH(int x,int y, int xp, int yp, const int tsteps, PDE *PDEfield)       
{
  int DH = 0;
  int i, sxy, sxyp;
  int neighsite;

  /* Compute energydifference *IF* the copying were to occur */
  sxy = sigma[x][y];
  sxyp = sigma[xp][yp];
    
  /* DH due to cell adhesion */
  for (i=1;i<=n_nb;i++) {
    int xp2,yp2;
    xp2=x+nx[i]; yp2=y+ny[i];
    if (par.periodic_boundaries) 
    {
      
      // since we are asynchronic, we cannot just copy the borders once 
      // every MCS
      
      if (xp2<=0)
	      xp2=sizex-2+xp2;
      if (yp2<=0)
	      yp2=sizey-2+yp2;
      if (xp2>=sizex-1)
	      xp2=xp2-sizex+2;
      if (yp2>=sizey-1)
	      yp2=yp2-sizey+2;
    
      neighsite=sigma[xp2][yp2];
	
    } 
    else 
    {
      if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1)
	      neighsite=-1;
      else
	      neighsite=sigma[xp2][yp2];
    }
    
    if (neighsite==-1) 
    { // border 
      DH += (sxyp==0?0:par.border_energy)-(sxy==0?0:par.border_energy);
    } 
    else 
    {
      // UP TO HERE!
      // DH += (*cell)[sxyp].CalculateJfromKeyLock((*cell)[neighsite].get_locks_bool(), (*cell)[neighsite].get_keys_bool()) 
      // - 
      if (tsteps < par.end_program)
        DH += (*cell)[sxyp].EnDif((*cell)[neighsite]) - (*cell)[sxy].EnDif((*cell)[neighsite]);
      else
        DH += (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) - (*cell)[sxy].EnergyDifference((*cell)[neighsite]);
      // debugging. 
      // cout << "COPYING: " << (*cell)[sxyp].getTau() << (*cell)[sxy].getTau() << std::endl;
      // cout << "sxyp is type: " << (*cell)[neighsite].getTau() << " with val: " << (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) 
      // << ". sxy is type:" << (*cell)[neighsite].getTau() << " with val: " << (*cell)[sxy].EnergyDifference((*cell)[neighsite]) << std::endl;
      
    }
  }

  
  // lambda is determined by chemical 0
  double lambda = (*cell)[sxy].get_lambda();
    
  //cerr << "[" << lambda << "]";
  if ( sxyp == MEDIUM ) {
    DH += (int)(lambda *  (1. - 2. *   
			       (double) ( (*cell)[sxy].Area() - (*cell)[sxy].TargetArea()) ));
  }
  else if ( sxy == MEDIUM ) {
    DH += (int)((lambda * (1. + 2. *  
			       (double) ( (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()) )));
  }
  else
    DH += (int)((lambda * (2.+  2.  * (double) 
			       (  (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()
			       - (*cell)[sxy].Area() + (*cell)[sxy].TargetArea() )) ));


  /* Chemotaxis */
  // if (PDEfield && (par.vecadherinknockout || (sxyp==0 || sxy==0))) {
    
  //   // copying from (xp, yp) into (x,y)
  //   // If par.extensiononly == true, apply CompuCell's method, i.e.
  //   // only chemotactic extensions contribute to energy change
  //   if (!( par.extensiononly && sxyp==0)) {
  //     int DDH=(int)(par.chemotaxis*(sat(PDEfield->Sigma(0,x,y))-sat(PDEfield->Sigma(0,xp,yp))));
    
  //     DH-=DDH;
  //   }
  // }
  double lambda2=0; 
  if (tsteps > par.end_program)
    lambda2 = (*cell)[sxyp].get_lambda_2(); // par.lambda2;

  // const double lambda_r=(*cell)[sxy].get_lambda_2();
  
  /* Length constraint */
  // sp is expanding cell, s is retracting cell  
  if ( sxyp == MEDIUM ) {
    DH -= (int)(lambda2*( DSQR((*cell)[sxy].Length()-(*cell)[sxy].TargetLength())
		       - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) - 
			      (*cell)[sxy].TargetLength()) ));
    
  }
  else if ( sxy == MEDIUM ) {
    DH -= (int)(lambda2*(DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength())
			 -DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y)-(*cell)[sxyp].TargetLength())));
    
  }
  else {
    DH -= (int)(lambda2*((DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength())
		     -DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y)-(*cell)[sxyp].TargetLength())) +
		    ( DSQR((*cell)[sxy].Length()-(*cell)[sxy].TargetLength())
		      - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) - 
			     (*cell)[sxy].TargetLength()) )) );
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

void CellularPotts::ConvertSpin(int x,int y,int xp,int yp)
{
  int tmpcell;
  if ( (tmpcell=sigma[x][y]) ) { // if tmpcell is not MEDIUM
    (*cell)[tmpcell].DecrementArea();
    (*cell)[tmpcell].RemoveSiteFromMoments(x,y);
        
    if (!(*cell)[tmpcell].Area()) {
      (*cell)[tmpcell].Apoptose();
      // cerr << "Cell " << tmpcell << " apoptosed\n";
    }
  }
  
  if ( (tmpcell=sigma[xp][yp]) ) {// if tmpcell is not MEDIUM
    (*cell)[tmpcell].IncrementArea();
    (*cell)[tmpcell].AddSiteToMoments(x,y);
    
  }
  sigma[x][y] = sigma[xp][yp];


}


/** PUBLIC **/
int CellularPotts::CopyvProb(int DH,  double stiff) {

  double dd; 
  int s;
  s=(int)stiff;
  if (DH<=-s) 
    return 2;
  
  // if DH becomes extremely large, calculate probability on-the-fly
  if (DH+s > BOLTZMANN-1)
    dd=exp( -( (double)(DH+s)/par.T ));
  else
    dd=copyprob[DH+s]; 

  if (RANDOM(s_val)<dd) 
    return 1; 
  else 
    return 0;
} 

void CellularPotts::CopyProb(double T) {
  int i;
  for ( i = 0; i < BOLTZMANN; i++ )
    copyprob[i] = exp( -( (double)(i)/T ) );
}

void CellularPotts::FreezeAmoebae(void) 
{
  if (frozen) 
    frozen=FALSE;
  else
    frozen=TRUE;
}

//! Monte Carlo Step. Returns summed energy change
int CellularPotts::AmoebaeMove(long tsteps, PDE *PDEfield)
{
  int loop,p;
  //int updated=0;
  thetime++;
  int SumDH=0;
  
  if (frozen) 
    return 0;

  loop=(sizex-2)*(sizey-2);
 
  for (int i=0;i<loop;i++) {  
    // take a random site
    int xy = (int)(RANDOM(s_val)*(sizex-2)*(sizey-2));
    int x = xy%(sizex-2)+1;
    int y = xy/(sizex-2)+1; 
    
    // take a random neighbour
    int xyp=(int)(n_nb*RANDOM(s_val)+1);
    int xp = nx[xyp]+x;
    int yp = ny[xyp]+y;
    
    int k=sigma[x][y];
    
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
        
        int D_H=DeltaH(x,y,xp,yp, tsteps, PDEfield);
        
        if ((p=CopyvProb(D_H,H_diss))>0) 
        {
          ConvertSpin( x,y,xp,yp );
          SumDH+=D_H;
        }
        // if (Probability(D_H)) 
        // {
        //   ConvertSpin( x,y,xp,yp );
        //   SumDH+=D_H;
        // }
      }
    } 
  }
  return SumDH;
  
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
      if (par.target_area) {
	c->SetTargetArea(par.target_area);
      } else	 {
	c->SetTargetArea(mean_area);
      }
    }
  }
}

void CellularPotts::MeasureCellSizes(void) {
  
  // Clean areas of all cells, including medium

  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->SetTargetArea(0);
    c->area = 0;
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
  
  // set the actual area to the target area
  {
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->SetAreaToTarget();

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


            if (par.polarity_on)
            {
              vector<double> &TF = daughterp->get_genes();
              for (int k=0;k<par.n_TF;++k)
                if (polarity[k])
                { 
                  TF[k+par.n_TF] = 0;
                }
            }

            if (par.gene_record)
            {
              daughterp->RecordDivision(t); // record division only in daughter cell.
              daughterp->reset_recordings(); 
            }

            
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

  // if (par.velocities)
  // {
  //   vector<Cell>::iterator c;
  //   for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  //   {
  //     if (which_cells[c->Sigma()])
  //     {
  //       vector<double>& xcens = c->get_xcens();
  //       vector<double>& ycens = c->get_ycens();

  //       int len = xcens.size();

        



  //     }
  //   }
  // }



  if (celldir) 
    delete[] (celldir);
  
  if (divflags)
    free(divflags);
}   


int CellularPotts::VerticalLine(int id)
{
  int llength{};
  int yval{};

  for (int y=0;y<sizey;y++)
  {
    int temp{};
    for (int x=0;x<sizex;x++) 
      if (sigma[x][y]==id) 
      {
         ++temp;
      }
    if (temp > llength)
    {
      llength = temp;
      yval = y;
    } 
  }

  int minl{sizex};
  int maxl{};

  for (int x=0;x<sizex;++x)
  {
    if (sigma[x][yval] == id)
    {
      if (x < minl)
        minl = x;
      if (x > maxl)
        maxl = x;
    }
  }
  int line = ceil((maxl + minl) / 2);

  return line;
}

int CellularPotts::HorizontalLine(int id)
{

  int llength{};
  int xval{};

  for (int x=0;x<sizex;x++)
  {
    int temp{};
    for (int y=0;y<sizey;y++) 
      if (sigma[x][y]==id) 
      {
         ++temp;
      }
    if (temp > llength)
    {
      llength = temp;
      xval=x;
    } 
  }

  int minl{sizey};
  int maxl{};

  for (int y=0;y<sizey;++y)
  {
    if (sigma[xval][y] == id)
    {
      if (y < minl)
        minl = y;
      if (y > maxl)
        maxl = y;
    }
  }
  int line = ceil((maxl + minl) / 2);

  return line;
}



void CellularPotts::xyCellDivision(int id, bool direction, int t)
{
  // if direction is true, draw vertical line. If direction is false, draw horizontal line. 
  int divisionline{};
  if (direction)
  {
    divisionline = VerticalLine(id);
  }
  else
  {
    divisionline = HorizontalLine(id);
  }

  bool found_division=false;

  
  /* division */
  {
    for (int i=0;i<sizex;i++)
      for (int j=0;j<sizey;j++) 
        if (sigma[i][j]==id) // i.e. not medium and not border state (-1)
        { 

          // Pointer to mother. Warning: Renew pointer after a new
          // cell is added (push_back). Then, the array *cell is relocated and
          // the pointer will be lost...
          Cell *motherp=&((*cell)[sigma[i][j]]);
          Cell *daughterp;
        
          /* Divide if NOT medium and if DIV bit set or divide_always is set */
          // if which_cells is given, divide only if the cell
          // is marked in which_cells.

          if (!found_division) 
          {
            found_division = true;

            // add daughter cell, copying states of mother
            daughterp=new Cell(*(motherp->owner));
            daughterp->CellBirth(*motherp);
            cell->push_back(*daughterp);
            
            // renew pointer to mother
            motherp=&((*cell)[sigma[i][j]]);

            delete daughterp;
            
            // array may be relocated after "push_back"
            
            // renew daughter pointers
            daughterp=&(cell->back());


            if (par.gene_record)
            {
              daughterp->RecordDivision(t); // record division only in daughter cell. 
            }

    
          } 
            
          /* Now the actual division takes place */
        
          if (direction && i>divisionline) // vertical division
          { 
            // cout << "position: " << i << "  " << j << " is changing" << endl;
            motherp->DecrementArea();
            motherp->DecrementTargetArea();
            motherp->RemoveSiteFromMoments(i,j);
            sigma[i][j]=daughterp->Sigma();
            daughterp->AddSiteToMoments(i,j);
            daughterp->IncrementArea();
            daughterp->IncrementTargetArea();

          }
          else if (!direction && j>divisionline) // horizontal
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

  
int CellularPotts::GrowInCells(int n_cells, int cell_size, double subfield) {

  
  int sx = (int)((sizex-2)/subfield);
  int sy = (int)((sizey-2)/subfield);
  
  int offset_x = (sizex-2-sx)/2;
  int offset_y = (sizey-2-sy)/2;
  
  if (n_cells==1) {
    return GrowInCells(1, cell_size, sizex/2, sizey/2, 0, 0);
  } else {
    return GrowInCells(n_cells, cell_size, sx, sy, offset_x, offset_y);
  }
}

int CellularPotts::GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y) {
  
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

  if (n_cells>1) {
    
    
    
    { for (int i=0;i<n_cells;i++) {
      
      sigma[RandomNumber(sx, s_val)+offset_x][RandomNumber(sy, s_val)+offset_y]=++cellnum;
      
    }}
  } else {
    sigma[sx][sy]=++cellnum;

  }

  // Do Eden growth for a number of time steps
  {for (int i=0;i<cell_size;i++) {
    for (int x=1;x<sizex-1;x++)
      for (int y=1;y<sizey-1;y++) {
	
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
	  {  for (int x=1;x<sizex-1;x++) {
      for (int y=1;y<sizey-1;y++) {
	sigma[x][y]=new_sigma[x][y];
      }
    }
  }}}
  free(new_sigma[0]);
  free(new_sigma);
  
  return cellnum;
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
  if (n_borders>2 && ( (stackp+1)>2 || one_of_neighbors_medium) ) {
    return false;
  }
  else 
    return true;

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

double CellularPotts::MeanCellArea(void) const {
  
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
  
  cerr << "Mean cell length is " << sum_length/((double)n) << endl;
  return (double)sum_area/(double)n;
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
    
      if (c->Area()>par.target_area) {
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
  if (CountCells() < par.max_cells)
  {
    vector<bool> which_cells(cell->size());
  
    int cell_division=0;
    
    vector<Cell>::iterator c;
    for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    {
      if (c->AliveP())
      {
        // if (c->Area()==1)
        // {
        //   int sig = c->Sigma();
        //   for (int i=0;i<sizex;i++)
        //     for (int j=0;j<sizey;j++)
        //     {
        //       if (sigma[i][j]==sig)
        //       {
        //         sigma[i][j]==0;
        //       }
        //     }
        //     c->Apoptose();
        //     cout << "apoptosing new way" << endl;
        // }

        int TA = c->TargetArea();
        int area = c->Area();
        int gthresh = par.gthresh;
        int sthresh;

        // if (c->get_stem())
        //   gthresh = par.stem_gthresh;
        // else
        //   gthresh = par.gthresh;

        if (c->get_shrink())
          sthresh = par.s_shrink;
        else
          sthresh = par.shrink;



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
          
          // int count=0;
          // if (!par.all_divide)
          // {
          //   vector<double> &genes = c->get_genes();

          //   int i = par.n_diffusers + par.n_MF;
          //   int k = i + par.n_stem;
          //   while (i < k)
          //   {
          //     if (genes[i] > 0.5)
          //       ++count;
          //     ++i;
          //   }
          // }
          // // only divide cells if they have stem cell genes on
          // // divide all cells on the rare case that there area is too large. 
          // if (count == par.n_stem || area > (double)(par.div_threshold) * 1.5) 
          // {
            which_cells[c->Sigma()]=true;
            cell_division++;

          // }
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
}



vector<bool> CellularPotts::divide_vector(void) 
{
  vector<bool> to_divide(cell->size());
  vector<Cell>::const_iterator i;
  for ( (i=cell->begin(),i++); i!=cell->end(); i++) 
  {
    if (i->AliveP()) 
    {
      to_divide[i->Sigma()]=true;
    } 
  }
  return to_divide;
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


// return the approxiamte location of the cell centre
vector<int> CellularPotts::MiddleOfCell(int sig)
{
  int lowx=sizex;
  int lowy=sizey;
  int highx=0;
  int highy=0;

  vector<int> pos;
  for (int i=0;i<sizex-1;++i)
    for (int j=0;j<sizey-1;++j)
    {
      if (sigma[i][j] == sig)
      {
        if (i<lowx)
          lowx=i;

        if (i>highx)
          highx=i;

        if (j<lowy)
          lowy=j;
        
        if (j>highy)
          highy=j;
      }
    }
  int midx = (lowx+highx)/2;
  int midy = (lowy+highy)/2;

  pos.push_back(sig);
  pos.push_back(midx);
  pos.push_back(midy);

  return pos;
}



void CellularPotts::set_MF(vector<vector<int>> middles, int gene)
{
  int xdif = middles.at(0).at(1) - middles.at(1).at(1);
  int ydif = middles.at(0).at(2) - middles.at(1).at(2);
  if (abs(xdif) > abs(ydif))
  {
    if (middles.at(0).at(1) > middles.at(1).at(1))
    {
      // furthest on right gets MF at 4 --> 1  
      // genes expression defualt is at 0 so only have to modify one
      std::vector<double>& g_list = cell->at(middles.at(0).at(0)).get_genes();
      g_list.at(gene) = 1;
      int val = g_list.at(2) * 4 + g_list.at(3)*3;
      cell->at(middles.at(0).at(0)).set_ctype(val);
    }
    else
    {
      std::vector<double>& g_list = cell->at(middles.at(1).at(0)).get_genes();
      g_list.at(gene) = 1;
      int val = g_list.at(2) * 4 + g_list.at(3)*3;
      cell->at(middles.at(1).at(0)).set_ctype(val);
    }
  }
  else 
  {
    if (middles.at(0).at(2) > middles.at(1).at(2))
    {
      // furthest on bottom(maybe decreasing y from top??) gets MF at 4 --> 1  
      std::vector<double>& g_list = cell->at(middles.at(0).at(0)).get_genes();
      g_list.at(gene) = 1;
      int val = g_list.at(2) * 4 + g_list.at(3)*3;
      cell->at(middles.at(0).at(0)).set_ctype(val);
    }
    else
    {
      std::vector<double>& g_list = cell->at(middles.at(1).at(0)).get_genes();
      g_list.at(gene) = 1;
      int val = g_list.at(2) * 4 + g_list.at(3)*3;
      cell->at(middles.at(1).at(0)).set_ctype(val);
    }
  }  
}




void CellularPotts::Programmed_Division(void)
{
  vector<bool> to_divide = divide_vector();
  int n_cells = CountCells();
  
  // set first maternal factor
  if (n_cells < 2)
  {
    // DivideCells(to_divide);

    int id{};
    vector<Cell>::const_iterator i;
    for ( (i=cell->begin(),i++); i!=cell->end(); i++)
    {
      if (i->AliveP()) 
      {
        id = i->Sigma();
      }
    }
    xyCellDivision(id, true);

    vector<vector<int>> middles;

    for ( (i=cell->begin(),i++); i!=cell->end(); i++)
    {
      if (i->AliveP()) 
      {
        middles.push_back(MiddleOfCell(i->Sigma()));
      }
    }
    set_MF(middles, 2);  
  }
  // set second maternal factor  
  else if (n_cells < 4)
  {
    // DivideCells(to_divide);  
    
    vector<int> id_list{};
    vector<Cell>::iterator i; 
    for ( (i=cell->begin(),i++); i!=cell->end(); i++)
    {
      if (i->AliveP()) 
      {
        id_list.push_back(i->Sigma());
      }
    }
    for (int id : id_list)
    {
      xyCellDivision(id, false);
    }    


    // determine which cells have first maternal factor on, and which have first maternal factor off
    vector<int> g4_on;
    vector<int> g4_off;
    
    for ( (i=cell->begin(),i++); i!=cell->end(); i++)
    {
      if (i->AliveP()) 
      {
        std::vector<double>& g = i->get_genes();
        if (g.at(2) > 0.5)
          g4_on.push_back(i->Sigma());
        else
          g4_off.push_back(i->Sigma());
      }
    }

    vector<vector<int>> middles1;
    for (int j : g4_on)
    {
      middles1.push_back(MiddleOfCell(cell->at(j).Sigma()));
    }
    set_MF(middles1, 3);
    
    vector<vector<int>> middles2;
    for (int j : g4_off)
    {
      middles2.push_back(MiddleOfCell(cell->at(j).Sigma()));
    }
    set_MF(middles2, 3);

    
  }
  else 
    DivideCells(to_divide);
} 



void CellularPotts::randomise_network()
{
  for (int i = 0; i < par.n_genes; ++i)
  {
    for (int j = 0; j < par.n_activators; ++j)
    {
      double val = RANDOM(s_val);
      // doing slight ON bias at the moment (about to remove)
      if (val < 0.03)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.16)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.78)
      {
        matrix[i][j] = 0;
      }
      else if (val < 0.94)
      {
        matrix[i][j] = 1;
      }
      else
      {
        matrix[i][j] = 2;
      }
    }
  }

  for (int i=0;i<par.n_TF;++i)
  {
    double val = RANDOM(s_val);
    if (val < 0.7)
      polarity[i]=0;
    else
      polarity[i]=1;
  }

  cout << "{ ";
  for (int i = 0; i < par.n_genes; ++i)
  {
    cout << " { ";
    for (int j = 0; j < par.n_activators; ++j)
    {
      cout << matrix[i][j] << ", "; 
    }
    cout << "}, ";
  }
  cout << " }" << endl;
  cout << "TF polarities { ";
  for (int k = 0; k < par.n_TF; ++k)
  {
    cout << polarity[k] << ", "; 
  }
  cout << "}" << endl;
}


void CellularPotts::start_network(vector<vector<int>> start_matrix, vector<bool> start_pol)
{
  matrix.clear();
  matrix.resize(par.n_genes);

  polarity.clear();
  polarity.resize(par.n_TF);

  for (vector<int> &i : matrix)
  {
      i.resize(par.n_activators);
  }
  

  if (par.randomise)
  {
    randomise_network();
  }
  else 
  {
    for (int i = 0; i < par.n_genes; ++i)
    {
      for (int j = 0; j < par.n_activators; ++j)
      { 
        matrix[i][j] = start_matrix[i][j]; 
      }
    }
    for (int i=0;i<par.n_TF;++i)
    {
      polarity[i] = start_pol[i];
    }      
  } //set matrix to input, either from simulation or par.set_matrix


  vector<double> new_g;
  vector<double> new_d;
  for (int i=0;i<par.gene_vector_size;++i)
  {
    if (i >= par.n_diffusers+par.n_MF && i < par.n_diffusers+par.n_MF+par.n_TF) 
      new_g.push_back(1.0);
    else
      new_g.push_back(0.0);

    if (i<par.n_diffusers)
      new_d.push_back(0.0);
  }

  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      c->set_genes(new_g);
      c->set_diffusers(new_d);
      c->set_lists();
    }
  }
}


double CellularPotts::discrete_decay(vector<double>& gene_list, double conc, int gene_n)
{

  double x_1 = 0;

  // if (12 > gene_n > 5)
  // {
  //   for (int j=0; j < par.n_activators; ++j)
  //   {
  //     x_1 = x_1 + (matrix[gene_n][j] * gene_list.at(j));
  //   }
  //   x_1 += par.theta;
  //   if (x_1 > 2)
  //     return 1;
  //   else if (x_1 < -2)
  //     return 0;
  //   else
  //     return (1 / (1 + exp(-20 * x_1)));
  // }
  // else
  // {
  for (int j=0; j < par.n_activators; ++j)
  {
    x_1 += (matrix[gene_n][j] * gene_list.at(j));
  }
  x_1 += par.theta;

  x_1 = (1 / (1 + exp(-20 * x_1))) * 0.25 + conc * par.d_rate;
    
  return x_1;

  // }
}



void CellularPotts::update_network(int tsteps)
{
  vector<Cell>::iterator c;

  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      vector<double>& genes = c->get_genes();

      vector<double>& diffusers = c->get_diffusers(); 

      vector<double>& locks = c->get_locks();
      vector<double>& keys = c->get_keys();
      vector<double>& meds = c->get_medp();

      vector<double> gene_copy = c->get_genes();

      if (genes.empty() == true)
        cout << "EMPTY" << endl;

      // iterate through genes and update them according to GRN
      // iteration is a bit weird through different vectors so dont screw it up. 
      int j=0;
      int k=0;
      int m=0;
      for (int i = 0; i < par.n_genes; ++i)
      {
        if (i < par.n_diffusers)
          diffusers[i] = discrete_decay(gene_copy, diffusers[i], i);
        else if (i < par.n_genes - par.n_lockandkey - par.n_mediums)
          genes[i] = discrete_decay(gene_copy, genes[i], i);
        else if (i < par.n_genes - par.n_locks - par.n_mediums)
        {
          locks[j] = discrete_decay(gene_copy, locks[j], i);
          ++j;
        }
        else if (i < par.n_genes - par.n_mediums)
        {
          keys[k] = discrete_decay(gene_copy, keys[k], i);
          ++k;
        }
        else 
        {
          meds[m] = discrete_decay(gene_copy, meds[m], i);
          ++m;
        }
      }

      // set to shrinker
      if (par.shrink_on && !c->get_shrink() && genes[par.shrink_loc] > 0.5 )
      {
        c->set_shrink(true);
      }
      else if (c->get_shrink())
      {
        c->set_shrink(false);  
      }


      // make boolean set. 
      vector<bool>& full_set = c->get_set();

      // for recording differentiation events to output.
      array<bool, 17> cp;
      if (par.gene_record && tsteps > par.end_program)
      {
        for (int i = 0;i<(int)full_set.size();i++)
        {
          cp[i] = full_set[i];
        }
      }


      //create bool based on lock and keys.
      vector<bool>& l_bool = c->get_locks_bool();
      vector<bool>& k_bool = c->get_keys_bool();
      vector<bool>& m_bool = c->get_medp_bool();

      for (int i=0; i < par.n_locks; ++i)
      {

        l_bool[i] = (locks[i]>0.5) ? true : false;
        full_set[i] = l_bool[i];

        k_bool[i] = (keys[i]>0.5) ? true : false;
        full_set[i+par.n_locks] = k_bool[i];
      }

      for (int i=0; i < par.n_mediums; ++i)
      {
        m_bool[i] = (meds[i]>0.5) ? true : false;
        full_set[i+par.n_lockandkey] = m_bool[i];
      }



      // hard coding the two target genes in for now. 
      full_set[15] = ((genes.at(par.tloc1)>0.5) ? true : false);
      full_set[16] = ((genes.at(par.tloc2)>0.5) ? true : false);


      if (par.gene_record && tsteps > par.end_program)
      {
        
        c->AddPhenotype();
        c->RecordLongSwitch(cp, RandomNumber(INT_MAX, s_val));
        if (tsteps > par.adult_begins)
        {
          c->AddType();
        }
        if (!par.potency_edges)
          c->RecordSwitch(cp, RandomNumber(INT_MAX, s_val));
        else if (par.potency_edges && tsteps > par.adult_begins)
        {
          c->RecordSwitch(cp, RandomNumber(INT_MAX, s_val));
        }
        // unfortunately need to add old and new phenotype just in case a node is missed. 
        c->Phenotype();
        // Have to record types before and after switch otherwise a node could be skipped.
        c->AddPhenotype();
        if (tsteps > par.adult_begins)
        {
          c->AddType();
        }
      }



      // set the type of the cell based on network arrangement.
      if (tsteps > par.end_program)
      {
        c->set_ctype(set_type(full_set) * c->getTau());
        if (par.set_colours)
        {
          SetColours();
        }
      }

      // c->add_to_cycle();

      
      /// change target length based on boolified values of two target genes at location 12 and 13
      if (genes.at(par.tloc1) > 0.5 && genes.at(par.tloc2) > 0.5)
      {
        c->SetTargetLength(round(c->Area() / par.tlength2)); 
        c->set_lambda_2(par.lambda2);
      }
      else if (genes.at(par.tloc1) > 0.5 || genes.at(par.tloc2) > 0.5)
      {
        c->SetTargetLength(round(c->Area() / par.tlength1));
        c->set_lambda_2(par.lambda2);
      }
      else
      {
        c->SetTargetLength(0.0);
        c->set_lambda_2(0);    
      }

    }

  }
}


bool CellularPotts::CycleCheck()
{
  int amount=0;
  int cycling=0;
  vector<Cell>::iterator i;
  for ( (i=cell->begin(),i++); i!=cell->end(); i++) 
  {
    if (i->AliveP()) 
    {
      amount++;
      cycling += i->limit_cycle();
    } 
  }
  if (cycling*4 > amount)
  {
    return true;
  }
  else
  {
    return false;
  } 

}




void CellularPotts::print_random_cell(void)
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
  int val = RandomNumber(amount, s_val);
  cell->at(val).print_genes();
  cout << "The random cell is: " << val << endl;
  cout << "There are: " << amount << " cells in the organism." << endl;
  // get_neighbours(val);
}



int CellularPotts::random_cell(void)
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
  int val = RandomNumber(amount, s_val);
  return val; 
}


void CellularPotts::swap_cells(void)
{
  bool acell1 = false;
  bool acell2 = false;
  int choose1 = 1;
  int choose2 = 1;
  while (!acell1)
  {
    choose1 = random_cell();
    if (cell->at(choose1).AliveP() && cell->at(choose1).Area() > 50)
      acell1 = true;
  }
  while (!acell2)
  {
    choose2 = random_cell();
    if (cell->at(choose2).AliveP() && cell->at(choose2).Area() > 50)
      acell2 = true;
  }
  

  int sig1 = cell->at(choose1).Sigma();
  int sig2 = cell->at(choose2).Sigma();

  // now swap the location of the two cells:
  for (int x=1; x<sizex; ++x)
    for (int y=1; y<sizey; ++y)
    {
      if (sigma[x][y] == sig1)
      {
        sigma[x][y] = sig2;

      }
      else if (sigma[x][y] == sig2)
      {
        sigma[x][y] = sig1;
      }
    }

  int T1 = (*cell)[sig1].TargetArea();
  int T2 = (*cell)[sig2].TargetArea();

  MeasureCellSize((*cell)[sig1]);
  MeasureCellSize((*cell)[sig2]);

  (*cell)[sig1].SetTargetArea(T2);
  (*cell)[sig1].SetTargetArea(T1);

  // now have to recount cell area etc. 

}




double CellularPotts::diffuser_check(int n, int x, int y)
{
  vector<double> &diff = (*cell)[sigma[x][y]].get_diffusers();
  if (diff[n] > 0.01)
  {
    // cout << "pos: " << x << "  " << y << ". Concentration: " << diff[n] << endl;
    return diff[n];
  }
  else
    return 0;
}

double CellularPotts::get_enzyme_conc(int n, int x, int y)
{
  double conc;
  if (n == 0) 
    conc = (*cell)[sigma[x][y]].return_enzyme(par.e1_loc);
  else
    conc = (*cell)[sigma[x][y]].return_enzyme(par.e2_loc);
  return conc;

}




int CellularPotts::set_type(vector<bool>& set)
{

  // iterate through type list to see if type is already there.
  auto it = find(type_list.begin(), type_list.end(), set);

  if (it != type_list.end())
  {
    return (it - type_list.begin() + 4);
  }
  else 
  {
    type_list.push_back(set);
    
    // print new cell type network
    // for (int i=0; i < new_list.size();++i)
    // {
    //   cout << new_list.at(i) << " ";
    // }
    // cout << endl;
    return static_cast<int>(type_list.size()) + 4;
  }
}


int CellularPotts::get_ntypes()
{
  std::vector<vector<bool>> cell_types{};
  vector<Cell>::iterator c;
  for ( (c=cell->begin(),c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      std::vector<bool>& set = c->get_set();

      auto it = find(cell_types.begin(), cell_types.end(), set);

      if (it == cell_types.end())
      {
        cell_types.push_back(set);
      }
    }
  }
  int n_types = static_cast<int>(cell_types.size());

  return static_cast<int>(cell_types.size());

}



int CellularPotts::hamming_distance(vector<bool> &str1, vector<bool> &str2)
{
  int dist=0;
  
  if (str1.size() != str2.size())
  {
    return -1;
  }

  int length = str1.size();

  for (int i=0; i < length; ++i)
  {
    dist += (str1[i] != str2[i]);
  }
  return dist;
}




int hamming(int t1, int t2)
{
  string bin1 = "";
  string bin2 = ""; 

  int count = 0;
  while (count < par.n_lockandkey + par.n_mediums + par.n_length_genes)
  {
    bin1 = to_string(t1 % 2) + bin1;
    t1 = t1 / 2;

    bin2 = to_string(t2 % 2) + bin2;
    t2 = t2 / 2;


    ++count;
  }

  int distance = 0;
  for (int i = 0; i < bin1.length(); i++) {
      if (bin1[i] != bin2[i]) {
          distance++;
      }
  }
  return distance;



}



void CellularPotts::CellHammingDifferences()
{

  // Need to get the cell phenotype.
  unordered_map<int, int> phenotypes{};

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
      int p = c->GetPhenotype();
      phenotypes[p] += 1;

    }
  }
  
  ofstream outfile;
  string out = data_file + "/hammings.dat";
  outfile.open(out, ios::app);


  vector<int> types;

  for (auto type : phenotypes)
  {
    if (type.second > 4)
      types.push_back(type.first);

  }

  int max_hamming_phenotype = 0;
  // compare hamming distance for all types
  for (int i=0;i<types.size();++i)
  {
    for (int j = i + 1; j < types.size(); ++j)
    {
      int t1 = types[i];
      int t2 = types[j];
      int distance = hamming(t1, t2);
      outfile << t1 << " " << t2 << " " << distance << endl;
      if (distance > max_hamming_phenotype)
        max_hamming_phenotype = distance;

    }
  }
  outfile.close();

  // now do it for regulation:

  unordered_map<int, int> reg_types{};


  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      int p = c->RegPhenotype();
      reg_types[p] += 1;

    }
  }
  
  out = data_file + "/regulatory-hammings.dat";
  outfile.open(out, ios::app);


  vector<int> regs{};

  for (auto type : reg_types)
  {
    if (type.second > 4)
      regs.push_back(type.first);

  }  

  int max_hamming_reg=0;
  // compare hamming distance for all types
  for (int i=0;i<regs.size();++i)
  {
    for (int j = i + 1; j < regs.size(); ++j)
    {
      int t1 = regs[i];
      int t2 = regs[j];
      int distance = hamming(t1, t2);
      outfile << t1 << " " << t2 << " " << distance << endl;
      if (distance > max_hamming_reg)
        max_hamming_reg = distance;
    }
  }
  outfile.close();


  // data to central file. I want: number of regulatory states, number of phenotype states, max hamming distance for both regulatory and phenotype
  // only need to do it at final time point. 



  out = "hamming_data.dat";
  outfile.open(out, ios::app);


  outfile << phenotypes.size() << "\t" << max_hamming_phenotype << "\t" << regs.size() << "\t" << max_hamming_reg << endl;

  outfile.close();









} 









void CellularPotts::TypeFitness()
{
  vector<vector<bool>> cell_types{};
  vector<Cell>::iterator c;
  for ( (c=cell->begin(),c++); c!=cell->end(); c++) 
  {
    if (c->AliveP() && c->checkforcycles(par.cycle_threshold) == false)
    {

      std::vector<bool>& set = c->get_set();
      
      auto it = find(cell_types.begin(), cell_types.end(), set);

      if (it == cell_types.end())
      {
        cell_types.push_back(set);
      }
    }
  }

  int n_types = static_cast<int>(cell_types.size());

  
  // static cast to double so can divide by 10.
  // 10 is arbirtrary number so fitness bonus from cells isn't too large. n_som > 10 ---> log(n_som/10) > 0 (otherwise negative !bad!). 
  // should i change from log to linear?
    
  // if (n_types > 1)
  // {
  //   // find average distance between all types. will use log this for now. 
  //   double average_hamming{};
  //   // find smallest hamming distance between types. Not using currently.   
  //   // int min_hamming{100};


  //   int i=0;
  //   int j=1;
  //   for (; i < n_types -1; ++i)
  //   {
  //     for (int k=j; k < n_types; ++k)
  //     {
  //       int h = hamming_distance(cell_types[i], cell_types[k]);
  //       average_hamming += h;
  //       // min_hamming = ((h < min_hamming) ? h : min_hamming);
  //     }
  //     ++j; 
  //   }
  //   int n_comparisons = (n_types * (n_types-1)) / 2;
    
  //   average_hamming = average_hamming / n_comparisons;

  //   // cout << "average hamming distance is: " << average_hamming << ". Number of types is: " << n_types << 
  //   // ". Number of somatic cells cells is: " << n_som << endl;
  //   type_fitness_list.push_back(static_cast<double>(n_types) + (log(average_hamming)));
  // }
  // else
  //   type_fitness_list.push_back(1);


  if (par.print_fitness)
    cout << "Latest type diversity is: " << type_fitness_list.back() << endl;
}





//new type fitness function. CURRENTLY DEPRACATED
int CellularPotts::TypeFitness2()
{


  // SECOND OPTION: NUMBER OF ATTRACTORS:
  // If nodes are MORE THAN A TWO CYCLE = ONLY 2 POINTS. NEED TO SOMEHOW COLLATE NODES TOGETHER. 
  // WHAT IM SAYING IS -> FOR EVERY STRONGLY CONNECTED COMPONENT = TWO POINTS (or one if its a fixed point)
  // ONLY IF THE NODE HAS A SIZE > 3 PERCENT (or the entire SCC???) !!!
  int n_types{};

  // create graph.
  // map<int, int> phens = get_phenotype_time();
  // map<int, int> diffs = get_AdultTypes();  

  // vector<int> edge_start{};
  // vector<int> edge_end{};
  // set_switches(edge_start, edge_end);    


  // Graph newgraph(phens.size());
  // newgraph.CreateDiGraph(phens, diffs, edge_start, edge_end);

  // vector<vector<int>> pruned_comps = newgraph.ReturnPrunedComponents();

  // for (vector<int> i : pruned_comps)
  // {
  //   if (i.size() > 1)
  //     n_types += 2;
  //   else
  //     ++n_types;
  //   for (int j : i)
  //     cout << j << " ";
  //   cout << endl;
  // }
  // cout << "n types is " << n_types << endl;


  return n_types;


 // We are going to do it this way: if a cell EVER switches back to itself, it is excluded. 
  // Only cells that go in one direction or stay in one spot are included. 

  // unordered_map<int, int> cell_types;
  // int total_time = 0;

  // vector<Cell>::iterator c;
  // for ((c=cell->begin(), c++); c!=cell->end(); c++)
  // {
  //   if (c->AliveP())
  //   {  
  //     vector<int> visited_types;
  //     unordered_map<int,int>& phent = c->PhenTime();
  //     bool cycles = false;
  //     for (auto m : phent)
  //     {
  //       //check if cell has visited this state before:
  //       if (find(visited_types.begin(), visited_types.end(), m.first) != visited_types.end())
  //       {
  //         // removed_cells[distance(cell->begin(), c)] = true;
  //         cycles = true;
  //         break;
  //       }
  //       else
  //       {
  //         visited_types.push_back(m.first);
  //       }
  //     }
  //     if (!cycles)
  //     {
  //       for (int i : visited_types)
  //       {
  //         cell_types[i] += 1;
  //         ++total_time;
  //       }
  //     }
  //   }

  // }
  // cout << "HERE" << endl;
  // double remove_threshold = 0.04 * total_time;
  
  // for (auto it = cell_types.begin() ;it!= cell_types.end();)
  // {
  //   if (it->second < remove_threshold)
  //     it = cell_types.erase(it);
  //   else
  //     ++it;

  // }

  // for (auto m : cell_types)
  // {
  //   cout << "starting..." << m.first << "  " << m.second << "  done."  << endl;
  // }

  // int n_types = 0;

  // // now look at current phenotypes:
  // for ((c=cell->begin(), c++); c!=cell->end(); c++)
  // {
  //   int cphen = c->GetPhenotype();
  //   if (cell_types.count(cphen) > 0)
  //   {
  //     ++n_types;
  //     cell_types.erase(cphen);
  //   }
  // }
  
  // return n_types;

}


// return the fitness: combination of n cell types + hamming distance between them + number of cells. 
void CellularPotts::update_fitness()
{
 
  // shape fitness
  if (ShapeMaintained)
  {
    // Fitness is a combination of Circle Deviation (deformation - force applied) and segments (coordinated cell movements) = complexity

    double dev = (DeviationFromCircle()) * 1.4; // The question is whether to add total mass into this equation. Sqrt cell mass maybe?
    double wspc = sqrt((WhiteSpace())) * 2;
    double asymmetry = 0;
    if (par.asymmetry_selection)
    {
      asymmetry = TraverseFitness();
    }

    if (par.print_fitness)
    {
      cout << "Circle Deviation fitness contribution: " << dev << endl;
      cout << "Segment contribution to fitness: " <<  wspc << endl;
    }     


    if (par.asymmetry_selection && par.asym_only)
    {
      shape_fitness_list.push_back(asymmetry);
    }
    else
    {
      asymmetry = asymmetry / 4.;
      shape_fitness_list.push_back(dev + wspc + asymmetry);
    }
  }
  else 
  {
    if (par.print_fitness)
      cout << "SHAPE NOT MAINTAINED IN ORG N " << org_num << endl;
  }

}




// called once at end of sim.
double CellularPotts::get_fitness()
{

  // Calculate the shape. All connected has been checked for every 2000 steps in "evolution.cpp".
  double shp{};
  // double circumference = 0;
  if (ShapeMaintained && shape_fitness_list.size() > 0)
  {
    for (double k : shape_fitness_list)
    {
      shp += k;
    }
    // take average to get rid of some variation
    shp = shp / static_cast<double>(shape_fitness_list.size());
  }
  else 
  {
    shp = 0;
  }
  
  // multiplier to ensure shape dominates fitness. 
  double fitness = shp; 

  if (par.print_fitness)
  {
    cout << "Shape: " << shp << "   Fitness: " << fitness << endl;
  }


  // function to select for growth
  if (par.growth_selection)
  {
    return CellDensity();
  }


     

  return fitness;

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



// void CellularPotts::RecordTypes()
// {
//   vector<Cell>::iterator c;
//   for ((c=cell->begin(), c++); c!=cell->end(); c++)
//   {
//     if (c->AliveP())
//     {
//       c->AddPhenotype();
//     }
//   }
// }

// void CellularPotts::RecordAdultTypes()
// {
//   vector<Cell>::iterator c;
//   for ((c=cell->begin(), c++); c!=cell->end(); c++)
//   {
//     if (c->AliveP())
//     {
//       c->AddType();
//     }
//   }
// }




void CellularPotts::record_GRN()
{
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
      c->add_to_vectors();
    }
  }
}




unordered_map<string, int> CellularPotts::transitions(bool cycling)
{
  vector<uint64_t> RNGers{};
  unordered_map<string, int> switch_tally{};
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      vector<tuple<int,int,uint64_t>> switches{};
      if (cycling)
      {
        switches = c->get_long_switches();
      }
      else
      {
        switches = c->get_switches();
      }

      for (auto i : switches)
      {
        // Ensure phylogenies are independent (aren't picking up the same differentiation event twice)
        uint64_t rng = get<2>(i);
        if (find(RNGers.begin(), RNGers.end(), rng) != RNGers.end())
        {
          continue;
        }
        else
        {
          string enums;

          enums = "('" + to_string(get<0>(i)) + "', '" + to_string(get<1>(i)) + "')";
          RNGers.push_back(get<2>(i));
          switch_tally[enums] += 1;
        }
      }
    }
  }

  if (par.gene_output)
  {
    if (mkdir(data_file.c_str(), 0777) == -1)
      cerr << "Error : " << strerror(errno) << endl;
    else
      cout << "Directory created." << endl;

    ofstream outfile;
    string switch_out = data_file + "/differentiation_tally.dat";
    outfile.open(switch_out, ios::app);
    for (auto i : switch_tally)
    {
      outfile << i.first << '\t' << i.second << endl;
      // outfile << i.first.substr(0,4) << '\t' << i.first.substr(4,7) << '\t' << i.second << endl;
    }
    outfile.close();

  }

  return switch_tally;
}

// Used in "potency tester." Add all the edges from different simulations of the same network together. 
void CellularPotts::set_switches(map<pair<int,int>,int>& tally)
{
  vector<uint64_t> RNGers{};

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      vector<tuple<int,int,uint64_t>> &switches = c->get_switches();
      for (auto i : switches)
      {
        
        // Ensure phylogenies are independent (aren't picking up the same differentiation event twice)
        uint64_t rng = get<2>(i);
        if (find(RNGers.begin(), RNGers.end(), rng) != RNGers.end())
        {
          continue;
        }
        else
        {
          int f1 = get<0>(i);
          int f2 = get<1>(i);
          // cout << f1 << "  " << f2 << endl;
          pair<int,int> newpair = {f1, f2};

          tally[newpair] += 1;

          // bool add = true;
          // int maxcon = edge_start.size();
          // for (int j = 0; j < maxcon; ++j)
          // {
          //   if (f1 == edge_start[j] && f2 == edge_end[j])
          //   {
          //     add = false;
          //     j = maxcon;
          //   }
          // }
          // if (add)
          // {
          //   edge_start.push_back(f1);
          //   edge_end.push_back(f2);
          //   // cout << "adding... " << f1 << "  " << f2 << endl;
          // }
          RNGers.push_back(get<2>(i));
        }
      }
    }
  }
}

// used when there is long limit cycles
void CellularPotts::set_long_switches(map<pair<int,int>,int>& tally)
{
  vector<uint64_t> RNGers{};

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      vector<tuple<int,int,uint64_t>> &switches = c->get_long_switches();
      for (auto i : switches)
      {
        
        // Ensure phylogenies are independent (aren't picking up the same differentiation event twice)
        uint64_t rng = get<2>(i);
        if (find(RNGers.begin(), RNGers.end(), rng) != RNGers.end())
        {
          continue;
        }
        else
        {
          int f1 = get<0>(i);
          int f2 = get<1>(i);
          // cout << f1 << "  " << f2 << endl;
          pair<int,int> newpair = {f1, f2};

          tally[newpair] += 1;
          RNGers.push_back(get<2>(i));
        }
      }
    }
  }
}




void CellularPotts::cell_concentrations()
{
  string fnamen = data_file + "/conc";

  if (mkdir(fnamen.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;  

  if (mkdir(data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      string var_name = data_file + "/conc/cell_conc_" + to_string(c->Sigma()) + ".dat";
      ofstream outfile;
      outfile.open(var_name, ios::app);
      
      vector<vector<double>>& gene_history = c->get_history();

      for (unsigned int i=0;i<gene_history.at(0).size();++i)
      {
        for (int j=0;j<(int)(gene_history.size()); ++j)
        {
          outfile << gene_history.at(j).at(i) << '\t';
        }
        outfile << endl;
      }
      outfile.close();
    }
  }
}




void CellularPotts::phenotype_time()
{

  unordered_map<int, int> phenotype_time{};
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {  
      unordered_map<int,int>& phent = c->PhenTime();
      for (auto m : phent)
      {
        phenotype_time[m.first] += m.second;
      } 
    }
  }

  if (par.gene_output)
  {
    if (mkdir(data_file.c_str(), 0777) == -1)
      cerr << "Error : " << strerror(errno) << endl;
    else
      cout << "Directory created." << endl;


    ofstream outfile;
    string fname = data_file + "/phenotype_time.dat";
    outfile.open(fname, ios::app);
    for (auto m : phenotype_time)
    {
      // how long cells have spent in each state
      outfile << m.first << '\t' << m.second << endl;
    }
    outfile.close();
  }
}


map<int, int> CellularPotts::get_phenotype_time()
{
  map<int, int> phenotype_time{};
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {  
      unordered_map<int,int>& phent = c->PhenTime();
      for (auto m : phent)
      {
        phenotype_time[m.first] += m.second;
      } 
    }
  }
  return phenotype_time;
}

map<int, int> CellularPotts::get_AdultTypes()
{
  map<int, int> diff_time{};
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {  
      unordered_map<int,int>& phent = c->AdultTime();
      for (auto m : phent)
      {
        diff_time[m.first] += m.second;
      } 
    }
  }
  return diff_time;

}










void CellularPotts::cell_divisions()
{
  string fnamen = data_file + "/types";
  if (mkdir(fnamen.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;



  vector<pair<int, int>> divisions{};

  unordered_map<int, int> divphentally{};


  vector<pair<int,int>> buffer{};


  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      vector<pair<int,int>>& cell_d = c->get_divisions();
      vector<pair<int,int>>& cell_p = c->DivisionPhenotype();
      for (pair<int,int> &t : cell_d)
      {
        if (find(divisions.begin(), divisions.end(), t) != divisions.end())
          continue;
        else
          divisions.push_back(t);
      }
      for (auto &t : cell_p)
      {
        if (find(buffer.begin(), buffer.end(), t) != buffer.end())
        {
          continue;
        }
        else
        {
          buffer.push_back(t);
          divphentally[t.first] += 1;
        }
      }


    }
  }
  string fname = data_file + "/divisions.dat";
  ofstream outfile;
  outfile.open(fname, ios::app);
  for (pair<int,int> &i : divisions)
  {
    // first is cell identity, second is time step. 
    outfile << i.first << '\t' << i.second << endl;
  }
  outfile.close();

  string dpt = data_file + "/div_phen_tally.dat";
  outfile.open(dpt, ios::app);
  auto it = divphentally.begin();
  for (; it != divphentally.end();++it)
  {
    // divisions by cell number
    outfile << it->first << '\t' << it->second << endl;
  }

  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      string var_name = data_file + "/types/cell_types_" + to_string(c->Sigma()) + ".dat";
      outfile.open(var_name, ios::app);
      
      vector<int>& hist = c->TypeHistory();

      for (int p : hist)
      {
        outfile << p <<  '\t';

        int q = p;
        int m = par.n_functional-1;
        while (m>=0)
        {
          int mltp =  pow(2,m);
          if (q >= mltp)
          {
            outfile << 1 << '\t';
            q = q - mltp;
          }
          else 
          {
            outfile << 0 << '\t';
          }
          --m;
        }
        outfile << endl;
      }

      outfile.close();
    }
  }      
}



void CellularPotts::print_cell_GRN()
{
  #ifdef FILESYSTEM
  // need to add  -lstdc++fs -std=c++17 to CXXFLAGS in makefile in order for filesystem to compile.
  std::uintmax_t n = std::filesystem::remove_all(data_file);
  if (n)
    cout << "Removed " << n << " folders!!" << endl;
  #endif

  cell_concentrations();


  cell_divisions();

  phenotype_time();

  transitions(false);

  PrintTypesTime(true);

  CellHammingDifferences();

  BindingBetweenCells();

}




void CellularPotts::PrintPhenotypes()
{
  unordered_map<int, int> phenotypes{};

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
      int p = c->GetPhenotype();
      phenotypes[p] += 1;

    }
  }
  for (auto i : phenotypes)
  {
    cout << i.first << "\t" << i.second << endl;
  }
}


void CellularPotts::PrintColours()
{
  unordered_map<int, int> colours{};

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
      int p = c->GetPhenotype();
      colours[p] = c->Colour();

    }
  }
  for (auto i : colours)
  {
    cout << "Colour of cell number: " << i.first << "  is: " << i.second << endl;
  }

}

void CellularPotts::SetColours()
{
  // map<int, int> colours = {{20233, 102}, {79363, 77}, {20235, 10}, {78463, 9}}; // this is for asym2

  // map<int, int> colours = {{61443, 81}, {28675, 91}, {129663, 10}, {64003, 9}, {64087, 28}}; // this is for pluri6  have to do 


  //this is for pluri46
  // map<int, int> colours = {{7427, 106}, {6527, 52}, {7675, 95}, {7531, 12}, {7659, 11}, {6651, 47} , {75137, 95} };
  
  // map<int, int> colours = {{63559, 109}, {116099, 91}, {11619, 10}, {59523, 88}, {59579, 137}, {63491, 55}, {67522, 5}, {118147, 85}, {59555, 73},
  // {126339, 56}, {118147, 56}, {83906, 47}}; // this is for arrow
  // 1160099 - connects layers 2 and 3, 63599 - cells at bottom, 118147 = cells center of top



  // map<int, int> colours = {{51331, 103}, {103235, 91}, {103171, 12}, {51699, 88}, 
  // {105227, 29}, {18443, 37}, {103234, 27}, {51339, 81}}; // this is for jellyfish

  // map<int,int> colours = {{93571, 81}, {23919, 6}, {23811, 73}, {19843, 102}, {19911, 4}}; // this is for megamind

  // map<int,int> colours = {{102403, 27}, {111357, 106}, {110595, 91}, {111105, 6}, {111107, 107}}; // this is for star??

  // map<int,int> colours = {{119675, 45}, {110595, 143},{115579, 84},{45059, 71}, {127867, 88}}; // this is for shield


<<<<<<< HEAD
  // map<int,int> colours = {{25600, 107}, {28160, 88}, {32256, 25}, {91136, 56}, 
  // {15914, 66}, {15874, 85}, {11947, 22}, {11907, 103}, {16130, 66}, {16186, 118}}; // this is for mushroom

  map<int,int> colours = {{108034, 252}, {123107, 254}, {123011, 253}, {123043, 255}, 
  {107010, 249}, {115075, 250}, {107651, 251}}; // fungi
=======
  // map<int,int> colours = {{25600, 107}, {28160, 88}, {32256, 25}, {91136, 56}, {15914, 66}, {15874, 85}, {11947, 22}, {11907, 103}, {16130, 66}, {16186, 118}}; // this is for mushroom-old


  map<int,int> colours = {{91137, 256}, {25601, 257}, {27651, 258}, {31747, 259}, {25603, 260}, {25602, 261}, {11266, 262}, {16130, 263}, {16186, 264}, {16314, 265}, {11962, 266}, {11906, 267}}; // this is for mushroom
>>>>>>> 254dd4f910c4256fe5c6b329d4f5b30a9d111066


  // map<int,int> colours = {{129287, 49}, {129295, 102}, {129293, 81}, {129039, 160}, 
  // {64771, 53}, {93951, 22}, {60547, 91}, {92927, 47}, {27779, 73}}; // this is for patterns

  // map<int,int> colours = {{59395, 13}, {65027, 17}, {64515, 14}, {51327, 5}, {63491, 12}, {59519, 7}, {63615, 16}, {59511, 22} }; //this is for pluri13

  // map<int, int> colours = {{45059, 122}, {31347, 118}, {31345, 52}}; // this is for inducer2

  // map<int,int> colours = {{16707, 91}, {11963, 72}, {11907, 25}};  // this is for asym4

  // map<int,int> colours = {{92811, 117}, {64892, 47}, {31747, 106}, {76427, 119}};  // this is for asym8



  // 45059, 3137, 31345

  // 103 is deep red, 13 is dark purle, 102 is fluoro purle, 81 is fluro yellow/green, 55 is grey, 85 is light green, 91 is forest green
  // 88 is light purple, 9 is beige, 73 is pinky red, 107 is camo, 108 is super light green, 109 is bluey green, 12 is boring, 81 also super light green
  // 8 is light black, 6 is light blue, 106 is light yellow, 10 is a normal green, 9 is dark grey, 22 is yellow/orange
  // 66 is dark purple, 72 is deep red, 95 is fluoro purple, 56 is light black, 25 is light grey, 26 is dark black,, 44 is white, 45 is light green
  // 46 is dark grey, 47 is sky blue, 48, 49, 50 is grey, 51 is black, 52 is red/pink, 53 is black, 112 is fluor green
  // 115 grey-red, 116 is dark grey, 117 is purple-blue, 118 is fluoro blue, 119 is deepy-light blue

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
      int p = c->GetPhenotype();

      map<int,int>::iterator it;
      it = colours.find(p);
      if (it != colours.end())
      {
        c->set_ctype(it->second);
      }
    }
  }
}




  void CellularPotts::CountTypesTime(void)
  {
    // iterate through cells, add to vector. Call this after "update network"


    map<int,int> new_list{};

    vector<Cell>::iterator c;
    for ((c=cell->begin(), c++); c!=cell->end(); c++)
    {
      if (c->AliveP())
      {
        c->Phenotype();
        int p = c->GetPhenotype();
        new_list[p] += 1;
      }
    }
    
    TypeCounts.push_back(new_list);
  }



  void CellularPotts::PrintTypesTime(bool prune)
  {

    // We are going to want to order everything correctly.
    // Need to create a new list where things are now ordered by state->time instead of time->state
    // then we can remove the transition types.

    map<int,int> ordered_types{};


    vector<vector<int>> timestate{};

    for (auto i : TypeCounts)
    {
      for (auto j : i)
      {
        ordered_types[j.first] += j.second;
        // if (find(ordered_types.begin(), ordered_types.end(), j.first) == ordered_types.end())
        // {
        //   ordered_types[]
        // }
      }
    }
    

    vector<vector<int>> stemdiff;


    // create new vector of small states, remove from ordered types
    vector<int> to_remove{};

    if (prune)
      for (auto it = ordered_types.begin(); it != ordered_types.end();)
      {
        if (it->second < 300)
        {
          to_remove.push_back(it->first);
          it = ordered_types.erase(it);
        }
        else
        {
          ++it;
        }



      }

    // now print - convert typecounts to timestate
    
    string var_name = data_file + "/type-counts.dat";
    ofstream outfile;
    outfile.open(var_name, ios::app);

    for (auto type : ordered_types)
    {
      outfile << type.first << '\t';
    }
    outfile << endl;



    for (auto time : TypeCounts)
    {
      for (auto type : ordered_types)
      {
        outfile << time[type.first] << '\t';

        // else
        // {
        //   outfile << 0 << '\t';
        // }
      }
      outfile << endl;
    }
    outfile.close();


    if (par.stem_counts)
    {
      var_name = data_file + "/stem-counts.dat";
      outfile.open(var_name, ios::app);


      // now do it just stem and differentiated
      map<int,int> stemtypes{};
      
      for (auto i : TypeCounts)
      {
        for (auto j : i)
        {
          stemtypes[j.first] += j.second;
          // if (find(ordered_types.begin(), ordered_types.end(), j.first) == ordered_types.end())
          // {
          //   ordered_types[]
          // }
        }
      }


      for (auto time : TypeCounts)
      {
        int stem{};
        int diff{};

        for (auto type : stemtypes)
        {
          if (type.first > 20000)
            stem += time[type.first];
          else
            diff += time[type.first];

        }
        outfile << stem << '\t' << diff << endl;
      }

    outfile.close();
    }




  }







void CellularPotts::DestroyCellsByPhenotype(int type, bool save, int type2, int type3, int type4)
{
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
    }
  }

  for (int x=1; x<sizex; ++x)
    for (int y=1; y<sizey; ++y)
    {
      if (sigma[x][y] > 0)
      {

        int p = (*cell)[sigma[x][y]].GetPhenotype();
        if (save && p != type && p != type2 && p != type3 && p != type4)
        {
          sigma[x][y] = 0;
        }
        else if (!save && (p == type || p == type2 || p == type3 || p == type4))
        {
          sigma[x][y] = 0;
        }
        
      }
    }

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
    }
  }
  cout << "Total cells killed: " << deadcells << endl;
  
}


void CellularPotts::DestroyCellsByMorphogen(int morph, double conc)
{
  // Kill all cells that are expressing the gene for a specific morphogen. 
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      vector<double>& gene_list = c->get_genes();
      if (gene_list[morph] > conc)
      {
        cout << "iterating morph... " <<  gene_list[morph] << endl;
        c->set_death_tag(true);
      }
      else 
      {
        c->set_death_tag(false);
      }
    }
  }

  for (int x=1; x<sizex; ++x)
    for (int y=1; y<sizey; ++y)
    {
      if (sigma[x][y] > 0)
      {

        if ((*cell)[sigma[x][y]].get_death_tag())
          sigma[x][y] = 0;
      }
    }

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
    }
  }
  cout << "Total cells killed: " << deadcells << endl;
    
}


void CellularPotts::DestroyCellsByRadius(double rad)
{
  double center[] = {0., 0., 0.,};
  get_center(center);

  for (int x=1;x<sizex;++x)
    for (int y = 1;y<sizey;++y)
    {
      double xdist = center[0] - x;
      double ydist = center[1] - y;
      double v = sqrt(pow(xdist, 2) + pow(ydist, 2));
     
      if (v > rad)
      {
        sigma[x][y] = 0;
      }
      // else if ((*cell)[sigma[x][y]].GetPhenotype() != 65027)
      // {
      //   sigma[x][y] = 0;
      // }
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
    }
  }
  cout << "Total cells killed: " << deadcells << endl;

}



void CellularPotts::ConvertToStem(int xloc, int yloc, int rad, int type, PDE *field, bool clear, int clear_rad)
{
  if (clear_rad == 0)
    clear_rad = rad;
  // x and y is center, rad is going to be the length (we will just make a square)

  int stem_sig{};
  vector<double> genes;
  vector<double> diffusers;
  vector<double> locks;
  vector<double> keys;
  vector<double> meds;


  // need to grab the identity of a stem cell
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
      int phen = c->GetPhenotype();
      cout << phen << endl;
      if (phen == type)
      {
        // get morphogen concentration at that point. 
        stem_sig = c->Sigma();
        genes = c->get_genes();
        diffusers = c->get_diffusers(); 
        locks = c->get_locks();
        keys = c->get_keys();
        meds = c->get_medp();

        break;
      }
    }
  }
  if (stem_sig==0)
  {
    cout << "ERROR - NO CELL WITH THAT IDENTITY FOUND" << endl;
  }

  vector<double> morph_conc(3);

  //get morphogen concentration at site of cell
  for (int x=1;x<sizex;x++)
    for (int y=1;y<sizey;y++)
    {
      if (sigma[x][y] == stem_sig)
      {
        for (int n=0;n<par.n_diffusers;++n)
        {
          morph_conc[n] = field->Sigma(n,x,y);
        }
      }
    }



  // set all cells in the differentiated area to the stem cell identity (NEED TO INCLDUE X AND Y)

  vector<int> cell_stack{};

  for (int x = xloc-rad;x < xloc + rad;++x)
    for (int y = yloc-rad;y<yloc +rad;++y)
    {
      if (sigma[x][y])
      {
        auto it = find(cell_stack.begin(), cell_stack.end(), sigma[x][y]);
        if (it == cell_stack.end())
        {
          (*cell)[sigma[x][y]].set_genes(genes);
          (*cell)[sigma[x][y]].set_diffusers(diffusers);
          (*cell)[sigma[x][y]].set_locks(locks);
          (*cell)[sigma[x][y]].set_keys(keys);
          (*cell)[sigma[x][y]].set_genes(genes);

          cell_stack.push_back(sigma[x][y]);
        }
      }

    }
  // if clearing, wash the morphogen concentration on that part of the grid
  if (clear)
    for (int x = xloc-clear_rad;x < xloc + clear_rad;++x)
      for (int y = yloc-clear_rad;y<yloc +clear_rad;++y)
        for (int n=0;n<par.n_diffusers;++n)
        {
          field->setValue(n,x,y, morph_conc[n]);
          morph_conc[n] = field->Sigma(n,x,y);
        }



  // for ((c=cell->begin(), c++); c!=cell->end(); c++)
  // {
  //   if (c->AliveP())
  //   {
  //     c->set_genes(genes);
  //     c->set_diffusers(diffusers);
  //     c->set_locks(locks);
  //     c->set_keys(keys);
  //     c->set_meds(meds);
    
  //   }    
  // }

  
  


}









// check if there are any lonely cells. Currently not in use. 
bool CellularPotts::SoloCheck()
{
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    if (c->AliveP())
    {
      bool lonely=true; 
      int id = c->Sigma();
      for (int x=1;x<sizex-2;++x)
        for (int y=1;y<sizey-2;++y)
        {
          if (sigma[x][y] == id)
          {
            if (sigma[x][y-1] > 0 && sigma[x][y-1] != id)
            {
              lonely = false;
              x = sizex;
              y = sizey;
              break;
            }
            
            if (sigma[x][y+1] > 0 && sigma[x][y+1] != id)
            {
              lonely = false;
              x = sizex;
              y = sizey;
              break;
            }
              
          }
        }
      if (lonely)
        return false;
    }
  return true;  
}


//IMPORTANT METHOD:  Function to ensure all cells are connected indirectly to all other cells on lattice.  
bool CellularPotts::CheckAllConnected()
{
  vector<unordered_set<int>> all_connections{};

  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      unordered_set<int> tempcon{};
      int id = c->Sigma();
      tempcon.emplace(id);
      int xp, yp;
      for (int x=1;x<sizex-2;++x)
        for (int y=1;y<sizey-2;++y)
        {
          if (sigma[x][y] == id)
          {
            for (int i = 1;i<=nbh_level[1];++i)
            {
              xp = x + nx[i];
              yp = y + ny[i];
              if (sigma[xp][yp] > 0 && sigma[xp][yp] != id)
                tempcon.emplace(sigma[xp][yp]);
            }
          }
        }
      all_connections.push_back(tempcon);
    }
  }

  if (all_connections.size() < 2)
    return false;


  unordered_set<int> MaxConnections{};
  for (int n : all_connections[0])
  {
    MaxConnections.emplace(n);
  }

  all_connections.erase(all_connections.begin());
  // auto it = all_connections.begin();
  // *it = move(all_connections.back());
  // all_connections.pop_back();

  unsigned int x = all_connections.size();
  for (unsigned int i=0; i<x;++i)
  {
    bool fbreak = false;
    unordered_set<int> nbrs = all_connections[i];
    for (int connection : nbrs)
    {
      unordered_set<int>::iterator f = MaxConnections.find(connection);
      if (f == MaxConnections.end())
        continue;
      else
      {
        for (int connection : nbrs)
        {
          MaxConnections.emplace(connection);
        }
        all_connections.erase(all_connections.begin() + i);
        
        // auto it = all_connections.begin() + i;
        // *it = move(all_connections.back());
        // all_connections.pop_back();
        fbreak = true;
        break;   
      }
    }
    if (fbreak == true)
    {
      x = all_connections.size();
      if (x > 0)
        i = 0;
      else
        break;
    }
  }

  if ((int)MaxConnections.size() == CountCells())
    return true;
  else
    return false;
}


bool CellularPotts::CheckShape()
{
  if (ShapeMaintained)
  {

    bool connected = CheckAllConnected();
    // double avg_bind = AverageBinding();
    // double key_check = AvgMedsOn();

    if (!connected)
    {
      ShapeMaintained = false;
    }

    // if (avg_bind < 4)
    // {
    //   ShapeMaintained = false;
    // }

    // if (key_check < 1.6)
    // {
    //   ShapeMaintained = false;
    // }

    if (TotalArea() < 3000)
    {
      ShapeMaintained = false;
    }
  }
  return ShapeMaintained;

}



void CellularPotts::SetCellCenters()
{
  map<int, double> xvals;
  map <int, double> yvals;


  // get center of mass for all cells. Probably faster than searching through the vector every time! Turn this into function
  for (int x=1;x<sizex;++x)
    for (int y=1;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
      {
        xvals[sigma[x][y]] += x;
        yvals[sigma[x][y]] += y;
      }
    }
  int count=1;
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      xvals[count] = xvals[count] / double(c->Area());
      yvals[count] = yvals[count] / double(c->Area());

    }
    ++count;
  }

  // give each cell its x and y coords
  map<int, double>::iterator iter;
  for (iter = xvals.begin(); iter != xvals.end(); iter++)
  {
    (*cell)[iter->first].set_xcen(iter->second);
    // cout << iter->first << " : " << iter->second << endl;
  }
  for (iter = yvals.begin(); iter != yvals.end(); iter++)
  {
    (*cell)[iter->first].set_ycen(iter->second);
    //cout << iter->first << " : " << iter->second << endl;
  }

}

void CellularPotts::RecordMasses()
{
  SetCellCenters();
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      c->RecordMass();
    }
  }
}

void CellularPotts::CellVelocities()
{
  // find the difference between every recording. 
  // need to do simple vector length calculation.


  string fnamen = data_file + "/velocities";

  if (mkdir(fnamen.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;  

  if (mkdir(data_file.c_str(), 0777) == -1)
    cerr << "Error : " << strerror(errno) << endl;
  else
    cout << "Directory created." << endl;



  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
  {
    if (c->AliveP())
    {
      string var_name = data_file + "/velocities/cell_" + to_string(c->Sigma()) + ".dat";
      ofstream outfile;
      outfile.open(var_name, ios::app);


      vector<double>& xm = c->get_xcens();
      vector<double>& ym = c->get_ycens();
      int s = xm.size();

      for (int i = 10; i < s; ++i)
      {
        // we want displacement from a while ago to account for back and forth motion
        double x = xm[i-10];
        double y = ym[i-10];
        double x1 = xm[i];
        double y1 = ym[i];

        double len = sqrt(pow(x1-x,2) + pow(y1-y,2));

        // cout << "length is: " << len << endl;
        outfile << len << endl;

      }

      outfile.close();


    }
  }
}











void CellularPotts::DrawListofCAC(Graphics *g, vector<array<int,2>> cac)
{
  for (unsigned int i=0;i<cac.size();++i)
  {
    g->Point(12, cac[i][0]*2, cac[i][1]*2);
    g->Point(12, cac[i][0]*2+1, cac[i][1]*2+1);
    g->Point(12, cac[i][0]*2+1, cac[i][1]*2-1);
    g->Point(12, cac[i][0]*2+1, cac[i][1]*2);
    g->Point(12, cac[i][0]*2-1, cac[i][1]*2+1);
    g->Point(12, cac[i][0]*2-1, cac[i][1]*2);
    g->Point(12, cac[i][0]*2-1, cac[i][1]*2+1);
    g->Point(12, cac[i][0]*2, cac[i][1]*2-1);
    g->Point(12, cac[i][0]*2, cac[i][1]*2+1);
  }
}




void CellularPotts::DrawPerimeter(Graphics *g, vector<int> pcells)
{

  for (unsigned int i=0;i<pcells.size()-1;++i)
  {
    int x0 = 2*(*cell)[pcells.at(i)].xcen;
    int y0 = 2*(*cell)[pcells.at(i)].ycen;
    int x1 = 2*(*cell)[pcells.at(i+1)].xcen;
    int y1 = 2*(*cell)[pcells.at(i+1)].ycen;
    g->Line(x0, y0, x1, y1, 1);
  }
}



void CellularPotts::PerimeterGrid()
{
  // clear the grid
  for (int i=0;i<sizex*sizey;++i)
  {
    outside[0][i] = 0;
  }

  // search with neighbour distance two:
  vector<array<int,2>> new_neighbours; // store the neighbours that are outside medium. MEDIUm = 2, PERIMETER = 1
  array<int,2> corner = {1,1};
  new_neighbours.push_back(corner);
  outside[1][1] = 2;

  long count = 0;
  while (new_neighbours.size() > 0)
  {
    
    // number of old neighbours to delete after
    int vecsize = new_neighbours.size();
    for (int i = 0; i < vecsize;++i)
    {
      array<int,2> p = new_neighbours.at(i);
      for (int j=1;j<=nbh_level[2];++j)
      {
        int xp2,yp2;
        xp2=p[0]+nx[j]; 
        yp2=p[1]+ny[j];
        // cout << xp2 << " " << yp2 << endl;
        if (sigma[xp2][yp2] == 0 && outside[xp2][yp2] != 2)
        {
          array<int,2> nbh = {xp2, yp2};
          new_neighbours.push_back(nbh);
          outside[xp2][yp2] = 2;
          ++count;
        }
      }
    }
    // remove old neighbours from new_neighbours for the next neighbour search
    new_neighbours.erase(new_neighbours.begin(), new_neighbours.begin()+vecsize);
  }
}


vector<array<int, 2>> CellularPotts::PerimeterCAC()
{
  vector<array<int, 2>> gridp{};
  // cant simply iterate through grid because there are concealed areas...


  PerimeterGrid();

  for (int x=1;x<sizex-1;++x)
    for (int y=1;y<sizey-1;++y)
    {
      if (sigma[x][y] > 0)
      {
        int xp2, yp2;
        for (int j=1;j<=nbh_level[2];++j)
        {
          xp2=x+nx[j]; 
          yp2=y+ny[j];
          if (outside[xp2][yp2] == 2)
          {
            outside[x][y] = 1;
            array<int,2> newp = {x, y};
            gridp.push_back(newp);
            break;
          } 
        } 
      }
    } 
  return gridp;
}



vector<int> CellularPotts::CellsFromCAC(vector<array<int,2>> cac)
{

  // Now we have all perimeter grid CA cells in a vector. Need to find the cells to which they belong.
  vector<int> pcells;

  for (array<int,2> s : cac)
  {
    int c = sigma[s[0]][s[1]];
    auto it = find(pcells.begin(), pcells.end(), c);
    if (it == pcells.end())
    {
      pcells.push_back(c);
    }
  }
  
  // now we have an unordered list. Need to make them all connected. 

  // Find neighbours --> Check for neighbours in list --> add to chain (How do I know chain direction?)
  return pcells;
}



vector<vector<int>> CellularPotts::CellNeighbours(vector<int> cell_list) 
{
  vector<vector<int>> nbh_vector{};
  nbh_vector.resize(cell_list.size());

  for (int x = 1; x<sizex - 1; ++x)
    for (int y = 1; y<sizey - 1; ++y)
    {
      if (outside[x][y] == 1)
      {
        for (int i = 1; i <= nbh_level[1];++i)
        {
          int xp = x + nx[i];
          int yp = y + ny[i];

          if (sigma[x][y] != sigma[xp][yp] && outside[xp][yp] == 1)
          {
            int nb = sigma[xp][yp];
            int nc = sigma[x][y];
            auto it = find(cell_list.begin(), cell_list.end(), nc);
            int val = it - cell_list.begin();
            auto iter = find(nbh_vector.at(val).begin(), nbh_vector.at(val).end(), nb);
            if (iter == nbh_vector.at(val).end())
            {
              nbh_vector.at(val).push_back(nb);
            }
          }
        }        
      }
    }

  return nbh_vector;

}  


vector<int> CellularPotts::LinkPerimeter()
{

  // ensure cells have correct mass centers. 
  SetCellCenters();

  // How about this.. Iterate around the perimeter to link it, and then iterate around those cells.
  vector<int> cell_list = CellsFromCAC(PerimeterCAC());  

  vector<int> ordered_list{};

  vector<vector<int>> nbs = CellNeighbours(cell_list);
  // the iterator of nbs is equal to the key in cell_list. 

  ordered_list.push_back(cell_list.front());
  // cell_list.erase(cell_list.begin());


  auto it = cell_list.begin();
  int veclen = it - cell_list.begin();
  int cn = cell_list.front();

  bool finished = false;

  while (!finished)
  {
    int next{};
    double minvec = sizex;


    int xcen = (*cell)[cn].xcen;
    int ycen = (*cell)[cn].ycen;

    vector<int> &nbh = nbs.at(veclen); 
    
    for (int &i : nbh)
    {
      auto ifused = find(ordered_list.begin(), ordered_list.end(), i);
      if (ifused == ordered_list.end())
      {
        auto newit = find(cell_list.begin(), cell_list.end(), i);
        int vlen = newit - cell_list.begin();
        if (nbs.at(vlen).size() > 1)
        {
          int xdist = (*cell)[i].xcen - xcen;
          int ydist = (*cell)[i].ycen - ycen;
          
          double vec = sqrt(pow(xdist, 2) + pow(ydist, 2));
          if (vec < minvec)
          {
            minvec = vec;
            next = i;
          }
        }
      }
    }
    if (next)
    {
      it = find(cell_list.begin(), cell_list.end(), next);
      veclen = it - cell_list.begin();
      cn = next;
      ordered_list.push_back(cn);
    }
    else
    {
      break;
    }
  }
  ordered_list.push_back(ordered_list.front());

  return ordered_list;

}




void CellularPotts::CellExposure()
{
  // reset cells
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    if (c->AliveP())
    {
      c->resetCAC();
    }


  for (int x=1;x<sizex;++x)
    for (int y=1;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
      {
        int count = 0;
        bool border = false;
        for (int i=0;i<nbh_level[1];++i)
        {
          int xp2,yp2;
          xp2=x+nx[i]; 
          yp2=y+ny[i];

          if (!sigma[xp2][yp2])
          {
            if (outside[x][y] == 1)
            {
              (*cell)[sigma[x][y]].cellperim();
              break;
            }  
          }
          else if (sigma[xp2][yp2] != sigma[x][y])
          {
            border = true;
            ++count;
          }
          else 
          {
            ++count;
          }
        }
        if (count == nbh_level[1] && border)
          (*cell)[sigma[x][y]].cellcell();
      }
    }


  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    if (c->AliveP())
    {
      // cout << "perim CAC: " << c->cell_perim << "  cell-cell CAC: " << c->cell_contact << endl;
      if ((double)c->cell_perim * 0.5 > c->cell_contact)
        c->exposed = true;
      else
        c->exposed = false;
    }
    
}



long CellularPotts::TotalArea()
{
  int count=0;
  for (int x=0;x<sizex;++x)
    for (int y=0;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
        ++count;
    }
  return count;
}




void CellularPotts::get_center(double* center)
{

  // changing algorithm to compute exact center
  int xtotal{};
  int ytotal{};
  int mass{};

  for (int x=1;x<sizex;++x)
    for (int y=1;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
      {
        xtotal += x;
        ytotal += y;
        ++mass;
      }
    }
  
  center[0] = (double)xtotal / mass;
  center[1] = (double)ytotal / mass;


  // int minx = sizex;
  // int miny = sizey;
  // int maxx = 0;
  // int maxy = 0;

  // for (int i=0;i<sizex-1;++i)
  //   for (int j=0;j<sizey-1;++j)
  //   {
  //     if (sigma[i][j] > 0)
  //     {
  //       if (i < minx)
  //         minx = i;
  //       else if (i > maxx)
  //         maxx = i;

  //       if (j < miny)
  //         miny = j;
  //       else if (j > maxy)
  //         maxy = j;
  //     }
  //   }
  // calculate approximate center:
  // center[0] = double(maxx + minx) / 2.;
  // center[1] = double(maxy + miny) / 2.;

  // get average distance to get stdev for gradient (/2 = radius, /4 = both x and y, / 12 to get 1/8 of distance)
  // center[2] = double((maxx - minx) + (maxy - miny)) / 8.;

}



double CellularPotts::BresenhamForCircle(int x1, int y1, int x2, int y2, int dx, int dy, int decide, double *center, double rad)
{
	// pk is initial decision making parameter
	// Note:x1&y1,x2&y2, dx&dy values are interchanged
	// and passed in plotPixel function so
	// it can handle both cases when m>1 & m<1
	int pk = 2 * dy - dx;

  bool encounter{false};
  int maxy{0};
  int miny{sizey};
  int xmin{};
  int xmax{};

  double dev{};

	for (int i = 0; i <= dx; i++) 
  {
    if (x1 > rows || y1 > cols)
        break;

    if (decide == 0)
    {
      // cout << x1 << "," << y1 << endl;
      if (sigma[x1][y1] > 0)
      {
        if (y1>maxy)
        {
          maxy = y1;
          xmax = x1;
          encounter = true;
        }

        if (y1 < miny)
        {
          miny = y1;
          xmin = x1;
        }
      }
    }
    else
    {
      // cout << y1 << "," << x1 << endl;
      if (sigma[y1][x1] > 0)
      {
        if (x1>maxy)
        {
          maxy = x1;
          xmax = y1;
          encounter = true;
        }
        
        if (x1 < miny)
        {
          miny = x1;
          xmin = y1;
        }
      }
        
    }

		// checking either to decrement or increment the
		// value if we have to plot from (0,100) to (100,0)
		x1 < x2 ? x1++ : x1--;
		if (pk < 0) 
    {
			// decision value will decide to plot
			// either x1 or y1 in x's position
			if (decide == 0) 
      {
				// putpixel(x1, y1, RED);
				pk = pk + 2 * dy;
			}
			else 
      {
				//(y1,x1) is passed in xt
				// putpixel(y1, x1, YELLOW);
				pk = pk + 2 * dy;
			}
		}
		else 
    {
			y1 < y2 ? y1++ : y1--;
			if (decide == 0) {

				// putpixel(x1, y1, RED);
			}
			else {
				// putpixel(y1, x1, YELLOW);
			}
			pk = pk + 2 * dy - 2 * dx;
        
		}
	}
  if (encounter)
  {
    double vec1 = sqrt(pow((double)xmax-center[0], 2) + pow((double)maxy-center[1], 2));
    double vec2 = sqrt(pow((double)xmin-center[0], 2) + pow((double)miny-center[1], 2));

    dev += abs(vec1 - rad); // / rad;
    dev += abs(vec2 - rad); // / rad;
    ++rad_count;
    // cout << "xmax: " << xmax << "  xmin: " << xmin << "  maxy: " << maxy << "  miny: " << miny << "  center[0]: " << center[0] << "  center[1]" << center[1] << endl;
    return dev;
  }
  return 0;
}





double CellularPotts::CircleNegGrad(double m, double *center, double rad)
{
  double dev{};

  int x1, y1, x2, y2, dx, dy;

  double skipx=0.;
  double skipy=0.;

  if (m == 0)
      skipx = rows;
  else if (m > -1)
      skipx = (abs(1 / m));
  else if (m < -1)
      skipy = (abs(m));
  
  double error = 0;


  for (int line = 1; line <= rows+cols;line++)
  {
      int skip{};

      // only skipx or skipy can be greater than 0 at a time. 
      // Need to calculate the error on non-integer gradients to calculate correct jumps
      if (line < rows && skipx > 0)
      {
          skip = floor(skipx);
          error += skipx - skip;

          while (error >= 1)
          {
              ++skip;
              --error;
          }
          line += skip - 1;
          
      }
      else if (line > rows && skipy > 0)
      {
          skip = floor(skipy);
          error += skipy - skip;

          while (error >= 1)
          {
              ++skip;
              --error;
          }
          line += skip - 1;
      }

      y1 = nmax(0, line-rows);
      x1 = minu(rows, line);
      x2 = 0; 
      y2 = round(m * (x2 - x1) + y1);

      dx = abs(x2 - x1);
      dy = abs(y2 - y1);

      // If slope is less than one
      if (dx > dy) 
      {
        // passing argument as 0 to plot(x,y)
        dev += BresenhamForCircle(x1, y1, x2, y2, dx, dy, 0, center, rad);

        // plotPixel(x1, y1, x2, y2, dx, dy, 0);
      }

      // if slope is greater than or equal to 1
      else 
      {
        // passing argument as 1 to plot (y,x)
        dev += BresenhamForCircle(y1, x1, y2, x2, dy, dx, 1, center, rad);

        // plotPixel(y1, x1, y2, x2, dy, dx, 1);
      }
  }
  // cout << " GRADIENT = " << m << "  returned circle deviation of: " << dev << endl;
  // griditcount=0;
  return dev;
}





// positive gradient
double CellularPotts::CirclePosGrad(double m, double *center, double rad)
{
  double dev{};
  int x1, y1, x2, y2, dx, dy;

  double skipx=0;
  double skipy=0;

  if (m == 0)
      skipx = cols;
  else if (m > 1)
      skipy = abs(m); 
  else if (m < 1)
      skipx = abs(1 / m);
  
  double error{};

  // only skipx or skipy can be greater than 0 at a time. 
  // Need to calculate the error on non-integer gradients to calculate correct jumps
  for (int line = rows + cols; line >= 0;line--)
  {
    int skip{};

    if (line > rows && skipy > 0)
    {
        skip = floor(skipy);
        error += skipy - skip;

        while (error >= 1)
        {
            ++skip;
            --error;
        }
        line -= skip - 1;
    }
    else if (line < rows && skipx > 0)
    {
        skip = floor(skipx);
        error += skipx - skip;

        while (error >= 1)
        {
            ++skip;
            --error;
        }
        line -= skip - 1;
    }

    y1 = nmax(0,line-rows);
    x1 = nmax(0, cols-line);
    y2 = rows;
    x2 = round((y2-y1)/m + x1);


    dx = abs(x2 - x1);
    dy = abs(y2 - y1);

    // If slope is less than one
    if (dx > dy) 
    {
      // passing argument as 0 to plot(x,y)
      dev += BresenhamForCircle(x1, y1, x2, y2, dx, dy, 0, center, rad);
    }

    // if slope is greater than or equal to 1
    else 
    {
      // passing argument as 1 to plot (y,x)
      dev += BresenhamForCircle(y1, x1, y2, x2, dy, dx, 1, center, rad);
    } 

  }
  // cout << "GRADIENT = " << m << "  returned circle deviation of: " << dev << endl;
  // griditcount = 0;
  return dev;
}




double CellularPotts::DeviationFromCircle()
{
  long cellmass = TotalArea();
  
  rad_count = 0;

  //calculate radius if perfect circle
  double rad = sqrt((double)(cellmass) / M_PI);

  // calculate approximate centroid based only on x and y
  double center[] = {0., 0., 0.,};
  get_center(center);
 
 
  // for (int i=0;i<3;++i)
  //   cout << "center val: " << center[i] << endl;


  // deviation from circle, calculated by summing radius differnce
  double dev{};
  double xdev{};
  double ydev{};
  

  for (double i : n_grads)
  {
    dev += CircleNegGrad(i, center, rad);
  }

  for (double i : p_grads)
  {
    dev += CirclePosGrad(i, center, rad);
  }

  // algorithm is wrong, need to draw lines to radius based on degrees. 
  // angle -> line, multiplied by radius to get x and y val, get actual x and y val (can draw bresenhaum line) and then compute squared vector difference

  for (int x=0;x<sizex;++x)
  {
    int maxy=0;
    int miny=sizey;
    
    bool if_sig = false;
    for (int y=0;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
      {
        if (y > maxy)
        {
          maxy = y;
          if_sig = true;
        }
        
        if (y < miny)
          miny = y;
      }        
    }

    if (if_sig)
    {
      // calculate vector based on max y for every x
      double vec1 = sqrt(pow((double)x-center[0], 2) + pow((double)maxy-center[1], 2));
      double vec2 = sqrt(pow((double)x-center[0], 2) + pow((double)miny-center[1], 2));

      // Divide by radius to normalise by size (dividing by sqrt(rad) increases values for large radius)
      xdev += abs(vec1 - rad);// rad;
      xdev += abs(vec2 - rad);// rad;
      ++rad_count;
    }
  }
  // and do the same for x across y
  for (int y=0;y<sizey;++y)
  {
    int maxx=0;
    int minx=sizex;
    
    bool if_sig = false;
    for (int x=0;x<sizex;++x)
    {
      if (sigma[x][y] > 0)
      {
        if (x > maxx)
        {
          maxx = x;
          if_sig = true;
        }
        
        if (x < minx)
          minx = x;
      }        
    }

    if (if_sig)
    {
      // calculate vector based on max and min x for every y
      double vec1 = sqrt(pow((double)y-center[0], 2) + pow((double)maxx-center[1], 2));
      double vec2 = sqrt(pow((double)y-center[0], 2) + pow((double)minx-center[1], 2));

      ydev += abs(vec1 - rad);// rad;
      ydev += abs(vec2 - rad);// rad;
      ++rad_count;
    }
  }
  // divide by the number of iterations, but multiply by cellmass to select for larger organisms. 
  // dev = double(cellmass) * (dev / (double)n);
  dev += xdev;
  dev += ydev;

  // cout << "GRADIENT = " << 0 << "  returned circle deviation of: " << xdev << endl;
  // cout << "GRADIENT = inf" << "  returned circle deviation of: " << ydev << endl;

  dev = (dev / rad_count); // * sqrt((double)cellmass);

  // if (par.print_fitness)
  // {
  //   cout << "Deviation from circle: " << dev << "  Total mass: " << cellmass << endl;
  // }

  return dev;  
}


void CellularPotts::DeviationCheck()
{
  EarlyWhiteSpace();

  if (par.print_fitness)
  {
    if (early_contig)
    {
      cout << "Organism Rule Violated!! Contig found!" << endl;
    }
  }
}






void CellularPotts::plotPixel(int x1, int y1, int x2, int y2, int dx,
			int dy, int decide)
{
	// pk is initial decision making parameter
	// Note:x1&y1,x2&y2, dx&dy values are interchanged
	// and passed in plotPixel function so
	// it can handle both cases when m>1 & m<1
	int pk = 2 * dy - dx;
	for (int i = 0; i <= dx; i++) 
    {
        if (x1 > rows || y1 > cols)
            break;
        // if (decide == 0)
		//     cout << x1 << "," << y1 << endl;
        // else
        //     cout << y1 << "," << x1 << endl;
        ++griditcount;
		// checking either to decrement or increment the
		// value if we have to plot from (0,100) to (100,0)
		x1 < x2 ? x1++ : x1--;
		if (pk < 0) {
			// decision value will decide to plot
			// either x1 or y1 in x's position
			if (decide == 0) {
				// putpixel(x1, y1, RED);
				pk = pk + 2 * dy;
			}
			else {
				//(y1,x1) is passed in xt
				// putpixel(y1, x1, YELLOW);
				pk = pk + 2 * dy;
			}
		}
		else {
			y1 < y2 ? y1++ : y1--;
			if (decide == 0) {

				// putpixel(x1, y1, RED);
			}
			else {
				// putpixel(y1, x1, YELLOW);
			}
			pk = pk + 2 * dy - 2 * dx;
        
		}
	}
}




double CellularPotts::BresenhamLine(int x1, int y1, int x2, int y2, int dx, int dy, int decide, int m_contig)
{
	// pk is initial decision making parameter
	// Note:x1&y1,x2&y2, dx&dy values are interchanged
	// and passed in plotPixel function so
	// it can handle both cases when m>1 & m<1
	int pk = 2 * dy - dx;

  bool encounter{false};
  int contig{};
  double space{};

	for (int i = 0; i <= dx; i++) 
  {
    if (x1 > rows || y1 > cols)
        break;

    if (decide == 0)
    {
      // cout << x1 << "," << y1 << endl;
      if (sigma[x1][y1] > 0)
      {
        encounter = true;
        if (contig >= m_contig)
        {
          space += sqrt((double)contig);
          ++n_segments;
          // cout << "CONTIG FOUND ON 1st DIAGONAL OF SIZE: " << contig << endl;
        }
        contig = 0;
      }
      else if (encounter == true && sigma[x1][y1] == 0)
      {
        ++contig;
      }
      
    }
    else
    {
      // cout << y1 << "," << x1 << endl;
      if (sigma[y1][x1] > 0)
      {
        encounter = true;
        if (contig >= m_contig)
        {
          space += sqrt((double)contig);
          ++n_segments;
          // cout << "CONTIG FOUND ON DIAGONAL OF GRADIENT ?? of size: " << contig << endl;
        }
        contig = 0;
      }
      else if (encounter == true && sigma[y1][x1] == 0)
      {
        ++contig;
      }          
    }
    // ++griditcount; //  COUNTING HOW MUCH OF THE GRID WE GET THROUGH


		// checking either to decrement or increment the
		// value if we have to plot from (0,100) to (100,0)
		x1 < x2 ? x1++ : x1--;
		if (pk < 0) 
    {
			// decision value will decide to plot
			// either x1 or y1 in x's position
			if (decide == 0) 
      {
				// putpixel(x1, y1, RED);
				pk = pk + 2 * dy;
			}
			else 
      {
				//(y1,x1) is passed in xt
				// putpixel(y1, x1, YELLOW);
				pk = pk + 2 * dy;
			}
		}
		else 
    {
			y1 < y2 ? y1++ : y1--;
			if (decide == 0) {

				// putpixel(x1, y1, RED);
			}
			else {
				// putpixel(y1, x1, YELLOW);
			}
			pk = pk + 2 * dy - 2 * dx;
        
		}
	}
  return space;
}






double CellularPotts::neg_grad(double m, int m_contig)
{
  double whspace{};

  int x1, y1, x2, y2, dx, dy;

  double skipx=0.;
  double skipy=0.;

  if (m == 0)
      skipx = rows;
  else if (m > -1)
      skipx = (abs(1 / m));
  else if (m < -1)
      skipy = (abs(m));
  
  double error = 0;


  for (int line = 1; line <= rows+cols;line++)
  {
      int skip{};

      // only skipx or skipy can be greater than 0 at a time. 
      // Need to calculate the error on non-integer gradients to calculate correct jumps
      if (line < rows && skipx > 0)
      {
          skip = floor(skipx);
          error += skipx - skip;

          while (error >= 1)
          {
              ++skip;
              --error;
          }
          line += skip - 1;
          
      }
      else if (line > rows && skipy > 0)
      {
          skip = floor(skipy);
          error += skipy - skip;

          while (error >= 1)
          {
              ++skip;
              --error;
          }
          line += skip - 1;
      }

      y1 = nmax(0, line-rows);
      x1 = minu(rows, line);
      x2 = 0; 
      y2 = round(m * (x2 - x1) + y1);

      dx = abs(x2 - x1);
      dy = abs(y2 - y1);

      // If slope is less than one
      if (dx > dy) 
      {
        // passing argument as 0 to plot(x,y)
        whspace += BresenhamLine(x1, y1, x2, y2, dx, dy, 0, m_contig);

        // plotPixel(x1, y1, x2, y2, dx, dy, 0);
      }

      // if slope is greater than or equal to 1
      else 
      {
        // passing argument as 1 to plot (y,x)
        whspace += BresenhamLine(y1, x1, y2, x2, dy, dx, 1, m_contig);

        // plotPixel(y1, x1, y2, x2, dy, dx, 1);
      }
  }
  // cout << " GRADIENT = " << m << "  returned white space of: " << whspace << endl;
  // griditcount=0;
  return whspace;
}





double CellularPotts::pos_grad(double m, int m_contig)
{
  double whspace{};
  int x1, y1, x2, y2, dx, dy;

  double skipx=0;
  double skipy=0;

  if (m == 0)
      skipx = cols;
  else if (m > 1)
      skipy = abs(m); 
  else if (m < 1)
      skipx = abs(1 / m);
  
  double error{};

  // only skipx or skipy can be greater than 0 at a time. 
  // Need to calculate the error on non-integer gradients to calculate correct jumps
  for (int line = rows + cols; line >= 0;line--)
  {
    int skip{};

    if (line > rows && skipy > 0)
    {
        skip = floor(skipy);
        error += skipy - skip;

        while (error >= 1)
        {
            ++skip;
            --error;
        }
        line -= skip - 1;
    }
    else if (line < rows && skipx > 0)
    {
        skip = floor(skipx);
        error += skipx - skip;

        while (error >= 1)
        {
            ++skip;
            --error;
        }
        line -= skip - 1;
    }

    y1 = nmax(0,line-rows);
    x1 = nmax(0, cols-line);
    y2 = rows;
    x2 = round((y2-y1)/m + x1);


    dx = abs(x2 - x1);
    dy = abs(y2 - y1);

    // If slope is less than one
    if (dx > dy) 
    {

      // passing argument as 0 to plot(x,y)
      whspace += BresenhamLine(x1, y1, x2, y2, dx, dy, 0, m_contig);


      // plotPixel(x1, y1, x2, y2, dx, dy, 0);
    }

    // if slope is greater than or equal to 1
    else 
    {

      // passing argument as 1 to plot (y,x)
      whspace += BresenhamLine(y1, x1, y2, x2, dy, dx, 1, m_contig);

      //plotPixel(y1, x1, y2, x2, dy, dx, 1);
    } 

  }
  // cout << "GRADIENT = " << m << "  returned white space of: " << whspace << endl;
  // griditcount = 0;
  return whspace;
}







double CellularPotts::WhiteSpace()
{
  // This function iterates across x,y, and both diagonals to determine medium to get a gauge on curvature. 
  double space{};

  n_segments = 0;

  for (double i : n_grads)
  {
    space += neg_grad(i, par.min_contig);
  }

  for (double i : p_grads)
  {
    space += pos_grad(i, par.min_contig);
  }
  

  // do iteration normal way for 0 degrees and 90 degrees because drawing lines is more computationally expensive. 
  double xspace{};
  double yspace{};


  for (int x=1;x<rows;++x)
  {
    int contig{};
    bool enc{false};
    for (int y=1;y<cols;++y)
    {
      if (sigma[x][y] > 0)
      {
        enc = true;
        if (contig >= par.min_contig)
        {
          xspace += sqrt((double)contig);
          ++n_segments;
          // cout << "CONTIG FOUND ON y straight with length: " << contig << endl;
        }
        contig = 0;
      }
      else if (enc == true && !sigma[x][y])
      {
        ++contig;
      }
    }
  }
  // Repeat process for x iteration across y.
  for (int y=1;y<cols;++y)
  {
    int contig{};
    bool enc{false};
    for (int x=1;x<rows;++x)
    {
      if (sigma[x][y] > 0)
      {
        enc = true;
        if (contig >= par.min_contig)
        {
          yspace += sqrt((double)contig);
          ++n_segments;
          // cout << "CONTIG FOUND ON x straight with length: " << contig << endl;
        }
        contig = 0;
      }
      else if (enc == true && !sigma[x][y])
      {
        ++contig;
      }
    }
  }
  // cout << "GRADIENT = 0" << " returned white space of: " << xspace << endl;
  // cout << "GRADIENT = inf" <<  " returned white space of: " << yspace << endl;

  space += xspace;
  space += yspace;

  // if (par.print_fitness)
  // {
  //   cout << "Total Whitespace: " << space << endl;
  //   cout << "Number of segments: " << n_segments << endl;
  // }


  return n_segments;


  // // Diagonal iteration
  // for (int i = 2;i <= (rows + cols - 1);i++)
  // {
  //   if (i - rows == 1)
  //     continue;

  //   int start_col = nmax(1,i - rows);
  //   int count = nmin(i-1, (cols - start_col), rows);

  //   int contig{};
  //   bool enc{false};

  //   for (int j = 0; j < count; j++)
  //   {
  //     int x = minu(rows, i) - j-1;
  //     int y = start_col + j;

  //     if (sigma[x][y] > 0)
  //     {
  //       enc = true;
  //       if (contig >= par.min_contig)
  //       {
  //         space += contig;
  //         // cout << "CONTIG FOUND ON 1st DIAGONAL OF SIZE: " << contig << endl;
  //       }
  //       contig = 0;
  //     }
  //     else if (enc == true && !sigma[x][y])
  //     {
  //       ++contig;
  //     }
  //     // cout << minu(rows, i) - j-1 << "-" << start_col + j << endl;
  //   }
  // }

  // // Rotate 90 degrees for this diagonal process
  // for (int i = rows + cols -1;i > 1; i--)
  // {
  //   if (i-cols == 1)
  //     continue;

  //   int start_row = nmax(1, i-cols);

  //   int count = nmin(i-1, (rows-start_row), cols);
    
  //   int contig{};
  //   bool enc{false};
    
  //   for (int j = 0; j < count; j++)
  //   {
  //     int x = start_row + j;
  //     int y = nmax(0, rows-i) + j + 1;
      
  //     if (sigma[x][y] > 0)
  //     {
  //       enc = true;
  //       if (contig >= par.min_contig)
  //       {
  //         space += contig;
  //         // cout << "CONTIG FOUND ON 2nd DIAGONAL OF SIZE: " << contig << endl;
  //       }
  //       contig = 0;
  //     }
  //     else if (enc == true && !sigma[x][y])
  //     {
  //       ++contig;
  //     }
  //     // cout << start_row + j << "-" << nmax(0, rows-i) + j + 1 << endl;
  //   }
  // } 
  // if (par.print_fitness)
  //   cout << "Whitespace: " << space << endl;
  // return space;
}




void CellularPotts::EarlyWhiteSpace()
{
  // This function iterates across x,y, and both diagonals to determine medium to get a gauge on curvature. 
  double space{};
  int e_cont = 15;

  n_segments = 0;

  for (double i : n_grads)
  {
    space += neg_grad(i, e_cont);
  }

  for (double i : p_grads)
  {
    space += pos_grad(i, e_cont);
  }
  
  for (int x=1;x<rows;++x)
  {
    int contig{};
    bool enc{false};
    for (int y=1;y<cols;++y)
    {
      if (sigma[x][y] > 0)
      {
        enc = true;
        if (contig >= e_cont)
        {
          space += sqrt((double)contig);
          ++n_segments;
        }
        contig = 0;
      }
      else if (enc == true && !sigma[x][y])
      {
        ++contig;
      }
    }
  }
  // Repeat process for x iteration across y.
  for (int y=1;y<cols;++y)
  {
    int contig{};
    bool enc{false};
    for (int x=1;x<rows;++x)
    {
      if (sigma[x][y] > 0)
      {
        enc = true;
        if (contig >= e_cont)
        {
          space += sqrt((double)contig);
          ++n_segments;
        }
        contig = 0;
      }
      else if (enc == true && !sigma[x][y])
      {
        ++contig;
      }
    }
  }

  if (space > 10)
  {
    early_contig = true;
  }

}



double CellularPotts::TraverseFitness()
{


  // Another option is just to test for center of mass and see how far it is shifted. 

  double center[] = {0., 0.};
  get_center(center);

  // could I just try for any asymmetry??

  double midx = ((double)sizex - 1) / 2;
  double midy = ((double)sizey - 1) / 2;

  double total_movement{};

  double dist = sqrt(pow((center[0] - midx), 2) + pow((center[1] - midy), 2));

  total_movement += pow(dist, 1.5);

  // total_movement += pow(abs(center[0] - midx), 1.5);
  // total_movement += pow(abs(center[1] - midy), 1.5);


  if (par.print_fitness)
    cout << "Asymmetry metric: " << total_movement << endl;
  

  return total_movement;

}

// higher J means less binding with medium
int MedBinding(vector<bool>& med)
{
  int Jval = 0;
  for (int i = 0; i < par.n_mediums; ++i)
  {
    
    Jval += med[i]*par.med_table[i]; // medp_bool[i]*4;
  }
  Jval += 6; //  += 6 offset so interaction with medium is not 0     
  return Jval;
}


int LKScore(vector<bool>& l1, vector<bool>& k1, vector<bool>& l2, vector<bool>& k2)
{
  int score{};
  for (int i =0; i < par.n_locks; ++i)
  {
    score += ( k1[i] != l1[i] )?1:0; // (( keys_bool[i] == lock2[i] )?1:0) * par.med_table[i];
    score += ( k2[i] != l2[i] )?1:0; // (( key2[i] == locks_bool[i] )?1:0) * par.med_table[i];
  }

  // perfect score is 10 (all locks and keys match). 
  int J = 4 + (int)( 20. - 20 * ((double)score) / (par.n_lockandkey)); //4 10 10     20-16 

  return J; 
}


void CellularPotts::BindingBetweenCells()
{
  vector<int> phenotypes{};
  map<int, vector<bool>> keys{};
  map<int, vector<bool>> locks{};
  map<int, vector<bool>> meds{};

  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->Phenotype();
      int p = c->GetPhenotype();
      if (find(phenotypes.begin(), phenotypes.end(), p) == phenotypes.end())
      {
        phenotypes.push_back(p);
        vector<bool> lock = c->get_locks_bool();
        vector<bool> key = c->get_keys_bool();
        vector<bool> med = c->get_medp_bool();
        keys[p] = key;
        locks[p] = lock;
        meds[p] = med;
      }
    }
  }

  ofstream outfile;
  string out = data_file + "/cell-bindings.dat";
  outfile.open(out, ios::app);

  int sizeL = phenotypes.size();
  for (int i = 0; i < sizeL; ++i)
  {
    for (int j = i; j < sizeL; ++j)
    {
      //compare the two with binding function
      int pi = phenotypes[i];
      int pj = phenotypes[j];
      int score = LKScore(locks[pi], keys[pi], locks[pj], keys[pj]);

      // higher score = higher binding = lower J
      outfile << pi << '\t' << pj << '\t' << score << endl;
    }
  }

  // should output mediums as well
  outfile << endl;
  for (int i = 0; i < sizeL; ++i)
  {
    int p = phenotypes[i];
    int binding = MedBinding(meds[p]);
    outfile << p << '\t' << binding << endl;
  }



}





double CellularPotts::AverageBinding()
{
  // int cell_count{};
  // int bad_cells{};
  double avg_score{};
  int count=0;

  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++)
  {
    if (c->AliveP())
    {
      vector<Cell>::iterator c2=c;
      for (;c2!=cell->end();c2++)
      {
        if (c->AliveP())
        {
          vector<bool>& locks = c2->get_locks_bool();
          vector<bool>& keys = c2->get_keys_bool();

          int score = c->LocksKeysScore(locks, keys);
          avg_score += score;
          ++count;
        }
      }
    }
  }
  avg_score = avg_score / count;
  // cout << "average score is: " << avg_score << endl;
  return avg_score;
}


double CellularPotts::AvgMedsOn()
{
  double avg{};
  int count=0;
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++)
  {
    if (c->AliveP())
    { 
      avg += c->CheckMedsOn();
      ++count;
    }
  }
  return avg / count;
}



// This is going to be part of the mother of all alogithms. 
double CellularPotts::CompareGrid(int **grid2)
{
  int overlap{};
  int outer{};
  double proportion{};

  for (int x = 0; x<sizex*sizey;++x)
  {
    if (!sigma[0][x] && !grid2[0][x])
    {
      continue;
    }
    else if (sigma[0][x] > 0 && grid2[0][x] > 0)
    {
      ++overlap;
    }
    else if (sigma[0][x] != -1)
    {
      ++outer;
    }
  }
  cout << "OVERLAP: " << overlap << "   OUTER: " << outer << endl;
  proportion = (double)overlap / (double)(overlap + outer);
  return proportion;
}

int** CellularPotts::ReturnGrid()
{
  return sigma;
}















































// ALL OF IS DEPRACATED


double CellularPotts::AngleCurvature()
{
  // find perimeter cells
  // find center of mass for perimeter of cells
  // if we assume cells are equally spaced around the border
  // calculate change in angle from one center of mass to the next around the circle.

  vector<int> pcells = LinkPerimeter();
  if (pcells.size() < 5)
  {
    cout << "Angle curvature incorrectly called." << endl;
    return 0.;
  }

  // Check which cells are too exposed (we exclude these from angle calculations. )
  CellExposure();


  // MUST CALL CELL MASS TO ENSURE CELLS HAVE CORRECT CENTER OF MASS VARIABLES ( this is called in link perimeter or derivatives )

  // storage
  vector<double> angles;
  vector<double> arc_lengths;

  // how far around the unit circle we hav gone
  double sumtheta = 0;

  //total perimeter (we don't want a hyper-accurate perimeter because prone to fluctuations.)
  double perimeter = 0;

  for (unsigned int i=0; i<pcells.size()-2;i++)
  {

    // Only calculate angle if both cells are approporiately positioned.
    if (!(*cell)[pcells.at(i+2)].exposed && !(*cell)[pcells.at(i)].exposed)
    {

      double dx = (*cell)[pcells.at(i+2)].xcen - (*cell)[pcells.at(i)].xcen;
      double dy = (*cell)[pcells.at(i+2)].ycen - (*cell)[pcells.at(i)].ycen;

      // add arc lengths and to perimeter
      double length = sqrt(pow(dx, 2) + pow(dy, 2));
      arc_lengths.push_back(length);
      perimeter += length;

      // formula is: angle = how far around the unit circle - how far went round in previous theta
      double angle{};

      if (dx > 0 && dy > 0)
      {
        angle = abs(atan(dx/dy));
      }
      else if (dx > 0)
      {
        angle = M_PI_2 + abs(atan(dy/dx));
      }
      else if (dy < 0)
      {
        angle = M_PI + abs(atan(dx/dy));
      }
      else
      {
        angle = M_PI + M_PI_2 + abs(atan(dy/dx));      
      }
      // cout << "UNIT CIRCLE: " << angle << endl;

      double chg = angle - sumtheta;
      if (chg > M_PI)
      {
        chg = -2. * M_PI + angle - sumtheta;
      }
      else if (chg < -M_PI)
      {
        chg = 2. * M_PI + angle - sumtheta;
      }
      angles.push_back(chg);
      sumtheta = angle;
    }
    else
    {
      angles.push_back(100.);
    }

  }

  // unfortunately need to remove first and redo last to go around the unit circle correctly. 
  angles.erase(angles.begin());

  if (!(*cell)[pcells.at(2)].exposed && !(*cell)[pcells.at(0)].exposed)
  {
    double dx = (*cell)[pcells.at(2)].xcen - (*cell)[pcells.at(0)].xcen;
    double dy = (*cell)[pcells.at(2)].ycen - (*cell)[pcells.at(0)].ycen;
    double length = sqrt(pow(dx, 2) + pow(dy, 2));
    double angle{};

    if (dx > 0 && dy > 0)
    {
      angle = abs(atan(dx/dy));
    }
    else if (dx > 0)
    {
      angle = M_PI_2 + abs(atan(dy/dx));
    }
    else if (dy < 0)
    {
      angle = M_PI + abs(atan(dx/dy));
    }
    else
    {
      angle = M_PI + M_PI_2 + abs(atan(dy/dx));      
    }
    double chg = angle - sumtheta;
    if (chg > M_PI)
    {
      chg = -2. * M_PI + angle - sumtheta;
    }
    else if (chg < -M_PI)
    {
      chg = 2. * M_PI + angle - sumtheta;
    }
    angles.insert(angles.begin(), chg);
    sumtheta = angle;
  }
  else
  {
    angles.insert(angles.begin(), 100.);
  }




  /// Now we compute first and second derivative with respect to arc lengths. Compute first and then second 
  vector<double> firstd{};

  double curve=0;

  for (unsigned int i=0; i < arc_lengths.size()-1;++i)
  {
    if (angles[i+1] > 90 || angles[i] > 90)
    {
      // exposed cell detected
      // cout << "skipped a cell.  " << angles[i+1] << "  " << angles[i] << endl;
      firstd.push_back(100.);
      continue;
    }
    else
    {
      double dydx = (angles[i+1] - angles[i]) / (arc_lengths[i] + arc_lengths[i+1]);
      firstd.push_back(dydx);
      // cout << "Angle: " << angles[i+1] << "  " << angles[i] << endl;
    }
  }
  
  int skip{};
  int noskip{};

  vector<double> secondd;
  for (unsigned int i=0; i < firstd.size()-1;++i)
  {
    if (firstd[i+1] > 90 || firstd[i] > 90)
    {
      // exposed cell detected
      ++skip;
      continue;
    }
    else
    {
      double dydx = (firstd[i+1] - firstd[i]) / (arc_lengths[i] + arc_lengths[i+1] + arc_lengths[i+2]);
      secondd.push_back(dydx);
      curve += abs(dydx);
      ++noskip;
    }
  }
  curve = curve / firstd.size();
  double scalar = perimeter * curve;
  // cout << "Shape value: " << scalar << endl;

  // cout << "Skipped: " << skip << "   Not skipped: " << noskip << endl;


  // curve = curve / static_cast<double>(secondd.size());

  // 


  // for (double i : angles)
  // {
  //   cout << "Angle: " << i << endl;
  // }

  // cout << endl << "NOW DERIVATIVES!" << endl;
  // for (double i : secondd)
  // {
  //   cout << "second d: " << i << endl;
  // }
  
  return scalar;

}



vector<vector<bool>> CellularPotts::ReturnGridBad()
{
  vector<vector<bool>> grid;
  grid.resize(sizex); 
  for (vector<bool>& i : grid)
  {
    i.resize(sizey);
  }

  for (int x=0;x<sizex;++x)
    for (int y=0;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
        grid[x][y] = true;
      else
        grid[x][y] = false;
    }
  return grid;
}



double CellularPotts::VecDoubleDeriv(vector<double> &vex)
{
  // calculate first derivative
  vector<double> firstder{};

  for (unsigned int i=0;i<vex.size()-1;++i)
  {
    firstder.push_back(vex[i+1]-vex[i]);
  }

  double sum=0;
  // calculate second derivative and sum absolute values
  for (unsigned int i=0;i<firstder.size()-1;++i)
  {
    sum += abs(firstder[i+1]-firstder[i]);
    cout << "double dir: " << abs(firstder[i+1]-firstder[i]) << endl;
  }

  return sum;
}



// Shape stuff is not working AND IS LEAKING!! 
// void CellularPotts::InitShape(int n)
// {
//   Shape = new int**[n];
//   for (int l=0;l<n;l++)
//   {
//     Shape[l] = new int*[sizex];
//     for (int i=0;i<sizex;i++)
//     {
//       Shape[l][i] = new int[sizey];
//     }
//   }
//   for (int l=0;l<n;l++)
//     for (int x=0;x<sizex;x++)
//       for (int y=0;y<sizey;y++)
//         Shape[l][x][y] = 0;




//   // Shape = (int ***)malloc(n*sizeof(int **));

//   // if (Shape == NULL)
//   //   cerr << "MEMORY FUCKERY\n";
  
  
//   // Shape[0]=(int **)malloc(n*sizex*sizeof(int *));
//   // if (Shape[0]==NULL)  
//   //     cerr << "MEMORY FUCKERY\n";
  
//   // for (int i=1;i<n;i++) 
//   //   Shape[i]=Shape[i-1]+sizex;
  
//   // Shape[0][0]=(int *)malloc(n*sizex*sizey*sizeof(int));
//   // if (Shape[0][0]==NULL)  
//   //   cerr << "MEMORY FUCKERY\n";

//   // for (int i=1;i<n*sizex;i++) 
//   //   Shape[0][i]=Shape[0][i-1]+sizey;

  
//   // for (int i=0;i<n*sizex*sizey;i++) 
//   //   Shape[0][0][i]=0;
// }


// void CellularPotts::AddNewShape()
// {
//   bool con = CheckAllConnected();
//   if (con && ShapeMaintained)
//   {
//     for (int x=0;x<sizex;++x)
//       for (int y=0;y<sizey;++y)
//       {
//         if (sigma[x][y] > 0)
//           Shape[scount][x][y] = 1;
//         else
//           Shape[scount][x][y] = 0;
//       }
//   }
//   else 
//   {
//     ShapeMaintained = false;
//   }
//   ++scount;
// }


// double CellularPotts::ChangeInShape()
// {
//   double count=0;

//   bool filled[sizex][sizey];

//   for (int n=1;n<scount;++n)
//   {
//     for (int x=0;x<sizex;++x)
//       for (int y=0;y<sizey;++y)
//       {
//         if (Shape[n][x][y] != Shape[n-1][x][y])
//         {
//           if (filled[x][y])
//             count += 0.2;
//           else
//           {
//             count += 1;
//             filled[x][y] = true;
//           } 
//         }
//       }
//    }
//   return count; 
// }




double getcurve(vector<int> &vec)
{
  int samples = vec.size()-2;
  int* fd = new int[samples];
  for (int i=0;i<samples;++i)
  {
    fd[i] = (vec[i+2] - vec[i+1]) - (vec[i+1] - vec[i]);
  }

  // get absolute value of curvature
  double var=0;
  for (int i=0;i<samples;++i)
  {
    var += abs(fd[i]);
  }
  var = var / (double)samples;

  delete[] fd;

  return var;

}



double CellularPotts::Curvature()
{
  vector<int> miny;
  vector<int> maxy;
  vector<int> maxx;
  vector<int> minx;

  for (int x=0;x<sizex;)
  {
    int max=0;
    int min=1000;
    bool if_sig = false;
    
    for (int y=0;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
      {
        if_sig = true;
        if (y > max)
          max = y;

        if (y < min)
          min = y;
      }        
    }
    if (if_sig)
    {
      miny.push_back(min);
      maxy.push_back(max); 
    }
    x += 5; // needs to be greater than 1 so that we can get a good look at global curvature and not random fluctuations. 
  }
  // iterate opposite direction
  for (int y=0;y<sizey;)
  {
    int max=0;
    int min=1000;
    bool if_sig = false;
    
    for (int x=0;x<sizex;++x)
    {
      if (sigma[x][y] > 0)
      {
        if_sig = true;
        if (x > max)
          max = x;

        if (x < min)
          min = x;
      }        
    }
    if (if_sig)
    {
      minx.push_back(min);
      maxx.push_back(max); 
    }
    y += 5;
  }


  double curvature = 0;
  curvature += getcurve(minx);
  curvature += getcurve(miny);
  curvature += getcurve(maxx);
  curvature += getcurve(maxy);

  return curvature;


}
















// void CellularPotts::morphogenWave() // FUNCTION IS SCREWED UP BECAUSE OF CENTER ATM
// {
// // get organism center, assumming it is an approximate circle at time of calling morphogen wave:

//   double center[] = {0., 0., 0.,};
//   get_center(center);

//   for (int i=0;i<sizex-1;++i)
//     for (int j=0;j<sizey-1;++j)
//     {
//       //find vector to value (function is symmetric so only need absolute value)

//       if (sigma[i][j] > 0)
//       {
//         int lx = i - center[0];
//         int ly = j - center[1];
//         double vec_length = sqrt(pow(lx, 2) + pow(ly,2));
//         // cout << vec_length << endl;

//         // We are doing a gaussian distribution for the morphogen gradient (e^(-(x)^2)).
//         double val = pow((-vec_length / center[2]), 2);
//         double amount = exp(-val);
//         // cout << amount << endl;

//         (*cell)[sigma[i][j]].add_morphogen(amount);
//       } 
//     }
  
//   vector<Cell>::iterator c;
//   for ( (c=cell->begin(), c++);c!=cell->end();c++) 
//     if (c->AliveP())
//     {
//       c->calc_morphogen();
//     }
// }

// void CellularPotts::decay_morph(double& morph)
// {
//   morph = morph * par.morphdecay;
// }



int CellularPotts::CountSomaticCells(void) 
{
  int amount=0;
  vector<Cell>::iterator c;
  for ( (c=cell->begin(),c++); c!=cell->end(); c++) 
  {
    if (c->AliveP() && c->checkforcycles(1) == false) 
    {
      vector<double> &genes = c->get_genes();
      int count=0;
      int i = par.n_diffusers + par.n_MF;
      int k = i + par.n_stem;
      while (i < k)
      {
        if (genes[i] > 0.5)
          ++count;
        ++i;
      }
      if (count == par.n_stem)
        amount++;
    }
    // Old way of counting somatic cells. 
    // {
    //   vector<bool> &set = c->get_set();
    //   int count=0;
    //   for (int i=0;i<par.n_lockandkey + par.n_length_genes;i++)
    //   {
    //     if (set[i])
    //       ++count;
    //   }
    //   if (count > par.max_on)
    //     amount++;

    // }
  }
  return amount;
}

void CellularPotts::som_fitness()
{
  som_cell_list.push_back(CountSomaticCells());
}
