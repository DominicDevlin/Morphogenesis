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
#include <set>
#include <bits/stdc++.h>

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
  CopyProb(par.eT);
  cell=cells;
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else 
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4]).";
  
}

CellularPotts::CellularPotts(vector<Cell> *cells,
			     const int sx, const int sy ) {
  
  sigma=0;
  frozen=false;
  thetime=0;
  zygote_area=0;

  
  BaseInitialisation(cells);
  sizex=sx;
  sizey=sy;

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

  if (Shape) {
    free(Shape[0][0]);
    free(Shape[0]);
    free(Shape);
    Shape=0;
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
  if (PDEfield && (par.vecadherinknockout || (sxyp==0 || sxy==0))) {
    
    // copying from (xp, yp) into (x,y)
    // If par.extensiononly == true, apply CompuCell's method, i.e.
    // only chemotactic extensions contribute to energy change
    if (!( par.extensiononly && sxyp==0)) {
      int DDH=(int)(par.chemotaxis*(sat(PDEfield->Sigma(0,x,y))-sat(PDEfield->Sigma(0,xp,yp))));
    
      DH-=DDH;
    }
  }
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

#include <fstream>
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
    if (par.periodic_boundaries) {
      
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
      
      
    } else {
      
      if (xp<=0 || yp<=0 
	  || xp>=sizex-1 || yp>=sizey-1)
	kp=-1;
      else
	kp=sigma[xp][yp];

    }

    
    
    // test for border state (relevant only if we do not use 
    // periodic boundaries)
    if (kp!=-1) {  
      // Don't even think of copying the special border state into you!
    
     
      if ( k  != kp ) {

	/* Try to copy if sites do not belong to the same cell */
	
	// connectivity dissipation:
	int H_diss=0;
	if (!ConnectivityPreservedP(x,y)) 
    H_diss=par.conn_diss;
	
	int D_H=DeltaH(x,y,xp,yp, tsteps, PDEfield);
	
	if ((p=CopyvProb(D_H,H_diss))>0) {
	  ConvertSpin ( x,y,xp,yp );
	  SumDH+=D_H;
	}
	
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

void CellularPotts::DivideCells(vector<bool> which_cells)
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

  if (celldir) 
    delete[] (celldir);
  
  if (divflags)
    free(divflags);
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

  static int stack[8]; // stack to count number of different surrounding cells
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
    if (s_nb) {
      for (j=stackp;j>=0;j--) {
	if (s_nb==stack[j]) {
	  on_stack_p=true;
	  break;
	}
      }
      if (!on_stack_p) {
	if (stackp>6) {
	  cerr << "Stack overflow, stackp=" << stackp << "\n";
	}
	stack[++stackp]=s_nb;
      }
    } else {
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

void CellularPotts::CellGrowthAndDivision(void) 
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
        int gthresh;
        if (c->get_stem())
          gthresh = par.stem_gthresh;
        else
          gthresh = par.gthresh;


        if ( (area-TA)>gthresh) // && area <= (double)(par.div_threshold) * 1.1) //  
        {
          int count= area-TA; //area-TA;
          while (count>0)
          {
            c->IncrementTargetArea();
            --count;
          }
        }
        else if ( (area-TA)<par.shrink ) 
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
      DivideCells(which_cells);
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
    DivideCells(to_divide);

    vector<vector<int>> middles;
    vector<Cell>::const_iterator i; 
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
    DivideCells(to_divide);

    // determine which cells have first maternal factor on, and which have first maternal factor off
    vector<int> g4_on;
    vector<int> g4_off;
    
    vector<Cell>::iterator i; 
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
      if (val < 0.02)
      {
        matrix[i][j] = -2;
      }
      else if (val < 0.2)
      {
        matrix[i][j] = -1;
      }
      else if (val < 0.7)
      {
        matrix[i][j] = 0;
      }
      else if (val < 0.98)
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
  for (int i=0;i<par.n_diffusers+par.n_MF+par.n_TF+par.n_length_genes;++i)
  {
    if (i >= par.n_diffusers+par.n_MF && i < par.n_diffusers+par.n_MF+par.n_TF) // this is wrong. Should be i=4,5 so i>=4, i<6
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


double CellularPotts::discrete_decay(vector<double>& gene_list, double conc, int gene_n, double morph)
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
    if (par.morphogen && j == par.n_activators-1)
      x_1 += (matrix[gene_n][j] * morph);
    else
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

  int i = 0;
  for ( (c=cell->begin(), c++); c!=cell->end(); c++) 
  {
    if (c->AliveP())
    {
      vector<double>& genes = c->get_genes();

      vector<double>& diffusers = c->get_diffusers(); 

      vector<double>& locks = c->get_locks();
      vector<double>& keys = c->get_keys();

      vector<double> gene_copy = c->get_genes();

      
      double& morph = c->get_morphogen();
      // cout << morph << endl;

      if (genes.empty() == true)
        cout << "EMPTY" << endl;

      // iterate through genes and update them according to GRN
      // iteration is a bit weird through different vectors so dont screw it up. 
      int j=0;
      int k=0;
      for (int i = 0; i < par.n_genes; ++i)
      {
        if (i < par.n_diffusers)
          diffusers.at(i) = discrete_decay(gene_copy, diffusers.at(i), i, morph);
        else if (i < par.n_genes - par.n_lockandkey)
          genes.at(i) = discrete_decay(gene_copy, genes.at(i), i, morph);
        else if (i < par.n_genes - par.n_locks)
        {
          locks.at(j) = discrete_decay(gene_copy, locks.at(j), i, morph);
          ++j;
        }
        else
        {
          keys.at(k) = discrete_decay(gene_copy, keys.at(k), i, morph);
          ++k;
        }
      }
      if (par.morphogen)
        decay_morph(morph);

      // make boolean set. 
      vector<bool>& full_set = c->get_set();


      //create bool based on lock and keys.
      vector<bool>& l_bool = c->get_locks_bool();
      vector<bool>& k_bool = c->get_keys_bool();

      for (int i=0; i < par.n_locks; ++i)
      {

        l_bool[i] = (locks[i]>0.5) ? true : false;
        full_set[i] = l_bool[i];

        k_bool[i] = (keys[i]>0.5) ? true : false;
        full_set[i+par.n_locks] = k_bool[i];
      }

      // hard coding the two target genes in for now. 
      full_set.at(10) = ((genes.at(par.tloc1)>0.5) ? true : false);
      full_set.at(11) = ((genes.at(par.tloc2)>0.5) ? true : false);

      // set the type of the cell based on network arrangement.
      if (tsteps > par.end_program)
        c->set_ctype(set_type(full_set) * c->getTau());

      c->add_to_cycle();


      // get neighbours and set internal gene value based on surrounding cells 
      // THIS IS GOING TO BE DEPRECATED 
      // for (int i = 0; i < par.n_diffusers; ++i)
      // {
      //   get_neighbours(i, genes);
      // }
      
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

      // change the cell fitness (growth rate) based on gene network state of TF1 and TF2 (stem cell signature). 
      if (genes.at(par.n_diffusers + par.n_MF) > 0.5 && genes.at(par.n_diffusers + par.n_MF + 1) > 0.5)
      {
        // if (!c->get_stem())
        //   cout << "Cell number: " << c->Sigma() << " DIFFERENTIATED SOMATIC -> STEM! " << endl;
        c->set_stem(true);
        // c->set_ctype(160);
        // cout << "cell number: " << c->Sigma() << "  is a stem cell!" << endl;
      }
      else
      {
        // if (c->get_stem())
        //   cout << "Cell number: " << c->Sigma() << " DIFFERENTIATED STEM -> SOMATIC! " << endl;
        c->set_stem(false);
      }

    }

    ++i;
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

// DEPRECATED FUNCTION!! USING PDE FIELD INSTEAD FOR DIFFUSION !
// for now I will put the connectivity proteins in here. Remember that 0 refers to medium. Called after decay. 
void CellularPotts::get_neighbours(int cell_n, std::vector<double>& genes) 
{
  std::set<int> neighbour_list;

  for ( int i = 0; i < sizex-1; i++ )
  {
    for ( int j = 0; j < sizey-1; j++ ) 
    {
      if (sigma[i][j] == cell_n)
      {
        int neighbour_1 = sigma[i][j-1];
        int neighbour_2 = sigma[i][j+1];
        if (neighbour_1 != cell_n)
        {
          if (neighbour_list.count(neighbour_1) == 0)
          {
            // add neighbour to the list do neighbour stuff
            neighbour_list.insert(neighbour_1);
            // cout << "neighbour site is: " << i << "-" << j+1 << ", with cell number:" << neighbour_1 << endl;
          }
        }

        if (neighbour_2 != cell_n)
        {
          if (neighbour_list.count(neighbour_2) == 0)
          {
            // do neighbour stuff
            neighbour_list.insert(neighbour_2);
            
          }
        }
      } 
    }
  }
  // update neighbour genes based on neighbours. 
  for (int j=0;j<par.n_diffusers;++j)
  {
    double aved=0;
    for (const int &i : neighbour_list)
    {
      if (i > 0)
      {
        vector<double>& diffusers = cell->at(i).get_diffusers();
        aved += diffusers.at(j);
      }
    }
    genes.at(j) = aved / neighbour_list.size();
  }  
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
  for (int i=0; i < par.n_lockandkey + par.n_length_genes; ++i)
  {
    dist += (str1[i] != str2[i]);
  }
  return dist;
}



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


// return the fitness: combination of n cell types + hamming distance between them + number of cells. 
void CellularPotts::type_fitness()
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
    
  if (n_types > 1)
  {
    // find average distance between all types. will use log this for now. 
    double average_hamming{};
    // find smallest hamming distance between types. Not using currently.   
    // int min_hamming{100};

    int i=0;
    int j=1;
    for (; i < n_types -1; ++i)
    {
      for (unsigned int k=j; k < n_types; ++k)
      {
        int h = hamming_distance(cell_types[i], cell_types[k]);
        average_hamming += h;
        // min_hamming = ((h < min_hamming) ? h : min_hamming);
      }
      ++j; 
    }
    int n_comparisons = (n_types * (n_types-1)) / 2;
    
    average_hamming = average_hamming / n_comparisons;

    // cout << "average hamming distance is: " << average_hamming << ". Number of types is: " << n_types << 
    // ". Number of somatic cells cells is: " << n_som << endl;
    type_fitness_list.push_back(static_cast<double>(n_types) + (log(average_hamming)));
  }
  else
    type_fitness_list.push_back(1);

  if (par.print_fitness)
    cout << "Latest type diversity is: " << type_fitness_list.back() << endl;
}




// called once at end of sim.
double CellularPotts::get_fitness()
{
  // Cell type diversity and expression distance
  double min_diversity=1000;
  
  for (unsigned int j=0; j < type_fitness_list.size(); ++j)
    if (type_fitness_list.at(j) < min_diversity)
    {
      min_diversity = type_fitness_list.at(j);
    }
  // log the number of cell types so we aren't selecting for a ridiculous number + doesn't take over fitness
  min_diversity = log2(min_diversity) / log(1.5);


  // Calculate the change in shape. 
  double Dshape =  0;
  double curvature = 0;
  
  if (ShapeMaintained) 
  {
    Dshape = ChangeInShape();
    curvature = Curvature();
  }

  if (par.print_fitness)
    cout << "Type diversity: " << min_diversity << "  Change in Shape: " << Dshape << "  Curvature: " << curvature << endl;

  return min_diversity + ((Dshape * curvature)  / 1000.);



  // double somatic_rate{};

  // double avg_shape{};

  // for (unsigned int i=0;i<som_cell_list.size()-1;++i)
  // {
  //   // What if I go to negative fitness because losing somatic cells? Is that okay?
  //   // No its not... Need to come up with a way to ensure that its just "production"
  //   //  Also need to log this to get straight line instead of exp. 
    
  //   somatic_rate += (double)(som_cell_list[i+1]-som_cell_list[i]+1);
  //   if (par.print_fitness)
  //     cout << "Som cells in sector: " << i << " is " << som_cell_list[i] << endl;
  // }

  // if (som_cell_list.size())
  //   somatic_rate = somatic_rate / ((double)(som_cell_list.size()-1)*3.);
  // else
  //   cerr << "fitness error" << '\n';




  // double n_new_cells=0;
  // if (n_new_cells > 64)
  //   n_new_cells = log(CountCells() - 64) + 1;
  // else 
  //   n_new_cells = 0;


  // for (unsigned int j=0; j < shape_fitness_list.size(); ++j) 
  // {
  //   avg_shape += shape_fitness_list.at(j);
  // }
  // avg_shape = avg_shape / (double)shape_fitness_list.size();



  // final fitness value to return, have turned off germ-soma thing or now. 
  // return min_diversity + somatic_rate + n_new_cells; 

  
  
  // if (CheckAllConnected())
  //   return CalculateShape() + min_diversity / 4;
  // else 
  //   return 0;
}





// Multithreading method
void CellularPotts::set_num(int in)
{
  org_num = in;
}
// Multithreading method
void CellularPotts::set_seed()
{
  s_val[0] = Seed(org_num); 
  cout << "Seed is: " << s_val[0] << endl;
}




void CellularPotts::record_GRN()
{
  vector<Cell>::iterator c;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      c->add_to_vectors();
    }
  }
}

void CellularPotts::print_cell_GRN()
{

  vector<Cell>::iterator c;
  int i=0;
  for ((c=cell->begin(), c++); c!=cell->end(); c++)
  {
    if (c->AliveP())
    {
      string var_name = "cell_conc_" + to_string(i) + ".dat";
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
      // break;
      ++i;
    }
  }

}





// check if there are any lonely cells. Currently deprecated. 
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


//IMPORTANT METHOD:  Function to ensure all cells are connected indirectly to all cells on lattice.  
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
      for (int x=1;x<sizex-2;++x)
        for (int y=1;y<sizey-2;++y)
        {
          if (sigma[x][y] == id)
          {
            if (sigma[x][y-1] > 0 && sigma[x][y-1] != id)
              tempcon.emplace(sigma[x][y-1]);

            if (sigma[x][y+1] > 0 && sigma[x][y+1] != id)
              tempcon.emplace(sigma[x][y+1]);
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

void CellularPotts::InitShape(int n)
{
  int ***mem;
  mem = (int ***)malloc(n*sizeof(int **));

  if (mem == NULL)
    cerr << "MEMORY FUCKERY\n";
  
  
  mem[0]=(int **)malloc(n*sizex*sizeof(int *));
  if (mem[0]==NULL)  
      cerr << "MEMORY FUCKERY\n";
  
  for (int i=1;i<n;i++) 
    mem[i]=mem[i-1]+sizex;
  
  mem[0][0]=(int *)malloc(n*sizex*sizey*sizeof(int));
  if (mem[0][0]==NULL)  
    cerr << "MEMORY FUCKERY\n";

  for (int i=1;i<n*sizex;i++) 
    mem[0][i]=mem[0][i-1]+sizey;

  
  /* Clear PDE plane */
  for (int i=0;i<n*sizex*sizey;i++) 
    mem[0][0][i]=0;

  Shape = mem;

}


void CellularPotts::AddNewShape()
{
  bool con = CheckAllConnected();
  if (con && ShapeMaintained)
  {
    for (int x=0;x<sizex;++x)
      for (int y=0;y<sizey;++y)
      {
        if (sigma[x][y] > 0)
          Shape[scount][x][y] = 1;
        else
          Shape[scount][x][y] = 0;
      }
  }
  else 
  {
    ShapeMaintained = false;
  }
  ++scount;
}


double CellularPotts::ChangeInShape()
{
  double count=0;

  bool filled[sizex][sizey];

  for (int n=1;n<scount;++n)
  {
    for (int x=0;x<sizex;++x)
      for (int y=0;y<sizey;++y)
      {
        if (Shape[n][x][y] != Shape[n-1][x][y])
        {
          if (filled[x][y])
            count += 0.2;
          else
          {
            count += 1;
            filled[x][y] = true;
          } 
        }
      }
   }
  return count; 
}












int* CellularPotts::get_center()
{
// get organism center, assumming it is an approximate circle at time of calling morphogen wave:
  int minx = 1000;
  int miny = 1000;
  int maxx = 0;
  int maxy = 0;

  for (int i=0;i<sizex-1;++i)
    for (int j=0;j<sizey-1;++j)
    {
      if (sigma[i][j] > 0)
      {
        if (i < minx)
          minx = i;
        else if (i > maxx)
          maxx = i;

        if (j < miny)
          miny = j;
        else if (j > maxy)
          maxy = j;
      }
    }
  // calculate approximate center:
  int* center = new int[3];
  center[0] = round((maxx + minx) / 2);
  center[1] = round((maxy + miny) / 2);

  // get average distance to get stdev for gradient (/2 = radius, /4 = both x and y, / 12 to get 1/8 of distance)
  center[2] = ((maxx - minx) + (maxy - miny)) / 8;

  return center;
}


// Currently Not using
void CellularPotts::ShapeFitness()
{
  if (CheckAllConnected())
  {
    shape_fitness_list.push_back(CalculateShape());
  }
  else
  {
    shape_fitness_list.push_back(0);
  }
  if (par.print_fitness)
    cout << "Latest shape fitness is: " << shape_fitness_list.back() << endl;
}






double getcurve(vector<int> &vec)
{
  int samples = vec.size()-2;
  int* fd = new int[samples];
  for (unsigned int i=0;i<samples;++i)
  {
    fd[i] = (vec[i+2] - vec[i+1]) - (vec[i+1] - vec[i]);
    cout << "Curvature: " << fd[i] << endl;
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

  cout << "Curvature is: " << curvature << endl;

  return curvature;


}




// Currently Not using
double CellularPotts::CalculateShape()
{

  //calculate area
  int count=0;
  for (int x=0;x<sizex;++x)
    for (int y=0;y<sizey;++y)
    {
      if (sigma[x][y] > 0)
        ++count;
    }

  //calculate radius if perfect circle
  double rad = sqrt((double)(count) / M_PI);

  // calculate approximate centroid based only on x and y
  int *center = get_center();

  // deviation from circle, calculated by summing radius differnce
  double dev=0;
  int whitespacex{};
  int whitespacey{};

  //vector of math vectors, need one for each clock hand
  // vector<double> up_vectors{};
  // vector<double> down_vectors{};
  

  for (int x=0;x<sizex;++x)
  {
    int count=0;
    int maxy=0;
    int miny=1000;
    
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
    // going to add in some code that calculates the amount of medium in between the max and min y
    for (int y=miny+1;y<maxy;++y)
    {
      if (!sigma[x][y])
        ++whitespacex;
    }

    if (if_sig)
    {
      // calculate vector based on max y for every x
      double vec1 = sqrt(pow((double)(x-center[0]), 2) + pow((double)(maxy-center[1]), 2));
      double vec2 = sqrt(pow((double)(x-center[0]), 2) + pow((double)(miny-center[1]), 2));

      
      // Function doesn't work that well. 
      // up_vectors.push_back(vec1);
      // down_vectors.push_back(vec2);

      // Want to select for size, so divide by sqrqt of radius instead of radius ( dividing by radius normalises by size)
      vec1 = abs(vec1 - rad) / sqrt(rad);
      vec2 = abs(vec2 - rad) / sqrt(rad);
      dev += vec1;
      dev += vec2; 
      ++count;
    }
  }
  // and do the same for x across y
  for (int y=0;y<sizey;++y)
  {
    int count=0;
    int maxx=0;
    int minx=1000;
    
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
    // going to add in some code that calculates the amount of medium in between the max and min y
    for (int x=minx+1;x<maxx;++x)
    {
      if (!sigma[x][y])
        ++whitespacey;
    }

    if (if_sig)
    {
      // calculate vector based on max and min x for every y
      double vec1 = sqrt(pow((double)(y-center[0]), 2) + pow((double)(maxx-center[1]), 2));
      double vec2 = sqrt(pow((double)(y-center[0]), 2) + pow((double)(minx-center[1]), 2));

      // Want to select for size, so divide by sqrqt of radius instead of radius ( dividing by radius normalises by size)
      vec1 = pow(vec1 - rad,2) / sqrt(rad);
      vec2 = pow(vec2 - rad,2) / sqrt(rad);
      dev += vec1;
      dev += vec2; 
      ++count;
    }
  }

  dev = dev / (double)count;
  // delete pointer to array
  delete[] center;

  //get second derivative
  // double deviation = VecDoubleDeriv(up_vectors) + VecDoubleDeriv(down_vectors);

  double shape = (dev * 10.) + sqrt(whitespacex + whitespacey);

  if (par.print_fitness)
  {
    cout << "Deviation from circle: " << dev << "  White space x: " << whitespacex << "  White space y: " << whitespacey 
    << "  Multiplying: " << shape << endl;
  }


  return shape;  
}

// depracated function. 
double CellularPotts::VecDoubleDeriv(vector<double> vex)
{
  // calculate first derivative
  vector<double> firstder{};

  for (int i=0;i<vex.size()-1;++i)
  {
    firstder.push_back(vex[i+1]-vex[i]);
  }

  vector<double> secondder{};
  double mean{};
  // calculate second derivative
  for (int i=0;i<firstder.size()-1;++i)
  {
    double val = firstder[i+1]-firstder[i];
    secondder.push_back(val);
    mean += val;
  }
  mean = mean / static_cast<double>(secondder.size());
  double variance{};

  // calculate variance in second derivative. 
  for (double i : secondder)
  {
    variance += pow(i - mean, 2);
  }
  variance = variance / static_cast<double>(secondder.size());
  return variance;
}













void CellularPotts::morphogenWave() // FUNCTION IS SCREWED UP BECAUSE OF CENTER ATM
{
// get organism center, assumming it is an approximate circle at time of calling morphogen wave:

  int* center = get_center();

  for (int i=0;i<sizex-1;++i)
    for (int j=0;j<sizey-1;++j)
    {
      //find vector to value (function is symmetric so only need absolute value)

      if (sigma[i][j] > 0)
      {
        int lx = i - center[0];
        int ly = j - center[1];
        double vec_length = sqrt(pow(lx, 2) + pow(ly,2));
        // cout << vec_length << endl;

        // We are doing a gaussian distribution for the morphogen gradient (e^(-(x)^2)).
        double val = pow((-vec_length / center[2]), 2);
        double amount = exp(-val);
        // cout << amount << endl;

        (*cell)[sigma[i][j]].add_morphogen(amount);
      } 
    }
  
  delete[] center;

  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    if (c->AliveP())
    {
      c->calc_morphogen();
    }
}

void CellularPotts::decay_morph(double& morph)
{
  morph = morph * par.morphdecay;
}


void CellularPotts::colourbymorph(void)
{
  vector<Cell>::iterator c;
  for ( (c=cell->begin(), c++);c!=cell->end();c++) 
    if (c->AliveP())
    {
      int val = round(c->get_morphogen() * 30);
      c->set_ctype(val);
    }  
}
