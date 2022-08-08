/* 

Copyright 1996-2006 Roeland Merks, Paulien Hogeweg

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
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "dish.h"
#include "sticky.h"
#include "parameter.h"
#include "info.h"
#include "crash.h"
#include "pde.h"

#define EXTERNAL_OFF

extern Parameter par;

using namespace std;

Dish::Dish()
{
  ConstructorBody();
  
  CPM=new CellularPotts(&cell, par.sizex, par.sizey);
  if (par.n_diffusers)
    PDEfield=new PDE(par.n_diffusers,par.sizex, par.sizey);
  
  // Initial cell distribution is defined by user in INIT {} block
  // Init();
      
  if (par.target_area>0)
    for (std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) {
      c->SetTargetArea(par.target_area);
    } 
}



Dish::~Dish() {
    cell.clear();
    
    delete CPM;
    delete PDEfield;
	
 }

void Dish::Plot(Graphics *g) {
    if (CPM)
      CPM->Plot(g);
 }


void Dish::ConstructorBody() {
  
  maxsigma=0;
  // amount=0;
  
  // Allocate the first "cell": this is the medium (tau=0)
  cell.push_back(*(new Cell(*this,0)));
  
  // indicate that the first cell is the medium
  cell.front().sigma=0; 
  cell.front().tau=0;
  
  CPM=0;
  PDEfield=0;

}


bool Dish::CellLonelyP(const Cell &c, int **neighbours) const {

  int i;

  for (i=0;i<(int)cell.size();i++) {
    if (neighbours[c.sigma][i]==EMPTY) 
      break;
    else
      if (neighbours[c.sigma][i]>0)
	return false;
  }
  
  return true;
  
}


int Dish::CountCells(void) const {
  
  int a=0;
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++); i!=cell.end(); i++) {
    if (i->AliveP()) {
      a++;
    } else {
      cerr << "Dead cell\n";
    }
  }
  return a;
}

 

int Dish::Area(void) const {
  
  int total_area=0;
  
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++);
	i!=cell.end();
	++i) {
    
    total_area+=i->Area();
    
  }
  return total_area;
}

int Dish::TargetArea(void) const {
  
  int total_area=0;
  
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++);
	i!=cell.end();
	++i) {
    
    if (i->AliveP()) 
      total_area+=i->TargetArea();
    
  }
  return total_area;
}



void Dish::SetCellOwner(Cell &which_cell) {
  which_cell.owner=this;
}



void Dish::ClearGrads(void) {

  vector<Cell>::iterator i;
  for ( (i=cell.begin(), i++); i!=cell.end(); i++) {
    i->ClearGrad();
  }
}


int Dish::ZygoteArea(void) const {
    return CPM->ZygoteArea();
}

int Dish::Time(void) const {
    return CPM->Time();
}


void Dish::MeasureChemConcentrations(void) {
 
  // clear chemical concentrations
  for (vector<Cell>::iterator c=cell.begin();
       c!=cell.end();
       c++) {
    for (int ch=0;ch<par.n_chem;ch++) 
      c->chem[ch]=0.;
  }

  // calculate current ones
  for (int ch=0;ch<par.n_chem;ch++)
    for (int i=0;i<SizeX()*SizeY();i++) {
      
      int cn=CPM->Sigma(0,i);
      if (cn>=0) 
	cell[cn].chem[ch]+=PDEfield->Sigma(ch,0,i);
	
    }

    for (vector<Cell>::iterator c=cell.begin();
       c!=cell.end();
       c++) {
      for (int ch=0;ch<par.n_chem;ch++) 
	c->chem[ch]/=(double)c->Area();
    }

}

int Dish::SizeX(void) { return CPM->SizeX(); }
int Dish::SizeY(void) { return CPM->SizeY(); }	


void Dish::add_maxsigma(void)
{
  ++maxsigma;
}

int Dish::get_maxsigma(void)
{
  return maxsigma;
}

void Dish::set_maxsigma(int max)
{
  maxsigma = max;
}



void Dish::AverageChemCell() // d is number of diffusers (2?)
{

  const int sizex = par.sizex;
  const int sizey = par.sizey;

  for (vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) 
  {
    for (int ch=0;ch<par.n_diffusers;ch++) 
      c->diffs[ch]=0.;
  }

  for (int i=0; i<par.n_diffusers; ++i)
    for (int x=0; x<sizex; ++x)
      for (int y=0; y<sizey; ++y)
      {
        int cn = CPM->Sigma(x,y);
        if (cn > 0)
        {
          (cell)[cn].diffs[i] += PDEfield->Sigma(i,x,y);
        }
      }
  




  // printing and debugging is rife below

  // double max_conc=0;
  // double max_conc1=0;
  // double average_conc;
  // double average_conc1;

  vector<Cell>::iterator c;
  for ( (c=cell.begin(),c++); c!=cell.end(); c++) 
    if (c->AliveP())
    {
      c->average_chem();


      // if (c->chem_conc(0) > max_conc)
      //   max_conc = c->chem_conc(0);
      // if (c->chem_conc(1) > max_conc1)
      //   max_conc1 = c->chem_conc(1);

      // average_conc += c->chem_conc(0);
      // average_conc1 += c->chem_conc(1);

    }

  // int tc = CPM->CountCells();

  // cout << "Max for diff 1: " << max_conc << ". Average: " << average_conc/(double)(tc) 
  // << ". Max for diff 2: " << max_conc1 << ". Average: " << average_conc1/(double)(tc) << endl;


}