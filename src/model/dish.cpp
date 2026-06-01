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

#define EXTERNAL_OFF

extern Parameter par;

using namespace std;

Dish::Dish()
{
  ConstructorBody();
  
  CPM=new CellularPotts(&cell, par.sizex, par.sizey);

  
  // Initial cell distribution is defined by user in INIT {} block
  // Init();
      
  if (par.init_area>0)
    for (std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) {
      c->SetTargetArea(par.init_area);
    } 
}



Dish::~Dish() {
    // for (auto& c : cell) {
    //     delete &c;
    // }
    cell.clear();
    
    delete CPM;
	
 }

void Dish::Plot(Graphics *g) {
    if (CPM)
      CPM->Plot(g);
 }


void Dish::ConstructorBody() {
  
  maxsigma=0;
  // amount=0;
  
  // Allocate the first "cell": this is the medium (tau=0)
  cell.emplace_back(*this, 0);
  // cell.push_back(*(new Cell(*this,0)));
  
  // indicate that the first cell is the medium
  cell.front().sigma=0; 
  cell.front().tau=0;
  
  CPM=0;

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




int Dish::ZygoteArea(void) const {
    return CPM->ZygoteArea();
}

int Dish::Time(void) const {
    return CPM->Time();
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

