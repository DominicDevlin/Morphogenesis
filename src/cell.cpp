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
#include <list>
#include <vector>
#include <stdio.h>
#include <fstream>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include "cell.h"
#include "sticky.h"
#include "parameter.h"
#include "dish.h"
#include <unordered_map>

#define HASHCOLNUM 255

extern Parameter par;

// int **Cell::J=0;
// int Cell::amount=0;
// int Cell::capacity=0;

// int Cell::maxsigma=0;

// int Cell::maxtau=0;

//Cell::Cell(const Dish &who) : Cytoplasm(who);
// Note: g++ wants to have the body of this constructor in cell.hh
// body is defined in "ConstructorBody" below
class Dish;

using namespace std;

Cell::~Cell(void) 
{
  delete[] chem;
  delete[] diffs;
}

void Cell::CellBirth(Cell &mother_cell) {

  colour = mother_cell.colour;
  alive = mother_cell.alive;
  v[0] = mother_cell.v[0];
  v[1] = mother_cell.v[1];
  
  // Administrate ancestry
  mother_cell.daughter=this->sigma;
  mother=mother_cell.sigma;
  times_divided=++mother_cell.times_divided;
  owner=mother_cell.owner;
  
  date_of_birth=owner->Time();

  colour_of_birth=mother_cell.colour;
  colour=mother_cell.colour;
  
  alive=mother_cell.alive;
  
  tau=mother_cell.tau;
  target_length = mother_cell.target_length;

  fitness=mother_cell.fitness;
  genes=mother_cell.genes;
  diff_genes=mother_cell.diff_genes;
  lambda_2 = mother_cell.lambda_2;
  lambda = mother_cell.lambda;

  locks=mother_cell.locks;
  locks_bool=mother_cell.locks_bool;
  keys=mother_cell.keys;
  keys_bool=mother_cell.keys_bool;
  full_set=mother_cell.full_set;
  cycles=mother_cell.cycles;
  gene_recordings=mother_cell.gene_recordings;
  morph=mother_cell.morph;
  stemness=mother_cell.stemness;



  for (int i=0;i<par.n_diffusers;i++)
  {
    diffs[i]=mother_cell.diffs[i];
  }


  c_type=mother_cell.c_type;

  
  for (int ch=0;ch<par.n_chem;ch++)
    chem[ch]=mother_cell.chem[ch];
  
  n_copies=0;

  grad[0]=mother_cell.grad[0];
  grad[1]=mother_cell.grad[1];
  
}


void Cell::ConstructorBody(int settau) {
  
  // Note: Constructor of Cytoplasm will be called first
  alive=true;
  colour=1; // undifferentiated

  c_type=1;
  
  colour_of_birth=1;
  date_of_birth=0;
  times_divided=0;
  mother=0;
  daughter=0;
    
  // add new elements to each of the dimensions of "J"
  // J deprecated here

  sigma = owner->get_maxsigma();
  owner->add_maxsigma();
  
  // if (!J) {
  //   ReadStaticJTable(par.Jtable);
  // }
  
  tau=settau;
  area=0;
  target_area=0;
  length=0;
  target_length=par.target_length;
  sum_x=0;
  sum_y=0;
  sum_xx=0;
  sum_yy=0;
  sum_xy=0;

  lambda = par.lambda;

  diffs = new double[par.n_diffusers];

  v[0]=0.; v[1]=0.;
  n_copies=0;

  chem = new double[par.n_chem];

}


/*! \brief Read a table of static Js.
 First line: number of types (including medium)
 Next lines: diagonal matrix, starting with 1 element (0 0)
 ending with n elements */
// void Cell::ReadStaticJTable(const char *fname) {

//   cerr << "Reading J's...\n";
//   ifstream jtab(fname);
//     if (!jtab)  {
//         perror(fname);
//         exit(1);
//     }
  
//   int n; // number of taus
//   jtab >> n;
//   cerr << "Number of celltypes:" <<  n << endl; 
//   maxtau=n-1;
  
//   // Allocate
//   if (J) { free(J[0]); free(J); }
//   J=(int **)malloc(n*sizeof(int *));
//   J[0]=(int *)malloc(n*n*sizeof(int));
//   for (int i=1;i<n;i++) {
//     J[i]=J[i-1]+n;
//   }
  
//   capacity = n;
//   {for (int i=0;i<n;i++) {
//     for (int j=0;j<=i;j++) {
//       jtab >> J[i][j];
//       // symmetric...
//       J[j][i]=J[i][j];
//     }
  
//   }}
// }


// return energies for programmed stage of development
int Cell::EnDif(Cell &cell2)
{ 
  
  if (sigma==cell2.sigma) 
    return 0;

  if (tau == 0 && cell2.tau == 0)
    return 0;
  else if (tau == 1 && cell2.tau == 1)
    return 30;
  else
    return 40;
  
}


// return energies by calculating lock & key products switched on by cells. 
int Cell::EnergyDifference(Cell &cell2)
{ 
  if (sigma==cell2.sigma) 
    return 0;
  else if (sigma==0)
    return CalculateJfromMed(cell2.get_keys_bool());
  else if (cell2.sigma == 0)
    return CalculateJwithMed();
  else
    return CalculateJfromKeyLock(cell2.get_keys_bool(), cell2.get_locks_bool());


  // return J[tau][cell2.tau];
  
}

// void Cell::ClearJ(void) {

//   for (int i=0;i<capacity*capacity;i++) {
//     J[0][i]=EMPTY;
//   }
// }


int Cell::CalculateJfromMed(vector<bool>& key2 )
{
  int Jval = 0;
  for (int i = 0; i < par.n_locks; ++i)
  {
    Jval += key2[i]*2; // key2[i] * par.med_table[i]; // third option (not in use)  key2[i] * 2; 
  }
  
  Jval += 2; // += 8 offset so interaction with medium is not 0
  return Jval;
}


int Cell::CalculateJwithMed()
{
  int Jval = 0;
  for (int i = 0; i < par.n_locks; ++i)
  {
    
    Jval += keys_bool[i]*2; // keys_bool[i] * par.med_table[i]; 
  }
  Jval +=2; //  += 8 offset so interaction with medium is not 0
  return Jval;
}


int Cell::CalculateJfromKeyLock(vector<bool>& key2, vector<bool>& lock2 )
{
  int score=0;

  for (int i =0; i < par.n_locks; ++i)
  {
    score += ( keys_bool[i] != lock2[i] )?1:0; // (( keys_bool[i] == lock2[i] )?1:0) * par.med_table[i];
    score += ( key2[i] != locks_bool[i] )?1:0; // (( key2[i] == locks_bool[i] )?1:0) * par.med_table[i];
  }
  // perfect score is 10 (all locks and keys match). 
  // This is fucked need to change it back. 
  int J = 2 + (int)( 10. - 10 * ((double)score) / (2.*par.n_locks)); // 20-16 

  return J; 
}

bool Cell::checkforcycles(int max)
{ 
  auto it = cycles.begin();

  unordered_map<vector<bool>, int> mapIt;

  for (vector<bool> i : cycles)
    ++mapIt[i];

  if (mapIt.size() > max)
    return true;
  else
    return false;

}

void Cell::add_to_vectors()
{
  int j=0;
  int k=0;

  for (int i=0;i<par.n_genes;++i)
  {
    // push back the concentration from respective vector:
    if (i < par.n_diffusers)
    {
      gene_recordings.at(i).push_back(diff_genes.at(i));
    }
    else if (i < par.n_genes - par.n_lockandkey)
      gene_recordings.at(i).push_back(genes.at(i));
    else if (i < par.n_genes - par.n_locks)
    {
      gene_recordings.at(i).push_back(locks.at(j));
      ++j; 
    }
    else 
    {
      gene_recordings.at(i).push_back(keys.at(k));
      ++k;
    }      
  }
  for (int i=0;i<par.n_diffusers;++i)
  {
    gene_recordings.at(par.n_genes + i).push_back(genes.at(i));
  }

}