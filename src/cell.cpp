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
#include <array>
#include <utility>
#include <random>
#include <chrono>


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
  medp=mother_cell.medp;
  medp_bool=mother_cell.medp_bool;
  full_set=mother_cell.full_set;
  cycles=mother_cell.cycles;
  gene_recordings=mother_cell.gene_recordings;
  stemness=mother_cell.stemness;
  shrinker = mother_cell.shrinker;

  div_time = mother_cell.div_time;
  div_phen = mother_cell.div_phen;
  phenotype_history = mother_cell.phenotype_history;
  phentime = mother_cell.phentime;
  adulttime = mother_cell.adulttime;
  phenotype = mother_cell.phenotype;
  switches = mother_cell.switches;
  long_switches = mother_cell.long_switches;

  xcen = mother_cell.xcen;
  ycen = mother_cell.ycen;
  xcens = mother_cell.xcens;
  ycens = mother_cell.ycens;
  vel_phens = mother_cell.vel_phens;

  gamma_list = mother_cell.gamma_list;
  mass_list = mother_cell.mass_list;
  time_created = mother_cell.time_created;


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
    return CalculateJfromMed(cell2.get_medp_bool()); // (cell2.get_medp_bool()); && (cell2.get_keys_bool());
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


int Cell::CalculateJfromMed(vector<bool>& medp2)
{
  int Jval = 0;
  for (int i = 0; i < par.n_mediums; ++i)
  {
    Jval += medp2[i]*par.med_table[i]; // medp2[i] * 4
  }
  
  Jval += par.minM; // += 6 offset so interaction with medium is not 0
  return Jval;
}

// higher J means less binding with medium
int Cell::CalculateJwithMed()
{
  int Jval = 0;
  for (int i = 0; i < par.n_mediums; ++i)
  {
    
    Jval += medp_bool[i]*par.med_table[i]; // medp_bool[i]*4;
  }
  Jval += par.minM; //  += 6 offset so interaction with medium is not 0     
  return Jval;
}


int Cell::CalculateJfromKeyLock(vector<bool>& key2, vector<bool>& lock2 )
{
  int score=0;

  for (int i=0; i < par.n_locks; ++i)
  {
    score += ( keys_bool[i] != lock2[i] )?1:0; // (( keys_bool[i] == lock2[i] )?1:0) * par.med_table[i];
    score += ( key2[i] != locks_bool[i] )?1:0; // (( key2[i] == locks_bool[i] )?1:0) * par.med_table[i];
  }
  int J = par.maxJ - par.interval2 * score; 
  // perfect score is 10 (all locks and keys match). 
  //  int J = 4 + (int)( 8. - 8 * ((double)score / par.n_lockandkey)); //4 10 10     20-16 
  return J; 
}

bool Cell::checkforcycles(int max)
{ 
  auto it = cycles.begin();

  unordered_map<vector<bool>, int> mapIt{};

  for (vector<bool> i : cycles)
    ++mapIt[i];

  if (mapIt.size() > max)
    return true;
  else
    return false;

}

//set a specific phenotype code. 
void Cell::Phenotype()
{
  int pcode{};
  int tot = full_set.size();
  for (int i=0;i<tot;++i)
  {
    int x = tot - 1 - i;
    pcode += full_set[i] * pow(2,x);
  }
  phenotype = pcode;
}

int Cell::RegPhenotype()
{

  int code{};
  for (int i=0;i<par.n_activators;++i)
  {
    bool val;
    if (genes[i] > 0.5)
      val = 1;
    else
      val = 0;
    int x = par.n_activators - 1 - i;
    code += val * pow(2,x);
  }
  return code;
}




void Cell::RecordLongSwitch(vector<bool> &v1, uint64_t rndm)
{
  int p1{};
  int p2{};
  int tot = full_set.size();

  for (int i=0;i<tot;++i)
  {
    int x = tot - 1 - i;
    p1 += v1[i] * pow(2,x);
    p2 += full_set[i] * pow(2,x);
  }

  if (p1 != p2)
  {
    //
    tuple<int,int, int> sch = {p1, p2, rndm};
    long_switches.push_back(sch);
    
  }
}


void Cell::RecordSwitch(vector<bool> &v1, uint64_t rndm)
{
  int p1{};
  int p2{};
  int tot = full_set.size();

  for (int i=0;i<tot;++i)
  {
    int x = tot - 1 - i;
    p1 += v1[i] * pow(2,x);
    p2 += full_set[i] * pow(2,x);
  }

  if (p1 != p2)
  {
    //
    tuple<int,int, int> sch = {p1, p2, rndm};
    switches.push_back(sch);
  }
}



void Cell::add_to_vectors()
{
  int j=0;
  int k=0;
  int m=0;
  for (int i=0;i<par.n_genes;++i)
  {
    // push back the concentration from respective vector:
    if (i < par.n_diffusers)
    {
      gene_recordings.at(i).push_back(diff_genes.at(i));
    }
    else if (i < par.n_genes - par.n_lockandkey - par.n_mediums)
    {
      gene_recordings.at(i).push_back(genes.at(i));

    }
    else if (i < par.n_genes - par.n_locks - par.n_mediums)
    {
      gene_recordings.at(i).push_back(locks.at(j));
      ++j; 
    }
    else if (i < par.n_genes - par.n_mediums) 
    {
      gene_recordings.at(i).push_back(keys.at(k));
      ++k;
    }      
    else 
    {
      gene_recordings.at(i).push_back(medp[m]);
      ++m;
    }
  }
  for (int i=0;i<par.n_diffusers;++i)
  {
    gene_recordings.at(par.n_genes + i).push_back(genes.at(i));
  }
  // get a specific phenotype code
  phenotype_history.push_back(phenotype);

}

int Cell::LocksKeysScore(vector<bool>& locks, vector<bool>& keys)
{
  int score{};
  for (int i =0; i < par.n_locks; ++i)
  {
    score += ( keys_bool[i] != locks[i] )?1:0; // (( keys_bool[i] == lock2[i] )?1:0) * par.med_table[i];
    score += ( keys[i] != locks_bool[i] )?1:0; // (( key2[i] == locks_bool[i] )?1:0) * par.med_table[i];
  }
  return score;
}

int Cell::CheckMedsOn()
{
  int n{};
  for (int i=0; i < par.n_locks; ++i)
  {
    n += medp_bool[i]; // medp_bool[i]; && keys_bool[i];
  }
  return n;
}


  bool Cell::limit_cycle()
  {
    vector<int> temp{};
    temp.push_back(0);
    int n_cycles=0;
    for (auto &i : phenotype_history)
    {
      if (temp.back() == i)
      {
        continue;
      }
      else 
      {
        int count = 0;
        for (int x = temp.size()-1; x >= 0; x--)
        {
          // cout << i << " " << temp[x] << endl;
          if (i == temp[x])
          {
            if (count > 12)
            {
              ++n_cycles;
            }
            break;
          } 

          ++count;
        }
        temp.push_back(i);
      }
    }
    // cout << "TOTAL CYCLES: " << n_cycles << endl;
    if (n_cycles > 10)
      return true;
    else
      return false;
  }