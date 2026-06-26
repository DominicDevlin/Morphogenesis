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

  phase_protein_conc = mother_cell.phase_protein_conc;
  phase_state = mother_cell.phase_state;
  medium_protein_conc = mother_cell.medium_protein_conc;
  medium_state = mother_cell.medium_state;

  temp_hexes = mother_cell.temp_hexes;
  temp_shapes = mother_cell.temp_shapes;

  epithelial = mother_cell.epithelial;

  synNotch_bound = mother_cell.synNotch_bound;
  synNotch_unbound = mother_cell.synNotch_unbound;
  synNotch_intra = mother_cell.synNotch_intra;
  E_cadherin = mother_cell.E_cadherin;
  CD19 = mother_cell.CD19;
  opposing_CD19 = mother_cell.opposing_CD19;
  opposing_E_cadherin = mother_cell.opposing_E_cadherin;
  random_binding_proteins = mother_cell.random_binding_proteins;
  touching_med = mother_cell.touching_med;
  mCherry=mother_cell.mCherry;
  GFP=mother_cell.GFP;
  opposing_GFP=mother_cell.opposing_GFP;
  P_cadherin=mother_cell.P_cadherin;
  N_cadherin=mother_cell.N_cadherin;
  spheroid_cell=mother_cell.spheroid_cell;

  opposing_CD19 = mother_cell.opposing_CD19;
  opposing_E_cadherin = mother_cell.opposing_E_cadherin;
  opposing_mCherry=mother_cell.opposing_mCherry;
  opposing_GFP=mother_cell.opposing_GFP;
  opposing_N_cadherin=mother_cell.opposing_N_cadherin;
  opposing_P_cadherin=mother_cell.opposing_P_cadherin;

  f_opposing_GFP = mother_cell.f_opposing_GFP;
  f_opposing_CD19 = mother_cell.f_opposing_CD19;
  f_opposing_mCherry=mother_cell.f_opposing_mCherry;
  f_opposing_E_cad = mother_cell.f_opposing_E_cad;
  f_opposing_N_cad = mother_cell.f_opposing_N_cad;
  f_opposing_P_cad = mother_cell.f_opposing_P_cad;

  constitutives=mother_cell.constitutives;
  GFP_induced=mother_cell.GFP_induced;
  mCherry_induced=mother_cell.mCherry_induced;
  CD19_induced=mother_cell.CD19_induced;
  
  cell_elastic_mod=mother_cell.cell_elastic_mod;
  motility_strength=mother_cell.motility_strength;
  leftover_area=mother_cell.leftover_area;


  for (int i=0;i<par.n_diffusers;i++)
  {
    diffs[i]=mother_cell.diffs[i];
  }

  velocity_histories_x = mother_cell.velocity_histories_x;
  velocity_histories_y = mother_cell.velocity_histories_y;
  prev_com_x = mother_cell.prev_com_x;
  prev_com_y = mother_cell.prev_com_y;
  avg_vx= mother_cell.avg_vx;
  avg_vy= mother_cell.avg_vy;

  c_type=mother_cell.c_type;
  
  perimeter = mother_cell.perimeter;
  target_perimeter = mother_cell.target_perimeter;
  
  for (int ch=0;ch<par.n_chem;ch++)
    chem[ch]=mother_cell.chem[ch];
  
  n_copies=0;

  grad[0]=mother_cell.grad[0];
  grad[1]=mother_cell.grad[1];

  centerx = mother_cell.centerx;
  centery=mother_cell.centery;
  
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

  perimeter = 0;
  target_perimeter = 0;

  lambda = par.lambda;

  diffs = new double[par.n_diffusers];

  v[0]=0.; v[1]=0.;
  n_copies=0;

  chem = new double[par.n_chem];

  if (par.active_motion)
  {
    velocity_histories_x.assign(par.persistence_time, 0.);
    velocity_histories_y.assign(par.persistence_time, 0.);
  }
  centerx = double(par.sizex)/2 - 1;
  centery = double(par.sizey)/2 - 1;

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
double Cell::EnDif(Cell &cell2)
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


double Cell::SheetDif(Cell &cell2, double &sJ, double &sheetmixJ)
{
  if (sigma==cell2.sigma) 
    return 0;
  else if (sheet_type != cell2.GetSheetType())
  {
    return sheetmixJ;
  }
  return sJ;

}


double Cell::SyntheticEnergy(Cell &cell2)
{
  if (sigma==cell2.sigma)
    return 0;
  else if (sigma == 0)
    return (par.synthetic_Jm)/par.neigh_multiplier;// + par.Jmed_scaling * cell2.getE_cadherin()) / par.neigh_multiplier;
  else if (cell2.sigma==0)
    return (par.synthetic_Jm)/par.neigh_multiplier;// + par.Jmed_scaling * E_cadherin) / par.neigh_multiplier;
  else
  {
    return (par.synthetic_Jcell_baseline 
     - par.JEcadherin_scaling * (E_cadherin * cell2.getE_cadherin())
     - par.Jrandom_scaling_E * ((E_cadherin*cell2.getRandomBindingProteins()) + (random_binding_proteins*cell2.getE_cadherin()))
     - par.JPcadherin_scaling * (P_cadherin * cell2.getP_cadherin()) 
     - par.JNcadherin_scaling * (N_cadherin * cell2.getN_cadherin())
     - par.Jrandom_scaling_N * (N_cadherin * cell2.getRandomBindingProteins() + random_binding_proteins*cell2.getN_cadherin())
     - par.Jrandom_scaling_P * (P_cadherin * cell2.getRandomBindingProteins() + random_binding_proteins*cell2.getP_cadherin())
    ) * 2 / par.neigh_multiplier;
  }
}




double Cell::Melt(Cell &cell2, int x)
{
  if (sigma==cell2.sigma)
  {
    return 0;
  }
  else if (sigma==0)
  {
    // cout << "1: " << phaseJfromMed(cell2.getmJ()) << endl;
    return phaseJfromMed();
  }
  else if (cell2.sigma==0)
  {
    // cout << "2: " <<  PhaseJwithMed() << endl;
    return PhaseJwithMed();
  }
  else
  {
    // cout << "3: " << PhaseJ(cell2.GetPhase(), Jstemdiff) << endl;
    return J_equation(x);
  }
}


double Cell::EnergyDifference(Cell &cell2, bool phase, double Jstemdiff)
{
  if (sigma==cell2.sigma)
  {
    return 0;
  }
  else if (sigma==0)
  {
    // cout << "1: " << phaseJfromMed(cell2.getmJ()) << endl;
    return phaseJfromMed();
  }
  else if (cell2.sigma==0)
  {
    // cout << "2: " <<  PhaseJwithMed() << endl;
    return PhaseJwithMed();
  }
  else
  {
    // cout << "3: " << PhaseJ(cell2.GetPhase(), Jstemdiff) << endl;
    return PhaseJ(cell2.GetPhase(), Jstemdiff, cell2.epithelial);
  }
}

double Cell::PhaseJ(bool &phase, double &Jstemdiff, bool &epith)
{

  if (epithelial && epith)
  {
    return par.epiJ;
  }

  if (epithelial || epith)
  {
    return par.epiJelse;
  }

  if (phase && phase_state)
  {
    return par.J_L;
  }
  else if (!phase && !phase_state)
  {
    return par.J_S;
  }
  else
  {
    return Jstemdiff;
  }
}

double Cell::PhaseJwithMed()
{
  if (epithelial)
    return par.epiM;
  else if (phase_state)
    return par.J_med;
  else
    return par.J_med2;
}

double Cell::phaseJfromMed()
{
  if (epithelial)
    return par.epiM;
  else if (phase_state)
    return par.J_med;
  else
    return par.J_med2;
}


double Cell::J_equation(int x)
{
  double eq = (double(par.J_S-par.J_L) / (1.+exp((double(par.tip_max - par.melt - x))/par.slope))) + par.J_L;
  // double first = double(par.xtip - par.melt - x)/par.slope;
  // double second = 1.0 + exp(first);
  // cout << par.xtip << '\t' << par.melt << '\t' << x << '\t' << par.slope << '\t' << first << '\t' << second << endl;
  return eq; 
}


// return energies by calculating lock & key products switched on by cells. 
double Cell::EnergyDifference(Cell &cell2)
{ 
  if (sigma==cell2.sigma) 
    return 0;
  else if (sigma==0)
    return CalculateJfromMed(cell2.get_medp_bool()) / par.neigh_multiplier; // (cell2.get_medp_bool()); && (cell2.get_keys_bool());
  else if (cell2.sigma == 0)
    return CalculateJwithMed() / par.neigh_multiplier;
  else
    return CalculateJfromKeyLock(cell2.get_keys_bool(), cell2.get_locks_bool()) * 2 / par.neigh_multiplier;


  // return J[tau][cell2.tau];
  
}

// void Cell::ClearJ(void) {

//   for (int i=0;i<capacity*capacity;i++) {
//     J[0][i]=EMPTY;
//   }
// }


double Cell::CalculateJfromMed(vector<bool>& medp2)
{
  double Jval = 0;
  for (int i = 0; i < par.n_mediums; ++i)
  {
    Jval += medp2[i]*par.med_table[i]; // medp2[i] * 4
  }
  
  Jval += par.minM; // += 6 offset so interaction with medium is not 0
  return Jval;
}

// higher J means less binding with medium
double Cell::CalculateJwithMed()
{
  double Jval = 0;
  for (int i = 0; i < par.n_mediums; ++i)
  {
    
    Jval += medp_bool[i]*par.med_table[i]; // medp_bool[i]*4;
  }
  Jval += par.minM; //  += 6 offset so interaction with medium is not 0     
  return Jval;
}


double Cell::CalculateJfromKeyLock(vector<bool>& key2, vector<bool>& lock2 )
{
  int score=0;

  for (int i=0; i < par.n_locks; ++i)
  {
    score += ( keys_bool[i] != lock2[i] )?1:0; // (( keys_bool[i] == lock2[i] )?1:0) * par.med_table[i];
    score += ( key2[i] != locks_bool[i] )?1:0; // (( key2[i] == locks_bool[i] )?1:0) * par.med_table[i];
  }
  double J = par.maxJ - par.interval2 * score; 
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