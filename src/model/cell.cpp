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

  xcen = mother_cell.xcen;
  ycen = mother_cell.ycen;
  xcens = mother_cell.xcens;
  ycens = mother_cell.ycens;

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




// Smooth stand-in for the old hard >0.2 concentration cutoff, so Sox2/Sox17
// -driven adhesion effects ramp up gradually around par.sox_threshold instead
// of jumping at it. This runs on every candidate copy attempt in the CPM's
// inner loop, so it uses a cubic Hermite smoothstep (clamp + multiplies)
// rather than a logistic sigmoid, which would need a exp() call per pixel.
static inline double SoxCommitment(double concentration)
{
  double t = (concentration - (par.sox_threshold - par.sox_threshold_width))
             / (2.0 * par.sox_threshold_width);
  if (t <= 0.0) return 0.0;
  if (t >= 1.0) return 1.0;
  return t * t * (3.0 - 2.0 * t);
}

double Cell::EmbryoEnergy(Cell &cell2, int zona_sigma)
{
  if (sigma==cell2.sigma)
    return 0;
  else if (sigma==0)
  {
    // Undifferentiated (comparable Sox2/Sox17) cells get an extra pull
    // towards the medium, on top of the usual Sox17+ (hypoblast) one, so
    // that unsorted cells are gradually sorted out of the tissue.
    return par.Jblasto - cell2.sox2_internal_adhesion * par.sox2_blasto_adhesion
                        - cell2.sox17_internal_adhesion * par.sox17_blasto_adhesion;
                        //- cell2.IsUndifferentiated() * par.undifferentiated_blasto_adhesion;
  }
  else if (cell2.sigma==0)
  {
    return par.Jblasto - sox2_internal_adhesion * par.sox2_blasto_adhesion;
                       - sox17_internal_adhesion * par.sox17_blasto_adhesion;
                        // - IsUndifferentiated() * par.undifferentiated_blasto_adhesion;
  }
  else if (cell2.sigma==zona_sigma) // 1 is zona pellucida
  {
    // cout << "here" << endl;
    return par.J_cell_zona;
  }
  else
  {
    double t2 = sox2_internal_adhesion;
    double t17 = sox17_internal_adhesion;
    double cell2_t2 = cell2.sox2_internal_adhesion;
    double cell2_t17 = cell2.sox17_internal_adhesion;

    // "Looser" = undifferentiated: t2 and t17 are comparable (both big or
    // both small), so neither lineage clearly dominates; 0.5 if intermediate.
    double is_looser = max(t2 * t17, (1. - t2) * (1. - t17));
    double cell2_is_looser = max(cell2_t2 * cell2_t17, (1. - cell2_t2) * (1. - cell2_t17));
    double one_of_both_loosers = max(is_looser, cell2_is_looser);

    // lambda_epi_epi and lambda_hypo_hypo may be equal. Same for
    // lambda_epi_hypo and lambda_hypo_epi.
    return par.J_cell_baseline - t2 * cell2_t2 * par.sox2binding
                              - t17 * cell2_t17 * par.sox17binding;
                              - t17 * cell2_t2 * par.sox2vs17binding;
                              - t2 * cell2_t17 * par.sox2vs17binding;
    
    // one_of_both_loosers * par.lambda_both_loosers
    //      - (1. - one_of_both_loosers) * (
    //            t2 * cell2_t2 * par.lambda_epi_epi
    //          + t17 * cell2_t17 * par.lambda_hypo_hypo
    //          + t2 * cell2_t17 * par.lambda_epi_hypo
    //          + t17 * cell2_t2 * par.lambda_hypo_epi
    //        );
  }


  // bool this_has_chem = false, other_has_chem = false;
  // for (int ch = 0; ch < par.n_cell_chem; ch++)
  //   if (par.getInitChem(tau, ch) >= 0.0) { this_has_chem = true; break; }
  // for (int ch = 0; ch < par.n_cell_chem; ch++)
  //   if (par.getInitChem(cell2.tau, ch) >= 0.0) { other_has_chem = true; break; }

  // if (this_has_chem && other_has_chem) {
  //   double weighted_dist2 = 0.0;
  //   for (int ch = 0; ch < par.n_cell_chem; ch++) {
  //     double d = chem[ch] - cell2.chem[ch];
  //     weighted_dist2 += par.chem_adhesion_weights[ch] * d * d;
  //   }
  //   return J[tau][cell2.tau] + (int)weighted_dist2;
  // }

  // return J[tau][cell2.tau];
}
