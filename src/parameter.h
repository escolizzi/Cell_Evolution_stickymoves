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
#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <iostream>
#include <vector>
using namespace std;
class Parameter {
 public:
  Parameter();
  ~Parameter();
  void CleanUp(void);
  int ReadArguments(int argc, char *argv[]);
  void PrintWelcomeStatement(void);
  void Read(const char *filename);
  void CreateRule(const char * Jmed_rule_input);
  void Read_KeyLock_list_fromfile(const char *filename);
  void Write(ostream &os) const;
  double T;
  int target_area;
  int half_div_area;
  int half_div_area_2;
  int target_length;
  double lambda;
  double lambda2;
  char * keylock_list_filename;
  char * Jtable;
  char * Jmed_rule_input;//see below for rule for J with med
  int conn_diss;
  bool vecadherinknockout;
  bool extensiononly;
  int chemotaxis;
  int border_energy;
  int neighbours;
  int min_area_for_life;
  int key_lock_length;
  bool periodic_boundaries;
  int n_chem;
  double * diff_coeff;
  double * decay_rate;
  double * secr_rate;
  double saturation;
  double dt;
  double dx;
  int pde_its;
  int n_init_cells;
  int size_init_cells;
  int sizex;
  int sizey;
  int divisions;
  int mcs;
  int rseed;
  double subfield;
  int relaxation;
  int storage_stride;
  bool graphics;
  bool store;
  char * datadir;
  char * datafile;
  int save_text_file_period;
  char *food_influx_location;
  int initial_food_amount;
  double foodinflux;
  double eatprob;
  double growth;
  double ardecay;
  double gradnoise;
  double gradscale;
  int min_contact_duration_for_preying;
  double frac_contlen_eaten;
  double metabolic_conversion;
  double chancemediumcopied;
  bool readcolortable;
  char *colortable_filename;
  int maxtau;
  double mut_rate;
  double startmu;
  double init_chemmu;
  int persduration;
  int scaling_cell_to_ca_time;
  int howmany_makeit_for_nextgen;
  int popsize;
  char * backupdir;
  char * backupfile;
  int starttime;
  int save_backup_period;
  double init_maintenance_fraction;
  double init_k_mf_0;
  double init_k_mf_A;
  double init_k_mf_P;
  double init_k_mf_C;
  double init_k_ext_0;
  double init_k_ext_A;
  double init_k_ext_P;
  double init_k_ext_C;
  double init_k_chem_0;
  double init_k_chem_A;
  double init_k_chem_P;
  double init_k_chem_C;

  int the_line;
  double init_weight_for_chemotaxis;
  int evolsim;  //will the simulation end after the first time cells arrive at target?
  bool evolreg; //do regulation parameters evolve?
  bool is_there_food; // is there food?
  bool zero_persistence_past_theline; // set persdur to zero after line is crossed
  
  bool season_experiment;
  int season_duration;
  
  struct key_lock_pair{
    int tau;
    vector<int> key;
    vector<int> lock;
  };
  vector< key_lock_pair > keylock_list;

  struct Jmed_rule{
    int offset;
    //char rule; //for now only o //'s' or 'p'
    int keypos_formedium;
    vector<int> lookup_table;
  } Jmedr;

 private:
};

ostream &operator<<(ostream &os, Parameter &p);
const char *sbool(const bool &p);


#endif
