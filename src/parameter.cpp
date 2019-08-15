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


#include "parameter.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "output.h"
#include "parse.h"


//parameter constructor - initialise default?
Parameter::Parameter() {

  T = 50.;
  target_area = 50;
  half_div_area = 50 ;
  half_div_area_2 = -1;
  target_length = 60;
  lambda = 50;
  lambda2 = 5.0;
  //Jtable = strdup("J.dat");
  keylock_list_filename = strdup("KLcellevol.dat");
  conn_diss = 2000;
  vecadherinknockout = false;
  extensiononly = false;
  chemotaxis = 1000;
  border_energy = 100;
  neighbours = 2;
  min_area_for_life = 5;
  key_lock_length = 10;
  periodic_boundaries = false;
  n_chem = 1;
  diff_coeff = new double[1];
  diff_coeff[0] = 1e-13;
  decay_rate = new double[1];
  decay_rate[0] = 1.8e-4;
  secr_rate = new double[1];
  secr_rate[0] = 1.8e-4;
  saturation = 0;
  dt = 2.0;
  dx = 2.0e-6;
  pde_its = 15;
  n_init_cells = 100;
  size_init_cells = 10;
  sizex = 200;
  sizey = 200;
  divisions = 0;
  mcs = 10000;
  rseed = -1;
  subfield = 1.0;
  relaxation = 0;
  storage_stride = 10;
  graphics = true;
  store = false;
  datadir = strdup("data_film");
  datafile = strdup("data_cellcount.txt");
  save_text_file_period = 100;
  food_influx_location = strdup("nowhere");
  foodinflux=0.;
  eatprob=0.;
  growth=0;
  ardecay=0.;
  gradnoise=0.1;
  gradscale=1.0;
  min_contact_duration_for_preying = 10;
  frac_contlen_eaten = 1.;
  metabolic_conversion = 0.5;
  readcolortable = false;
  colortable_filename = "default.ctb";
  mut_rate=0.01;
  startmu=0.0;
  persduration=0;
  scaling_cell_to_ca_time = 1;
  backupdir=strdup("backup");
  save_backup_period=0;
  init_maintenance_fraction = 0.85;
  init_chemmu=0.;
  backupfile=strdup("");
  starttime=0;
  init_k_mf_0 = 0.5;
  init_k_mf_A = 0.;
  init_k_mf_P = 0.;
  init_k_mf_C = 0.;

  init_k_ext_0 = 0.5;
  init_k_ext_A = 0.;
  init_k_ext_P = 0.;
  init_k_ext_C = 0.;

  init_weight_for_chemotaxis=0.;
  init_k_chem_0 = 0.5;
  init_k_chem_A = 0.;
  init_k_chem_P = 0.;
  init_k_chem_C = 0.;

  howmany_makeit_for_nextgen=30;
  popsize=100;
  the_line=50;
  evolsim=0;
  is_there_food=false;
}

Parameter::~Parameter() {

  // destruct parameter object

  // free string parameter

  CleanUp();

}

void Parameter::CleanUp(void) {
  if (Jtable)
     free(Jtable);
  if (diff_coeff)
     free(diff_coeff);
  if (decay_rate)
     free(decay_rate);
  if (secr_rate)
     free(secr_rate);
  if (datadir)
     free(datadir);
  if (backupdir)
     free(backupdir);
  if (backupfile)
        free(backupfile);

}
void Parameter::PrintWelcomeStatement(void)
{
  cerr<<"CellEvol: v0.something (very much a prototype)"<<endl;
  cerr<<"Usage is: "<<endl;
  cerr<<"./cell_evolution path/to/data [optional arguments]"<<endl;
  cerr<<"Arguments: "<<endl;
  cerr<<" -datafile path/to/datafile"<<endl;
  cerr<<" -datadir path/to/datadir"<<endl;
  cerr<<" -backupdir path/to/backupdir"<<endl;
  cerr<<" -keylockfilename path/to/keylockfilename"<<endl;
  cerr<<" -seed INT_NUMBER"<<endl;
  cerr<<" -maxtime INT_NUMBER"<<endl;
  cerr<<" -halfdiv_area_predator INT_NUMBER"<<endl;
  cerr<<" -persduration INT_NUMBER"<<endl;
  cerr<<" -mutrate FLOAT_NUMBER [0,1)"<<endl;
  cerr<<endl<<"Will not execute if datafile and datadir already exist"<<endl;
  cerr<<"Also, parameter file and Jtable should be in the same directory (unless you used option -keylockfilename)"<<endl;
  cerr<<"Have fun!"<<endl;
}

int Parameter::ReadArguments(int argc, char *argv[])
{
  cerr<<endl<<"Reading arguments from command line"<<endl;
  //starts from 2 because 0 is filename, 1 is parameter file path
  for(int i = 2;i<argc;i++){
    if( 0==strcmp( argv[i],"-datafile") ){
      i++; if(i==argc){
        cerr<<"Something odd in datafile?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      //strcpy(datafile, argv[i]); //this can be buggy because it copies over a dynamically allocated char* (datafile) that can be a lot shorter
      free(datafile);
      datafile = (char *)malloc( 5+strlen(argv[i])*sizeof(char) )  ; //strlen(argv[i]) is ok because argv[i] is null terminated
      datafile = strdup(argv[i]);
      cerr<<"New value for datafile: "<<datafile<<endl;
//       exit(1);
    }else if( 0==strcmp(argv[i],"-datadir") ){
      i++; if(i==argc) {
        cerr<<"Something odd in datadir?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      //strcpy(datadir, argv[i]);
      free(datadir);
      datadir = (char *)malloc( 5+strlen(argv[i])*sizeof(char) )  ; //strlen(argv[i]) is ok because argv[i] is null terminated
      datadir = strdup(argv[i]);

      cerr<<"New value for datadir: "<<datadir<<endl;

    }else if( 0==strcmp(argv[i],"-backupdir") ){
      i++; if(i==argc) {
        cerr<<"Something odd in backupdir?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      //strcpy(datadir, argv[i]);
      free(backupdir);
      backupdir = (char *)malloc( 5+strlen(argv[i])*sizeof(char) )  ; //strlen(argv[i]) is ok because argv[i] is null terminated
      backupdir = strdup(argv[i]);

      cerr<<"New value for backupdir: "<<backupdir<<endl;

    }else if( 0==strcmp(argv[i],"-keylockfilename") ){
      i++; if(i==argc) {
        cerr<<"Something odd in keylockfilename?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      //strcpy(keylock_list_filename, argv[i]);
      free(keylock_list_filename);
      keylock_list_filename = (char *)malloc( 5+strlen(argv[i])*sizeof(char) )  ; //strlen(argv[i]) is ok because argv[i] is null terminated
      keylock_list_filename = strdup(argv[i]);

      cerr<<"New value for keylock_list_filename: "<<keylock_list_filename<<endl;

    }else if( 0==strcmp(argv[i],"-seed") ){
      i++; if(i==argc){
        cerr<<"Something odd in seed?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      rseed = atoi( argv[i] );
      cerr<<"New value for rseed: "<<rseed<<endl;
    }else if( 0==strcmp(argv[i],"-maxtime") ){
      i++; if(i==argc){
        cerr<<"Something odd in maxtime?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      mcs = atoi( argv[i] );
      cerr<<"New value for maxtime (mcs in the code): "<<mcs<<endl;
    }else if( 0==strcmp(argv[i],"-halfdiv_area_predator") ){
      i++; if(i==argc){
        cerr<<"Something odd in halfdiv_area_predator?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      half_div_area_2 = atoi( argv[i] );
      cerr<<"New value for half division area of predators (half_div_area_2 in the code): "<<half_div_area_2<<endl;
    }else if( 0==strcmp(argv[i],"-mutrate") ){
      i++; if(i==argc){
        cerr<<"Something odd in mutrate?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      mut_rate = atof( argv[i] );
      cerr<<"New value for mutation rate: "<<mut_rate<<endl;
    }else if( 0==strcmp(argv[i],"-persduration") ){
      i++; if(i==argc){
        cerr<<"Something odd in persduration?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      persduration = atoi( argv[i] );
      cerr<<"New value for persistence of movement: "<<persduration<<endl;
    }else if( 0==strcmp(argv[i],"-backupfile") ){
      i++; if(i==argc){
        cerr<<"Something odd in backupfile?"<<endl;
        return 1;  //check if end of arguments, exit with error in case
      }
      free(backupfile);
      backupfile = (char *)malloc( 5+strlen(argv[i])*sizeof(char) )  ; //strlen(argv[i]) is ok because argv[i] is null terminated
      backupfile = strdup(argv[i]);

      cerr<<"New value for backupfile: "<<backupfile<<endl;
    }else
      return 1;
  }
  return 0;
}

void Parameter::Read(const char *filename) {

  static bool ReadP=false;

  if (ReadP) {

    //throw "Run Time Error in parameter.cpp: Please Read parameter file only once!!";
    CleanUp();

  } else
    ReadP=true;

  FILE *fp=OpenReadFile(filename);


  T = fgetpar(fp, "T", 50., true);
  target_area = igetpar(fp, "target_area", 100, true);
  half_div_area = igetpar(fp, "half_div_area", 100, true);
  half_div_area_2 = igetpar(fp, "half_div_area_2", 100, true);
  target_length = igetpar(fp, "target_length", 60, true);
  lambda = fgetpar(fp, "lambda", 50, true);
  lambda2 = fgetpar(fp, "lambda2", 5.0, true);
  //Jtable = sgetpar(fp, "Jtable", "J.dat", true);
  keylock_list_filename = sgetpar(fp, "keylock_list_filename", "KLcellevol.dat", true);
  conn_diss = igetpar(fp, "conn_diss", 2000, true);
  vecadherinknockout = bgetpar(fp, "vecadherinknockout", false, true);
  extensiononly = bgetpar(fp, "extensiononly", false, true);
  chemotaxis = igetpar(fp, "chemotaxis", 1000, true);
  border_energy = igetpar(fp, "border_energy", 100, true);
  neighbours = igetpar(fp, "neighbours", 2, true);
  min_area_for_life = igetpar(fp, "min_area_for_life", 5, true);
  key_lock_length = igetpar(fp, "key_lock_length", 10, true);
  periodic_boundaries = bgetpar(fp, "periodic_boundaries", false, true);
  n_chem = igetpar(fp, "n_chem", 0, true);
  if(n_chem){
    diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, true);
    decay_rate = dgetparlist(fp, "decay_rate", n_chem, true);
    secr_rate = dgetparlist(fp, "secr_rate", n_chem, true);
    saturation = fgetpar(fp, "saturation", 0, true);
    dt = fgetpar(fp, "dt", 2.0, true);
    dx = fgetpar(fp, "dx", 2.0e-6, true);
    pde_its = igetpar(fp, "pde_its", 15, true);
  }
  n_init_cells = igetpar(fp, "n_init_cells", 100, true);
  size_init_cells = igetpar(fp, "size_init_cells", 10, true);
  sizex = igetpar(fp, "sizex", 200, true);
  sizey = igetpar(fp, "sizey", 200, true);
  divisions = igetpar(fp, "divisions", 0, true);
  mcs = igetpar(fp, "mcs", 10000, true);
  rseed = igetpar(fp, "rseed", -1, true);
  subfield = fgetpar(fp, "subfield", 1.0, true);
  relaxation = igetpar(fp, "relaxation", 0, true);
  storage_stride = igetpar(fp, "storage_stride", 10, true);
  graphics = bgetpar(fp, "graphics", true, true);
  store = bgetpar(fp, "store", false, true);
  datadir = sgetpar(fp, "datadir", "data_film", true);
  datafile = sgetpar(fp,"datafile" , "data_cellcount.txt",true);
  save_text_file_period = igetpar(fp, "save_text_file_period", 100, true);
  food_influx_location = sgetpar(fp,"food_influx_location" , "nowhere",true);
  initial_food_amount = fgetpar(fp, "initial_food_amount", 0, true);
  foodinflux = fgetpar(fp, "foodinflux", 0., true);
  eatprob = fgetpar(fp, "eatprob", 0., true);
  ardecay = fgetpar(fp, "ardecay", 0., true);
  growth = fgetpar(fp, "growth", 0., true);
  gradnoise = fgetpar(fp, "gradnoise", 0.1, true); //did I put these in?
  gradscale = igetpar(fp, "gradscale", 1, true);
  min_contact_duration_for_preying = fgetpar(fp, "min_contact_duration_for_preying", 1., true);
  frac_contlen_eaten = fgetpar(fp, "frac_contlen_eaten", 1., true);
  metabolic_conversion = fgetpar(fp, "metabolic_conversion", 0.5, true);
  chancemediumcopied = fgetpar(fp, "chancemediumcopied", 0.0001, true);
  readcolortable = bgetpar(fp, "readcolortable", false, true);
  colortable_filename = sgetpar(fp,"colortable_filename" , "default.ctb",true);
  evolsim = igetpar(fp, "evolsim", 0, true);
  mut_rate = fgetpar(fp, "mut_rate", 0.01, true);
  persduration = igetpar(fp, "persduration", 0, true);
  startmu = fgetpar(fp, "startmu", 0.0, true);
  init_chemmu = fgetpar(fp, "init_chemmu", 0.0, true);
  Jmed_rule_input = sgetpar(fp, "Jmed_rule_input", "0a0", true);
  scaling_cell_to_ca_time = igetpar(fp, "scaling_cell_to_ca_time", 1, true);
  backupdir = sgetpar(fp, "backupdir", "backup", true);
  save_backup_period = igetpar(fp, "save_backup_period", 0, true);
  init_maintenance_fraction = fgetpar(fp, "init_maintenance_fraction", 1.0, true);
  init_k_mf_0 = fgetpar(fp, "init_k_mf_0", 0., true);
  init_k_mf_A = fgetpar(fp, "init_k_mf_A", 0., true);
  init_k_mf_P = fgetpar(fp, "init_k_mf_P", 0., true);
  init_k_mf_C = fgetpar(fp, "init_k_mf_C", 0., true);
  init_k_ext_0 = fgetpar(fp, "init_k_ext_0", 0., true);
  init_k_ext_A = fgetpar(fp, "init_k_ext_A", 0., true);
  init_k_ext_P = fgetpar(fp, "init_k_ext_P", 0., true);
  init_k_ext_C = fgetpar(fp, "init_k_ext_C", 0., true);
  init_weight_for_chemotaxis = fgetpar(fp, "init_weight_for_chemotaxis", 0., true);
  init_k_chem_0 = fgetpar(fp, "init_k_chem_0", 0., true);
  init_k_chem_A = fgetpar(fp, "init_k_chem_A", 0., true);
  init_k_chem_P = fgetpar(fp, "init_k_chem_P", 0., true);
  init_k_chem_C = fgetpar(fp, "init_k_chem_C", 0., true);
  howmany_makeit_for_nextgen = igetpar(fp, "howmany_makeit_for_nextgen", 1, true);
  popsize = igetpar(fp, "popsize", 1, true);
  the_line = igetpar(fp, "the_line", 1, true);
  is_there_food = bgetpar(fp,"is_there_food",false, true);
  evolreg = bgetpar(fp,"evolreg",false, true);
}

//creates a rule for lookup table, by setting values,
// and setting a pointer to a function that sums values or multiplies them 8O
// typical input looks like 10o5_4_3_2_1
// 10 is the offset, rest is lookup_table, of which we also need to get the length
void Parameter::CreateRule(const char * Jmed_rule_input)
{
  int i=0;
  //char crule='n'; //means not set yet
  char coffset[10]={0};
  char cvalue[10]={0};
  int cvalcounter=0;

  while(Jmed_rule_input[i]!='o'){
    coffset[i]=Jmed_rule_input[i];
    i++;
  }
  Jmedr.offset=atoi(coffset);
  i++;
  while(Jmed_rule_input[i] != '\0'){
    if(Jmed_rule_input[i]=='_'){
      int value = atoi(cvalue);
      Jmedr.lookup_table.push_back(value);
      for(int j=0;j<10;j++) cvalue[j]='\0';
      i++;
      cvalcounter=0;
    }
    cvalue[cvalcounter]=Jmed_rule_input[i];
    i++;
    cvalcounter++;
  }
  if(cvalue[0]!='\0'){
    int value = atoi(cvalue);
    Jmedr.lookup_table.push_back(value);
  }

  Jmedr.keypos_formedium = static_cast<int>(Jmedr.lookup_table.size()); // this forces conversion from size to INT
                                                           // I can't imagine a way for this to overflow, but you know...
  cerr<<"Offset : "<< Jmedr.offset<<", Table:";
  for(auto x: Jmedr.lookup_table) cerr<<" "<<x;
  cerr<<endl;
  cerr<<"Length = "<< Jmedr.keypos_formedium<<endl;
  //exit(1);

}

// In the future the parser for the rules for key to J val tau,medium
// will be more developed, maybe even evolvable 8O
// int Parameter::SumLookupTableValue(int *lookup_table){
//   return -1;
// }
// int Parameter::MultiplyLookupTableValue(int *lookup_table){
//   return -1;
// }

// key lock pair files start with the initial number of taus (incl medium)
// then there is key, then lock, then key then lock etc...
void Parameter::Read_KeyLock_list_fromfile(const char *filename)
{
  string line;
  int current_tau=0;
  ifstream file(filename);

  cerr<<"Reading Key Lock list file"<<endl;

  //getline gets next line in file
  while( getline(file, line) ){
    istringstream iss(line);   // turn this into array of int and put it into an array of arrays...
    int a;

    key_lock_pair this_kl;
    vector<int> key;
    vector<int> lock;

    if(current_tau==0){
      if(!(iss >> a)) {
        cerr<<"Read_KeyLock_list_fromfile(): Error, the file for initial keys and locks seems empty."<<endl;
        exit(1);
      }else{
        //since keylock_list is indicised with tau and tau=0 is medium
        // we let the first element of keylock_list be a "mock" entry
        this_kl.tau = 0;
        this_kl.key = vector<int>(1, -1);
        this_kl.lock= vector<int>(1, -1);
        keylock_list.push_back(this_kl);
        // end of mockery :P

        maxtau=a-1;
        cerr<<"Got maxtau = "<<maxtau<<endl;
        //increase tau
        current_tau++;
      }
    }else{
      //cerr<<"Into reading key-lock pairs"<<endl;
      //get key for this tau
      while(iss >> a) key.push_back(a);

      //we read key, next line is lock
      if( !getline(file, line) ){
        cerr<<"Read_KeyLock_list_fromfile(): Error, odd number of lines?"<<endl;
        exit(1);
      }
      //iss.str("");
      //iss.clear();

      istringstream iss( (line) );
      //get lock for this tau
      while(iss >> a) lock.push_back(a);

      //assign these values to key lock pair
      this_kl.tau = current_tau;
      this_kl.key = key;
      this_kl.lock = lock;

      keylock_list.push_back(this_kl);

      cerr<<"New key-lock pairs for tau = "<<this_kl.tau<<endl;
      for (auto i: this_kl.key)
          cerr << i << " ";
      cerr<<endl;
      for (auto i: this_kl.lock)
          cerr << i << " ";
      cerr<<endl;

      key.clear();
      lock.clear();
      //increase tau
      current_tau++;

    }


    // process pair (a,b)
  }

}

const char *sbool(const bool &p) {

  const char *true_str="true";
  const char *false_str="false";
  if (p)
    return true_str;
  else
    return false_str;
}

void Parameter::Write(ostream &os) const {
  setlocale(LC_NUMERIC, "C");

  os << " T = " << T << endl;
  os << " target_area = " << target_area << endl;
  os << " div_area = " << half_div_area << endl;
  os << " div_area_2 = " << half_div_area_2 << endl;
  os << " target_length = " << target_length << endl;
  os << " lambda = " << lambda << endl;
  os << " lambda2 = " << lambda2 << endl;
  if(keylock_list_filename)
    os << " keylock_list_filename = " << keylock_list_filename << endl;
  //if (Jtable)
  //  os << " Jtable = " << Jtable << endl;
  os << " conn_diss = " << conn_diss << endl;
  os << " vecadherinknockout = " << sbool(vecadherinknockout) << endl;
  os << " extensiononly = " << sbool(extensiononly) << endl;
  os << " chemotaxis = " << chemotaxis << endl;
  os << " border_energy = " << border_energy << endl;
  os << " neighbours = " << neighbours << endl;
  os << " min_area_for_life = " << min_area_for_life << endl;
  os << " key_lock_length = " << key_lock_length << endl;
  os << " periodic_boundaries = " << sbool(periodic_boundaries) << endl;
  os << " n_chem = " << n_chem << endl;
  os << " diff_coeff = "<< diff_coeff[0] << endl;
  os << " decay_rate = "<< decay_rate[0] << endl;
  os << " secr_rate = "<< secr_rate[0] << endl;
  os << " saturation = " << saturation << endl;
  os << " dt = " << dt << endl;
  os << " dx = " << dx << endl;
  os << " pde_its = " << pde_its << endl;
  os << " n_init_cells = " << n_init_cells << endl;
  os << " size_init_cells = " << size_init_cells << endl;
  os << " sizex = " << sizex << endl;
  os << " sizey = " << sizey << endl;
  os << " divisions = " << divisions << endl;
  os << " mcs = " << mcs << endl;
  os << " rseed = " << rseed << endl;
  os << " subfield = " << subfield << endl;
  os << " relaxation = " << relaxation << endl;
  os << " storage_stride = " << storage_stride << endl;
  os << " graphics = " << sbool(graphics) << endl;
  os << " store = " << sbool(store) << endl;
  os << " initial_food_amount = "<< initial_food_amount << endl;
  os << " food_influx_location = "<< food_influx_location <<endl;
  os << " foodinflux = " << foodinflux << endl;
  os << " eatprob = " << eatprob << endl;
  os << " ardecay = " << ardecay << endl;
  os << " growth = " << growth << endl;
  os << " gradnoise = " << gradnoise << endl;
  os << " gradscale = " << gradscale << endl;
  os << " min_contact_duration_for_preying = " << min_contact_duration_for_preying;
  os << " frac_contlen_eaten = " << frac_contlen_eaten << endl;
  os << " metabolic_conversion = " << metabolic_conversion << endl;
  os << " chancemediumcopied = " << chancemediumcopied << endl;
  os << " datafile = " << datafile << endl;
  os << " save_text_file_period = " << save_text_file_period << endl;
  os << " readcolortable = " << readcolortable <<endl;
  os << " colortable_filename = " << colortable_filename <<endl;
  os << " mut_rate = " << mut_rate <<endl;
  os << " evolsim = " << evolsim <<endl;
  os << " persduration = " << persduration <<endl;
  os << " startmu = " << startmu <<endl;
  os << " scaling_cell_to_ca_time = " << scaling_cell_to_ca_time <<endl;
  os << " backupdir = " << backupdir <<endl;
  os << " save_backup_period = " << save_backup_period <<endl;
  if (datadir)
    os << " datadir = " << datadir << endl;
  //os << " init_maintenance_fraction = " << init_maintenance_fraction << endl;
  os << " init_k_mf_0 = " << init_k_mf_0 << endl;
  os << " init_k_mf_A = " << init_k_mf_A << endl;
  os << " init_k_mf_P = " << init_k_mf_P << endl;
  os << " init_k_mf_C = " << init_k_mf_C << endl;
  os << " init_k_ext_0 = " << init_k_ext_0 << endl;
  os << " init_k_ext_A = " << init_k_ext_A << endl;
  os << " init_k_ext_P = " << init_k_ext_P << endl;
  os << " init_k_ext_C = " << init_k_ext_C << endl;
  os << " init_k_chem_0 = " << init_k_chem_0 << endl;
  os << " init_k_chem_A = " << init_k_chem_A << endl;
  os << " init_k_chem_P = " << init_k_chem_P << endl;
  os << " init_k_chem_C = " << init_k_chem_C << endl;
  //os << " init_weight_for_chemotaxis = " << init_weight_for_chemotaxis << endl;
  os << " howmany_makeit_for_nextgen = " <<  howmany_makeit_for_nextgen << endl;
  os << " popsize = " << popsize << endl;
  os << " the_line = " << the_line <<endl;
  os << " is_there_food = " << is_there_food << endl;
  os << " evolreg = " << evolreg <<endl;
}

ostream &operator<<(ostream &os, Parameter &p) {
  p.Write(os);
  return os;
}

Parameter par;
