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
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <functional> // for bind, see InitIncreaseVal()
#include "crash.h"
#include "parameter.h"
#include "ca.h"
#include "intplane.h"
#include "conrec.h"

using std::placeholders::_1;

/* STATIC DATA MEMBER INITIALISATION */
const int IntPlane::nx[9] = {0, 1, 1, 1, 0,-1,-1,-1, 0 };
const int IntPlane::ny[9] = {0, 1, 0,-1,-1,-1, 0, 1, 1 };

extern Parameter par;


/** PRIVATE **/

IntPlane::IntPlane(const int sx, const int sy) {
  
  sigma=0;
  thetime=0;
  sizex=sx;
  sizey=sy;
  
  sigma=AllocateSigma(sx,sy);
}


IntPlane::IntPlane(void) {

  sigma=0;
  sizex=0; sizey=0;
  thetime=0;
  
}

// destructor (virtual)
IntPlane::~IntPlane(void) {
  if (sigma) {
    free(sigma[0]);
    free(sigma);
    sigma=0;
  }
}

int **IntPlane::AllocateSigma(const int sx, const int sy) {
  
  int **mem;
  sizex=sx; sizey=sy;
  
  mem=(int **)malloc(sizex*sizeof(int *));
  
  if (mem==NULL)
    MemoryWarning();
  
  mem[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (mem[0]==NULL)  
      MemoryWarning();
  
  {  for (int i=1;i<sizex;i++) 
    mem[i]=mem[i-1]+sizey;}
  
  /* Clear IntPlane plane */
  { for (int i=0;i<sizex*sizey;i++) 
    mem[0][i]=0.; }

   return mem;
}

void IntPlane::Plot(Graphics *g2) {
  // l=layer: default layer is 0
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
      // Make the pixel four times as large
      // to fit with the CPM plane
      g2->Point(sigma[x][y],2*x,2*y);
      g2->Point(sigma[x][y],2*x+1,2*y);
      g2->Point(sigma[x][y],2*x,2*y+1);
      g2->Point(sigma[x][y],2*x+1,2*y+1);
    } 
  
}

// Plot the value of the intplane only in the medium of the CPM
void IntPlane::Plot(Graphics *g2, CellularPotts *cpm) {
  
  // this has to take into account stuff from cpm plane (maybe x,y info should give a tau of the cpm plane)
  
  // cpm->sigma[x][y] returns sigma, which I can use to indicise the vector of cells... can I? 
  // not really, food doesn't know about cell
  
  // maybe this function should really be in dish...
  
  // suspend=true suspends calling of DrawScene
  for(int x=1;x<sizex-1;x++)
    for(int y=1;y<sizey-1;y++) 
      //if (cpm->Sigma(x,y)==0) {
      if(sigma[x][y]!=0){
	    // Make the pixel four times as large
	    // to fit with the CPM plane
	    g2->Point(10+sigma[x][y],2*x,2*y);
        g2->Point(10+sigma[x][y],2*x+1,2*y);
        g2->Point(10+sigma[x][y],2*x,2*y+1);
        g2->Point(10+sigma[x][y],2*x+1,2*y+1);
      }
}

// private
void IntPlane::NoFluxBoundaries(void) {

  // all gradients at the edges become zero, 
  // so nothing flows out
  // Note that four corners points are not defined (0.)
  // but they aren't used in the calculations
  
  
    for (int x=0;x<sizex;x++) {
      sigma[x][0]=sigma[x][1];
      sigma[x][sizey-1]=sigma[x][sizey-2];
    }
  
    for (int y=0;y<sizey;y++) {
      sigma[0][y]=sigma[1][y];
      sigma[sizex-1][y]=sigma[sizex-2][y];
    }
  
}


// private
void IntPlane::AbsorbingBoundaries(void) {

  // all boundaries are sinks, 
  
    for (int x=0;x<sizex;x++) {
      sigma[x][0]=0.;
      sigma[x][sizey-1]=0.;
    }
  
    for (int y=0;y<sizey;y++) {
      sigma[0][y]=0.;
      sigma[sizex-1][y]=0.;
    }
  
}

// private
void IntPlane::PeriodicBoundaries(void) {

  // periodic...
  
    for (int x=0;x<sizex;x++) {
      sigma[x][0]=sigma[x][sizey-2];
      sigma[x][sizey-1]=sigma[x][1];
    }
    for (int y=0;y<sizey;y++) {
      sigma[0][y]=sigma[sizex-2][y];
      sigma[sizex-1][y]=sigma[1][y];
    }
  
}

//copy of this function in ca.cpp
int IntPlane::SetNextVal(int pos){
  //the plane has a 1 px boundary on all size, therefore we place the pixels
  //within that boundary 
  static int xcount=1, ycount=1;
  
  if(xcount>=sizex-1 ||ycount>=sizey-1){
    return 1;
  }
  
  sigma[xcount][ycount]=pos;
  ycount++;
  if(ycount==sizey-1){
    ycount=1;
    xcount++;
  }
  return 0;
}

// This function initialises a functional (from <function>) 
// depending on parameters. The function is responsible for updating the field
// If 'nowhere' option is used, parameter food_influx should be set
void IntPlane::InitIncreaseVal(CellularPotts *cpm) {
  
  // change values of initial_food_amount and others 
  // if food_influx_location== "nowhere"
  if( strcmp(par.food_influx_location,"everywhere") == 0 ){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValEverywhere, this);
    
    //exit(1);
  }else if(strcmp(par.food_influx_location,"nowhere") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    par.foodinflux=0.; //this is not necessary, because NotIncreaseVal is an empty function, bu well...
    IncreaseVal = std::bind(&IntPlane::NotIncreaseVal, this);
//     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"notonprey") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValIfEmpty, this, cpm);
//     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"patchy") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValPatchy, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"somewhere") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSomewhere, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"somewhere_notonprey") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSomewhereIfEmpty, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"patchy_random") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValPatchyRandom, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"patchy_random_persistence") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValPatchyRandomPersistence, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"food_growth") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSelfGrowth, this, cpm);
    //     exit(1);
    ;
  }
  else{
    cerr<<"INIT: Error. Got unidentified food influx location: "<<par.food_influx_location<<endl;
    exit(1);
  }
  
  //exit(1);
}

// we check beforeahand how many grid points we update and then radnomly generate coordinates.
// Typically, scaling_cell_to_ca_time>1, so foodinflux/scaling_cell_to_ca_time < 1/2,
// so we update less than half of the pixels, so this method is better
void IntPlane::IncreaseValEverywhere(void)
{
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
    int totalpixels= (sizex-1)*(sizey-1);
    double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
    for(int i=0; i<howmany_pixels_updated; i++){
      int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1[ 
      int posj = RandomNumber(sizey-2);
      if(sigma[posi][posj]<10) sigma[posi][posj]++;
    }  
//   }else{
//     for (int i=1; i<sizex-1; i++)
//       for (int j=1; j<sizey-1;j++){
//         if(RANDOM()<prob_food_influx && sigma[i][j]<10) //replace random() with proper function
//           sigma[i][j]++;
//       }
//   }
  
  
}

void IntPlane::NotIncreaseVal(void)
{
  return;
}

void IntPlane::IncreaseValIfEmpty(CellularPotts *cpm)
{
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1[
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10 && !cpm->Sigma(posi,posj)) sigma[posi][posj]++;
  }  
  
//     for (int i=1; i<sizex-1; i++)
//         for (int j=1; j<sizey-1;j++)
//         {
//             if(RANDOM()<par.foodinflux && sigma[i][j]<10 && !cpm->Sigma(i,j)) //replace random() with proper function
//                 sigma[i][j]++;
// //                 cerr<<"Hello food increased"<<endl;
//         }
}

void IntPlane::IncreaseValPatchy(CellularPotts *cpm)
{
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    for(int ii=-1;ii<=1;ii++)for(int jj=-1;jj<=1;jj++){
      int posii=posi+ii;
      int posjj=posj+jj;
      if(posii >= sizex-1 || posii < 1) continue;
      if(posjj >= sizey-1 || posjj < 1) continue;
      if(sigma[posi+ii][posj+jj]<10) sigma[posi+ii][posj+jj]=10;
    }
  }  
  //   }else{
  //     for (int i=1; i<sizex-1; i++)
  //       for (int j=1; j<sizey-1;j++){
  //         if(RANDOM()<prob_food_influx && sigma[i][j]<10) //replace random() with proper function
  //           sigma[i][j]++;
  //       }
  //   }
  
  
}

// food increases everwhere, but somewhere it's a lot, always there
// this is a 3x3 grid
void IntPlane::IncreaseValSomewhere(CellularPotts *cpm){

//   std::vector<int> meanposx = {sizex/8, sizex/4, 3*sizex/8, sizex/2, 5*sizex/8, 3*sizex/4, 7*sizex/8};
//   std::vector<int> meanposy = {sizey/8, sizey/4, 3*sizey/8, sizey/2, 5*sizey/8, 3*sizey/4, 7*sizey/8};
//   
  
  std::vector<int> meanposx = {sizex/4, sizex/2, 3*sizex/4};
  std::vector<int> meanposy = {sizey/4, sizey/2, 3*sizey/4};
  
//   std::vector<int> meanposx = {sizex/4,  3*sizex/4};
//   std::vector<int> meanposy = {sizey/4,  3*sizey/4};
  
//   std::vector<int> meanposx = {sizex/2};
//   std::vector<int> meanposy = {sizey/2};
  
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10) sigma[posi][posj]++;
  }  
  
  for(auto mpx: meanposx){
    for(int k=-10; k<10;k++){
      for(auto mpy: meanposy){
        for(int l=-10; l<10;l++){
          if( mpx+k >= sizex-1 || mpx+k < 1 || mpy+l >= sizey-1 || mpy+l < 1 ) continue;
          
//           if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<par.foodinflux/double(par.scaling_cell_to_ca_time) ) {
          if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<0.2/double(par.scaling_cell_to_ca_time) ) {
//             cerr<<"Updating: "<< meanposx[i]+k <<" "<< meanposy[i]+l <<endl;
            sigma[ mpx+k ][ mpy+l ]++;
          }
          
        }
      }
    }
    
  }
  
}

void IntPlane::IncreaseValPatchyRandom(CellularPotts *cpm){
  
  // standard trickling-in of food
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10) sigma[posi][posj]++;
  }  
  
  double food_patch_frequency = 0.025;
  
  if( RANDOM() < food_patch_frequency/double(par.scaling_cell_to_ca_time) ){
    int mpx = RandomNumber(sizex-2); //mean position x of patch
    int mpy = RandomNumber(sizey-2);
    for(int k=-100; k<100;k++){
      for(int l=-100; l<100;l++){
        int distsquare= k*k+l*l;
        if(distsquare > 100*100) continue; //makes a circle
        int px=mpx+k;
        int py=mpy+l;
        //boundary wrapping
        if( mpx+k >= sizex-1) px -= (sizex-1);
        else if(mpx+k < 1) px += (sizex-1);
        if (mpy+l >= sizey-1) py -= (sizey-1);
        else if( mpy+l < 1 ) py += (sizey-1);
        
        if( distsquare<25*25) 
          sigma[ px ][ py ] = 10;
        else if( distsquare< 50*50  && RANDOM()< 0.25 ) 
          sigma[ px ][ py ] = 10;
        else if( distsquare< 75*75  && RANDOM()< 0.5 ) 
          sigma[ px ][ py ] = 10;
        else if(RANDOM()< 0.25)
          sigma[ px ][ py ] =  10;
        
        // if( distsquare>75*75  && RANDOM()< 0.25 ) 
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>50*50  && RANDOM()< 0.5 ) 
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>25*25  && RANDOM()< 0.75 ) 
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // 
        // else if( distsquare>10*10  && RANDOM()< 0.75 ) 
        // sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else 
        //   sigma[ px ][ py ] =  10;
        
      }
    }
  }
  
}

//like above, but lower increase rate, and longer time
void IntPlane::IncreaseValPatchyRandomPersistence(CellularPotts *cpm){
  int food_update_persistence=10000; //every 1000 time steps food will update somewhere else
  static int time_since_update=0;
  // double food_patch_frequency = 0.02;
  double food_patch_prob=0.2;
  int food_radius =25;
  int food_radius_square =food_radius*food_radius;
  
  // static vector<int> v_time_since_update={0,0,0};
  // static vector<int> meanposx = {0,0,0};
  // static vector<int> meanposy = {0,0,0};
  // 
  
  // int mpx,mpy;
  static vector<vector<int>> v_t_since_update_mposxy(3, vector<int>(3,0) );
  
  
  // 
  // static int mpx=0,mpy=0;
  
  // standard trickling-in of food
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10) sigma[posi][posj]++;
  }  
  
  for(vector<vector<int>>::iterator it = v_t_since_update_mposxy.begin() ; it != v_t_since_update_mposxy.end(); ++it){
    if( (*it)[0] % food_update_persistence == 0){
      (*it)[0]= food_update_persistence/2 + RANDOM()*food_update_persistence/2; //5000 + random(0,5000) timesteps
      (*it)[1]= RandomNumber(sizex-2);
      (*it)[2]= RandomNumber(sizey-2);
      
    } 
  // }
  
  
  // if(time_since_update%food_update_persistence == 0){
  //   time_since_update=0;
  //   mpx = RandomNumber(sizex-2); //mean position x of patch
  //   mpy = RandomNumber(sizey-2); 
  // 
  // }
  
  
  // if( RANDOM() < food_patch_frequency/double(par.scaling_cell_to_ca_time) ){
    int mpx = (*it)[1];
    int mpy = (*it)[2];
    
    
    for(int k=-food_radius; k<food_radius;k++){
      for(int l=-food_radius; l<food_radius;l++){
        int distsquare= k*k+l*l;
        if(distsquare > food_radius*food_radius) continue; //makes a circle
        int px=mpx+k;
        int py=mpy+l;
        //boundary wrapping
        if( mpx+k >= sizex-1) px -= (sizex-1);
        else if(mpx+k < 1) px += (sizex-1);
        if (mpy+l >= sizey-1) py -= (sizey-1);
        else if( mpy+l < 1 ) py += (sizey-1);
        
        if( distsquare<food_radius_square/(4*4) && RANDOM()< food_patch_prob/double(par.scaling_cell_to_ca_time) ){ 
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}
        else if( distsquare< food_radius_square/(2*2)  && RANDOM()< 0.25*food_patch_prob/double(par.scaling_cell_to_ca_time) ) {
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}
        else if( distsquare< food_radius_square/((4./3.) * (4./3.) )  && RANDOM()< 0.5*food_patch_prob/double(par.scaling_cell_to_ca_time) ) {
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}
        else if(RANDOM()< 0.25*food_patch_prob/double(par.scaling_cell_to_ca_time)){
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}
        
        // if( distsquare>75*75  && RANDOM()< 0.25 ) 
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>50*50  && RANDOM()< 0.5 ) 
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>25*25  && RANDOM()< 0.75 ) 
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // 
        // else if( distsquare>10*10  && RANDOM()< 0.75 ) 
        // sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else 
        //   sigma[ px ][ py ] =  10;
        
      }
    }
    
    (*it)[0]++; // increase time in the counter
  }
  
  
  
}

//the idea is a population of bacteria like things that selfreplicate locally
void IntPlane::IncreaseValSelfGrowth(CellularPotts *cpm)
{
  double bacteria_birthrate=0.01;
  for(int i=1;i<sizex;i++)for(int j=1;j<sizey;j++){
    if(sigma[i][j]!=0){ 
      // there are so many bacteria, each with their birth rate, so:
      int new_borns = BinomialDeviate( sigma[i][j] , bacteria_birthrate/double(par.scaling_cell_to_ca_time) );
      if(new_borns==0) continue;
      
      while(new_borns>0){
        //where does it replicate? randomly in the 9 neighbourhood (including self)
        int xpos = -1 + (int)( 3.*RANDOM() ); // int number in [-1,1]
        int ypos = -1 + (int)( 3.*RANDOM() ); // int number in [-1,1]
        xpos+=i;
        ypos+=j;
        if(xpos>=sizex-1) xpos -= sizex-2;
        if(xpos<=0) xpos += sizex-2;
        if(ypos>=sizey-1) ypos -= sizey-2;
        if(ypos<=0) ypos += sizey-2;
        //std::cerr << "Hello3" << '\n';
        if(sigma[xpos][ypos]<2) sigma[xpos][ypos]++;
        //std::cerr << "Hello4" << '\n';
        
        new_borns--;
      }
    }
    
  }
  // return;
}

void IntPlane::IncreaseValSomewhereIfEmpty(CellularPotts *cpm){
  
  //   std::vector<int> meanposx = {sizex/8, sizex/4, 3*sizex/8, sizex/2, 5*sizex/8, 3*sizex/4, 7*sizex/8};
  //   std::vector<int> meanposy = {sizey/8, sizey/4, 3*sizey/8, sizey/2, 5*sizey/8, 3*sizey/4, 7*sizey/8};
  //   
  
  std::vector<int> meanposx = {sizex/4, sizex/2, 3*sizex/4};
  std::vector<int> meanposy = {sizey/4, sizey/2, 3*sizey/4};
  
  //   std::vector<int> meanposx = {sizex/4,  3*sizex/4};
  //   std::vector<int> meanposy = {sizey/4,  3*sizey/4};
  
  //   std::vector<int> meanposx = {sizex/2};
  //   std::vector<int> meanposy = {sizey/2};
  
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10 && !cpm->Sigma(posi,posj)) sigma[posi][posj]++;
  }  
  
  for(auto mpx: meanposx){
    for(int k=-10; k<10;k++){
      for(auto mpy: meanposy){
        for(int l=-10; l<10;l++){
          if( mpx+k >= sizex-1 || mpx+k < 1 || mpy+l >= sizey-1 || mpy+l < 1 ) continue;
          
          //           if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<par.foodinflux/double(par.scaling_cell_to_ca_time) ) {
          if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<0.2/double(par.scaling_cell_to_ca_time) && !cpm->Sigma(mpx+k,mpy+l)) {
            //             cerr<<"Updating: "<< meanposx[i]+k <<" "<< meanposy[i]+l <<endl;
            sigma[ mpx+k ][ mpy+l ]++;
          }
          
        }
      }
    }
    
  }
  
}



