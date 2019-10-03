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
// #include <iostream>
// #include <fstream>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <sstream>
#include <errno.h>
#include <math.h>
#include "dish.h"
#include "sticky.h"
#include "parameter.h"
#include "info.h"
#include "crash.h"
#include "pde.h"
#include "intplane.h"

#define EXTERNAL_OFF

extern Parameter par;

using namespace std;

Dish::Dish(void) {

  ConstructorBody();

  CPM=new CellularPotts(&cell, par.sizex, par.sizey);
  Food= new IntPlane(par.sizex, par.sizey);

  if (par.n_chem)
    PDEfield=new PDE(par.n_chem,par.sizex, par.sizey);

  cout<<"Starting the dish. Initialising..."<<endl;
  // Initial cell distribution is defined by user in INIT {} block
  Init();



  //cout<<cell[1].neighbours[0].first<<endl;
  //cout<<cell[1].TargetArea()<<endl;
}



Dish::~Dish() {
  cell.clear();

  delete CPM;
  delete Food;
}

void Dish::InitTargetArea(void)
{
  if (par.target_area>0)
    for(std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) {
      c->SetTargetArea(par.target_area);
    }
}

void Dish::InitKeyLock(void)
{
  size_t key_lock_length=(size_t)par.key_lock_length; //take size of key and lock
  for(auto kl_pair: par.keylock_list){
    if( kl_pair.tau == 0 ) continue;
    if(kl_pair.key.size() != key_lock_length || kl_pair.lock.size() != key_lock_length){
      cerr<<"Dish::InitKeyLock(): error. Initial key and lock vectors are not of size par.key_lock_length = "<<par.key_lock_length<<endl;
      exit(1);
    }
  }
  //cout << "here 2"<<endl;
  /*

  if(  init_keyprey.size()     != key_lock_length ||
       init_lockprey.size()    != key_lock_length ||
       init_keypredator.size() != key_lock_length ||
       init_lockpredator.size()!= key_lock_length
    ){
    cerr<<"Dish::InitKeyLock(): error. Initial key and lock vectors are not of size par.key_lock_length = "<<par.key_lock_length<<endl;
    exit(1);
  }*/

  //cerr<<"Checking that prey and predator works: "<<PREY<<" "<<PREDATOR<<endl;
  //exit(1);
//   for(auto c:cell)
//     cerr<<c.tau<<endl;
//   exit(1);
//
  //cout << "here 1"<<endl;
  vector<Cell>::iterator c;
  for(c=cell.begin(), ++c; c!=cell.end(); ++c){
    int current_tau=c->getTau();
    if(current_tau == PREY){
      //cerr<<"Hello, got prey"<<endl;
      //cerr<<"Now: "<<PREY<<endl;

      c->setJkey( par.keylock_list[ PREY ].key );
      c->setJlock( par.keylock_list[ PREY ].lock );

      //cerr<<"Check "<< par.keylock_list[ PREY ].key[0] <<endl;


    }else if(current_tau == PREDATOR){
      //cerr<<"Hello, got predator"<<endl;
      //cerr<<"Now: "<<PREDATOR<<endl;

      c->setJkey( par.keylock_list[ PREDATOR ].key );
      c->setJlock( par.keylock_list[ PREDATOR ].lock );
      //cerr<<"Check "<< par.keylock_list[ PREDATOR ].key[0] <<endl;

    }else{
      cerr<<"Dish::InitKeyLock(): error. Are there more than two types? got tau = " << c->getTau()<<endl;
      exit(1);
    }
//     for(auto bla: c->getJkey())
//       cerr<<bla<<" ";
//     cerr<<endl;
//     for(auto bla: c->getJlock())
//       cerr<<bla<<" ";
//     cerr<<endl;
//
  }

  //exit(1);
}

// In this version we can let user define the lookup table for J with medium
//notice we use Paulien's method of summing powers
// Also, is this half J value? think so... - YES
int Dish::CalculateJwithMedium( vector<int> key )
{
  int keypos_formedium = par.Jmedr.keypos_formedium;
  vector <int> lookup_table = par.Jmedr.lookup_table;
  int offset = par.Jmedr.offset;
  int Jval= 0;


  // Per renske suggestion, I will try a different range, more contained
  // if I just add 10 to the range (0->15)
  //I get range (10->25)
  for(int i=0;i<keypos_formedium;i++)
    //Jval += key[i]*pow(2.0,keypos_formedium-i-1); //so that zeroth bit is most significant
    Jval += key[i]*lookup_table[i];
  Jval += offset;//so that interaction with medium can't be 0

  return (int)Jval;
}

int Dish::CalculateJfromKeyLock( vector<int> key1, vector<int> lock1, vector<int> key2, vector<int> lock2 )
{
  int score=0;

  for(int i=0;i<par.key_lock_length;i++){
    score += ( key1[i] != lock2[i] )?1:0;
    score += ( key2[i] != lock1[i] )?1:0;
  }
  //now perfect score is 20, make... a sigmoid response?
  // with 20 you should get very low J val (high adhesion)
  // with 0 you should get high J val (low adh)
  //This is a arbitrary function 3+40*exp(-0.01*x^2)
  //add 0.5 before truncation to int makes it apprx to closest integer

  int Jfromkeylock =(int)( 52. - 48. * ((double)score) /(2.*par.key_lock_length) );
  // int Jfromkeylock = 3 + (int)(0.5+ 40.*exp( -pow( (score/double(par.key_lock_length)) , 2.) ));


  /*
  cout<<"CalculateJfromKeyLock: I got:"<<endl<< "key1: ";
  for (auto i: key1)
    cout << i << " ";
  cout<<"lok1: ";
  for (auto i: lock1)
    cout << i << " ";
  cout<<endl<<"lok2: ";
  for (auto i: lock2)
    cout << i << " ";
  cout<<"key2: ";
  for (auto i: key2)
    cout << i << " ";
  cout<<endl<<"score = "<<score<<" Jval = "<<Jfromkeylock<<endl;
  */

  return Jfromkeylock;
}

void Dish::InitVectorJ(void) //Initialise vector of J values for each cell
{
  std::vector<Cell>::iterator ci;
  std::vector<Cell>::iterator cj;

  //cerr<<"InitVectorJ begin"<<endl;

  int thisval;
  for(ci=cell.begin(); ci!=cell.end(); ++ci){
    for(cj=ci, ++cj; cj!=cell.end(); ++cj){
      if(ci->Sigma() == MEDIUM){
        if(cj->Sigma() == MEDIUM) thisval=0;
        else thisval = CalculateJwithMedium( cj->getJkey() );  //J with medium is calculated from part of cell's
        ci->setVJ_singleval(cj->Sigma(), thisval);
        cj->setVJ_singleval( MEDIUM , thisval); //this cj.sigma=0, we update its J values as well

      }else{
        thisval = CalculateJfromKeyLock( ci->getJkey(), ci->getJlock(), cj->getJkey(), cj->getJlock() );
        ci->setVJ_singleval(cj->Sigma(), thisval);
        cj->setVJ_singleval(ci->Sigma(), thisval);

//         if(thisval<=1){
//           cerr<<"ci="<<ci->Sigma()<<", cj->sigma="<<cj->Sigma()<<", J=" <<thisval<<endl;
//           cerr<<"keylock"<<endl;
//           for (auto i = ci->getJkey().begin(); i != ci->getJkey().end(); ++i)
//             std::cout << *i ;
//           std::cout << ' ';
//           for (auto i = ci->getJlock().begin(); i != ci->getJlock().end(); ++i)
//             std::cout << *i ;
//           std::cout<<endl;
//           for (auto i = cj->getJkey().begin(); i != cj->getJkey().end(); ++i)
//             std::cout << *i ;
//           std::cout << ' ';
//           for (auto i = cj->getJlock().begin(); i != cj->getJlock().end(); ++i)
//             std::cout << *i ;
//           std::cout<<endl;
//
//           exit(1);
//         }

      }
    }
  }

  cerr<<"InitVectorJ done"<<endl;
//   for(auto c: cell){
//     cerr<<"Cell "<<c.Sigma()<<", tau: "<<c.getTau()<<": " ;
//     for(auto jval: c.getVJ())
//       cerr<<jval<<" ";
//     cerr<<endl;
//   }

}

//This function and the one above could be merged
// sigma_to_update is a vector of int, there are a lot of zeros,
// but the numbers are the new sigmas (see DivideCells in ca.cpp)
// Note: some optimisation is possible here because we are updating twice the new cells (all the upd_sigma)
void Dish::UpdateVectorJ(vector<int> sigma_to_update)
{
  //cerr<<"UpdateVectorJ begin, cell vector size: "<<cell.size()<<endl;

  vector<Cell>::iterator c;

  for(auto upd_sigma: sigma_to_update){
    if(upd_sigma != 0){
      for(c=cell.begin(); c!=cell.end(); ++c){
        if( 0 == c->Area() && 0 == c->AliveP() ) continue; //used to be if !c->AliveP(), but some cells are dead but not disappeared yet.
        if(c->Sigma() != MEDIUM){
          if(c->Sigma()==upd_sigma) continue;
          //update for each cell, their interactions with cell at position sigma to update,
          int jval = CalculateJfromKeyLock( c->getJkey(), c->getJlock(), cell[upd_sigma].getJkey(), cell[upd_sigma].getJlock() );  //k1,l1,k2,l2
//           if(jval<=1){
//             cerr<<"upd_sigma="<<upd_sigma<<", c->sigma="<<c->Sigma()<<", J=" <<jval<<endl;
//             cerr<<"keylock"<<endl;
//             for (auto i = cell[upd_sigma].getJkey().begin(); i != cell[upd_sigma].getJkey().end(); ++i)
//               std::cout << *i ;
//             std::cout << ' ';
//             for (auto i = cell[upd_sigma].getJlock().begin(); i != cell[upd_sigma].getJlock().end(); ++i)
//               std::cout << *i ;
//             std::cout<<endl;
//             for (auto i = c->getJkey().begin(); i != c->getJkey().end(); ++i)
//               std::cout << *i ;
//             std::cout << ' ';
//             for (auto i = c->getJlock().begin(); i != c->getJlock().end(); ++i)
//               std::cout << *i ;
//             std::cout<<endl;
//
//             exit(1);
//           }
          c->setVJ_singleval(upd_sigma, jval);
          //the same number can be used by cell at which sigma is being updated
          cell[upd_sigma].setVJ_singleval( c->Sigma() , jval);
        }else{
          int jval = CalculateJwithMedium( cell[upd_sigma].getJkey() ); // needs only key

          c->setVJ_singleval( upd_sigma , jval );
          cell[upd_sigma].setVJ_singleval( c->Sigma() , jval ); //update nrg of medium with cell c

        }
      }
    }
  }
  //cerr<<"UpdateVectorJ ends"<<endl;

//   for(auto c: cell){
//     cerr<<"Cell "<<c.Sigma()<<", tau: "<<c.getTau()<<": " ;
//     for(auto jval: c.getVJ())
//       cerr<<jval<<" ";
//     cerr<<endl;
//   }
//   exit(1);
}

//Initialise init_maintenance_fraction from parameters.h
void Dish::InitMaintenanceFraction()
{
  vector<Cell>::iterator c;
  for(c=cell.begin(); c!=cell.end(); ++c){
    c->maintenance_fraction = par.init_maintenance_fraction;
  }
}

// sigma_newcells is an int vector as lognas there are cells,
// it is zero everywhere, except at the positions of a mother cell's sigma,
// where it contains as value the sigma of the daughter cell
void Dish::MutateCells(vector<int> sigma_to_update)
{
  for(auto upd_sigma: sigma_to_update){
    if(upd_sigma != 0){
      cell[upd_sigma].MutateKeyAndLock();
      //cell[upd_sigma].MutateMu();
      //cerr<<"hello from before mutation"<<endl;
      if(par.evolreg == true){
        //cell[upd_sigma].MutateMaintenanceFractionParameters();
        cell[upd_sigma].MutateExtProtFractionParameters();
        //cell[upd_sigma].MutateChemotaxisParameters();
      }
    }
  }
}

// Notice that at this stage cells are completely isolated,
// thus completely surrounded by medium
// Also, they should be far enough from borders,
//nevertheless, we are checking that
void Dish::InitContactLength()
{
   int k;
   int celltype, celltypeneigh;
   int sigma, sigmaneigh;
   std::vector<Cell>::iterator c;

   for(c=cell.begin(); c!=cell.end(); ++c)
   {
     c->clearNeighbours();
   }

   for(int x=1; x<par.sizex-1; x++){
     for(int y=1; y<par.sizey-1; y++){
       sigma=CPM->Sigma(x,y); //focus is on a cell - provided it is not medium
       if( sigma ){
         for(k=1; k<=CPM->n_nb; k++)//go through neighbourhood of the pixel
         {

           int neix = x + CPM->nx[k];
           int neiy = y + CPM->ny[k];
           if(neix<=0 || neix>=par.sizex-1 || neiy<=0 || neiy>=par.sizey-1){
             cerr<<"InitContactLength(): warning. Neighbour is beyond borders"<<endl;
             if( par.periodic_boundaries ){
               cerr<<"Wrapped boundary condition applies"<<endl;
               if(neix<=0) neix=par.sizex-2+neix;
               if(neix>=par.sizex-1) neix=neix-par.sizex+2;
               if(neiy<=0) neiy=par.sizey-2+neiy;
               if(neiy>=par.sizey-1) neiy=neiy-par.sizey+2;
             }else{
               cerr<<"Fixed boundary condition applies: neighbour contact discarded."<<endl;
               continue;
             }
           }

           sigmaneigh = CPM->Sigma( neix , neiy );
           if( sigmaneigh != sigma  )//medium can also be a neighbour!
           {
             cell[sigma].updateNeighbourBoundary(sigmaneigh,1);
             //cout<<"We updated cell "<<sigma<<" to have boundary with "<<sigmaneigh<<endl;
           }
         }
       }
     }
   }
   //PrintContactList();
   // PRINTING INITIALISED CONTACTS - THIS SHOULD PRINT ONLY MEDIUM -- TRUE
   //   cout<<"cell 1 has "<<cell[1].neighbours[0].first<<" contacts with cell 0"<<endl;
}

 //call this function after every MCS to update the contact duration between neighbours and to erase neighbours
 //with whom contact has been lost
void Dish::UpdateNeighDuration()
{
   std::vector<Cell>::iterator c;
   std::map<int, pair<int,int> >::iterator n, prev;

//    cerr<<"Printing neighbor list"<<endl;
//    for(c=cell.begin(),c++; c!=cell.end(); ++c){
//      cerr << c->Sigma() <<": ";
//      n=c->neighbours.begin();
//      while( n != c->neighbours.end() ){
//        if(n->second.first!=0 && n->second.second!=0)
//          cerr << n->first <<", ";
//        n++;
//     }
//     cerr << endl;
//    }

   //cerr<<"Hello UpdateNeighDuration begin"<<endl;
   for(c=cell.begin(),c++; c!=cell.end(); ++c)
   {
     if(! c->AliveP() ) continue;
     n=c->neighbours.begin();
     //cerr<<"Hello UpdateNeighDuration 0"<<endl;
     while( n != c->neighbours.end() ){
       if(n->first == -1 ){
         cerr << "We got a cell that has boundary as neighbour, n=-1:"<<endl;
         exit(1);
       }
       if(n->second.first==0){
        prev=n;
        n++;
        //cerr<<"Hello UpdateNeighDuration 1"<<endl;
        c->neighbours.erase(prev);
        //cerr<<"Hello UpdateNeighDuration 2"<<endl;
        //cout<<"erasing contact of cell "<<c->Sigma()<<" with cell "<<n->first<<endl;
       }
       else{
         n->second.second++;
         n++;
       }

     }

   }
   //cerr<<"Hello UpdateNeighDuration end"<<endl;
}


//Prints how much medium is neigh of a given cell (which is sigma of the cell you want)
void Dish::PrintReality(int which)
{
  int counter=0;
  bool crossborder=false;
  if(!which){
    cerr<<"PrintReality(): Error use me properly"<<endl;
    exit(1);
  }
  cerr<<"Printing reality"<<endl;
  for(int x=1; x<par.sizex-1; x++){
    for(int y=1; y<par.sizey-1; y++){
      int sigma=CPM->Sigma(x,y); //focus is on a cell - provided it is not medium
      if( sigma == which ){
        for(int k=1; k<=CPM->n_nb; k++)//go through neighbourhood of the pixel
        {

          int neix = x + CPM->nx[k];
          int neiy = y + CPM->ny[k];
          if(neix<=0 || neix>=par.sizex-1 || neiy<=0 || neiy>=par.sizey-1){
            crossborder=true;
            //cerr<<"InitContactLength(): warning. Neighbour is beyond borders"<<endl;
            if( par.periodic_boundaries ){
              //cerr<<"Wrapped boundary condition applies"<<endl;
              if(neix<=0) neix=par.sizex-2+neix;
              if(neix>=par.sizex-1) neix=neix-par.sizex+2;
              if(neiy<=0) neiy=par.sizey-2+neiy;
              if(neiy>=par.sizey-1) neiy=neiy-par.sizey+2;
            }
          }

          int sigmaneigh =  CPM->Sigma(neix,neiy);
          if(sigmaneigh == MEDIUM)
            counter++;
        }
      }
    }
  }
  cerr<< "Reality, medium in neighbourhood = "<<counter<<endl;
  if(crossborder) cerr<<"cross border"<<endl;
  else cerr<<"NEVER cross border"<<endl;
}

void Dish::PrintContactList(int which)
{
  int flag=0;
  std::vector<Cell>::iterator c;
  std::map<int, pair<int,int> >::iterator n, prev;

  if(which!=-1){
    cerr << "Cell id: " << cell[which].Sigma() << ", contact with:" <<endl;
    n=cell[which].neighbours.begin();
    while( n != cell[which].neighbours.end() ){
      //if(n->second.first!=0 && n->second.second!=0)
      cerr << "\t" << n->first <<", contact length: " << n->second.first << ", duration: " << n->second.second << endl;
      if(n->first == -1 ){
        cerr << "We got a cell that has boundary as neighbour, n=-1:"<<endl;
        flag=1;
      }
      n++;
    }

  }
  else{
  cerr<<"Printing contact list - BEGIN"<<endl;
  for(c=cell.begin(); c!=cell.end(); ++c){
    if(! c->AliveP() ) continue;

    cerr << "Cell id: " << c->Sigma() << ", contact with:" <<endl;
    n=c->neighbours.begin();
    while( n != c->neighbours.end() ){
      //if(n->second.first!=0 && n->second.second!=0)
      cerr << "\t" << n->first <<", contact length: " << n->second.first << ", duration: " << n->second.second << endl;
      if(n->first == -1 ){
        cerr << "We got a cell that has boundary as neighbour, n=-1:"<<endl;
        flag=1;
      }
      n++;
    }
    //cerr << endl;
  }
  if(flag==1){
    cerr<<"Someone has -1 as neighbour"<<endl;
    exit(1);
  }
  cerr<<"Printing contact list - END"<<endl;
  }
}

void Dish::PrintCellParticles(void)
{
  for(auto c:cell){
    cerr<<" sigma: "<<c.sigma<<", particles: "<<c.particles<<endl;
  }
}

void Dish::FoodPlot(Graphics *g)
{
  // cpm->sigma[x][y] returns sigma, which I can use to indicise the vector of cells... can I?

  // suspend=true suspends calling of DrawScene
  for(int x=1;x<par.sizex-1;x++)
    for(int y=1;y<par.sizey-1;y++)

  //for(int x=0;x<par.sizex;x++)
  //  for(int y=0;y<par.sizey;y++)
      //if (cpm->Sigma(x,y)==0) {
      if(Food->Sigma(x,y) != 0){
        if(CPM->Sigma(x,y)==0){
          // Make the pixel four times as large
          // to fit with the CPM plane
          g->Point(10+Food->Sigma(x,y),2*x,2*y);
          g->Point(10+Food->Sigma(x,y),2*x+1,2*y);
          g->Point(10+Food->Sigma(x,y),2*x,2*y+1);
          g->Point(10+Food->Sigma(x,y),2*x+1,2*y+1);
        }else{
          ;
          // it's getting a bit cumbersome to look at this, for now I'll do without
          // g->Point(60+Food->Sigma(x,y),2*x,2*y);
          // g->Point(60+Food->Sigma(x,y),2*x+1,2*y);
          // g->Point(60+Food->Sigma(x,y),2*x,2*y+1);
          // g->Point(60+Food->Sigma(x,y),2*x+1,2*y+1);
        }
      }
}

void Dish::Plot(Graphics *g, int colour) {
    if (CPM){
      CPM->Plot(g, colour);
    }

    // return;


    //here food plotting, with info from cpm and cell
    FoodPlot(g);
    // return;
    //Plot direction arrows, with line function from X11?
    if(par.startmu>0){
      for(auto c: cell){
        if(c.sigma==0 or !c.alive) continue;
        int x1=2*c.meanx;
        int y1=2*c.meany;
        int x2= 2*(c.meanx+5*c.tvecx) ;
        if(x2>=2*par.sizex) x2=2*par.sizex; //if too large or too small, truncate it
        else if(x2<0) x2=0;
        int y2=2*(c.meany+5*c.tvecy);
        if(y2>=2*par.sizey) y2=2*par.sizey;
        else if(y2<0) y2=0;
        //now we have to wrap this
        // do we really? we could just truncate vectors up to the max size..
        g->Line(x1,y1,x2, y2, 1); //notice that Line just calls Point for drawing,
                                  // so it does not intrinsically suffer from x,y inversion
      }
    }
    // return;

    //get info where the peak is and draw a line for box where who_made_it should register stuff
    int peakx = Food->GetPeakx();
    int peaky = Food->GetPeaky();
    //std::cerr << "peak is at "<<peakx<<" "<<peaky << '\n';
    int minx,maxx,miny,maxy;
    // notice that the box is now radial
    if(strcmp(par.food_influx_location,"specified_experiment") == 0){


      if(peakx==1) {
        //then peaky = sizey/2 and
        minx = 1;
        maxx = par.the_line;
        miny = par.sizey/2 - par.the_line; //circle
        maxy = par.sizey/2 + par.the_line;
        // miny = 1;
        // maxy = par.sizey-1;
        // g->Line(2*par.the_line, 1 , 2*maxx , 2*par.sizey-1 , 1);
      }else if(peakx==par.sizex-1){
        minx = par.sizex-par.the_line;
        maxx = par.sizex-1;
        miny = par.sizey/2 - par.the_line; //circle
        maxy = par.sizey/2 + par.the_line;
        // g->Line(2*minx,1,2*minx, 2*par.sizey-1, 1);
      }else if(peaky==1){
        minx = par.sizex/2 - par.the_line; //circle
        maxx = par.sizex/2 + par.the_line;
        miny = 1;
        maxy = par.the_line;
        // g->Line(1,2*par.the_line,2*par.sizex-1, 2*par.the_line, 1);
      }else if(peaky==par.sizey-1){
        minx = par.sizex/2 - par.the_line; //circle
        maxx = par.sizex/2 + par.the_line;
        miny = par.sizey-par.the_line;
        maxy = par.sizey-1;
        // g->Line(1, 2*(par.sizey-par.the_line) , 2*par.sizex-1  , 2*(par.sizey-par.the_line), 1);
      }else{
        cerr<<"Plot(): Error. Got weird peakx and peaky position: peakx, peaky = "<<peakx<<", "<<peaky<<endl;
        std::cerr << "Don't know what to do with this, program exits now." << '\n';
        exit(1);
      }

      // return;
      // Draw the_line, but only if it's inside the field
      int dx,dy;
      double dist;
      for(int i=minx-1; i<=maxx+1;++i) for(int j=miny-1; j<=maxy+1;++j){
        if(i>=par.sizex -1 || i<=1 || j>=par.sizey-1 || j<=1) continue; //don't draw on the borders, or beyond them
        dx = peakx - i;
        dy = peaky - j;
        dist = hypot(dx,dy);
        if(par.the_line - 0.5 < dist && dist < par.the_line+0.5) g->Point(1,2*i,2*j);
      }
    }

    else if (strcmp(par.food_influx_location,"boundarygradient") == 0){

      if(peakx==1) {
        //then peaky = sizey/2 and
        minx = 1;
        maxx = par.the_line;
        miny = 1;
        maxy = par.sizey-1;
      }else if(peakx==par.sizex-1){
        minx = par.sizex-par.the_line;
        maxx = par.sizex-1;
        miny = 1;
        maxy = par.sizey-1;
      }else if(peaky==1){
        minx = 1;
        maxx = par.sizex-1;
        miny = 1;
        maxy = par.the_line;
      }else if(peaky==par.sizey-1){
        minx = 1;
        maxx = par.sizex-1;
        miny = par.sizey-par.the_line;
        maxy = par.sizey-1;
        // g->Line(1, 2*(par.sizey-par.the_line) , 2*par.sizex-1  , 2*(par.sizey-par.the_line), 1);
      }else{
        cerr<<"Plot(): Error. Got weird peakx and peaky position: peakx, peaky = "<<peakx<<", "<<peaky<<endl;
        std::cerr << "Don't know what to do with this, program exits now." << '\n';
        exit(1);
      }
      // return;
      // Draw the_line, but only if it's inside the field
      int dx, dy;
      for(int i=minx-0; i<=maxx+1;++i) for(int j=miny-1; j<=maxy+1;++j){
        if(i>=par.sizex -1 || i<=1 || j>=par.sizey-1 || j<=1) continue; //don't draw on the borders, or beyond them
        //dx = peakx - i;
        //dy = peaky - j;

        if(peakx==1 || peakx==par.sizex-1){
          if(par.the_line - 0.5 < i && i < par.the_line+0.5){
            g->Point(1,2*i,2*j);
            //cerr<<"yes!";
          }
        }
        else if(peaky==1 || peaky==par.sizex-1){
          if(par.the_line - 0.5 < j && j < par.the_line+0.5){
            g->Point(1,2*i,2*j);
            //cerr<<"yes2!";
          }
        }
      }

    }
    else{
      return;
    }



 }


void Dish::ConstructorBody() {

  Cell::maxsigma=0;

  // Allocate the first "cell": this is the medium (tau=0)
  cell.push_back(*(new Cell(*this,0)));

  // indicate that the first cell is the medium
  cell.front().sigma=0;
  cell.front().tau=0;

  CPM=0;
  PDEfield=0;
  Food=0; // SORT OF NULL INITIALISATION
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


// DEPRECATED
// Based on code by Paulien Hogeweg.
void Dish::CellGrowthAndDivision(void) {

  cerr<<"Do not use this function Dish::CellGrowthAndDivision(), it is outdated"<<endl<<"The program will terminate now"<<endl;
  exit(1);

  vector<bool> which_cells(cell.size());

  static int mem_area=0;

  // if called for the first time: calculate mem_area
  if (!mem_area) {
    mem_area=TargetArea()/CountCells();
  }

  int cell_division=0;

  vector<Cell>::iterator c;
  for ( (c=cell.begin(), c++);
	c!=cell.end();
	c++) {

    if ( (c->Area()-c->TargetArea())>c->GrowthThreshold() ) {
      c->IncrementTargetArea();

    }

    if ( (c->TargetArea() > 2 * mem_area ) ) {
      which_cells[c->Sigma()]=true;
      cell_division++;
    }
  }

  // Divide scheduled cells
  if (cell_division) {
    CPM->DivideCells(which_cells);
  }

}


// This function does not change contacts, only target areas!
void Dish::Predate(void)
{
  //vector<int> scrambled_cells(cell.size()); //order of cells eating what from whom - this has to be scrambled
  vector<int> list_cell_id; //order of cells eating what from whom - this has to be scrambled -
                                               // this is also not ok, because some cells have died but they have not been removed from cell vector????
  //go through pop size and see who is alive, and add their ids to list_cell_id
  vector<Cell>::iterator c;

  for (c=cell.begin(), c++ ; c!=cell.end() ; c++){
    if( c->AliveP() && c->getTau()==PREDATOR) list_cell_id.push_back( c->Sigma() );
  }

  int i=1;

  //Initialise array
  vector<int>::iterator it;

  // first shuffle
  int sizeOfVector = list_cell_id.size();
  for (int k = sizeOfVector-1; k > 0; k--) {
    int r = (int)( (k+1)*RANDOM() );  //r from 0 to k included
    swap(list_cell_id[k], list_cell_id[r]);
  }
  //now for each PREDATOR in random order we go take its contact
  for(it=list_cell_id.begin() ; it!=list_cell_id.end() ; it++) {
    if(cell[*it].getTau() != PREDATOR) {
      cerr<<"Predate(): warning. Got sigma of prey in list of predators list_cell_id"<<endl;
      continue;    // Tau, NOT sigma !!! this is just a check,
    }

    //print its neighbours
    map<int, pair<int,int> >::iterator m; //iterates over neigh of a cell
    //iterates over neighbours
    //THIS WILL HAVE TO BE SCRAMBLED... maybe... maybe not - it eats from eveybody
    for(m=cell[*it].neighbours.begin(); m!=cell[*it].neighbours.end(); m++ ){
      //cerr << m->first << ':' << (m->second).first << ',' << (m->second).second << endl;
      int idnei = m->first;
      int contlength = (m->second).first;
      int contduration = (m->second).second;

      if(idnei == 0 || cell[idnei].getTau() == PREDATOR){
        //cout << "Predate(): Warning. A neighbour has id = 0" << endl;
        continue;
      }

      //*******************************************************************
      // In this model predators eat fast, all the cell in one time step
      //

      // check that the prey target area >0 and cell is alive, so that it has not already been eaten by another predator
      if(contduration >= par.min_contact_duration_for_preying && cell[idnei].TargetArea() && cell[idnei].AliveP()) {
        int shrinkage_prey = cell[idnei].Area(); //predator eats the whole cell
        //PREY SHRINKS
        cell[idnei].SetTargetArea(0); //so
        cell[idnei].Apoptose();
        //PREDATOR GROWS
        cell[*it].particles += par.metabolic_conversion * shrinkage_prey;
        //*******************************************************************

      }


//         //*******************************************************************
//         // In this model predators eat slowly, proportional to
//         //contact length with prey
//         //
//         //if(contduration >= par.min_contact_duration_for_preying) {
//
//         //this makes sure that preys can be eaten down to the very last pixel
//         double dshrinkage_prey = par.frac_contlen_eaten*contlength;
//         int shrinkage_prey = par.frac_contlen_eaten*contlength;
//         shrinkage_prey += (  RANDOM() <  dshrinkage_prey - shrinkage_prey )? 1 : 0;
//
//         //PREY SHRINKS
//         //this could make predators eat more than there is prey, no.
//         int new_area_nei = cell[idnei].TargetArea() - shrinkage_prey;
//         if(new_area_nei > 0) cell[idnei].SetTargetArea(new_area_nei);
//         else cell[idnei].SetTargetArea(0);
//         //PREDATOR GROWS
//         cell[*it].particles += par.metabolic_conversion * shrinkage_prey;
//         //*******************************************************************
//
//      }
    }
  }
}

void Dish::CellsEat(void)
{

  int MAX_PARTICLES=1000000; //max particles that a cell can have inside it

  int foodload;
  for(int x=1; x<par.sizex-1;x++){
      for(int y=1; y<par.sizey-1;y++){
        if(CPM->Sigma(x,y) && cell[CPM->Sigma(x,y)].AliveP() && cell[CPM->Sigma(x,y)].getTau() == PREY && Food->Sigma(x,y)){
          int howmuchfood = BinomialDeviate( Food->Sigma(x,y) , cell[CPM->Sigma(x,y)].GetEatProb()/(double)par.scaling_cell_to_ca_time );
          Food->addtoValue(x,y,-howmuchfood);

          // we cannot add all the particles endlessly, otherwise it overflows
          int current_particles = cell[CPM->Sigma(x,y)].particles;
          int added_particles = ( current_particles <= (MAX_PARTICLES-howmuchfood) )?howmuchfood:(MAX_PARTICLES-current_particles);
          cell[CPM->Sigma(x,y)].particles += added_particles;

        }
      }
    }
}

//changes cells direction vector based on where more food is
void Dish::CellsEat2(void)
{
  // if(par.periodic_boundaries)
  // {
  //   std::cerr << "CellsEat2: does not work with wrapped boundaries" << '\n';
  //   exit(1);
  // }

  int MAX_PARTICLES=1000000; //max particles that a cell can have inside it
  std::vector<int> fsumx(cell.size(),0), fsumy(cell.size(),0),ftotal(cell.size(),0);
  int foodload;
  for(int x=1; x<par.sizex-1;x++){
      for(int y=1; y<par.sizey-1;y++){
        if(CPM->Sigma(x,y) && cell[CPM->Sigma(x,y)].AliveP() && cell[CPM->Sigma(x,y)].getTau() == PREY && Food->Sigma(x,y)){
          int cell_sigma=CPM->Sigma(x,y);
          if(Food->Sigma(x,y) > 0){
            //determine the mean position of the food that the cell sees
            int fx=x, fy=y;

            if(par.periodic_boundaries){
              double meanx = cell[cell_sigma].getXpos();
              double meany = cell[cell_sigma].getYpos();
              if( (fx-meanx)>0 && (fx-meanx)>(meanx-(fx-(par.sizex-2))) ) fx-=(par.sizex-2);
              else if( (meanx-fx)>0 && (meanx-fx)>(fx+(par.sizex-2)-meanx)) fx+=(par.sizex-2);
              if( (fy-meany)>0 && (fy-meany)> (meany - (fy - (par.sizey-2))) ) fy-=(par.sizey-2);
              else if( (meany-fy>0) && (meany-fy)>(fy+(par.sizey-2)-meany)) fy+=(par.sizey-2);
            }

            fsumx[cell_sigma]+=fx*Food->Sigma(x,y);
            fsumy[cell_sigma]+=fy*Food->Sigma(x,y);
            ftotal[cell_sigma]+=Food->Sigma(x,y);
          }else{
            int howmuchfood = BinomialDeviate( -1*Food->Sigma(x,y) , cell[cell_sigma].GetEatProb()/(double)par.scaling_cell_to_ca_time );
            Food->addtoValue(x,y,-1*-howmuchfood);

            // we cannot add all the particles endlessly, otherwise it overflows
            int current_particles = cell[cell_sigma].particles;
            int added_particles = ( current_particles <= (MAX_PARTICLES-howmuchfood) )?howmuchfood:(MAX_PARTICLES-current_particles);
            // std::cerr << "cell eats so many particles: "<< added_particles;
            cell[cell_sigma].particles += added_particles;
            // std::cerr << " therefore it now has: "<< cell[cell_sigma].particles << '\n';
            // ADD SIGNAL for regulation
          }

        }
      }
    }

    //update the cell's movement vector with respect to the location of food
    int count=0;

    double hyphyp;
    double xv,yv;
    for(auto &c: cell){
      if(c.sigma && ftotal[c.sigma]){
        //calculate "food" vector with respect to cell mean pos
        xv=fsumx[c.sigma]/(double)ftotal[c.sigma]-c.meanx;
        yv=fsumy[c.sigma]/(double)ftotal[c.sigma]-c.meany;

        hyphyp=hypot(xv,yv);

        // in a homogeneous medium, gradient is zero
        // we then pick a random direction
        //if(hyphyp > 0.0001){
        //  xv/=hyphyp;
        //  yv/=hyphyp;
        c.setChemVec(xv/hyphyp,yv/hyphyp);

        c.setChemVecLength(hyphyp);
        //}else{
        //  double theta = 2.*M_PI*RANDOM();
        //  c.setChemVec( cos(theta) , sin(theta) );
        //}
      }
      else if (c.sigma){
        c.setChemVec(0.0,0.0);
        c.setChemVecLength(0.0);
      }
      //if(c.chemvecx>1 || c.chemvecy>1){
      //  std::cerr << ", vector: "<< c.chemvecx <<" "<< c.chemvecy  << '\n';
      //  exit(1);
      //}
    }

}

//to initialise cells' mu, perstime and persdur
void Dish::InitCellMigration(void)
{
  auto icell = std::begin(cell);
  ++icell;  //discard first element of vector cell, which is medium

  //when the initialisation period has passed: start up the vectors and migration
  for(auto end=std::end(cell); icell != end; ++icell)
  {
    icell->setMu(par.startmu);
    icell->startTarVec();
    if (par.persduration<par.mcs){
      icell->setPersTime(int(par.persduration*RANDOM())); //so that cells don't all start turning at the same time...
    }
    else{
      icell->setPersTime(0); //special type of experiment
      icell->tvecx=1.;
      icell->tvecy=0.;
    }
    icell->setPersDur(par.persduration);

    icell->setChemMu(par.init_chemmu);
    icell->startChemVec();
    //cerr<<"Cell "<<icell->sigma<<" vecx="<<icell->tvecx<<" vecy="<<icell->tvecy<<endl;
    //cerr<<"Cell persdur "<<icell->persdur<<" perstime "<<icell->perstime<<endl;
  }
cerr<<"init chemmu is"<<par.init_chemmu<<endl;
}

//function to have cells update their persistence time (perstime);
//In the future perhaps also their persistence duration (persdur), or how long they remember their preferred direction;
//and force of migration (mu)
void Dish::CellMigration(void)
{
  auto icell = std::begin(cell);
  ++icell;  //discard first element of vector cell, which is medium
  for(auto end=std::end(cell); icell != end; ++icell)
  {
    if( ! icell->AliveP() ) continue; //if cell is not alive, continue

    icell->updatePersTime();

  }
}

//Function for cell growth and division based on uptake of food particles
void Dish::CellGrowthAndDivision2(void)
{
   //cout<<"Hello beginning CellGrowthAndDivision2 "<<endl;
  vector<bool> which_cells(cell.size()); //which cells will divide
  vector<int> sigma_newcells;
    // DIVISION VOLUME HAS BECOME INTERNAL CELL VARIABLE
    //static int mem_area=0;
    // if called for the first time: calculate mem_area
    //if (!mem_area) {
    //  mem_area=TargetArea()/CountCells();
    //}

    vector<Cell>::iterator c; //iterator to go over all Cells
    int area;
    double newar;
    int newarint;
    int celldivisions=0;

    //cells grow due to food intake, and shrink constantly due to metabolism
    for( c=cell.begin(), ++c; c!=cell.end(); ++c){
      if( c->AliveP() ){
        //cerr<<"This cell's sigma: "<<c->Sigma()<<", Alive? "<< c->AliveP() << ", target area: "<<c->TargetArea()<<endl;

        //growth and shrinkage due to particles
        area=c->TargetArea();
        newar=double(area);

        double particles_metabolised = 0.;
        double particles_for_movement = 0.;
        //in this version, maintenance_fraction is a function of evolvable parameters:
        // m = k0 + Area * kA/Division area + particles * kP/50 + <J contact>/length
        // and bounded inside [0,1]
        c->maintenance_fraction =1.;//c->CalculateMaintenance_or_ExtProtExpr_Fraction(c->k_mf_0,c->k_mf_A,c->k_mf_P,c->k_mf_C);

        //Next, also Js should be a function of that, bounded between [reasonably high, actual J values]
        // this is just a number that is going to be multiplied to the J values of this cell...
        // ... how? maybe directly in amoebamove? (it'd be easier than updating all J values for this guy
        // and for all those in contact with this guy)

        //same function, be careful which parameters you pass
        c->extprotexpress_fraction = c-> CalculateMaintenance_or_ExtProtExpr_Fraction(c->k_ext_0,c->k_ext_A,c->k_ext_P,c->k_ext_C);
        // if (c->extprotexpress_fraction<0. || c->extprotexpress_fraction>1.) {
        //  std::cerr <<"Sigma: "<<c->Sigma()<< ", extprotexpress_fraction: "<< c->extprotexpress_fraction << '\n';
        // }

        //same function for regulation of chemotaxis
        c->weight_for_chemotaxis =  c-> CalculateMaintenance_or_ExtProtExpr_Fraction(c->k_chem_0,c->k_chem_A,c->k_chem_P,c->k_chem_C);
        // if (c->weight_for_chemotaxis<0. || c->weight_for_chemotaxis>1.) {
          // std::cerr <<"Sigma: "<<c->Sigma()<< ", weight_for_chemotaxis: "<< c->weight_for_chemotaxis << '\n';

        // }

        if(area){
          // particles_metabolised are those particles that are used for maintenance
          // only a fraction should be used for this
          //the other should be used for movement
          particles_metabolised = c->maintenance_fraction* c->particles/(double)par.scaling_cell_to_ca_time;
          particles_for_movement = (1. - c->maintenance_fraction) * c->particles/(double)par.scaling_cell_to_ca_time;

          //newar+= c->growth*particles_metabolised/double(area);
          newar+= c->growth*particles_metabolised;

          if(RANDOM() < par.ardecay/(double)par.scaling_cell_to_ca_time )
            newar -= 1; //area decays of the same amount per time step, on average.
            //newar=double(area) - double(area)*par.ardecay + c->growth*double(c->particles)/double(area);
        }
        newarint=int(newar);
        if( RANDOM() < fabs(newar-double(newarint)) )
            newarint+=ceil(newar-newarint);

        //cerr<<"area growth?"<<endl;
        //exit(1);

        //if(newarint > 2*c->getHalfDivArea()) newarint=area;

        //OUTCOMMENT ME WHEN YOU ARE DONE-DONE
        if(newarint > 10*c->getHalfDivArea()) newarint=area; //if target area too large,
                                                             //it does not grow anymore!
                                                             //But food IS consumed

        c->SetTargetArea(newarint);


        //Setting mu:
        // particles_for_movement maps to mu as follows:
        // mu = particles_for_movement* 10 (i.e. )
        // if mu > 10 -> mu = 10 and some particles_for_movement are not consumed, and should be recycled
        //  i.e. you use at most one particle per time step for movement

        // FIRST: check how much movement this would be
        //cerr<<"particles for movement: "<< particles_for_movement <<",thus mu: ";
        double newmu;
        if(particles_for_movement > 1. ){
          newmu=3.;
          particles_for_movement=1.;
        }else{
          // newmu = 10.*particles_for_movement;
          newmu = 3.*particles_for_movement;
        }
        c->mu = newmu;

        c->mu = par.startmu; particles_for_movement=0; // specified experiments
        c->chemmu = par.init_chemmu;
        //cerr<<newmu<<endl;

        //if(c->growth<1.)
          //c->particles=int(double(c->particles)*c->growth);
        //else
        c->particles-= (particles_metabolised + particles_for_movement);
        //cout<<"cell nr "<<c->Sigma()<<", cell type "<<c->getTau()<<endl;
        //cout <<"area is "<<area<<", new area is "<<newarint<<endl;


        //OUTCOMMENT ME WHEN YOU ARE DONE - DONE
        //is cell big enough to divide?
        if( c->Area() > 2*c->half_div_area){
          which_cells[c->Sigma()]=true;
          celldivisions++;
        }



        //make sure target_area is not negative
        if(c->TargetArea()<0)
          c->SetTargetArea(0);
      }
      //check area:if cell is too small (whether alive or not) we remove its sigma
      // notice that this keeps the cell in the cell array, it only removes its sigma from the field
      if(c->Area()< par.min_area_for_life){
         c->SetTargetArea(0);
         c->Apoptose(); //set alive to false
         CPM->RemoveCell(&*c,par.min_area_for_life,c->meanx,c->meany);
      }


    }

//     cerr<<"Cell vector size:"<<cell.size()<<endl;

    if(celldivisions){
      // cerr<<"Hello7"<<endl;

      // All this can be deleted ***************
//        int i=0;
//        cout << "Following cells will divide: ";
//        for(vector<bool>::iterator d = which_cells.begin(); d != which_cells.end() ;d++){
//          if((*d)) cout << i <<" ";
//          i++;
//        }
//        cout << endl;
      //
      // till here   ***************************

//         cerr<<"Hello CellGrowthAndDivision2 0.1"<<endl;
        sigma_newcells = CPM->DivideCells(which_cells);
        // sigma_newcells is an int vector as lognas there are cells,
        // it is zero everywhere, except at the positions of a mother cell's sigma,
        // where it contains as value the sigma of the daughter cell

//         cerr<<"Hello CellGrowthAndDivision2 0.2"<<endl;


        //Mutate Cells:
        MutateCells(sigma_newcells);
        //cerr<<"Hello"<<endl;
        //for(auto asigma: sigma_newcells) cerr<<asigma<<" ";
        //cerr<<endl;
        // Update sigmas of everybody with the new born cells
        UpdateVectorJ(sigma_newcells);
//         cerr<<"Hello CellGrowthAndDivision2 0.3"<<endl;
        //cerr<<"Hello7.1"<<endl;

//         for(auto s267: sigma_newcells) if(s267==267){
//           cerr<<"Printing cell[55] jvalue with 267 as soon as 267 is born: "<<cell[55].vJ[267]<<endl;
//           cerr<<"Is cell 55 alive? "<<cell[55].alive<<endl;
//
//         }
    }
}


//Function that checks and changes cell parameters
void Dish::UpdateCellParameters(void)
{
   //cout<<"Hello beginning CellGrowthAndDivision2 "<<endl;
  vector<bool> which_cells(cell.size()); //which cells will divide
  vector<int> sigma_newcells;

    vector<Cell>::iterator c; //iterator to go over all Cells
    int area;
    double newar;
    int newarint;
    int celldivisions=0;

    for( c=cell.begin(), ++c; c!=cell.end(); ++c){
      if( c->AliveP() ){

        //in this version, maintenance_fraction is a function of evolvable parameters:
        // m = k0 + Area * kA/Division area + particles * kP/50 + <J contact>/length
        // and bounded inside [0,1]
        c->maintenance_fraction =1.;//c->CalculateMaintenance_or_ExtProtExpr_Fraction(c->k_mf_0,c->k_mf_A,c->k_mf_P,c->k_mf_C);

        //Next, also Js should be a function of that, bounded between [reasonably high, actual J values]
        // this is just a number that is going to be multiplied to the J values of this cell...
        // ... how? maybe directly in amoebamove? (it'd be easier than updating all J values for this guy
        // and for all those in contact with this guy)

        //same function, be careful which parameters you pass
        c->extprotexpress_fraction = c-> CalculateMaintenance_or_ExtProtExpr_Fraction(c->k_ext_0,c->k_ext_A,c->k_ext_P,c->k_ext_C);

        //same function for regulation of chemotaxis
        c->weight_for_chemotaxis =  c-> CalculateMaintenance_or_ExtProtExpr_Fraction(c->k_chem_0,c->k_chem_A,c->k_chem_P,c->k_chem_C);

        //c->mu = par.startmu;
        //c->chemmu = par.init_chemmu;

        //std::cerr << "Cell with sigma = "<<c->Sigma()<<" has particles: "<<c->particles << '\n';
      }
      //check area:if cell is too small (whether alive or not) we remove its sigma
      // notice that this keeps the cell in the cell array, it only removes its sigma from the field
      if(c->Area()< par.min_area_for_life){
         c->SetTargetArea(0);
         c->Apoptose(); //set alive to false
         CPM->RemoveCell(&*c,par.min_area_for_life,c->meanx,c->meany);
      }


    }

}

// As the one after, except that the reproduction zone is a semicircle around the peak of the gradient
int Dish::CheckWhoMadeitRadial(void){
  unsigned int howmany_makeit_for_nextgen = par.howmany_makeit_for_nextgen; //we do this to cast the par, which is int, to unsigned int

  //get info where the peak is
  int peakx = Food->GetPeakx();
  int peaky = Food->GetPeaky();

  //find if cells in some area around the peak are already in the list
  for(auto &c:cell){
    if( ! c.AliveP() ) continue;
    double distx = c.meanx - peakx;
    double disty = c.meany - peaky;
    double dist_from_peak = hypot(distx,disty);

    //cerr<< "Sigma: "<< c.Sigma()<<". Dist = "<<dist_from_peak<< ", mu = "<<c.mu<< endl;

    //cerr<< "Sigma: "<< c.Sigma()<<". Dist = "<<dist_from_peak<< ", persdur = "<<c.persdur<<", perstime = "<<c.perstime<< endl;

    if( dist_from_peak <= par.the_line ){
      who_made_it.insert( c.Sigma() ); //if already there it will not be duplicated in the set
      if( par.zero_persistence_past_theline ) {
        c.setPersDur(1); //zero persistence - hence moving - when you cross the line
                         // don't forget that this has to be reset to a value >0 after replication
        c.setPersTime(0);
      }
    }
  }

  // cerr<< "px,py: "<< peakx<<" "<<peaky <<" box mx,My my,My: "<< minx<<" "<<maxx<<" "<<miny<<" "<<maxy<<endl;
  // cerr<< "who_made_it has so many members: "<< who_made_it.size() << endl;

  //if list is large enough return 1, else 0
  if( who_made_it.size() >= howmany_makeit_for_nextgen ) {
    // std::cerr << "Many made it !" << '\n';
    //who_made_it.clear();
    return 1;
  }
  // return 1;
  return 0;
}


//checks if sufficiently many cells made it into the reproduction zone
//it keeps a list of sigmas of cells that are there
//if length of list is large enough -> returns 1
// should create a new list of cells based on some fitness function,
// kill everybody, place these new cells, change gradient direction + add new food
int Dish::CheckWhoMadeitLinear(void){
  unsigned int howmany_makeit_for_nextgen = par.howmany_makeit_for_nextgen; //we do this to cast the par, which is int, to unsigned int
  // int par.the_line = 41;

  //as gradients are now, there is always a coordinate that is either 1 or size_x_or_y,
  // while the other is size_y/2_or_x/2
  //we can predetermine where a minx maxx miny maxy rectangle should be based on peakx and peaky
  // and run a for loop only within these numbers, and check if CPM->sigma[][] != 0
  static int current_peakx=-1,current_peaky=-1;
  static int minx,maxx,miny,maxy;

  //get info where the peak is
  int peakx = Food->GetPeakx();
  int peaky = Food->GetPeaky();


  if(peakx != current_peakx || peaky != current_peaky){
    current_peakx=peakx;
    current_peaky=peaky;

    if(peakx==1) {
      //then peaky = sizey/2 and
      minx = 1;
      maxx = par.the_line;
      miny = 1;
      maxy = par.sizey-1;
    }else if(peakx==par.sizex-1){
      minx = par.sizex-par.the_line;
      maxx = par.sizex-1;
      miny = 1;
      maxy = par.sizey-1;
    }else if(peaky==1){
      minx = 1;
      maxx = par.sizex-1;
      miny = 1;
      maxy = par.the_line;
    }else if(peaky==par.sizey-1){
      minx = 1;
      maxx = par.sizex-1;
      miny = par.sizey-par.the_line;
      maxy = par.sizey-1;
    }else{
      cerr<<"CheckWhoMadeit(): Error. Got weird peakx and peaky position: peakx, peaky = "<<peakx<<", "<<peaky<<endl;
      std::cerr << "Don't know what to do with this, program exits now." << '\n';
      exit(1);
    }
  }

  //find if cells in some area around the peak are already in the list
  //easy: go in order through CPM and check the sigmas
  for(int i=minx;i<maxx;i++)for(int j=miny;j<maxy;j++){
    if( CPM->Sigma(i,j)!=0 ){
      who_made_it.insert( CPM->Sigma(i,j) ); //if already there it will not be duplicated in the set
    }
  }

  // cerr<< "px,py: "<< peakx<<" "<<peaky <<" box mx,My my,My: "<< minx<<" "<<maxx<<" "<<miny<<" "<<maxy<<endl;
  // cerr<< "who_made_it has so many members: "<< who_made_it.size() << endl;

  //if list is large enough return 1, else 0
  if( who_made_it.size() > howmany_makeit_for_nextgen ) {
    // std::cerr << "Many made it !" << '\n';
    //who_made_it.clear();
    return 1;
  }
  // return 1;
  return 0;
}

//remove cells from dish and CPM based on indexes in who_made_it
void Dish::RemoveWhoDidNotMakeIt(void)
{
  // std::cerr << "Who made it: ";
  // for(auto sig:who_made_it) std::cerr << sig<<" ";
  // std::cerr << '\n';

  vector<Cell>::iterator c; //iterator to go over all Cells
  for( c=cell.begin(), ++c; c!=cell.end(); ++c){

    // cerr<<"This sigma: "<<c->Sigma();

    if( who_made_it.count( c->Sigma() ) == 0 ) {
      // cerr<<" will be removed"<<endl;
      c->SetTargetArea(0);
      c->Apoptose(); //set alive to false
      CPM->RemoveCell(&*c,par.min_area_for_life,c->meanx,c->meany);
    }
    // else{
    //   cerr<<" will not be removed"<<endl;
    // }
  }
}

double Dish::FitnessFunction(int particles, double meanx, double meany)
{
  //get info where the peak is
  int peakx = Food->GetPeakx();
  int peaky = Food->GetPeaky();

  double dx = meanx - peakx;
  double dy = meany - peaky;
  double dist = hypot( dx,dy );

  double epsilon = 0.05;
  double h_food = 10;
  double h_dist = par.the_line;

  double fitness_food = ( particles + epsilon*h_food/(1.-2.*epsilon)  )/(particles + h_food/(1.-2.*epsilon)); //looks weird, it's not (if you plot it)
  double fitness_distance = 1. / ( 1. + pow( dist/h_dist , 2.) );
  // std::cerr << "particles: "<< particles <<", fitness food = " <<fitness_food<<", distance: "<<dist<<", fitness distance" << fitness_distance<<", tot = "<<fitness_food*fitness_distance<<'\n';
  return fitness_food*fitness_distance;

}

//Based on ReproduceWhoMadeIt3, this reproduces all cells,
// with fitness dependent on food and distance from target
void Dish::ReproduceEndOfSeason(void)
{
  vector<bool> which_cells(cell.size()); //which cells will divide
  vector<double> particles_of_those_whomadeit(cell.size());
  vector<int> sigma_of_cells_that_will_divide;
  vector<int> sigma_newcells;
  std::vector<double> fitness(cell.size(), 0.); //as many as there are cells (dead or alive), all with value 0.
  double tot_fitness = 0.;

  for(auto &c : cell){
    if(c.AliveP() && c.Sigma()>0){
      double cell_fitness = FitnessFunction( c.particles, c.getXpos(), c.getYpos());
      fitness[ c.Sigma() ] = cell_fitness;
      tot_fitness += cell_fitness;

      c.mu = 0.;  //also the other mu? yes-
      c.chemmu = 0.;
    }
  }
  // std::cerr << "tot fitness: "<<tot_fitness << '\n';
  //for(auto &f:fitness) f/=tot_fitness;  <- what the hell?
  for(int i = 0; i<par.howmany_makeit_for_nextgen; ++i){
    double sum_fitness=0.;
    double rn = tot_fitness*RANDOM(); //random choice of who replicates proportional to fitness
    // std::cerr << "rn = " <<rn << '\n';
    int which_sig=-1; //guaranteed to be initialised by the end of the loop
    int fitness_vector_size=(int)fitness.size(); // so that compiler stops whining about this
    for(int sig =1; sig< fitness_vector_size; ++sig ){
      //if(fitness[sig] < 0.) continue; //fitness vector is as big as cell, and contains -1 where sigma is not alive in cell
      sum_fitness+= fitness[sig];
      // std::cerr << "sum so far = " <<sum_fitness << '\n';
      if(sum_fitness>rn){
        which_sig=sig; //we found that sigma that replicates
        break;
      }
    }
    if(which_sig==-1){
      std::cerr << "ReproduceEndOfSeason(): Error. which_sig is still uninitialised" << '\n';
      exit(1);
    }
    sigma_of_cells_that_will_divide.push_back(which_sig);
  }

  // std::cerr << "These are the sigmas of the cells alive now, before cell division" << '\n';
  // for(auto c: cell)
  //   if(c.AliveP()) std::cerr << c.Sigma()<<" ";
  // std::cerr << " " << '\n';
  // std::cerr << "These are the sigmas of the cells dead now, before cell division" << '\n';
  // for(auto c: cell)
  //   if(!c.AliveP()) std::cerr << c.Sigma()<<" ";
  // std::cerr << " " << '\n';
  //
  // std::cerr << "These are the cells that will divide" << '\n';
  // for(auto sig:sigma_of_cells_that_will_divide) std::cerr << sig<<" ";
  // std::cerr << " " << '\n';

  //At this point we should orchestrate actual cell division:
  // some cells replicate 10 times, other zero: the idea is that we let cells divide
  // if they are bigger than a certain amount, if not, we run amoeabeamove until they have expanded enough
  std::vector<int> new_sigma_of_cells_that_will_divide;
  while( ! sigma_of_cells_that_will_divide.empty() ){
    for(auto sig: sigma_of_cells_that_will_divide){
      //if cell[sig] is large enough
      //if which_cells[sig] is not already marked
      // std::cerr << "Cell: "<<sig <<" has area "<< cell[sig].area << " and which_cells[sig] = "<<which_cells[sig] ;
      // put in which_cells
      if(cell[sig].Area() > 30 && which_cells[sig]==false){
        which_cells[sig]=true;
        // std::cerr << " therefore yes" << '\n';
      }
      // if not, put in new_sigma_of_cells_that_will_divide
      else{
        new_sigma_of_cells_that_will_divide.push_back(sig);
        // std::cerr << " therefore no" << '\n';
      }
    }
    //make cell division
    // std::cerr << "\nThese are the sigmas that replicate this round"<<endl;
    // for(auto x: which_cells) std::cerr << x <<" ";
    // std::cerr << " " << '\n';
    // std::cerr <<"cell.size()="<<cell.size()<< ", vector size = " << which_cells.size() << '\n';
    //it can still be that which cells is empty because cells that should divide are still too small
    // it is VERY unlikely that this happens, but it can happen:
    if(!which_cells.empty()){
      sigma_newcells = CPM->DivideCells(which_cells); //replicate cells
      // sigma_newcells is an int vector as long as there are cells,
      // it is zero everywhere, except at the positions of a mother cell's sigma,
      // where it contains as value the sigma of the daughter cell
      MutateCells(sigma_newcells);
      UpdateVectorJ(sigma_newcells);

      // std::cerr << "\nThese are the sigmas of new born cells"<<endl;
      // for(auto x: sigma_newcells) std::cerr << x <<" ";
      // std::cerr << " " << '\n';
      // std::cerr << "There are now so many alive cells: "<<CountCells() << '\n';
      //zero the which_cells vector
      std::fill(which_cells.begin(), which_cells.end(), false);
      //vector has to be resized because new cells are born
      which_cells.resize(cell.size(),false);

      // std::cerr << "\nNow the vector which cells has to be resized to cell.size() and all zeroed. check:"<<endl;
      // for(auto x: which_cells) std::cerr << x <<" ";
      // std::cerr << " " << '\n';
      // if(which_cells.size() == cell.size()) std::cerr << "The two vectors have the same size" << '\n';
      // exit(1);
    }
    //reset target area

    //just for checking
    //vector<Cell>::iterator c; //iterator to go over all Cells
    // int counter=0;

    //reset target area, because cell division changes
    // for(auto sig: sigma_of_cells_that_will_divide){
    //     cell[sig].SetTargetArea(par.target_area); // for good measure
    //
    // }
    //return;


    //do a bunch of AmoebaeMove2
    for(int i=0;i<5;i++) CPM->AmoebaeMove2(PDEfield); // let them expand a little
    //copy new_sigma_of_cells_that_will_divide into old one
    // std::cerr << "After AmoebaeMove2 there are so many alive cells: "<<CountCells() << '\n';


    sigma_of_cells_that_will_divide = new_sigma_of_cells_that_will_divide;
    new_sigma_of_cells_that_will_divide.clear();

    //checks - PASSED
    // std::cerr << "Just a check that vector copying went fine\n sigma_of_bla... should not be empty"<<endl;
    // for(auto x: sigma_of_cells_that_will_divide) std::cerr << x <<" ";
    // std::cerr << " " << '\n';
    // std::cerr << "new_sigma_of_bla... should be empty "<<endl;
    // for(auto x: new_sigma_of_cells_that_will_divide) std::cerr << x <<" ";
    // std::cerr << " " << '\n';


    // return;
  }

  //also needed, to set all daugher cells
  vector<Cell>::iterator c; //iterator to go over all Cells
  int counter=0;
  for( c=cell.begin(), ++c; c!=cell.end(); ++c){
    if(c->AliveP()){
      // c->SetTargetArea(par.target_area);
      c->particles = 0;

      c->mu = par.startmu;
      c->chemmu = par.init_chemmu;

      if(par.zero_persistence_past_theline ) c->setPersDur(par.persduration); //de-zero persistence so they can move again
                                                                              // might have to randomzie this
      counter++;
    }
  }
  // std::cerr << "Counter sees so many cells: "<<counter << '\n';

}

//random death until popsize is restored
void Dish::RemoveCellsUntilPopIs(int popsize)
{
  int current_popsize=0;
  std::vector<int> alivesigma;
  for(const auto c: cell){
    if(c.Sigma()>0 && c.AliveP()){
      current_popsize++;
      alivesigma.push_back(c.Sigma());
    }
  }
  // std::cerr << "going to remove cells until I get popsize = "<<popsize << '\n';
  while(current_popsize>popsize){
    // std::cerr << "current_popsize (faulty) = " << current_popsize << '\n';

    int rn = current_popsize*RANDOM();
    int sigtorm = alivesigma[rn];
    // std::cerr << "sig to rm = " << sigtorm << '\n';
    // std::cerr << "sig.back() that takes its place = "<< alivesigma.back() << '\n';
    alivesigma[rn] = alivesigma[current_popsize-1];
    current_popsize--; //so that next time rn interval is reduced

    cell[sigtorm].SetTargetArea(0);
    cell[sigtorm].Apoptose(); //set alive to false
    CPM->RemoveCell(&cell[sigtorm] ,par.min_area_for_life,cell[sigtorm].meanx,cell[sigtorm].meany);
    // std::cerr << "There are now so many cells: "<< CountCells() << '\n';
  }
}

//give cells a target area proportional to how many particles they have (relative to others)
// so that cell division does not lead to cells with zero area?
// the trick is: we let cells divide if their area is large enough
void Dish::ReproduceWhoMadeIt3(void)
{
  unsigned int popsize = par.popsize;

  vector<bool> which_cells(cell.size()); //which cells will divide
  vector<double> particles_of_those_whomadeit(cell.size());
  vector<int> sigma_of_cells_that_will_divide;
  vector<int> sigma_newcells;

  double tot_eaten_particles=0.;
  unsigned int current_popsize = who_made_it.size();

  // checked
  // std::cerr << "Just a check" << '\n';
  // for(auto c:cell){
  //   if(c.AliveP() && c.Sigma()) std::cerr << "Sigma: "<< c.Sigma() << ", particles: "<< c.particles << '\n';
  // }

  for(auto sig: who_made_it){
    tot_eaten_particles+=cell[sig].particles;
    particles_of_those_whomadeit[sig] = cell[sig].particles; //then use this to iterate
                                                             // the reason is that daughter cells themselves don't replicate
                                                             // so after halving particles with them a mother cell decreases its own exponentially
                                                             // which screws up the probability of replication
    //cell[sig].SetTargetArea(4*par.target_area);
    cell[sig].mu = 0.;  //also the other mu??
    cell[sig].chemmu = 0.;
  }

  // for(int i=0;i<100;i++) CPM->AmoebaeMove2(PDEfield); // let them expand a little more

  // What happens if no-one ate anything?
  // Should we give a little chance of reprodution to cells that ate nothing?
  // yes: espilon = 0.1
  double epsilon=0.1;
  tot_eaten_particles+= epsilon*current_popsize;

  //At this point we have eveybody's fitness ( (particles+epsilon)/tot_eaten_particles )
  // we first generate the sigmas that divide, and then we do the actual divisions
  while(current_popsize<popsize){
    double sum_particles=0.;
    double rn = tot_eaten_particles*RANDOM(); //random choice of who replicates proportional to particles
    int which_sig=-1; //guaranteed to be initialised by the end of the loop
    for(auto sig: who_made_it){
      sum_particles+= epsilon+particles_of_those_whomadeit[sig];
      if(sum_particles>rn){
        which_sig=sig; //we found that sigma that replicates
        break;
      }
    }
    if(which_sig==-1){
      std::cerr << "ReproduceWhoMadeIt3(): Error. which_sig is still uninitialised" << '\n';
      exit(1);
    }
    sigma_of_cells_that_will_divide.push_back(which_sig);
    current_popsize++;
  }
  // std::cerr << "These are the cells that will divide" << '\n';
  // for(auto sig:sigma_of_cells_that_will_divide) std::cerr << sig<<" ";
  // std::cerr << " " << '\n';
  //
  //At this point we should orchestrate actul cell division:
  // some cells replicate 10 times, other zero: the idea is that we let cells divide
  // if they are bigger than a certain amount, if not, we run amoeabeamove until they have expanded enough
  std::vector<int> new_sigma_of_cells_that_will_divide;
  while( ! sigma_of_cells_that_will_divide.empty() ){
    for(auto sig: sigma_of_cells_that_will_divide){
      //if cell[sig] is large enough
      //if which_cells[sig] is not already marked
      // std::cerr << "Cell: "<<sig <<" has area "<< cell[sig].area << " and which_cells[sig] = "<<which_cells[sig] ;
      // put in which_cells
      if(cell[sig].Area() > 30 && which_cells[sig]==false){
        which_cells[sig]=true;
        // std::cerr << " therefore yes" << '\n';
      }
      // if not, put in new_sigma_of_cells_that_will_divide
      else{
        new_sigma_of_cells_that_will_divide.push_back(sig);
        // std::cerr << " therefore no" << '\n';
      }
    }
    //make cell division
    // std::cerr << "\nThese are the sigmas that replicate this round"<<endl;
    // for(auto x: which_cells) std::cerr << x <<" ";
    // std::cerr << " " << '\n';
    sigma_newcells = CPM->DivideCells(which_cells); //replicate cells
    MutateCells(sigma_newcells);
    UpdateVectorJ(sigma_newcells);

    // std::cerr << "There are now "<< CountCells() <<" cells" << '\n';

    //zero the which_cells vector
    std::fill(which_cells.begin(), which_cells.end(), false);

    //reset target area

    //just for checking
    //vector<Cell>::iterator c; //iterator to go over all Cells
    // int counter=0;

    //reset target area, because cell division changes
    for(auto sig: sigma_of_cells_that_will_divide){
        cell[sig].SetTargetArea(par.target_area); // for good measure

    }
    //do a bunch of AmoebaeMove2
    for(int i=0;i<5;i++) CPM->AmoebaeMove2(PDEfield); // let them expand a little
    //copy new_sigma_of_cells_that_will_divide into old one
    sigma_of_cells_that_will_divide = new_sigma_of_cells_that_will_divide;
    new_sigma_of_cells_that_will_divide.clear();

    //checks - PASSED
    // std::cerr << "Just a check that vector copying went fine\n sigma_of_bla... should not be empty"<<endl;
    // for(auto x: sigma_of_cells_that_will_divide) std::cerr << x <<" ";
    // std::cerr << " " << '\n';
    // std::cerr << "new_sigma_of_bla... should be empty "<<endl;
    // for(auto x: new_sigma_of_cells_that_will_divide) std::cerr << x <<" ";
    // std::cerr << " " << '\n';
    //

    // return;
  }

  //also needed, to set all daugher cells
  vector<Cell>::iterator c; //iterator to go over all Cells
  int counter=0;
  for( c=cell.begin(), ++c; c!=cell.end(); ++c){
    if(c->AliveP()){
      c->SetTargetArea(par.target_area);
      if(par.zero_persistence_past_theline ) c->setPersDur(par.persduration); //de-zero persistence so they can move again
                                                                              // might have to randomzie this
      counter++;
    }
  }

  // cerr<<"At the end of ReproduceWhoMadeIt3 there are so many cells: "<<counter<<endl;
}

//food dependent replication (andif no food, neutral replication)
void Dish::ReproduceWhoMadeIt2(void)
{
  unsigned int popsize = par.popsize;

  vector<bool> which_cells(cell.size()); //which cells will divide
  vector<double> particles_of_those_whomadeit(cell.size());
  vector<int> sigma_newcells;

  double tot_eaten_particles=0.;
  unsigned int current_popsize = who_made_it.size();

  for(auto sig: who_made_it){
    tot_eaten_particles+=cell[sig].particles;
    particles_of_those_whomadeit[sig] = cell[sig].particles; //then use this to iterate
                                                             // the reason is that daughter cells themselves don't replicate
                                                             // so after halving particles with them a mother cell decreases its own exponentially
                                                             // which screws up the probability of replication
    cell[sig].SetTargetArea(4*par.target_area);
    cell[sig].mu = 0.;  //also the other mu??
    cell[sig].chemmu = 0.;
  }

  // for(auto sig: who_made_it){
  //   cell[sig].SetTargetArea(  );
  // }

  for(int i=0;i<100;i++) CPM->AmoebaeMove2(PDEfield); // let them expand a little more

  // What happens if no-one ate anything?
  // Should we give a little chance of reprodution to cells that ate nothing?
  // yes: espilon = 0.1
  double epsilon=0.1;
  tot_eaten_particles+= epsilon*current_popsize;

  while(true){
    //strategy is to fill which_cells until a repetition happens,
    // at which point we duplicate the cells in which_cells
    // then we zero which_cells
    // and start again to fill it
    // it's not the most optimal but better than replicating one cell at a time

    double sum_particles=0.;
    double rn = tot_eaten_particles*RANDOM(); //random choice of who replicates proportional to particles
    int which_sig=-1; //guaranteed to be initialised by the end of the loop
    for(auto sig: who_made_it){
      sum_particles+= epsilon+particles_of_those_whomadeit[sig];
      if(sum_particles>rn){
        which_sig=sig; //we found that sigma that replicates
        break;
      }
    }
    if(which_sig==-1){
      std::cerr << "ReproduceWhoMadeIt2(): Error. which_sig is still uninitialised" << '\n';
      exit(1);
    }
    current_popsize++;

    if(which_cells[which_sig] == false){
      which_cells[which_sig] = true;
      cell[which_sig].mu = 0.;        // set mu to zero here so that they do not migrate in amoeabeamove later
                                      // maybe also the other mu has to be set to zero?
      cell[which_sig].chemmu = 0.;
    }
    else{

      sigma_newcells = CPM->DivideCells(which_cells); //replicate cells
      MutateCells(sigma_newcells);
      UpdateVectorJ(sigma_newcells);
      std::cerr << "There are now "<< CountCells() <<" cells" << '\n';

      for(auto a: sigma_newcells){
        if(a) cell[a].SetTargetArea(par.target_area);
      }
      for(int i=0;i<200;i++) CPM->AmoebaeMove2(PDEfield); // let them expand a little <- this is not working super well
      // return;

      std::fill(which_cells.begin(), which_cells.end(), false); //zero the vector
      which_cells[which_sig] = true; //set the cell to replicate
      cell[which_sig].mu = 0.;  // set mu to zero here so that they do not migrate in amoeabeamove later
      cell[which_sig].chemmu = 0.;
    }

    //if we are done we should replicate the last cells in which_cells and we are done
    if(current_popsize==popsize) {
      if(which_cells.size()>0){
        sigma_newcells = CPM->DivideCells(which_cells);
        MutateCells(sigma_newcells);
        UpdateVectorJ(sigma_newcells);
        std::cerr << "Last repl done, are now "<< CountCells() <<" cells" << '\n';
      }
      for(auto a: sigma_newcells){
        if(a) cell[a].SetTargetArea(par.target_area);
      }

      for(int i=0;i<200;i++) CPM->AmoebaeMove2(PDEfield); // let them expand a little more
      std::cerr << "Done with repl, resulting cell sizes:" << '\n';
      for(auto c:cell) if(c.Sigma() && c.AliveP() ) std::cerr << c.area << '\n';

      break;
    }


  }

  vector<Cell>::iterator c; //iterator to go over all Cells
  int counter=0;
  for( c=cell.begin(), ++c; c!=cell.end(); ++c){
    if(c->AliveP()){
      c->SetTargetArea(par.target_area); // for good measure
      counter++;
    }
  }

  // cerr<<"At the end of ReproduceWhoMadeIt2 there are so many cells: "<<counter<<endl;
}




// previous version, without food-dependent fitness
void Dish::ReproduceWhoMadeIt(void)
{
  vector<bool> which_cells(cell.size()); //which cells will divide
  for(auto sig: who_made_it){
    which_cells[sig] = true;
  }
  vector<int> sigma_newcells;
  int howmany_rounds_of_division = 2;
  for(int i=0;i<howmany_rounds_of_division;i++){
    sigma_newcells = CPM->DivideCells(which_cells);
    MutateCells(sigma_newcells);
    UpdateVectorJ(sigma_newcells);
  }
  vector<Cell>::iterator c; //iterator to go over all Cells
  for( c=cell.begin(), ++c; c!=cell.end(); ++c){
    c->SetTargetArea(par.target_area); // for good measure
  }
}

int Dish::CountCells(void) const {

  int amount=0;
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),++i); i!=cell.end(); ++i) {
    if (i->AliveP()) {
      amount++;
    } //else {
      //cerr << "Dead cell\n";
      //}
  }
  return amount;
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
  for ( (i=cell.begin(),i++); i!=cell.end();++i) {

    if (i->AliveP())
      total_area+=i->TargetArea();

  }
  return total_area;
}

int Dish::CountPreys(void) {
  int howmany=0;
  vector<Cell>::iterator i;

  for ( (i=cell.begin(),i++); i!=cell.end(); ++i) {
    if( i->AliveP() && i->getTau() == PREY) howmany++;
  }
  //cout << "These many preys" << howmany<<endl;
  //exit(1);
  return howmany;
}
int Dish::CountPredators(void){
  int howmany=0;
  vector<Cell>::iterator i;

  for ( (i=cell.begin(),i++); i!=cell.end(); ++i) {
    if( i->AliveP() && i->getTau() == PREDATOR) howmany++;
  }
  //cout << "These many preys" << howmany<<endl;
  //exit(1);
  return howmany;
}

int Dish::SaveData(int Time)
{
  std::ofstream ofs;
  ofs.open (par.datafile, std::ofstream::out | std::ofstream::app);
  int pred=0,prey=0;

  //make file where evey so often you dump everybody - EASY!
  // save for each indiv
  // Time pred/prey key lock J with medium, neighbors
  auto icell = std::begin(cell);
  ++icell;  //discard first element of vector cell, which is medium
  for(auto end=std::end(cell); icell != end; ++icell){
    if( ! icell->AliveP() ) continue; //if cell is not alive, continue

    int itau=icell->getTau();
    int isigma=icell->Sigma();
    if( itau == PREDATOR ) pred++;
    else if(itau == PREY) prey++;
    else{
      cerr<<"SaveData(): Error. Got cell that is neither prey, nor predator"<<endl;
      exit(1);
    }
    ofs << Time << " "<<isigma<<" "<< itau << " "; // Time now, tau of me
    ofs<< icell->getXpos()<<" "<<icell->getYpos()<<" ";
    ofs<<icell->getXvec()<<" "<<icell->getYvec()<<" ";
    ofs<<icell->getChemXvec()<<" "<<icell->getChemYvec()<<" ";
    //date of birth is not a simple thing.
    // DateOfBirth() returns date_of_birth,
    // which is assigned at cell birth (i.e. in the function CellBirth() ) from owner->Time()
    // owner is dish, and in dish there is the function Time()
    // which returns the dish variable thetime
    // which is increased in the function AmoebaeMove2() in ca.cpp
    // ... and all this seems to work fine. except here...
    ofs << icell->DateOfBirth() << " "; // date of birth
    for( auto x: icell->getJkey() ) ofs<<x; //key
    ofs << " ";
    for( auto x: icell->getJlock() ) ofs<<x; //lock
    ofs << " ";

    ofs << icell->mu << " ";
    //ofs << icell->getXvec()<<" "<< icell->getYvec()<<" ";
    ofs << icell->particles << " ";

    //recalculated here because just-born cells do not have defined values
    // for maintenance_fraction or extprotexpress_fraction
    ofs << icell->CalculateMaintenance_or_ExtProtExpr_Fraction(icell->k_mf_0,
                                                               icell->k_mf_A,
                                                               icell->k_mf_P,
                                                               icell->k_mf_C) << " ";

    ofs << icell->k_mf_0 << " ";
    ofs << icell->k_mf_A << " ";
    ofs << icell->k_mf_P << " ";
    ofs << icell->k_mf_C << " ";

    ofs << icell->CalculateMaintenance_or_ExtProtExpr_Fraction(icell->k_ext_0,
                                                               icell->k_ext_A,
                                                               icell->k_ext_P,
                                                               icell->k_ext_C) << " ";
    ofs << icell->k_ext_0 << " ";
    ofs << icell->k_ext_A << " ";
    ofs << icell->k_ext_P << " ";
    ofs << icell->k_ext_C << " ";

    ofs << icell->CalculateMaintenance_or_ExtProtExpr_Fraction(icell->k_chem_0,
                                                               icell->k_chem_A,
                                                               icell->k_chem_P,
                                                               icell->k_chem_C) << " ";
    ofs << icell->k_chem_0 << " ";
    ofs << icell->k_chem_A << " ";
    ofs << icell->k_chem_P << " ";
    ofs << icell->k_chem_C << " ";


    // YOU SHOULD KEEP THIS AT THE LAST, because it's not constant
    for( auto i: icell->neighbours){
      int thisj=icell->getVJ()[ i.first ];
      ofs << cell[i.first].getTau() << " " << thisj << " "; // get J val with neighbors
    }

    ofs << endl;
  }

  //prey=CountPreys();
  //pred=CountPredators();
  //save
  //ofs << Time << " "<< prey << " " << pred << endl;
  ofs.flush();
  ofs.close();
  return (prey+pred);
}

void Dish::MakeBackup(int Time){
  std::ofstream ofs;

  //filename, c++11 strings are concatenated by summing them
  // but because we live in the middle age we are going to use sprintf;
  char filename[300];
  sprintf(filename,"%s/backup_t%09d.txt",par.backupdir,Time);

  ofs.open ( filename , std::ofstream::out | std::ofstream::app);
  ofs<<Time<<endl;
  ofs<<Food->GetPeakx()<<" "<<Food->GetPeaky()<<endl;
  //each cell's variables are written on a single line.
  //don't store area or neighbours: can be inferred from the stored plane
  for(auto c: cell){
    if(c.sigma==0) continue;
    ofs<<c.sigma<<" "<< c.tau<<" "<< c.alive<<" "<< c.tvecx<<" "<<c.tvecy<<" "<< c.prevx<<" "<< c.prevy<<" "
    << c.persdur<<" "<<c.perstime<<" "<<c.mu<<" "<<c.chemmu<<" "<<c.chemvecx<<" "<<c.chemvecy<<" "<< c.target_area<<" "<< c.half_div_area<<" "<< c.eatprob<<" "<<c.particles<<" "<< c.growth<<" "
    <<c.k_mf_0<<" "<<c.k_mf_A<<" "<<c.k_mf_P<<" "<<c.k_mf_C<<" "<<c.k_ext_0<<" "<<c.k_ext_A<<" "<<c.k_ext_P<<" "<<c.k_ext_C<<" ";
    for( auto x: c.jkey ) ofs<<x; //key
    ofs << " ";
    for( auto x: c.jlock ) ofs<<x; //lock
    ofs << endl;
  }

  // particle plane
  // ca plane
  // posix files can be at most 2048 characters long on this laptop (LINE_MAX)
  // so a big field cannot be represented in );an obvious way,
  // however we can save some lines if we print a few sigmas in the same line
  //(anyway, field size is already specified above)... but howmany?
  // well, conceivably I will never run anything larger than in the thousands by thousands pixels
  // and so population size should be less than 10^5 <- so 5 characters+1 space,
  // so I can print 100 sigmas, which should be max 6*100=600 << 2048
  // ok!
  int counter=0;
  int maxperline=100;
  ofs<<endl;
  for(int x=1; x<par.sizex-1; x++){
    for(int y=1; y<par.sizey-1; y++){
      int isigma=CPM->Sigma(x,y);
      ofs<<isigma<<" ";
      counter++;
      if(counter>=maxperline){
        ofs<<endl;
        counter=0;
      }
    }
  }
  ofs<<endl;

  counter=0;
  ofs<<endl;
  for(int x=1; x<par.sizex-1; x++){
    for(int y=1; y<par.sizey-1; y++){
      int isigma=Food->Sigma(x,y);
      ofs<<isigma<<" ";
      counter++;
      if(counter>=maxperline){
        ofs<<endl;
        counter=0;
      }
    }
  }
  ofs<<endl;
}

int Dish::ReadBackup(char *filename){
  //for file reading
  std::ifstream ifs;
  string line;
  Cell *rc;
  //return value: time at which the simulation restarts
  int starttime;

  string jkey, jlock; //read them first as a string, then put into vector?
  int pos;

  ifs.open( filename , std::ifstream::in );

 if (ifs.is_open()){
   //first read the time stamp
   getline(ifs, line);
   stringstream strstr(line);
   strstr >> starttime;

   //second read peakx and peaky
   getline(ifs, line);
   stringstream strstr2(line);
   int px,py;
   strstr2>> px>> py;
   Food->SetPeakx(px); Food->SetPeaky(py);

   //now read all the cell variables
   getline(ifs, line);
   while (line.length()){
     rc=new Cell(*this); //temporary pointer to cell object
     stringstream strstr(line);
     //read the straightforward cell variables from the line
     strstr>>rc->sigma>>rc->tau>>rc->alive>>rc->tvecx>>rc->tvecy>>rc->prevx>>rc->prevy>>rc->persdur
     >>rc->perstime>>rc->mu>>rc->chemmu>>rc->chemvecx>>rc->chemvecy>>rc->target_area>>rc->half_div_area>>rc->eatprob>>rc->particles>>rc->growth
     >>rc->k_mf_0>>rc->k_mf_A>>rc->k_mf_P>>rc->k_mf_C>>rc->k_ext_0>>rc->k_ext_A>>rc->k_ext_P>>rc->k_ext_C>>jkey>>jlock;
     //read the key and lock into the cell
     for (char& c : jkey){
       rc->jkey.push_back(c - '0');
     }
     for (char& c : jlock){
       rc->jlock.push_back(c - '0');
     }
     rc->meanx=0.5*par.sizex;
     rc->meany=0.5*par.sizey;
     cell.push_back(*rc);

     //read the next line
     getline(ifs, line);
   }

   //now read the planes
   //hopefully, plane allocation already happened during dish creation
   //this does assume size parameters match, so be careful!
   getline(ifs, line);
   while (line.length()){
    stringstream ss(line);
    while(ss>>pos){
      int wrong=CPM->SetNextSigma(pos);
      if (wrong){
        cerr<<"ReadBackup error in reading CA plane. more values than fit in plane."<<endl;
        exit(1);
      }
    }
    getline(ifs, line);
  }

  //read the food intplane
  getline(ifs, line);
  while (line.length()){
   stringstream ss(line);
   while(ss>>pos){
    int wrong=Food->SetNextVal(pos);
    if (wrong){
      cerr<<"ReadBackup error in reading Food plane. more values than fit in plane."<<endl;
      exit(1);
    }
   }
   getline(ifs, line);
  }
 }
 else{
   cerr<<"ReadBackup error: could not open file. exiting..."<<endl;
   exit(1);
 }


 return starttime;

}

void Dish::SetCellOwner(Cell &which_cell){
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

  cerr<<"Sandro's note:"<<endl;
  cerr<<"Function Dish::MeasureChemConcentrations looks tricky, and compiler complains"<<endl;
  cerr<<"Do not use this unless you know what is happening"<<endl;

  // clear chemical concentrations
  for (vector<Cell>::iterator c=cell.begin(); c!=cell.end(); c++) {
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


    for (vector<Cell>::iterator c=cell.begin(); c!=cell.end(); c++) {
      for (int ch=0;ch<par.n_chem;ch++)
	    c->chem[ch]/=(double)c->Area();
    }

}

int Dish::SizeX(void) { return CPM->SizeX(); }
int Dish::SizeY(void) { return CPM->SizeY(); }
