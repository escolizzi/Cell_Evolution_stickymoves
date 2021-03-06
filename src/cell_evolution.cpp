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
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <math.h>
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "info.h"
#include "parameter.h"
#include "sqr.h"
#include "output.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif


//NOTE: The bookkeeping for cell contacts is very extensive:
// When cells are initially placed (call dish->InitContactLength afterwards)
// When cells divide (in cpm->dividecells)
// When cells are killed and removed (cpm->removecells)
// During CPM updates (cpm->convertspin)
// I added pieces of code to take care of this in the various applicable functions
// We may want to add a parameter to make these parts optional in case we don't need it --it's a bit more costly

using namespace std;

INIT {
  
  
  
  try {
    
    // Define initial distribution of cells
    //CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);

    // THIS IS JUST FOR EXPERIMENTS
    //CPM->PlaceOneCellsAtXY(par.sizex/2,par.sizey/2, par.size_init_cells, 1);
    //CPM->PlaceOneCellsAtXY(par.sizex/4,par.sizey/4, par.size_init_cells, 2);
    //CPM->PlaceOneCellsAtXY(3*par.sizex/4,3*par.sizey/4, par.size_init_cells, 3);
    //CPM->PlaceOneCellsAtXY(par.sizex/4,3*par.sizey/4, par.size_init_cells, 4);
    //CPM->PlaceOneCellsAtXY(3*par.sizex/4,par.sizey/4, par.size_init_cells, 5);
    //CPM->PlaceOneCellsAtXY(par.sizex-3-(int)(sqrt(par.size_init_cells/3.14)),par.sizey/3, par.size_init_cells, 1);
    //CPM->PlaceOneCellsAtXY(par.sizey/3 , par.sizex-3-(int)(sqrt(par.size_init_cells/3.14)), par.size_init_cells , 2);
    //CPM->PlaceOneCellsAtXY((int)(sqrt(par.size_init_cells/3.14))+3,(2/3.)*par.sizey, par.size_init_cells, 3);
    //CPM->PlaceOneCellsAtXY((2/3.)*par.sizey, (int)(sqrt(par.size_init_cells/3.14))+3, par.size_init_cells, 4);
    
    
    if (! strlen(par.backupfile)) {
      //THIS IS TO USE FOR NORMAL INITIALISATION
      CPM->PlaceCellsRandomly(par.n_init_cells,par.size_init_cells);
      CPM->ConstructInitCells(*this); //within an object, 'this' is the object itself

      // Assign a random type to each of the cells, i.e. PREYS and PREDATORS
      CPM->SetRandomTypes();
      //cerr<<"Hello bla 0"<<endl;
      //Initialise key-lock pairs - we do it after types, because we need the information
      InitKeyLock();
     // cerr<<"Hello bla 1"<<endl;
      //Initialise vector of J values for each cell
      InitVectorJ();
      //cerr<<"Hello bla 2"<<endl;
      //Initialise the contactlength bookkeeping now that the cells are placed
      // at this stage, cells are only surrounded by medium
      InitContactLength();  // see dish.cpp - you don't need dish->InitContactLength because this part IS in dish
      //cerr<<"Hello bla 2.5"<<endl;
      InitMaintenanceFraction();
      
      // If we have only one big cell and divide it a few times
      // we start with a nice initial clump of cells.
      //
      //The behavior can be changed in the parameter file using
      //parameters n_init_cells, size_init_cells and divisions
      for(int howmanydivisions=0;howmanydivisions<par.divisions;howmanydivisions++){
        vector<int> sigma_newcells = CPM->DivideCells();
        UpdateVectorJ(sigma_newcells);
        cerr<<"dividing again: "<<howmanydivisions<<endl;
      }
      InitTargetArea();

      //PrintContactList();

      //Set function pointer for food update, depending on parameters
      Food->InitIncreaseVal(CPM); //a pointer to CPM is an argument to InitIncreaseVal
                                   // but NOT to IncreaseVal if it points to IncreaseValEverywhere

      // Initialises food plane
      for(int i=0;i<par.sizex;i++)
        for(int j=0;j<par.sizey;j++)
          Food->addtoValue(i,j,par.initial_food_amount);  //add initial amount of food for preys

       //cout<<"Hello bla 3"<<endl;
       par.starttime=0;
    }
    else {
      par.starttime=ReadBackup(par.backupfile);
      InitContactLength();
      InitVectorJ();
      Food->InitIncreaseVal(CPM); 

    }
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }  

}

TIMESTEP {

  try {
    static Dish *dish=new Dish(); //here ca planes and cells are constructed
    static Info *info=new Info(*dish, *this);
    static int i=par.starttime; //starttime is set in Dish. Not the prettiest solution, but let's hope it works.
    
    //cout << "running... "<< i<<endl;
    if( !(i%100000) ) cerr<<"TIME: "<<i<<endl;

//     cerr<<"target areas before step"<<endl;
//     for(auto c: cell){
//       cerr<<c.TargetArea()<<endl;
//     }


    //cerr<<"TIME: "<<i<<endl;
    //if(i==50) exit(1);
    //add food to plane
    //if(!(i%10))
    //if(!(i%par.scaling_cell_to_ca_time))
    //{
      //dish->Food->IncreaseValIfEmpty(dish->CPM);

    // TIME SCALING IS DONE INSIDE FUNCTIONS
  //  cout <<"hello1"<<endl;

    dish->Food->IncreaseVal(*(dish->Food)); // SCALED

    
//       // testing //
//
//       // Initialises food plane
//       for(int i=0;i<par.sizex;i++)
//         for(int j=0;j<par.sizey;j++)
//           if
//       exit(1);
//
//       //   testing    //
//
      dish->CellsEat(); // SCALED // HERE MAX PARTICLES IS DEFINED, should be a parameter
      
      dish->Predate(); //ALREADY SCALED //this does not changes neighs, only target areas!!!
          
      dish->CellGrowthAndDivision2(); // SCALED//this changes neighs (via DivideCells)
      //Recalculate the all vs. all J table.
      //this can be optimised by having some intelligent return flags from dish->CellGrowthAndDivision2();
      // for now it's every one vs everyone all the times.
   // }

     
    //deal with cell migration
    if(i==100){
      dish->InitCellMigration();
    }

    if(i>100){
     dish->CellMigration();//updates persistence time and targetvectors
    }
    
    

    //dish->CPM->AmoebaeMove(dish->PDEfield);  //this changes neighs
    dish->CPM->AmoebaeMove2(dish->PDEfield);  //this changes neighs
    //cout <<"hello2"<<endl;
    //cerr<<"Hello 1"<<endl;
    dish->UpdateNeighDuration();
  //  cout <<"hello3"<<endl;
    //BY THE WAY THIS IS HOW YOU CALLED CELL FROM HERE
    //cout<<i<<" "<<dish->getCell(1).getXpos()<<" "<<dish->getCell(1).getYpos()<<endl;
    
    
    // if(i%1000==0 ) {
    //   cerr<<"Time: "<<i<<endl;
    //   dish->PrintCellParticles();
    // }
    
    // TO SCREEN
    // UNUSED
    if (par.graphics && !(i%par.storage_stride)) {

      BeginScene();
      ClearImage();
      dish->Plot(this);
      dish->Food->Plot(this, dish->CPM);
      //char title[400];
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      EndScene();
      info->Menu();
    }
//cout <<"hello4"<<endl;
    // TO FILE FOR MOVIE
    if (par.store && !(i%par.storage_stride)) {
      char fname[300];
      sprintf(fname,"%s/ext%07d.png",par.datadir,i);
      BeginScene(); //this is an empty function for X11
      ClearImage(); //
      dish->Plot(this); //everything contained here
      //dish->Food->Plot(this,dish->CPM); //will this work?  YES !!!
      EndScene();
      Write(fname);
    }
  //  cout <<"hello5"<<endl;
    //exit(1);
    // TO FILE FOR TEXT
    if( !(i%par.save_text_file_period) ){
      int popsize=dish->SaveData(i); //saves data to text
      if( 0 == popsize ){
        cerr << "Global extinction after"<<i<<"time steps, simulation terminates now" << endl;
        exit(0);
      }
    }
  //  cout <<"hello2"<<endl;
    // TO FILE FOR BACKUP
    if( !(i%par.save_backup_period) ){
      dish->MakeBackup(i); //saves all permanent data
    }

    i++;
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
}

int PDE::MapColour(double val) {

  return (((int)((val/((val)+1.))*100))%100)+155;
}

//////////////////////////////
// ------------------------ //
// ---       MAIN       --- //
// ------------------------ //
//////////////////////////////
int main(int argc, char *argv[]) {


  try {

#ifdef QTGRAPHICS
    //QCoreApplication a(argc, argv);
    QApplication a(argc, argv);
    QTimer g;
    //QApplication a2(argc, argv);
#endif

    par.Read(argv[1]); // Read parameters from file

    //command line arguments overwrite whatever is in the parameter file
    if(argc>2){
      int exit_valarg = par.ReadArguments(argc,argv);
      if( 0 != exit_valarg ){
        par.PrintWelcomeStatement(); //see parameter.h
        exit(1);
      }
    }

    //Creates a rule for J val with medium that is going to be used in following function
    par.CreateRule(par.Jmed_rule_input);
    // Open file where initial key locks are specified, and assigns them to cells
    par.Read_KeyLock_list_fromfile(par.keylock_list_filename);

    cerr<<endl<<"WARNING, THIS VERSION IS ***NOT*** SUITABLE FOR PDE FIELD!!!"<<endl;
    //Depends on this: AddSiteToMoments (and Remove), FindCellDirections2, etc...
    cerr<<endl<<"WARNING, use wrapped boundaries if cells are A LOT smaller than sizex and sizey"<<endl<<endl;

    //check if directory for movies exists, create it if not, exit otherwise
    DoesDirExistsIfNotMakeit(par.datadir);  //see output.cpp
    DoesDirExistsIfNotMakeit(par.backupdir);  //see output.cpp

    //check if data file exists, if not exit
    if(FileExistsP(par.datafile)){
      cerr<<"File "<< par.datafile << " already exists, simulation not starting" << endl;
      exit(1);
    }

    if(par.periodic_boundaries && par.lambda2>0.){
      cerr<<"main(): Error. Cannot have wrapped periodic boundaries and lambda2>0"<<endl;
      cerr<<"(because I cannot calculate second moment for cells crossing boundaries)"<<endl;
      exit(1);
    }

    Seed(par.rseed);

    //QMainWindow mainwindow w;
#ifdef QTGRAPHICS
    cerr<<"wat is deze? "<<par.readcolortable<<endl;
    //exit(1);
    //QtGraphics g(par.sizex*2,par.sizey*2);

    QImage image(par.sizex*2,par.sizey*2, QImage::Format_ARGB32);
    QPainter painter(&image);
    QPaintDevice *device = painter.device();

    //QtGraphics g2(par.sizex*2,par.sizey*2);
    cerr<<"Hello 1"<<endl;
    //a->setMainWidget( &g );
    //a->connect(&g, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );
    g.connect (&g, SIGNAL(timeout()), &a, SLOT(quit()));
    cerr<<"Hello 1.1"<<endl;

    //a2.setMainWidget( &g2 );
    //a2.connect(&g2, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );

    if (par.graphics)
    {
      //g.show();
      //   // g2.show();
      cerr<<"Hello 2"<<endl;
    }
    a.exec();
    //a2.exec();
    cerr<<"Hello 3"<<endl;
#else
    cerr <<"Using X11 graphics (batch mode)..."<<endl;
    X11Graphics g(par.sizex*2,par.sizey*2);
    int t;

    for(t=0;t<par.mcs;t++){
      //cerr<<"Time: "<<t<<endl;
      g.TimeStep();

    }
#endif

  }catch(const char* error){
    std::cerr << error << "\n";
    exit(1);
  }
  catch(...) {
    std::cerr << "An unknown exception was caught\n";
  }
  return 0;
}
