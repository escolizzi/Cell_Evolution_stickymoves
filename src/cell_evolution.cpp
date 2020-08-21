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
    //CPM->PlaceOneCellsAtXY(par.sizex/2,par.sizey/2., par.size_init_cells, 1);
    //CPM->PlaceOneCellsAtXY(par.sizex/4,par.sizey/4, par.size_init_cells, 2);
    
    if (! strlen(par.backupfile)) {

      //THIS IS TO USE FOR NORMAL INITIALISATION
      //CPM->PlaceCellsRandomly(par.n_init_cells,par.size_init_cells);
      CPM->PlaceCellsOrderly(par.n_init_cells,par.size_init_cells);
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

      for(auto &c: cell) c.SetTargetArea(par.target_area); //sets target area because in dividecells the new target area = area

      // for(auto &c: cell) c.SetTargetArea(par.target_area); //sets target area because in dividecells the new target area = area

      //PrintContactList();

        //Set function pointer for food update, depending on parameters
        Food->InitIncreaseVal(CPM); //a pointer to CPM is an argument to InitIncreaseVal
                                     // but NOT to IncreaseVal if it points to IncreaseValEverywhere

      // Initialises food plane
      // for(int i=0;i<par.sizex;i++)
      //   for(int j=0;j<par.sizey;j++)
      //     Food->addtoValue(i,j,par.initial_food_amount);  //add initial amount of food for preys

      Food->IncreaseVal(*(Food));
      //cout<<"Hello bla 3"<<endl;
      // exit(1);
      for(int init_time=0;init_time<10;init_time++){
      //   // cerr<<"Init Time: "<<init_time<<endl;
      //   // for(auto c: cell){
      //   //   if(c.AliveP()){
      //   //     printf(" Sigma %d, weight_for_chemotaxis: %.15f\n", c.Sigma(), cell[c.Sigma()].weight_for_chemotaxis);
      //   //   }
      //   //   else
      //   //     printf(" Cell with sigma %d is dead\n", c.Sigma());
      //   // }

        CPM->AmoebaeMove2(PDEfield);  //this changes neighs
      }
      InitCellMigration();

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



  // std::cerr << "howmany cells? "<< cell.size() << '\n';
  // for(auto c: cell){
  //   if(c.AliveP()){
  //     printf("Sigma %d, weight_for_chemotaxis: %.15f\n", c.Sigma(), cell[c.Sigma()].weight_for_chemotaxis);
  //   }
  //   else
  //     printf("Cell with sigma %d is dead\n", c.Sigma());
  // }
  // std::cerr << "How can it be that weight_for_chemotaxis is different between here and inside the function?" << '\n';
  // it isn't, but there is something weird-
  // exit(1);
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

    // **************************************************** //
    // WE NOW CHANGE FOOD BELOW - SEE FUNCTION CheckWhoMadeit
    // dish->Food->IncreaseVal(*(dish->Food)); // SCALED
    // *************************************************** //


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
      // dish->CellsEat(); // SCALED // HERE MAX PARTICLES IS DEFINED, should be a parameter


      dish->CellsEat2();



      //dish->Predate(); //ALREADY SCALED //this does not changes neighs, only target areas!!!

      dish->UpdateCellParameters(i); // SCALED//this changes neighs (via DivideCells)
      //dish->CellGrowthAndDivision2(); // SCALED//this changes neighs (via DivideCells)

      //Recalculate the all vs. all J table.
      //this can be optimised by having some intelligent return flags from dish->CellGrowthAndDivision2();
      // for now it's every one vs everyone all the times.
   // }


    
    // if(i>100){
     dish->CellMigration();//updates persistence time and targetvectors
    // }

    //dish->CPM->AmoebaeMove(dish->PDEfield);  //this changes neighs
    dish->CPM->AmoebaeMove2(dish->PDEfield);  //this changes neighs
    //cout <<"hello2"<<endl;
    //cerr<<"Hello 1"<<endl;
    dish->UpdateNeighDuration();

    //dish->Food->IncreaseVal(*(dish->Food));

    if( i%25 == 0){
      if(par.evolsim){
        if(par.season_experiment){
          if(i>0 && i%par.season_duration==0){
            //reproduce people based on fitness criterion
            //remove random cells until popsize is back to normal
            //reset food and gradient
            std::cerr << "Time = "<<i << '\n';
            std::cerr << "End of season: there are "<< dish->CountCells() <<" cells" << '\n';
            dish->ReproduceEndOfSeason();
            std::cerr << "After reproduction there are "<< dish->CountCells() <<" cells" << '\n';
            dish->RemoveCellsUntilPopIs(par.popsize);
            std::cerr << "After remove there are "<< dish->CountCells() <<" cells" << '\n';

            dish->Food->IncreaseVal(*(dish->Food)); //this has to be last thing to do here
            std::cout << "End of season: Gradient switching at time (+/- 25 MCS) = "<< i << '\n';
          }
        }else{
          if( dish->CheckWhoMadeitRadial() ){
            //reset food
            // clone them with mutations
            // wipe out the previous pop
            // reseed
            //reset whomadeit vector
            dish->RemoveWhoDidNotMakeIt(); //remove those that did not makeit
            dish->ReproduceWhoMadeIt3(); //reproduction
            dish->ClearWhoMadeItSet(); //zeros the who_made_it set,
                                     // zero the particles eaten
            dish->Food->IncreaseVal(*(dish->Food)); //this has to be last thing to do here
                                                  // because we do some AmoebaeMove2 steps in
                                                  // ReproduceWhoMadeIt2 to let cells grow a little
                                                  // but we don't want this to go along the new gradient
                                                  // which would be unfair.
           std::cout << "Gradient switching at time (+/- 25 MCS) = "<< i << '\n';
         }
       }
      }else{
        //not evolutionary simulation
        if( ((strcmp(par.food_influx_location,"boundarygradient") == 0) && dish->CheckWhoMadeitLinear() ) || 
            ((strcmp(par.food_influx_location,"specified_experiment") == 0) && dish->CheckWhoMadeitRadial() )){
          //for printing switching times
          //write switching time to file
          static char timename[300];
          sprintf(timename,"%s/finaltime.txt",par.datadir);
          static ofstream myfile(timename, ios::out | ios::app);
          myfile << i << endl;
          myfile.close();
          exit(0);
        }
      }
    }
      
      
    //BY THE WAY THIS IS HOW YOU CALLED CELL FROM HERE
    //cout<<i<<" "<<dish->getCell(1).getXpos()<<" "<<dish->getCell(1).getYpos()<<endl;

    // if( i%25 == 0){
    //   cerr<<"by time: "<<i<<" there are so many cells: "<<dish->CountCells()<<endl;
    // }

    // if(i%1000==0 ) {
    //   cerr<<"Time: "<<i<<endl;
    //   dish->PrintCellParticles();
    // }

    // TO SCREEN
    // UNUSED
    if (par.graphics && !(i%par.storage_stride)) {

      BeginScene();
      ClearImage();
      if(par.readcolortable){
      //  dish->Plot(this,1);
        dish->Plot(this,2);
      }
      else{
        dish->Plot(this,0);
      }
      //dish->Food->Plot(this, dish->CPM);
      //char title[400];
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      EndScene();
      info->Menu();
    }
//cout <<"hello4"<<endl;
    // TO FILE FOR MOVIE
    if (par.store && !(i%par.storage_stride)) {
      if(par.readcolortable){
        char fname[300];
        sprintf(fname,"%s/angle%09d.png",par.datadir,i);
        BeginScene(); //this is an empty function for X11
        ClearImage(); //
        dish->Plot(this,1); //everything contained here
        EndScene();
        Write(fname);
        sprintf(fname,"%s/order%09d.png",par.datadir,i);
        BeginScene(); //this is an empty function for X11
        ClearImage(); //
        dish->Plot(this,2); //everything contained here
      //dish->Food->Plot(this,dish->CPM); //will this work?  YES !!!
        EndScene();
        Write(fname); //FIXED SO THAT CODE AND IMAGE MATCH!
      }else{
        char fname[300];
        sprintf(fname,"%s/tau%09d.png",par.datadir,i);
        // BeginScene(); //this is an empty function for X11
        ClearImage(); //
        
        //test
        // Point(1, 2*10,2*par.sizey/2);
        // Point(1, 2*10+1,2*par.sizey/2);
        // Point(1, 2*10,2*par.sizey/2+1);
        // Point(1, 2*10+1,2*par.sizey/2+1);
        dish->Plot(this,0); // this is g //everything contained here
        EndScene();
        Write(fname);
      }
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

  // exit(1);
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

    cerr<<endl<<"Warning, this version is ***NOT*** suitable for pde field!!!"<<endl;
    //Depends on this: AddSiteToMoments (and Remove), FindCellDirections2, etc...
    cerr<<endl<<"WARNING, use wrapped boundaries if cells are A LOT smaller than sizex and sizey"<<endl<<endl;
    cerr<<endl<<"WARNING: DO NOT EVOLVE CHEMMU, or if you do, change the replication function (where it is always reset to init_chemmu)"<<endl<<endl;
    
    // cerr<<endl<<"MEGA WARNING !!! \n THIS VERSION IS SUPER HACKED AND CRAMMED, \n DO NOT USE FOR ANYTHING ELSE THAN CAMLEY'S CHEMOTAX + CIL"<<endl<<endl;
    

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
    cerr <<"Using X11 graphics (batch mode). sizex and y are "<<par.sizex<<" "<< par.sizey <<endl;
    X11Graphics g(par.sizex*2,par.sizey*2);
    int t;

    for(t=0;t<=par.mcs;t++){
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
