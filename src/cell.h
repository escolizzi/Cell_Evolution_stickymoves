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
#ifndef _CELL_HH_
#define _CELL_HH_

#include "parameter.h"
//#define EMPTY -1
#include <math.h>
#include <map>
#include "random.h"

#define PREY 1
#define PREDATOR 2

extern Parameter par;
class Dish;

class Cell
{

  friend class Dish;
  friend class CellularPotts;
  friend class Info;

public:



  /*! \brief Constructor to insert a cell into Dish "who"

  Used to add a new Cell to the dish: new Cell(dish,
  celtype).
  */
  Cell(const Dish &who, int settau=1, int setrecycledsigma=-1) {

    owner=&who;
    ConstructorBody(settau,setrecycledsigma);
  }

  Cell(void) {
    if(par.n_chem){
      chem = new double[par.n_chem];
    }
  };

  ~Cell(void);


  //! Default copy constructor.
  // Copies a cell into a new one
  Cell(const Cell &src) {
    //cout<<"Tomato"<<endl;
    // make an exact copy (for internal use)
    sigma=src.sigma;
    amount++;
    area=src.area;
    target_area=src.target_area;
    half_div_area=src.half_div_area;
    length=src.length;
    target_area = src.target_area;
    target_length=src.target_length;
    growth_threshold=src.growth_threshold;
    mother=src.mother;
    daughter=src.daughter;
    times_divided=src.times_divided;
    date_of_birth=src.date_of_birth;
    //cerr<<"this?: "<<date_of_birth<<endl;
    colour_of_birth=src.colour_of_birth;
    tau=src.tau;
    alive=src.alive;
    v[0]=src.v[0];
    v[1]=src.v[1];
    n_copies=src.n_copies;
    sum_x=src.sum_x;
    sum_y=src.sum_y;
    sum_xx=src.sum_xx;
    sum_yy=src.sum_yy;
    sum_xy=src.sum_xy;

    meanx=src.meanx;
    meany=src.meany;
    tvecx=src.tvecx;
    tvecy=src.tvecy;
    prevx=src.prevx;
    prevy=src.prevy;
    persdur=src.persdur;
    perstime=src.perstime;
    mu=src.mu;

    chemmu=src.chemmu;
    chemvecx=src.chemvecx;
    chemvecy=src.chemvecy;
    chemveclength=src.chemveclength;

    owner=src.owner;
    particles=src.particles;
    eatprob=src.eatprob;
    growth=src.growth;

    //this is copied for if there is a first time step
    // where things depend on it and CellGrowthAndDivision2 is not called yet
    maintenance_fraction = src.maintenance_fraction;
    k_mf_0 = src.k_mf_0;
    k_mf_A = src.k_mf_A;
    k_mf_P = src.k_mf_P;
    k_mf_C = src.k_mf_C;

    extprotexpress_fraction = src.extprotexpress_fraction;
    k_ext_0 = src.k_ext_0;
    k_ext_A = src.k_ext_A;
    k_ext_P = src.k_ext_P;
    k_ext_C = src.k_ext_C;

    weight_for_chemotaxis = src.weight_for_chemotaxis;
    k_chem_0=src.k_chem_0;
    k_chem_A=src.k_chem_A;
    k_chem_P=src.k_chem_P;
    k_chem_C=src.k_chem_C;

    jlock = src.jlock;
    jkey = src.jkey;
    vJ = src.vJ;

    neighbours=src.neighbours;

    if(par.n_chem){
      chem = new double[par.n_chem];
      for (int ch=0;ch<par.n_chem;ch++)
        chem[ch]=src.chem[ch];
    }


  }

  /*! \brief Add a new cell to the dish.

     Call it as: new Cell(parent, true); mother will be modified for
     ancestry administration!

     \param settau
     Cell type of daughter cell.
  */
  void CellBirth(Cell &mother);

  /*! \brief Assignment operator.

  Called if one cell is assigned to another. Remember to change both
  assignment operator and copy constructor when adding new attributes
  to Cell.
  */
  inline Cell &operator=(const Cell &src) {
    //cout<<"Potato"<<endl;
    colour=src.colour;
    alive=src.alive;
    sigma=src.sigma;
    area=src.area;
    tau=src.tau;
    target_area=src.target_area;
    half_div_area=src.half_div_area;
    date_of_birth = src.date_of_birth;

    v[0]=src.v[0];
    v[1]=src.v[1];
    n_copies=src.n_copies;

    meanx=src.meanx;
    meany=src.meany;
    tvecx=src.tvecx;
    tvecy=src.tvecy;
    prevx=src.prevx;
    prevy=src.prevy;

    persdur=src.persdur;
    perstime=src.perstime;
    mu=src.mu;

    chemmu=src.chemmu;
    chemvecx=src.chemvecx;
    chemvecy=src.chemvecy;
    chemveclength=src.chemveclength;

    sum_x=src.sum_x;
    sum_y=src.sum_y;
    sum_xx=src.sum_xx;
    sum_yy=src.sum_yy;
    sum_xy=src.sum_xy;

    length=src.length;
    target_area = src.target_area;
    target_length=src.target_length;
    amount++;
    owner=src.owner;

    particles=src.particles;
    eatprob=src.eatprob;
    growth = src.growth;
    neighbours=src.neighbours;

    maintenance_fraction = src.maintenance_fraction;
    k_mf_0 = src.k_mf_0;
    k_mf_A = src.k_mf_A;
    k_mf_P = src.k_mf_P;
    k_mf_C = src.k_mf_C;

    extprotexpress_fraction = src.extprotexpress_fraction;
    k_ext_0 = src.k_ext_0;
    k_ext_A = src.k_ext_A;
    k_ext_P = src.k_ext_P;
    k_ext_C = src.k_ext_C;

    weight_for_chemotaxis = src.weight_for_chemotaxis;
    k_chem_0=src.k_chem_0;
    k_chem_A=src.k_chem_A;
    k_chem_P=src.k_chem_P;
    k_chem_C=src.k_chem_C;

    jlock = src.jlock;
    jkey = src.jkey;
    vJ = src.vJ;

    chem = new double[par.n_chem];
    for (int ch=0;ch<par.n_chem;ch++)
      chem[ch]=src.chem[ch];

    return *this;

  }
  //another overload, this is used for assigning a cell to another
  //but with specified sigma
  // except tghat this also doesn't work
//   inline Cell &operator=(const Cell &src, const int &recycled_sigma) {
//     cout<<"Potato with sigma"<<endl;
//     colour=src.colour;
//     alive=src.alive;
//     sigma=recycled_sigma;
//     area=src.area;
//     tau=src.tau;
//     target_area=src.target_area;
//     v[0]=src.v[0];
//     v[1]=src.v[1];
//     n_copies=src.n_copies;
//
//     sum_x=src.sum_x;
//     sum_y=src.sum_y;
//     sum_xx=src.sum_xx;
//     sum_yy=src.sum_yy;
//     sum_xy=src.sum_xy;
//
//     length=src.length;
//     target_length=src.target_length;
//     amount++;
//     owner=src.owner;
//
//     particles=src.particles;
//     growth = src.growth;
//     neighbours=src.neighbours;
//
//     chem = new double[par.n_chem];
//     for (int ch=0;ch<par.n_chem;ch++)
//       chem[ch]=src.chem[ch];
//
//     return *this;
//
//   }

  /*! \brief Returns false if Cell has apoptosed (vanished). */
  inline bool AliveP(void) const {
    return alive;
  }

  //! Returns the cell colour.
  inline int Colour(void) const {

    /* if (par.dynamicJ)
      return colour;
      else */
      return tau+1;

  };

//this function maps migration vector angle to a colour in radial_colour array (see misc.cpp)
  inline int AngleColour(void) const {
    double ang=atan(tvecy/tvecx);
    if(tvecx<0.000) ang+=M_PI;
    else if(tvecy<0.000) ang+=2*M_PI;
    ang/=2.*M_PI;

    return (int)(ang*254)+2;
  };

  //sets properties of cell
  // for now only growth   ---- THIS STUFF CREATES PROBLEMS!

  //predator is 1? No, prey is 1, predator is 2
  inline void SetCellTypeProperties(void)
  {
    if(tau==PREY){
      //growth = par.growth/2.;
      ;
      half_div_area = par.half_div_area;
    }else if(tau==PREDATOR){
      //growth = 2.*par.growth;

      if(par.half_div_area_2 > 0){
        half_div_area = par.half_div_area_2;
        //cout<<endl<<endl<<"Setting half_div_area to "<<par.half_div_area_2<<" because this guy has tau = "<<tau<<endl;
      }
      else half_div_area = par.half_div_area;

      // THIS IS NOT WHAT YOU WANT!!!
      // YOU WANT TO SET AREA AT WHICH DIVISION HAPPEN, NOT TARGET AREA (which is food dependent)
      //if(par.target_area_2!=-1) target_area = par.target_area_2;
      //else target_area = par.target_area;

    }
    //;
  }

  //Nonsensical object oriented stuff:
  inline vector<int> getJkey(void){
    return jkey;
  }

  inline vector<int> getJlock(void){
    return jlock;
  }

  inline void setJkey(vector<int> setjkey){
    jkey = setjkey;
  }

  inline void setJlock(vector<int> setjlock){
    jlock= setjlock;
  }

  inline vector<int> getVJ(void){
    return vJ;
  }

  inline void setVJ(vector<int> setvj){
    vJ = setvj;
  }

  inline void InitMeanX(double xpos){
    meanx=xpos;
  }
  inline void InitMeanY(double ypos){
    meany=ypos;
  }

  //Return values related to cell position and direction of migration
  inline double getXpos(void){
    return meanx;
  }
  inline double getYpos(void){
    return meany;
  }
  inline double getXvec(void){
    return tvecx;
  }
  inline double getYvec(void){
    return tvecy;
  }
  inline double getChemXvec(void){
    return chemvecx;
  }
  inline double getChemYvec(void){
    return chemvecy;
  }
  inline void setChemVecLength(double length){
    //cout<<"numu: "<<mu<<endl;
    chemveclength=length;
  }
  inline double getChemVecLength(void){
    //cout<<"numu: "<<mu<<endl;
    return chemveclength;
  }
  inline double getChemMu(void){
    //cout<<"numu: "<<mu<<endl;
    return chemmu;
  }
  inline double getMu(void){
    //cout<<"numu: "<<mu<<endl;
    return mu;
  }
  inline void setMu(double initmu){
    mu=initmu;
   // cout<<"initmu: "<<mu<<endl;
  }

  inline void setChemMu(double initweight){
    chemmu=initweight;
   // cout<<"initmu: "<<mu<<endl;
  }
  inline void setPersDur(int dur){
    persdur=dur;
  }
  inline void setPersTime(int time){
    perstime=time;
  }
  //initialise variables for target vector position
  inline void startTarVec(){
    prevx=meanx;
    prevy=meany;
    double pol=RANDOM()*2.*M_PI; //random angle

    //pol=0.;

    tvecx=cos(pol);  // try swapping these around
    tvecy=sin(pol);
  }
  inline void startChemVec(){ //to make sure hamiltonian has something to work with.

    double pol=RANDOM()*2.*M_PI; //random angle

    //pol=0.;

    chemvecx=cos(pol);  // try swapping these around
    chemvecy=sin(pol);
  }
  inline void setChemVec(double xx, double yy){ //to make sure hamiltonian has something to work with.
    chemvecx=xx;  // try swapping these around
    chemvecy=yy;
  }
  //update the target vector with the actual direction of motion
  inline void updateTarVec()
  {

    if((meanx-prevx)*(meanx-prevx)>0.0000001 || (meany-prevy)*(meany-prevy)>0.0000001)
    {
      tvecx=meanx-prevx;
      tvecy=meany-prevy;
      double hyphyp=hypot(tvecx,tvecy);
      tvecx/=hyphyp;
      tvecy/=hyphyp;
      prevx=meanx;
      prevy=meany;
    }
    else
    {
      //just keep moving in the same direction then
      prevx=meanx;
      prevy=meany;
    }
  }

  inline void updatePersTime(void){
   perstime++;
   if(perstime==persdur){
     updateTarVec();
     //cerr<<"persdur="<<persdur<<endl;
     //cerr<<"Cell "<<sigma<<" geupdate, "<<tvecx<<" "<<tvecy<<endl;
     perstime=0;
   }
  }

//   // it uses size_t instead of int to shut up warnings
  inline void setVJ_singleval(int pos, int val){
    size_t unsigned_pos = (size_t)pos;
    if( unsigned_pos >= vJ.size()){
      //cerr<<"pos larger than vJ size: pos = "<<unsigned_pos<<" size = "<< vJ.size() <<endl;
      for(size_t i=vJ.size(); i<unsigned_pos+1 ; i++)
        vJ.push_back(-1);
    }
    vJ[unsigned_pos] = val;
  }

  int MutateKeyAndLock(void);

  double MAXmu=30;
  inline void MutateMu(void){
    mu+= (RANDOM() -0.5)/1.;
    if(mu<0) mu= -mu;
    else if(mu>MAXmu) mu=MAXmu-mu;
  }

  //I am using really the same function twice,
  // no point in writing two functions
  double CalculateMaintenance_or_ExtProtExpr_Fraction(double k0, double kA , double kP, double kC);

  inline double GetExtProtExpress_Fraction(void){
    return extprotexpress_fraction;
  }
  // inline void MutateMaintenanceFraction(void){
  //   maintenance_fraction += (RANDOM() -0.5)/20.;
  //   if(maintenance_fraction<0) maintenance_fraction= -maintenance_fraction;
  //   else if(maintenance_fraction>1.) maintenance_fraction=2.-maintenance_fraction;
  // }

  inline void MutateMaintenanceFractionParameters(void){
    //if(RANDOM() < par.mut_rate){
      k_mf_0 += (RANDOM() -0.5)/10.;
      k_mf_A += (RANDOM() -0.5)/10.;
      k_mf_P += (RANDOM() -0.5)/10.;
      k_mf_C += (RANDOM() -0.5)/10.;
    // }
  }

  inline void MutateExtProtFractionParameters(void){
    // if(RANDOM() < par.mut_rate){
      k_ext_0 += (RANDOM() -0.5)/10.;
      // k_ext_A += (RANDOM() -0.5)/10.;
      k_ext_P += (RANDOM() -0.5)/10.;
      // k_ext_C += (RANDOM() -0.5)/10.;
    // }
  }

  inline void MutateChemotaxisParameters(void){
    // if(RANDOM() < par.mut_rate){
      k_chem_0 += (RANDOM() -0.5)/10.;
      k_chem_A += (RANDOM() -0.5)/10.;
      k_chem_P += (RANDOM() -0.5)/10.;
      k_chem_C += (RANDOM() -0.5)/10.;
    // }
  }

  //! Set cell type of this Cell.
  inline int getHalfDivArea(void) {
    return half_div_area;
  }

  //! Set cell type of this Cell.
  inline void setTau(int settau) {
    tau=settau;
  }

  //! Get cell type of this Cell.
  inline int getTau(void) {
    return tau;
  }
  //! Set color of this cell to new_colour, irrespective of type.
  inline int SetColour(const int new_colour) {
    return colour=new_colour;
  }

  /* \brief Returns the energy between this cell and cell2.

  Called from CellularPotts::DeltaH.
  **/
  int EnergyDifference(const Cell &cell2) const;

  //! Return Cell's actual area.
  inline int Area() const {
    return area;
  }

  //! Return Cell's target area.
  inline int TargetArea() const {
    return target_area;
  }

  /*! \brief Return Cell's target length.

  Length constraint is documented in Merks et al. 2006, Dev. Biol.
  */
  inline double TargetLength() const {
    return target_length;
  }

  //! Set the Cell's target length
  inline double SetTargetLength(double l) {
    return target_length=l;
  }


  //! Debugging function used to print the cell's current inertia tensor (as used for calculations of the length )
  inline void PrintInertia(void) {

    double ixx=(double)sum_xx-(double)sum_x*sum_x/(double)area;
    double iyy=(double)sum_yy-(double)sum_y*sum_y/(double)area;
    double ixy=(double)sum_xy-(double)sum_x*sum_y/(double)area;

    cerr << "ixx = " << ixx << "\n";
    cerr << "iyy = " << iyy << "\n";
    cerr << "ixy = " << ixy << "\n";

  }

  // return the current length
  inline double Length(void) {
    return length;
  }

/*! \brief Clears the table of J's.

This is only important for a
feature called "DynamicJ's", where J-values depend on internal states
of the cells (such as a genetic network; see e.g. Hogeweg et
al. 2000). The current version of TST does not include such functionality.
*/
  static void ClearJ(void);
  double polarvec[9];
  void RenormPolarVec(void);

  /*! \brief Returns the maximum cell identity number in the Dish.
    This would normally be the number of cells in the Dish, although
   the number includes apoptosed cells.
  */
  static inline int MaxSigma() {
    return maxsigma;
  }

  //! Returns the cell's cell identity number.
  inline int Sigma() const {
    return sigma;
  }

  // THIS REALLY DOES NOT WORK, because sigma is protected.
//   // set new sigma at birth
//   inline int SetSigmaIfRecycled(const int recycled_sigma) const {
//     return sigma=recycled_sigma;
//   }

  //! Sets the target area of the cell.
  inline int SetTargetArea(const int new_area) {
    return target_area=new_area;
  }

  //! Sends the current cell into apoptosis
  inline void Apoptose() {
    alive=false;
  }

  //! Decrement the cell's target area by one unit.
  inline int IncrementTargetArea() {
    return ++target_area;
  }

  //! Increment the cell's target area by one unit.
  inline int DecrementTargetArea() {
    return --target_area;
  }

  //! Cell lineage tracking: get the cell's parent
  inline int Mother(void) const { return mother; }

  //! Cell lineage tracking: get the cell's daughter
  inline int Daughter(void) const { return daughter; }

  //! Returns a counter keeping track of the number of divisions
  inline int TimesDivided(void) const { return times_divided; }

  //! Returns Monte Carlo Step (MCS) when this cell originated.
  inline int DateOfBirth(void) const { return date_of_birth; }

  //! Returns the cell type at the time of birth.
  inline int ColourOfBirth(void) const { return colour_of_birth; }

  //! Returns the bond energy J between this cell and cell c2.
  inline int GetJ(const Cell &c2) const {
    return J[sigma][c2.sigma];
  }


  //! Sets bond energy J between cell type t1 and t2 to val
  inline static int SetJ(int t1,int t2, int val) {
    return J[t2][t1]=J[t1][t2]=val;
  }


  // Deal with gradient measurements:

  //! Set the current gradient of the cell to g. Currently not in use.
  inline double* SetGrad(double *g) {
    grad[0]=g[0];
    grad[1]=g[1];
    return grad;
  }

  //! Returns the cell's measured gradient. Currently not in use.
  inline const double* GetGrad(void) const {
    return grad;
  }

  //! Returns the cell's measured gradient. Currently not in use.
  inline const double GradX() const {
    return grad[0];
  }

  //! Returns the cell's measured gradient. Currently not in use.
  inline const double GradY() const {
    return grad[1];
  }

  //! Currently not in use (remove?)
  inline double* AddToGrad(double *g) {
    grad[0]+=g[0];
    grad[1]+=g[1];
    return grad;
  }

  //! Currently not in use (remove?)
  inline void ClearGrad(void) {
    grad[0]=0.;
    grad[1]=0.;
  }

  inline double GetEatProb(void){
    return eatprob;
  }
  inline void SetEatProb(double par_eatprob){
    eatprob=par_eatprob;
  }

  /*! After introducing a new Cell (e.g. with GrowInCell)
    call this function to set the moments and areas right.
  */
  void MeasureCellSize(Cell &c);

  //This is for keeping track of who is neigh, how much contact and for how long
  std::map<int, pair<int,int> >neighbours; //stores neigh as {neighbouring cells(ID) <amount of contact,duration>
  //cell neighbours
  void setNeighbour(int neighbour, int boundarylength, int contactduration);
  void clearNeighbours();
  int returnBoundaryLength(int cell);
  int returnDuration(int cell);
  int SetNeighbourDurationFromMother(int cell, int motherduration);
  int updateNeighbourBoundary(int cell, int boundarymodification);
  int updateNeighbourDuration(int cell, int durationmodification);

private:
  //updated version: read key-lock pairs
  static void ReadKeyLockFromFile(const char *fname);
  /*! \brief Read a table of static Js.
    First line: number of types (including medium)
    Next lines: diagonal matrix, starting with 1 element (0 0)
    ending with n elements */
  static void ReadStaticJTable(const char *fname);

  // used internally by dish in "CellGrowthAndDivision"
  inline int GrowthThreshold(void) const { return growth_threshold; }

  // used internally by class CellularPotts
  inline void CleanMoments(void) {
    sum_x = sum_y = sum_xx = sum_xy = sum_yy = area = target_area  = 0;
  }

  // used internally by class CellularPotts
  // notice that it should be preceded always by IncrementArea()
  //
  // so far: problem is I am not changing sumx, only meanx


  inline double AddSiteToMoments(int x,int y, double new_l=-1. ) {

    //calculate higher moments ONLY if not wrapped boundaries
    length=0.; //set length to zero becaues we don't need it anyway

    // Add a site to the raw moments, then update and return the
    // length of the cell
    //cout<<"area now: "<<area<<"x: "<<x<<"y: "<<y<<endl;
    // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
    // Eventually this function may be used to carry
    // out all necessary adminstration at once
    if(par.periodic_boundaries){
      if(area>0){
        //also notice that this is very related to meanx...

        //double curr_avrg_x = sum_x/((double)(area-1)); //area-1 because we have just incremented area without counting this
        //double curr_avrg_y = sum_y/((double)(area-1));

        //if x is closer to running average when wrapped, we wrap it
        if( (x-meanx)>0 && (x-meanx)>(meanx-(x-(par.sizex-2))) ){
        //if( abs(meanx - x) > abs(meanx - (x - (par.sizex-2))) ){
          //print passing border
          //cerr<<"Astm: meanx ="<<meanx<<", pixel on the right: "<<x<<", we wrap it to "<<x-(par.sizex-2)<<endl;
          x-=(par.sizex-2); //par.sizex -2 because to compare floats with int we assume the pixel begins (on the left, with integer)
          //cerr<<"passb1"<<endl;
        //cerr<<"so sum_x will be: "<<sum_x+x<<endl;
        }
        else if( (meanx-x)>0 && (meanx-x)>(x+(par.sizex-2)-meanx) ){
        //else if( abs(meanx - x) > abs(meanx - (x + (par.sizex-2))) ) {
          //print passing border
          //cerr<<"Astm: meanx ="<<meanx<<", pixel on the left: "<<x<<", we wrap it to "<<x+(par.sizex-2)<<endl;
          x+=(par.sizex-2);
          //cerr<<"so sum_x will be: "<<sum_x+x<<endl;
          //cerr<<"passb2"<<endl;
        }
        //same for y
        if( (y-meany)>0 && (y-meany)> (meany - (y - (par.sizey-2))) ){
          y-=(par.sizey-2);
          //cerr<<"passb3 meany: "<<meany<<endl;
        }
        else if( (meany-y>0) &&  meany-y>  (y + (par.sizey-2)-meany) ){
          y+=(par.sizey-2);
          //cerr<<"passb4"<<endl;
        }


        sum_x+=x;
        //cerr<<"Astm sum_x= "<<sum_x<<endl;
        sum_y+=y;

        //cerr<<"Before add: meanx: "<<meanx<<", meany: "<<meany<<endl;

        //now if meanx or meany are outside borders we wrap them
        // and shift prevx or y accordingly
        meanx = sum_x/((double)(area));
        meany = sum_y/((double)(area));



        //if x is closer to running average when wrapped, we wrap it


        //notice that there is equal sign with if( meanx >= sizex-1 )
        // this means that the right edge of the CA does not exist
        // (it is the first point wrapped)
        if( meanx >= par.sizex-1 ){
          //wrap meanx
          //cerr<<"meanx>par.sizex-1 we wrap meanx, sum_x and prevx"<<endl;
          meanx-=(par.sizex-2);
          //change sumx and prevx
          sum_x-= area*(par.sizex-2); // sumx has to be shifted, as a whole (each pixel by sizex-2, so in total by area*(sizex-2))
          prevx-=(par.sizex-2);
          //cerr<<"passb5"<<endl;
        }else if(meanx<1){
//           cerr<<"meanx<1 we wrap it meanx, sum_x and prevx"<<endl;
          meanx+=(par.sizex-2);
          sum_x += area*(par.sizex-2);
          prevx+=(par.sizex-2);
          //cerr<<"passb6"<<endl;
        }
        //same for y
        if( meany >= par.sizey-1 ){
          meany-=(par.sizey-2);
          sum_y -= area*(par.sizey-2);
          prevy-=(par.sizey-2);
          //cerr<<"passb7"<<endl;
        }else if(meany<1){
          meany+=(par.sizey-2);
          sum_y += area*(par.sizey-2);
          prevy+=(par.sizey-2);
          //cerr<<"passb8 new meany: "<<meany<<endl;
        }

      }else{
        cerr<<"AddSiteToMoments(): Error. How can area be zero after incrementing it?"<<endl;
        exit(1);
      }
      //cerr<<"After add: meanx: "<<meanx<<", meany: "<<meany<<endl;
    }
    //else if not periodic_boundaries
    else{
      //update alll moments
      sum_x+=x;
      sum_y+=y;
      sum_xx+=x*x;
      sum_yy+=y*y;
      sum_xy+=x*y;

      // update length (see appendix. A, Zajac.jtb03), if length is not given
      // NB. 24 NOV 2004. Found mistake in Zajac's paper. See remarks in
      // method "Length(..)".
      if (new_l<0.) {
        length=Length(sum_x,sum_y,sum_xx,sum_yy,sum_xy,area);
      } else {
        length=new_l;
      }

      //mean position of cell
      if(area>0){
        meany=double(sum_y)/double(area);
        meanx=double(sum_x)/double(area);
      }
    }
    return length;
  }

  // used internally by class CellularPotts
  inline double RemoveSiteFromMoments(int x,int y, double new_l=-1.) {

    length=0.;

    if(par.periodic_boundaries){
      if(area>0){
        //if x is closer to running average when wrapped, we wrap it
        if( (x-meanx)>0 && (x-meanx)>(meanx-(x-(par.sizex-2))) ){
        //if( abs(meanx - x) > abs(meanx - (x - (par.sizex-2))) ){
          //print passing border
//           cerr<<"Rstm: meanx ="<<meanx<<", pixel on the right: "<<x<<", we wrap it to "<<x-(par.sizex-2)<<endl;
          x-=(par.sizex-2); //par.sizex -2 because to compare floats with int we assume the pixel begins (on the left, with integer)
          //cerr<<"so sum_x will be: "<<sum_x-x<<endl;
          //cerr<<"passb"<<endl;
        }
        else if( (meanx-x)>0 && (meanx-x)>(x+(par.sizex-2)-meanx) ){
        //else if( abs(meanx - x) > abs(meanx - (x + (par.sizex-2))) ) {
          //print passing border
//           cerr<<"Rstm: meanx ="<<meanx<<", pixel on the left: "<<x<<", we wrap it to "<<x+(par.sizex-2)<<endl;
          x+=(par.sizex-2);
          //cerr<<"so sum_x will be: "<<sum_x-x<<endl;
          //cerr<<"passb"<<endl;
        }
        //same for y
        if( (y-meany)>0 && (y-meany)> (meany - (y - (par.sizey-2))) ){
          y-=(par.sizey-2);
          //cerr<<"passb"<<endl;
        }
        else if( (meany-y>0) &&  meany-y>  (y + (par.sizey-2)-meany) ){
          y+=(par.sizey-2);
          //cerr<<"passb"<<endl;
        }

        sum_x-=x;
//         cerr<<"Rstm sum_x= "<<sum_x<<endl;
        sum_y-=y;

        //cerr<<"Before rem: meanx: "<<meanx<<", meany: "<<meany<<endl;

        //now if meanx or meany are outside borders we wrap them
        // and shift prevx or y accordingly
        meanx = sum_x/((double)(area));
        meany = sum_y/((double)(area));

        //if x is closer to running average when wrapped, we wrap it


        //this is still to correct for -1 or -2
        if( meanx >= par.sizex-1 ){
          //wrap meanx
//           cerr<<"meanx>par.sizex-1 we wrap meanx, sum_x and prevx"<<endl;
          meanx-=(par.sizex-2);
          //change sumx and prevx
          sum_x-= area*(par.sizex-2); // sumx has to be shifted, as a whole (each pixel by sizex-2, so in total by area*(sizex-2))
          prevx-=(par.sizex-2);
          //cerr<<"passb"<<endl;
        }else if(meanx<1){
//           cerr<<"meanx<1 we wrap it meanx, sum_x and prevx"<<endl;
          meanx+=(par.sizex-2);
          sum_x += area*(par.sizex-2);
          prevx+=(par.sizex-2);
          //cerr<<"passb"<<endl;
        }
        //same for y
        if( meany >= par.sizey-1 ){
          meany-=(par.sizey-2);
          sum_y -= area*(par.sizey-2);
          prevy-=(par.sizey-2);
          //cerr<<"passb"<<endl;
        }else if(meany<1){
          meany+=(par.sizey-2);
          sum_y += area*(par.sizey-2);
          prevy+=(par.sizey-2);
          //cerr<<"passb"<<endl;
        }

      }
      //cerr<<"After rem: meanx: "<<meanx<<", meany: "<<meany<<endl;
    }
    //else if not periodic_boundaries
    else{

      // Remove a site from the raw moments, then update and return the
      // length of the cell

      // sum_x, sum_y, sum_xx, sum_xy and sum_yy are adjusted
      // Eventually this function may be used to carry
      // out all necessary adminstration at once
      sum_x-=x;
      sum_y-=y;



      sum_xx-=x*x;
      sum_yy-=y*y;
      sum_xy-=x*y;

      // update length (see app. A, Zajac.jtb03), if length is not given
      if (new_l<0.) {
        length=Length(sum_x,sum_y,sum_xx,sum_yy,sum_xy,area);
      } else {
        length=new_l;
      }

      if(area>0){
        meany=double(sum_y)/double(area);
        meanx=double(sum_x)/double(area);
      }
    }
    return length;
  }



  //! \brief Calculates the length based on the given inertia tensor
  //components (used internally)
  inline double Length(long int s_x,long int s_y,long int s_xx,
				long int s_yy,long int s_xy,long int n) {

    // inertia tensor (constructed from the raw momenta, see notebook)
    double iyy=(double)s_xx-(double)s_x*s_x/(double)n;
    double ixx=(double)s_yy-(double)s_y*s_y/(double)n;
    double ixy=-(double)s_xy+(double)s_x*s_y/(double)n;

    double rhs1=(ixx+iyy)/2., rhs2=sqrt( (ixx-iyy)*(ixx-iyy)+4*ixy*ixy )/2.;

    double lambda_b=rhs1+rhs2;
    //double lambda_a=rhs1-rhs2;

    // according to Zajac et al. 2003:
    //return 2*sqrt(lambda_b);
    // Grumble, this is not right!!!
    // Must divide by mass!!!!!!

    // see: http://scienceworld.wolfram.com/physics/MomentofInertiaEllipse.html
    //    cerr << "n = " << n << "\n";
    return 4*sqrt(lambda_b/n);

    // 2*sqrt(lambda_b/n) give semimajor axis. We want the length.

  }

  // return the new length that the cell would have
  // if site (x,y) were added.
  // used internally by CellularPotts
  inline double GetNewLengthIfXYWereAdded(int x, int y) {
    double lengthifxywereadded = 0.;
    if(!par.periodic_boundaries)
      lengthifxywereadded = Length(sum_x+x,sum_y+y,sum_xx+x*x,sum_yy+y*y,sum_xy+x*y,area+1);
    return lengthifxywereadded;
  }

  // return the new length that the cell would have
  // if site (x,y) were removed
  // used internally by CellularPotts
  inline double GetNewLengthIfXYWereRemoved(int x, int y) {
    double lengthifxywereremoved = 0.;
    if(!par.periodic_boundaries)
      lengthifxywereremoved = Length(sum_x-x,sum_y-y,sum_xx-x*x,sum_yy-y*y,sum_xy-x*y,area-1);
    return lengthifxywereremoved;
  }


private:
//! Increments the cell's actual area by 1 unit.
  inline int IncrementArea() {
    return ++area;
  }

  //! Decrements the cell's actual area by 1 unit.
  inline int DecrementArea() {
    return --area;
  }

  inline void DecrementAreaBy(int change) {
    area-=change;
  }
  /*! \brief Sets target area to actual area, to remove "pressure".

  This is useful when reading an initial condition from an image.
  */
  inline int SetAreaToTarget(void) {
    return area=target_area;
  }

  //! Called whenever a cell is constructed, from constructor
  void ConstructorBody(int settau=1,int setrecycledsigma=-1);

  // returns the maximum cell type index
  // (depends on Jtable)
  static int MaxTau(void) {
    return maxtau;
  }

protected:
  int colour;
  bool alive;
  int sigma; // cell identity, 0 if medium
  int tau; // Cell type, when dynamicJ's are not used

  double meanx;
  double meany;
  //for Cell migration: target vector
  double tvecx;
  double tvecy;

  //stores old pos of cells for target vector calculations
  double prevx;
  double prevy;

  //store direction of chemokine gradient (int plane)
  double chemvecx;
  double chemvecy;
  double chemmu; //this is the max strength of chemotaxis
  double chemveclength;
  //migration parameters
  int persdur; //how long is this cell's persistent walk?
  int perstime; //counter for how long it has walked persistently
  double mu; //force of migration

  // Two dimensional (square) array of ints, containing the J's.
  double length; // length of the cell;
  double target_length;

  // key and lock pair for j values:
  vector<int> jkey; //= vector<int>(par.key_lock_length, -1);  // notice that part (half) of the key is used also for medium
  vector<int> jlock; //= vector<int>(par.key_lock_length, -1); //c++11 in-class declaration of vectors...
                                                            // so things have to be scaled a bit...

  //array of J values with every other possible sigma,
  // update it everytime a new cell is born
  vector<int> vJ;

  // this is not used anymore - vJ is dynamically increased
  // Dynamically increased when cells are added to the system
  // unless a static Jtable is used (currently this is the default situation)
  static int  **J;

  static int maxtau;

  // Amount: the number of Cell instantations, INCLUDING copies
  // For internal use only.
  // Reading amount is NOT the way to get the number of cells!!
  static int amount;
  static int capacity;
  static int maxsigma; // the last cell identity number given out


  // indices of mother and daughter
  // (Note: no pointers, cells may be relocated)
  int mother;
  int daughter;
  int times_divided;
  int date_of_birth;
  int colour_of_birth;

  int area;
  int target_area;
  int growth_threshold;
  int half_div_area;

  //food intake
  double eatprob;
  int particles;

  //food-conversion-to-growth rate
  double growth;
  //fraction of metabolised particles assigned to maintenance
  // 1 - maintenance_fraction is given to movement
  double maintenance_fraction;

  double k_mf_0; //intercept of maintenance_fraction function
  double k_mf_A; //evolvable factor for area component of maint.fract.
  double k_mf_P; // as above for particles
  double k_mf_C; // as above for contacts

  //fraction of surface proteins (i.e. adhesion) actually used
  double extprotexpress_fraction;
  double k_ext_0;
  double k_ext_A;
  double k_ext_P;
  double k_ext_C;

  double weight_for_chemotaxis;
  double k_chem_0;
  double k_chem_A;
  double k_chem_P;
  double k_chem_C;

  double v[2];
  int n_copies; // number of expansions of this cell
  // gradient of a chemical (to be extended to the total number chemicals)
  double grad[2];

  double *chem;
  // Raw moments of the cells
  // Are used to calculate minor and major axes
  // and center of mass
  // are locally adjusted, so axes are easily
  // and quickly calculated!

  // N.B: N is area!

  long int sum_x;
  long int sum_y;
  long int sum_xx;
  long int sum_yy;
  long int sum_xy;

  const Dish *owner; // pointer to owner of cell

};

#endif
