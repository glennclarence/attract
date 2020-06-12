#include <cmath>
#include "state.h"
#include <vector>
#include <iostream>

using namespace std;

extern CartState &cartstate_get(int cartstatehandle);


inline double get_stepVal(CartState & cs, int recType,int ligType,double dist){
    // double (&dists)[cs.step_nrDists][cs.step_nrTypes][cs.step_nrTypes]=cs.stepPot_dist;
    // double (&stepEnergies)[cs.step_nrDists][cs.step_nrTypes][cs.step_nrTypes]=cs.stepPot_energy;
    for (int i=0;i<cs.step_nrDists;i++){
        if (dist<cs.stepPot_dist[i]){     //[recType][ligType]) {
            //if (i==0) cerr<<"first step collision w/ "<<dist<<", atomstypes:" <<recType<<" "<<ligType<<endl;         
            return cs.stepPot_energy[i][recType][ligType]; 
            break;
        }
    } 
    return 0;
    //this is >10% faster, but fixed regarding stepping... 
    /*                if (dist>dists[1]) continue;
     *        else if (dist>dists[0]) enon+=stepEnergies[1][recType][ligType];
     *        else enon+=stepEnergies[0][recType][ligType];*/
}

extern "C" void nonbon_step_( double xl[],double xr[],  
                              int iacir[], int iacil[],int nonr[], int nonl[],  int& nonp,
                              double& enon, double& epote, int& cartstateHandle) {
    
    enon=0;
    epote=0;
    CartState & cs=cartstate_get(cartstateHandle);
    
    for (int pair=0;pair<nonp;pair++){
        int recAtom=nonr[pair];  //not with cpp-routine... -1; //offset from fortran...
        int ligAtom=nonl[pair];                          //-1;
        
        int recType=iacir[recAtom]-1; // atomtypes beginn at 1, index of params at 0...
        int ligType=iacil[ligAtom]-1;
        
        int recCoord=3*(recAtom); //test too ->passt
        int ligCoord=3*(ligAtom); 
        double dist=0;
        for (int i=0;i<3;i++){
            dist+=pow(xr[recCoord+i]-xl[ligCoord+i],2); //debug for lig-xori: cs.xori[cs.n3atom[0]+ligCoord]
        }
        dist=sqrt(dist);
        enon +=get_stepVal(cs, recType,ligType,dist);
    }
                              }
                              
                              
                              
                              double nonbon_step_cartstate(const int cartstateHandle, const int recType, const int ligType, const double dist){
                                  CartState & cs=cartstate_get(cartstateHandle);
                                  return get_stepVal(cs, recType,ligType,dist);
                              }
                              
                              
                              extern "C" bool   use_steppot_(int & cartstatehandle){ //fortran will not find upper-case function names...
                                  CartState & cs =cartstate_get(cartstatehandle);
                                  bool test=cs.useStepPot;
                                  return test;
                              }
