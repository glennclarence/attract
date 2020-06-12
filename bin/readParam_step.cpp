#include <string>
#include "state.h"
#include <cstdio>
#include <fstream>
#include <iostream>

using namespace std;
void readParam_step(string fileName, CartState & cs){
    
    
    ifstream in(fileName);
    
    if (!in) {
        cerr << "Cannot open file.\n";
        return;
    }
    //read step-energies
    for (int z=0; z < cs.step_nrDists; z++){
        for (int y = 0; y < cs.step_nrTypes; y++) {
            for (int x = 0; x < cs.step_nrTypes; x++) {
                
                in >> cs.stepPot_energy[z][x][y];
            }
        }
    }
    
    
    //cs.stepPot_dist.push_back(2.8);
    cs.stepPot_dist.push_back(1.5); //repulsion-step
    cs.stepPot_dist.push_back(4);
    cs.stepPot_dist.push_back(6);
    
    // //read step-len
    //   for (int z=0; z < cs.step_nrDists; z++){
    //   for (int y = 0; y < cs.step_nrTypes; y++) {
    //     for (int x = 0; x < cs.step_nrTypes; x++) {
    
    //       in >> cs.stepPot_dist[z][x][y];
    //     }
    //   }
    //   }
    //   in.close();
    
    // //scale first step-len
    // for (int y = 0; y < cs.step_nrTypes; y++) {
    //   for (int x = 0; x < cs.step_nrTypes; x++) {
    
    //      cs.stepPot_dist[0][x][y]= cs.stepPot_dist[0][x][y]*0.85 ;
    //   }
    // }
    
}
