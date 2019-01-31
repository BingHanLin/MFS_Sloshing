// main_MFS.cpp -- reading more than one string

// include header file
#include <iostream>
#include<iostream>
#include <gmm/gmm.h>
#include<vector>

// include header file in project
#include "GMMRectangular.h"
#include "GMMMFSSloshing2D.h"
#include "LRBF/GMMMQBasis2D.h"


using namespace std;
using namespace gmm;

int main()
{
    int m = 30;     // number of nodes in x direction
    int n = 20;     // number of nodes in z direction
    double Lx = 0.9;  // domain size in x direction
    double Lz = 0.6;  // domain size in z direction

    int StepN = 1000;
    double dt = 0.01;
    int SaveStepN = 5;
    int PrintStepN = 5;
    double TOLERANCEPHI = 10e-6;
    double TOLERANCEETA = 10e-6;
    
    double sloshing_amp = 0.002;
    double sloshing_fre = 5.5;

    string SaveFolderName = "OutputData";

    // Create Mesh
    GMMRECTANGLE MESH(m, n, Lx, Lz, 0, 0, 1.25);

    // Chose RBF Basis (for frees urface derivative) 
    GMMMQBasis2D RBFBasis(0.2);

    // Create Model 
    GMMMFSSLOSHING<GMMRECTANGLE, GMMMQBasis2D> MODEL(MESH, RBFBasis);

    // Setting Boundary Condition Type
    for(int i = 0; i < MESH.Nbou; i++)
    {
        if (i < MESH.Nsurf)
        {
            MODEL.SetNodeBCType(i, GMMMFSSLOSHING<GMMRECTANGLE, GMMMQBasis2D>::BoundaryType::FreeSurface);
        }
        
        else if (i >= MESH.Nsurf)
        {
            MODEL.SetNodeBCType(i, GMMMFSSLOSHING<GMMRECTANGLE, GMMMQBasis2D>::BoundaryType::Wall);
        }
    }


    // Start simulation
    MODEL.Start(sloshing_amp, sloshing_fre, StepN, dt, SaveStepN, PrintStepN, TOLERANCEPHI, TOLERANCEETA, SaveFolderName);


    cout<<"Simulation Completed"<<endl;   

    return 0;   

}




