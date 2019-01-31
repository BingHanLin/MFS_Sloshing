#ifndef GMMMFSSLOSHING2D_H
#define GMMMFSSLOSHING2D_H

// include header file
#include<stdio.h>
#include<assert.h>
#include<iostream>
#include<math.h>
#include<vector>
#include<string>
#include <iomanip> // setprecision
#include <time.h>
#include <stdlib.h>
#include <direct.h>
#include <sstream> // stringstream
#include <gmm/gmm.h>


// include header file in project
#include "LRBF/KDTreeTsaiAdaptor.h"
#include "LRBF/Collocation2D.h"

using namespace std;
using namespace gmm;


template <typename MeshType, typename RBFBasisType>
class GMMMFSSLOSHING
{

    public:

        enum class BoundaryType
            {
                Wall,
                FreeSurface
            };

        GMMMFSSLOSHING(MeshType& MMESH, RBFBasisType& RRBFBasis);
        ~GMMMFSSLOSHING() {};

        vector<double> SolveDerPhi(const vector<double> &Phi_temp);
        vector<double> SolveDerEta(const vector<double> &Eta_temp);
        vector<double> UpdateDFSBC(const vector<double> AccTank, const double AccG) const;
        vector<double> UpdateKFSBC() const;


        void SetNodeBCType( int i, BoundaryType bctype);
        void CheckNodeBCType();
        void SetNodeRHSValue( int i, double val);
        void CheckNodeRHSValue();
        void Start(double sloshing_amp, double sloshing_fre, int StepN, double ddt, int SaveStepN,
                   int PrintStepN, double TOLERANCEPHI, double TOLERANCEETA, string SaveFolderName);

    private:

        MeshType MESH;
        RBFBasisType RBFBasis;

        vector<vector<double>> BouNodes, SouNodes, SurfNodes;

        vector<BoundaryType> NodeBCType;
        int Nbou, Nsou, Nsurf;
        
        // why use 2 here?
        KDTreeTsaiAdaptor< vector<vector<double> >, double, 2> kdtree;
        
        vector<double> Eta_n1, Phi_n1, Eta_n, Phi_n;;
        vector<double> d_Eta_n1, d_Phi_n1, d_Eta_n, d_Phi_n;
        double dt;

        gmm::row_matrix<rsvector<double>> SystemPhiMatrix, SystemEtaMatrix; 
        gmm::row_matrix<rsvector<double>> Dx_PhiMatrix, Dy_PhiMatrix, Dx_EtaMatrix;

        void assembly();
        void DeleteOriginData(string);
        void SaveData(string, double);
};


template <typename MeshType, typename RBFBasisType >
GMMMFSSLOSHING<MeshType, RBFBasisType>::
GMMMFSSLOSHING(MeshType& MMESH, RBFBasisType& RRBFBasis)
:MESH(MMESH),RBFBasis(RRBFBasis)
{

    Nbou = MESH.Nbou;
    Nsou = MESH.Nsou;
    Nsurf = MESH.Nsurf;

    BouNodes.resize(Nbou);
    SouNodes.resize(Nsou);
    SurfNodes.resize(Nsurf);

    NodeBCType.resize(MESH.Nbou);
    // NodeRHSValue.resize(MESH.Nbou);

    gmm::resize(SystemEtaMatrix, Nsurf, Nsurf);
    gmm::resize(Dx_EtaMatrix, Nsurf, Nsurf);

    gmm::resize(SystemPhiMatrix, Nbou, Nsou);
    gmm::resize(Dx_PhiMatrix, Nbou, Nsou);
    gmm::resize(Dy_PhiMatrix, Nbou, Nsou);

    gmm::resize(Eta_n1, Nsurf);
    gmm::resize(Phi_n1, Nbou);
    gmm::resize(d_Eta_n1, Nsurf);
    gmm::resize(d_Phi_n1, 2*Nbou);
    gmm::resize(Eta_n, Nsurf);
    gmm::resize(Phi_n, Nbou);
    gmm::resize(d_Eta_n, Nsurf);
    gmm::resize(d_Phi_n, 2*Nbou);

    assembly();
    cout<<"MFS Sloshing Model is created"<<endl;

};


template <typename MeshType, typename RBFBasisType >
void GMMMFSSLOSHING<MeshType, RBFBasisType >::
assembly()
{
    
    // get position of each node
    vector<double> VX(2);

    for(int i = 0; i < Nbou; i++)
    {
        VX=MESH.GetBoundaryNode(i+1);
        BouNodes[i].resize(2);
        BouNodes[i][0]=VX[0];
        BouNodes[i][1]=VX[1];
    }

    for(int i = 0; i < Nsou; i++)
    {
        VX=MESH.GetSourceNode(i+1);
        SouNodes[i].resize(2);
        SouNodes[i][0]=VX[0];
        SouNodes[i][1]=VX[1];
    }

    for(int i = 0; i < Nsurf; i++)
    {
        VX=MESH.GetSurfaceNode(i+1);
        SurfNodes[i].resize(2);
        SurfNodes[i][0]=VX[0];
        SurfNodes[i][1]=VX[1];
    }

    // ==================================    
    // 1D local rbf obtain derivatives of eta
    // ==================================
    int near_num = 5;

    kdtree = KDTreeTsaiAdaptor<vector<vector<double>>, double, 2>(SurfNodes);

    dense_matrix<double> nodes_cloud(near_num,1);
    // using size_t instead of UINT for compatibility reason with nanoflann.hpp
    vector<size_t> neighbours(near_num); 
    vector<double> out_dists_sqr(near_num);
    vector<double> LocalVector;

    // ==================================
    // go through all surface nodes
    // ==================================
    gmm::clear(SystemEtaMatrix);
    gmm::clear(Dx_EtaMatrix);

    for(int i=0; i<Nsurf; i++)
    {
        // use kdtree find indexes of neighbor nodes
        kdtree.query(i, near_num, &neighbours[0],& out_dists_sqr[0]);

        // save nodes cloud in vector
        for(int j=0; j<near_num; j++)
        {
            nodes_cloud(j,0) = SurfNodes[neighbours[j]][0];
        }

        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::IdentityOperation);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for(int j=0; j<near_num; j++)
        {
            SystemEtaMatrix(i,neighbours[j]) = LocalVector[j];
        }

        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::Partial_D1);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for(int j=0; j<near_num; j++)
        {
            Dx_EtaMatrix(i,neighbours[j]) = LocalVector[j];
        }

    };


    // ==================================
    // go through all boundary nodes
    // ==================================
    gmm::clear(SystemPhiMatrix);
    gmm::clear(Dx_PhiMatrix);
    gmm::clear(Dy_PhiMatrix);

    double r1, r2, rs;
    vector<double> VNX(2);

    for(int i=0; i<Nbou; i++)
    {
        for(int j=0; j<Nsou; j++)
        {
        
            r1 = (BouNodes[i][0] - SouNodes[j][0]);
            r2 = (BouNodes[i][1] - SouNodes[j][1]);
            rs = ( r1*r1 + r2*r2 );

            if (NodeBCType[i] == BoundaryType::FreeSurface)
            {
                SystemPhiMatrix(i,j) = log(sqrt(rs));
            }
            else if (NodeBCType[i] == BoundaryType::Wall)
            {
                VNX = MESH.GetOutNormal(i+1);
                SystemPhiMatrix(i,j) = VNX[0]*r1/rs + VNX[1]*r2/rs;     
            }
            else
            {
                cout <<"Wrong NodeBCType Type"<<endl;
                assert(false);
            }

            // derivative(r1) operation for phi
            Dx_PhiMatrix(i,j) = r1/rs;
            // derivative(r2) operation for phi
            Dy_PhiMatrix(i,j) = r2/rs;

        }

    }
   
   cout << "Assemble system matrix successfully"<< endl;
};


template <typename MeshType, typename RBFBasisType >
vector<double> GMMMFSSLOSHING<MeshType, RBFBasisType >::
UpdateDFSBC(const vector<double> AccTank, const double AccG) const
{
    vector <double> Phi_temp(Nsurf);

    for(int i=0; i<Nsurf; i++)
    {

        Phi_temp[i] = ( Phi_n[i] + dt*( 
                        - 0.5*(  pow( (d_Phi_n[i]+d_Phi_n1[i])/2,2) 
                               + pow( (d_Phi_n[i+Nbou]+d_Phi_n1[i+Nbou])/2,2) ) 
                        - AccG* ( ( Eta_n[i] + Eta_n1[i] )/2 )
                        - AccTank[1]*( ( Eta_n[i] + Eta_n1[i] )/2 )
                        - AccTank[0]*SurfNodes[i][0] 
                    ) );
        
    }

    return Phi_temp;
}


template <typename MeshType, typename RBFBasisType >
vector<double> GMMMFSSLOSHING<MeshType, RBFBasisType >::
UpdateKFSBC() const
{
    vector <double> Eta_temp(Nsurf);

    for(int i=0; i<Nsurf; i++)
    {


        Eta_temp[i] = ( Eta_n[i] + dt*( 
                        + ( (d_Phi_n[i+Nbou] + d_Phi_n1[i+Nbou])/2 )
                        - ( (d_Phi_n[i] + d_Phi_n1[i])/2 )*( ( d_Eta_n[i] + d_Eta_n1[i] )/2 )
                       ));
    }

    return Eta_temp;
}


template <typename MeshType, typename RBFBasisType >
vector<double> GMMMFSSLOSHING<MeshType, RBFBasisType >::
SolveDerPhi(const vector<double> &Phi_temp)
{
    vector<double> d_Phi_temp(2*Nbou), Alpah_temp(Nbou), PhiRHS(Nbou);

    for(int i = 0; i < Nbou; i++)
    {
        if (i < Nsurf)
        {
            PhiRHS[i] = Phi_temp[i];
        }
        else if (i >= Nsurf)
        {
            PhiRHS[i] = 0.0;
        }

    }

    gmm::iteration iter(10e-4);// Iteration object with the max residu

    gmm::least_squares_cg(SystemPhiMatrix, Alpah_temp, PhiRHS, iter);

    gmm::mult(Dx_PhiMatrix, Alpah_temp, gmm::sub_vector(d_Phi_temp, gmm::sub_interval(0, Nbou)));
    gmm::mult(Dy_PhiMatrix, Alpah_temp, gmm::sub_vector(d_Phi_temp, gmm::sub_interval(Nbou, Nbou)));

    return d_Phi_temp;
}


template <typename MeshType, typename RBFBasisType >
vector<double> GMMMFSSLOSHING<MeshType, RBFBasisType >::
SolveDerEta(const vector<double> &Eta_temp)
{
    vector<double> d_Eta_temp(Nsurf), Alpah_temp(Nsurf);

    gmm::mult(Dx_EtaMatrix, Eta_temp, d_Eta_temp);
    
    return d_Eta_temp;
}

template <typename MeshType, typename RBFBasisType >
void GMMMFSSLOSHING<MeshType, RBFBasisType >::
Start(double sloshing_amp, double sloshing_fre, int StepN, double ddt, int SaveStepN,
      int PrintStepN, double TOLERANCEPHI, double TOLERANCEETA, string SaveFolderName)
{
    // parameters
    int CNIterMax = 20;
    double AccG = 9.81;
    dt = ddt;

    // initial value
    gmm::clear(Eta_n1);
    gmm::clear(Phi_n1);
    gmm::clear(d_Eta_n1);
    gmm::clear(d_Phi_n1);

    gmm::clear(Eta_n);
    gmm::clear(Phi_n);
    gmm::clear(d_Eta_n);
    gmm::clear(d_Phi_n);


    vector<double> VX(2);
    for(int i=0; i<Nsurf; i++)
    {
       Eta_n1[i] = SurfNodes[i][1];
       Eta_n[i] = SurfNodes[i][1];
    }


    // time loop star
    vector<double> AccTank(2);
    vector<double> Eta_pre(Nsurf), Phi_pre(Nbou), EtaDiff(Nsurf);

    gmm::clear(Eta_pre);
    gmm::clear(Phi_pre);
    gmm::clear(EtaDiff);

    DeleteOriginData(SaveFolderName);
    for(int Istep = 1; Istep<=StepN; Istep++)
    {
        AccTank[0] = -sloshing_amp*sloshing_fre*sloshing_fre*sin(sloshing_fre*dt*Istep);
        AccTank[1] = 0.0;

        for(int CNIter = 1; CNIter<=CNIterMax; CNIter++)
        {

            // Save Phi & Eta in previous iteration
            gmm::copy(Eta_n1 ,Eta_pre);
            gmm::copy(Phi_n1 ,Phi_pre);

            // Obtain derivative of Phi
            d_Phi_n  = SolveDerPhi( Phi_n );
            d_Phi_n1  = SolveDerPhi( Phi_n1 );

            // DFSBC update Phi on free surface
            gmm::copy( UpdateDFSBC(AccTank, AccG), gmm::sub_vector(Phi_n1, gmm::sub_interval(0, Nsurf)) );

            // Obtain derivative of Eta
            d_Eta_n  = SolveDerEta( Eta_n );
            d_Eta_n1  = SolveDerEta( Eta_n1 );

            // KFSBC updating Eta on free surface
            gmm::copy( UpdateKFSBC(), Eta_n1 );
            
            // Calculate error
            if( gmm::vect_dist2(Eta_n1, Eta_pre) <= TOLERANCEETA &&
                gmm::vect_dist2(Phi_n1, Phi_pre) <= TOLERANCEPHI )
            {
                break;
            }

        }
        
        // Save data
        if( Istep % SaveStepN == 0 || Istep == 0)
        {
            SaveData(SaveFolderName, (double)Istep*dt);
        }

        // Print data
        if( Istep % PrintStepN == 0 || Istep == 0)
        {
            cout <<"======================================= Time step  " << Istep <<"  finished."<< endl;
        } 

        // Update system matrix
        gmm::add(Eta_n1, gmm::scaled(Eta_n, -1.0), EtaDiff);
        MESH.UpateBoundaryNode(EtaDiff);
        assembly();

        // Update variables
        gmm::copy(Eta_n1 ,Eta_n);
        gmm::copy(Phi_n1 ,Phi_n);
    }

}



template <typename MeshType, typename RBFBasisType >
void GMMMFSSLOSHING<MeshType, RBFBasisType >::
SetNodeBCType(int i, BoundaryType bctype)
{
    NodeBCType[i] = bctype;
}

template <typename MeshType, typename RBFBasisType >
void GMMMFSSLOSHING<MeshType, RBFBasisType >::
CheckNodeBCType()
{
    for(int i = 0; i < MESH.Nbou; i++)
    {
        cout << static_cast<int>(NodeBCType[i]) << endl;
    }
}

template <typename MeshType, typename RBFBasisType >
void GMMMFSSLOSHING<MeshType, RBFBasisType >::
SaveData(string FolderName, double SaveTime)
{

    string FilePath, FolderPath;
    stringstream  SaveTimeStream;
    SaveTimeStream << fixed << setprecision(5) << SaveTime;
    string SaveTimeStr = SaveTimeStream.str();

    FilePath = ".\\"+ FolderName +"\\"+ SaveTimeStr + ".dat";
    FolderPath =".\\"+ FolderName;

    cout << "Save data at time = " + SaveTimeStr<<endl;

    char charFolderPath[FolderPath.size() + 1];
	strcpy(charFolderPath, FolderPath.c_str());

    mkdir(charFolderPath);

    // Save Boundadry Nodes
    ofstream fout(FilePath);

    for(int i=0; i<Nbou; i++)
    {
        fout<<BouNodes[i][0]<<", "<<BouNodes[i][1]<<", "<<SouNodes[i][0]<<", "<<SouNodes[i][1]<<endl;
    }

    fout.close();

}


template <typename MeshType, typename RBFBasisType >
void GMMMFSSLOSHING<MeshType, RBFBasisType >::
DeleteOriginData(string FolderName)
{

    string FolderPath;
    string DeleteString;

    FolderPath =".\\"+ FolderName;
    DeleteString = "del " + FolderPath + " /f /s /q";

    cout << "Delete Original Data in "<< FolderName <<endl;

    char charDeleteString[DeleteString.size() + 1];
	strcpy(charDeleteString, DeleteString.c_str());

    system(charDeleteString);

}

#endif