using namespace std;
using namespace gmm;


// collocation for LRBF
template<typename T, typename RbfBasisType>
vector<T> Collocation2D(const dense_matrix<T> & nodes_cloud, RbfBasisType& RbfBasis)
{
    int DIMENSION = mat_ncols(nodes_cloud);
    
    int Nearnum = mat_nrows(nodes_cloud);

    dense_matrix<T> Phi(Nearnum,Nearnum);
    vector<T> temp(Nearnum), LPhi(Nearnum);
    vector<T> VX1(DIMENSION), VX2(DIMENSION);

    for(int i=0; i<Nearnum; i++)
    {
        for(int j=0; j<Nearnum; j++)
        {
            for (int d = 0; d<DIMENSION; d++)
            {
                VX1[d] = nodes_cloud(i,d);
                VX2[d] = nodes_cloud(j,d);
            }

            Phi(j,i) /* transposed */ = RbfBasis.GetBasisValue(VX1, VX2, 1);

            if (i == 0)
                LPhi[j] = RbfBasis.GetBasisValue(VX1, VX2);
        }
    }

    // ! 考慮之後加上判斷 ill condition 
    lu_solve(Phi, temp, LPhi);
    
    return temp;
};