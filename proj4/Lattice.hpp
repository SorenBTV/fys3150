#include <armadillo>

class Lattice
{
    public:
    int L, N, J, k_b;
    double E, M, Z, Beta, T_in, 
    arma::mat S;


    //Constructor
    Lattice(int N)
    {
        arma::mat S = arma::mat(N,N);
    }

    //Function for energy
    double energy(arma::mat S, int N, int J);

    //Function for magnetization
    double magnetization(arma::mat S, int N);





};