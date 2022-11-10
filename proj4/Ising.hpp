#include <armadillo>

class Ising_model
{
    public:
    int L;
    double T;
    double N;
    arma::mat S;
    double Z;
    double exp_E;
    double exp_e;
    double exp_m;
    double spes_heat;
    double susceptibility;


    //Constructor
    Ising_model(double T, int L);

}