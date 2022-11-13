#ifndef __Lattice__
#define __Lattice__

#include <armadillo>


class Lattice
{
    private:

    public:
    int L;
    int N;
    double T;
    arma::Mat<int> spin_matrix;

    // A random number generator and the two basic 
    // probability distributions we need.
    /*
    std::mt19937 generator;
    std::uniform_real_distribution<double> draw_uniform_r;
    std::uniform_int_distribution<int> draw_uniform_index;
    */
    //Constructor
    Lattice(int L, double T);

    //Functions

    void fill_lattice(int seed);


};


#endif