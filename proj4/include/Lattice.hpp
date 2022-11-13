#ifndef __Lattice__
#define __Lattice__

#include <armadillo>
#include <random>
#include <map>
#include <vector>
#include <cmath>

class Lattice
{
    private:

    public:
    int L, N;
    double T, E, M, E_per, M_per;
    arma::mat spin_matrix;

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

    void fill_lattice(int seed=137);


};


#endif