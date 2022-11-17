#ifndef __Lattice__
#define __Lattice__

#include <iostream>
#include <fstream>
#include <armadillo>
#include <random>
#include <cmath>
#include <string>
#include <iomanip>
#include <chrono>

using namespace std;

class Lattice
{
    private:

    public:
    int L_;
    int N_;
    double T_;
    int MC_cycles_;
    double E_;
    double M_;
    arma::vec beta_;
    arma::Mat<int> spin_matrix_;

    double eps_;
    double eps2_;
    double m_abs;
    double m2;
    double c_v;
    double X;
    double *boltzmann_;
    arma::mat results_;

    mt19937 generator;
    uniform_int_distribution<int> int_01;
    uniform_int_distribution<int> dist;
    uniform_real_distribution<double> index;

    double expectedE_;
    double expectedE2_;
    double expectedM_;
    double expectedM2_;
    double epsilon_;
    double cv_;
    double chi_;
    //arma::vec boltzmann;

    // A random number generator and the two basic 
    // probability distributions we need.
    /*
    std::mt19937 generator;
    std::uniform_real_distribution<double> draw_uniform_r;
    std::uniform_int_distribution<int> draw_uniform_index;
    */

    

    //Functions
    void Initializer(int L, double T, int MC_cycles, bool ordered);
    int pbc(int i, int limit, int add);
    void Metropolis();
    void MCMC();
    void MCMC_burn_in_time(string filename_);
    void write_file_problem4();
    void write_file();
};


#endif