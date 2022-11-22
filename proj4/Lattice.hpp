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
    int64_t L_;
    int64_t N_;
    double T_;
    int64_t MC_cycles_;
    double E_;
    double M_;
    bool ordered_;
    arma::vec beta_;
    arma::mat spin_matrix_;
    arma::vec eps;
    arma::vec mag;

    double eps_;
    double eps2_;
    double m_abs;
    double m2;
    double c_v;
    double X;
    arma::vec boltzmann_;
    

    mt19937 generator;
    uniform_int_distribution<int> int_01;
    uniform_int_distribution<int> dist;
    uniform_real_distribution<double> index;

    arma::vec E_dist_;

    double expectedE_;
    double expectedE2_;
    double expectedM_;
    double expectedM2_;
    double epsilon_;
    double cv_;
    double chi_;
    arma::vec energy;
    //arma::vec boltzmann;

    // A random number generator and the two basic 
    // probability distributions we need.
    /*
    std::mt19937 generator;
    std::uniform_real_distribution<double> draw_uniform_r;
    std::uniform_int_distribution<int> draw_uniform_index;
    */

    

    //Functions
    void Initializer(int64_t L, double T, int64_t MC_cycles, bool ordered);
    int pbc(int64_t i, int64_t add);
    void Metropolis();
    void MCMC();
    void MCMC_burn_in_time_study(string filename);
    void write_file_problem4();
    void write_file();
};


#endif