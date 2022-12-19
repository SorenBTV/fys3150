#ifndef __schrodinger__
#define __schrodinger__

#include <armadillo>
#include <iostream>
#include <string>
#include <cstdlib>
#include <complex>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <fstream>

using namespace std;
using namespace std::complex_literals;


class schrodinger
{
    private:

    public:

    int M_;
    double h_;
    double yc_;
    double dt_;
    double T_;
    arma::cx_double r_;
    double x_c;
    double y_c;
    double sigma_x;
    double sigma_y;
    double p_x;
    double p_y;
    double v0_;
    bool save_prob_;


    int nslits_;

    arma::cx_vec a_;
    arma::cx_vec b_;
    arma::cx_vec u_;

    int n;
    arma::cx_cube res;


    arma::mat V_;
    arma::sp_cx_mat A_;
    arma::sp_cx_mat B_;


    void initialize(double h, double dt, double T, double xc, double yc, double sigmax, double sigmay, double px, 
                    double py, double v_0, int n_slits, bool save_prob);
                             
    void initilalize_u();
    int ij_k(int i, int j);
    void A_B();
    void update_u();
    void initialize_V();
    void simulation();

};

#endif