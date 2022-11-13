#include "Lattice.hpp"
#include <iostream>
#include <fstream>
#include <armadillo>
#include <random>
#include <cmath>
#include <string>

Lattice::Lattice(int L_in, double T_in)
{
  this->L = L;
  this->T = T;
  arma::Mat<int> spin_matrix = arma::Mat<int>(L, L, arma::fill::zeros);
  this->spin_matrix = spin_matrix;

}

void Lattice::fill_lattice(int seed)
{
  std::mt19937 generator (seed);
  //std::uniform_real_distribution<double> draw_uniform_r;
  std::uniform_int_distribution<int> dis(0,1);
  
  int num;
  for (int i = 0; i < L; i ++)
  {
    for (int j = 0; j < L; j ++)
    {
      num = dis(generator);
      num = num*2 - 1;
      this->spin_matrix(i, j) = num;
    }
  }

}
