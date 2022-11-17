#include "Lattice.hpp"
#include <iostream>
#include <fstream>
#include <armadillo>
#include <random>
#include <cmath>
#include <string>

Lattice::Lattice(int L, double T)
{
  this->L = L;
  this->T = T;
  this->J = 1;
  this->N = L*L;
  this->k_b=1;
  this->beta=1/T;
  arma::Mat<int> spin_matrix = arma::Mat<int>(L, L).fill(0);
  this->spin_matrix = spin_matrix;
  this->E = Lattice::tot_energy();

  this->eps=0;
  this->eps2=0;
  this->m_abs=0;
  this->m2=0;
  this->c_v=0;
  this->X=0;

}

void Lattice::fill_lattice(int seed)
{
  std::mt19937 generator (seed);
  std::uniform_int_distribution<int> int_01(0,1);
  
  int spin;
  for (int i = 0; i < L; i ++)
  {
    for (int j = 0; j < L; j ++)
    {
      spin = int_01(generator);
      spin = spin*2 - 1;
      this->spin_matrix(i, j) = spin;
    }
  }

}

// inline function for periodic boundary conditions
inline int pbc(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}


double Lattice::tot_energy()
{
  double E = 0;
  for(int i =0; i < L; i++)
  {
    for (int j= 0; j < L; j++)
    {
      E -= this->spin_matrix(i,j) * (this->spin_matrix(pbc(i,L,-1),j) + this->spin_matrix(i,pbc(j,L,-1)));
    }
  }
  return E;
}


double Lattice::tot_magnetization()
{
  double M = 0;
  for(int i =0; i < L; i++)
  {
    for (int j= 0; j < L; j++)
    {
      M += this->spin_matrix(i,j);
    }
  }
  return M;
}

void Lattice::flip_spin(int x, int y)
{
  this->spin_matrix(x,y) *= -1.0;
}



void Lattice::MCMC(int seed, int MC_cycles, arma::vec Expectation_values, double T, int L, int N, double E)
{
  std::mt19937 generator (seed);
  std::uniform_real_distribution<double> index(0,L);
  std::uniform_real_distribution<double> dist(0,1);

  arma::vec boltzmann = arma::vec(17, arma::fill::zeros);
  for (int i=-8; i<=8; i+=4)
  {
    boltzmann(i+8) = exp(-1 * beta * i);
  }

  //arma::vec dE = {-8, -4, 0, 4, 8};

  //std::cout << energy_change << std::endl;
  for(int cycles=1; cycles<=MC_cycles; cycles++)
  {
    // The sweep over the lattice, looping over all spin sites

        int x = (index(generator));
        int y = (index(generator));
        //std::cout << ix << iy << std::endl;
        int delta_E = 2 * this->spin_matrix(x,y)*
        (this->spin_matrix(x,pbc(y,L,-1))+
        this->spin_matrix(pbc(x,L,-1),y) +
        this->spin_matrix(x,pbc(y,L,1)) +
        this->spin_matrix(pbc(x,L,1),y));

      
        double w = boltzmann(delta_E/4+2);
        double r = dist(generator);
        if(w <= 0)
        {
          flip_spin(x, y);
          E += delta_E;
          M += 2*this->spin_matrix(x, y);
        }
        
        else {if (r <= w)
        {
          flip_spin(x, y);
          E += delta_E;
          M += 2*this->spin_matrix(x, y);
        }}
        
  
    Expectation_values(0) += E;
    Expectation_values(1) += E*E;
    Expectation_values(2) += abs(M);    
    Expectation_values(3) += M*M;
    //std::cout << Expectation_values(0) << std::endl;
  }
  
  //std::cout << Expectation_values(0) << std::endl;
  double c_v = (1/N)*(1/(k_b*T*T))*(Expectation_values(1) - Expectation_values(0)*Expectation_values(0));
  double X = (1/N)*(1/(k_b*T*T))*(Expectation_values(3) - Expectation_values(2)*Expectation_values(2));

  double eps = Expectation_values(0)/N;
  double eps2 = Expectation_values(1)/(N*N);
  double m_abs = Expectation_values(2)/N;
  double m2 = Expectation_values(3)/(N*N);
  std::cout << eps << std::endl;
}

void Lattice::analytical(int L, double T)
{
  double Z = 2*exp(beta*8) + 2*exp(-beta*8) + 12;
  double epsilon = (4*J/Z)*(exp(-8/T) - exp(8/T));
  double epsilon2 = (8*J/Z)*(exp(-8/T) + exp(8/T));
  double magnetization_abs= abs(2*J/Z * (exp(8/T) + 2));
  double magnetization2 = 2*J/Z * (exp(8/T) + 1);
  double cv = L*L/(T*T) * (epsilon2 - epsilon*epsilon);
  double chi = L*L/T * (magnetization2 - magnetization_abs*magnetization_abs);
  std::cout << epsilon << std::endl;
}