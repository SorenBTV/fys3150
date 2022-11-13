#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;

// output file
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &);
// prints to file the results of the calculations  
void output(int, int, double, vec);