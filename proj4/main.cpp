#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "omp.h"
#include <random>
// including the header and function files containing our classes
#include "Lattice.hpp"

using namespace std;

/*
Make:
g++ main.cpp src/Lattice.cpp -I include/ -o main.exe -fopenmp -larmadillo -Wall && ./main.exe
*/

int main()
{
    int L = 5;
    double T = 1;
    int seed = 137;
    Lattice mysystem = Lattice(L, T);
    mysystem.fill_lattice(seed);
    cout << mysystem.spin_matrix << endl;
    //cout << Lattice::spin_matrix << endl;


    return 0;
}