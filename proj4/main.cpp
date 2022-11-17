// including the header and function files containing our classes
#include "omp.h"
#include "Lattice.hpp"


/*
Make:
g++ main.cpp src/Lattice.cpp -I include/ -o main.exe -fopenmp -larmadillo -Wall && ./main.exe
*/

int main()
{
    
    int L = 2.0;
    double T = 1.0;
    int MC_cycles = 10000;
    bool ordered = true;
    int threads = 4;

    int seed = 123;
    int N = L * L;

    //double T_final = 1;
    //double dt = (T_final - T)/N;
    
    /*
    // Problem 4
    Lattice mysystem;
    mysystem.Initializer(L, 2.4, MC_cycles, false);
    mysystem.MCMC();
    //cout << mysystem.spin_matrix_<< endl;
    //mysystem.write_file_problem4();
    


*/
    //Problem 5
    Lattice mysystem2;
    mysystem2.Initializer(20, 1.0, 1000, false);
    mysystem2.MCMC_burn_in_time("problem_5_T1.0_cycles10000_not_ordered.txt");
/*
    Lattice mysystem3;
    mysystem2.Initializer(20, 1.0, 100000, true);
    mysystem2.MCMC_burn_in_time("problem_5_T1.0_cycles10000_ordered.txt");

    Lattice mysystem4;
    mysystem2.Initializer(20, 2.4, 100000, false);
    mysystem2.MCMC_burn_in_time("problem_5_T2.4_cycles10000_not_ordered.txt");

    Lattice mysystem5;
    mysystem2.Initializer(20, 2.4, 100000, true);
    mysystem2.MCMC_burn_in_time("problem_5_T2.4_cycles10000_ordered.txt");
    */


    //Problem 6




    //Problem 7,8


    //#pragma omp parallel{omp_set_num_threads(threads);}
    





    return 0;
}