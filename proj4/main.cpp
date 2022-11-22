// including the header and function files containing our classes
#include "omp.h"
#include "Lattice.hpp"


/*
Make:
g++ main.cpp src/Lattice.cpp -I include/ -o main.exe -fopenmp -larmadillo -Wall && ./main.exe
*/

void parallellization_temp_loop(int64_t L, int64_t cycles_per_thread)
{

    auto t1 = chrono::steady_clock::now();

    double T_init = 2.1;
    double T_end = 2.4;
    int64_t n = 6;
    arma::mat results(n, 5, arma::fill::zeros);
    arma::vec temp = arma::linspace(T_init, T_end, n);

    #pragma omp parallel
    {
        #pragma omp for
        for (int64_t i=0; i<n; i++)
        {
            double T_now = temp(i);
            Lattice mysystem;
            mysystem.Initializer(L, T_now, cycles_per_thread, false);
            mysystem.MCMC();

            results(i, 0) = temp(i);
            results(i, 1) = mysystem.expectedE_;
            results(i, 2) = mysystem.expectedM_;
            results(i, 3) = mysystem.cv_;
            results(i, 4) = mysystem.chi_;
        }
    }
    auto t2 = chrono::steady_clock::now();
    cout<< "Time used to run paralellization for L=" << L << " was "<< chrono::duration_cast<chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<endl;

    string output_file = "Paralellization_L" + to_string(L) + ".csv";
    ofstream ofile;
    ofile.open(output_file);
    int64_t width = 12;
    int64_t prec = 4;

    for (int64_t i=0; i<n; i++)
    {

        ofile << setw(width) << setprecision(prec) << scientific << results(i, 0)
        << setw(width) << setprecision(prec) << scientific << results(i, 1)
        << setw(width) << setprecision(prec) << scientific << results(i, 2)
        << setw(width) << setprecision(prec) << scientific << results(i, 3)
        << setw(width) << setprecision(prec) << scientific << results(i, 4) << endl;
    }
    ofile.close();
}


void write_analytical_vs_numerical_L2_T1(int64_t L, double T, int64_t MC_cycles, bool ordered)
{

  //Analytical results
  double Z = 2*exp(8/T) + 2*exp(-8/T) + 12;
  double epsilon = (4/Z)*(exp(-8/T) - exp(8/T));
  double epsilon2 = (8/Z)*(exp(-8/T) + exp(8/T));
  double magnetization_abs= 2/Z * (exp(8/T) + 2);
  double magnetization2 = 2/Z * (exp(8/T) + 1);
  double c_v = L*L/(1.0) * (epsilon2 - epsilon*epsilon);
  double _chi = L*L/(1.0) * (magnetization2 - magnetization_abs*magnetization_abs);

  Lattice mysystem;
  mysystem.Initializer(L, T, MC_cycles, ordered);
  mysystem.MCMC();


  string filename = "analytical_vs_numerical_L2_T1" + to_string(MC_cycles)+".csv";
  ofstream ofile;
  ofile.open(filename);
  int64_t width = 12;
  int64_t prec = 4;


  ofile << "Number of cycles = " << MC_cycles << endl;
  ofile << "Numerical results" << endl;
  ofile << "{epsilon}=" << setw(width) << setprecision(prec) << scientific << mysystem.expectedE_ << "  "
  << "{|m|}=" << setw(width) << setprecision(prec) << scientific << mysystem.expectedM_ << "  "
  << "{C_v}=" << setw(width) << setprecision(prec) << scientific << mysystem.cv_ << "  "
  << "{chi}=" << setw(width) << setprecision(prec) << scientific << mysystem.chi_ << endl;

  ofile << "Analytical results" << endl;
  ofile << "{epsilon}=" << setw(width) << setprecision(prec) << scientific << epsilon << "  "
  << "{|m|}=" << setw(width) << setprecision(prec) << scientific << magnetization_abs << "  "
  << "{C_v}=" << setw(width) << setprecision(prec) << scientific << c_v << "  "
  << "{chi}=" << setw(width) << setprecision(prec) << scientific << _chi << endl; 
  ofile.close();
}


void equilibration_time(int64_t L, double T, int64_t MC_cycles, bool ordered, string filename)
{

    Lattice mysystem;
    mysystem.Initializer(L, T, MC_cycles, ordered);
    mysystem.MCMC_burn_in_time_study(filename);
    
}

void probability_density(int64_t L, double T, int64_t MC_cycles, bool ordered, string filename)
{
    Lattice mysystem;
    mysystem.Initializer(L, T, MC_cycles, ordered);
    mysystem.MCMC();

      
    string fname = filename + ".csv";
    ofstream ofile;
    ofile.open(fname);
    int64_t width = 12;
    int64_t prec = 4;

    for (int64_t i=1; i<=MC_cycles; i++)
    {
        ofile << setw(width) << setprecision(prec) << i
        << setw(width) << setprecision(prec) << mysystem.energy(i-1) << endl;
    }
    ofile.close();
}

int main()
{
    cout << "Select each of these for the data needed:"<< endl;
    cout << "   [1] Write analytical vs numerical data for L=2" << endl;
    cout << "   [2] Equilibration time for L=20, T=1.0 and T=2.4" << endl;
    cout << "   [3] Probability density data for L=20, T=1.0 and T=2.4" << endl;
    cout << "   [4] Parallelliation for L=40" << endl;
    cout << "   [5] Parallelliation for L=60" << endl;
    cout << "   [6] Parallelliation for L=80" << endl;
    cout << "   [7] Parallelliation for L=100" << endl;
    
    cout << "Please select a number" << endl;

    int x;
    std::cin >> x;
    switch(x)
    {


        case 1:
        {
            write_analytical_vs_numerical_L2_T1(2, 1.0, 10000, false);
            write_analytical_vs_numerical_L2_T1(2, 1.0, 100000, false);
            write_analytical_vs_numerical_L2_T1(2, 2.4, 1000000, false);
            break;
        }

        case 2:
        {
            equilibration_time(20, 1.0, 100000, false, "Equilibrium_time_study_L20_T1.0_not_ordered");
            equilibration_time(20, 1.0, 100000, true, "Equilibrium_time_study_L20_T1.0_ordered");
            equilibration_time(20, 2.4, 100000, false, "Equilibrium_time_study_L20_T2.4_not_ordered");
            equilibration_time(20, 2.4, 100000, true, "Equilibrium_time_study_L20_T2.4_ordered");
            break;
        }

        case 3:
        {
            probability_density(20, 1.0, 30000, false, "probability_density_T1_L20");
            probability_density(20, 2.4, 30000, false, "probability_density_T2.4_L20");   
            break;  
        }


        case 4:
        {
            parallellization_temp_loop(40, 500000);
            break;
        }

         case 5:
        {
            parallellization_temp_loop(60, 500000);
            break;
        }

         case 6:
        {
            parallellization_temp_loop(80, 500000);
            break;
        }

         case 7:
        {
            parallellization_temp_loop(100, 500000);
            break;
        }
    }



    return 0;
}