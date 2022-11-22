#include "Lattice.hpp"


void Lattice::Initializer(int64_t L, double T, int64_t MC_cycles, bool ordered)
{
  L_ = L;
  T_ = T;
  N_ = L_*L_;
  MC_cycles_ = MC_cycles;
  ordered_ = ordered;
  energy = arma::vec (MC_cycles_).fill(0.);



  // Initial matrix
  arma::mat spin_matrix = arma::mat(L_, L_, arma::fill::ones);
  spin_matrix_ = spin_matrix;

  mt19937 generator (123);
  uniform_int_distribution<int64_t> int_01(0.0,1.0);

  if(ordered_ == false)
  {
    int spin;
    for (int64_t i = 0; i < L_; i ++)
    {
      for (int64_t j = 0; j < L_; j ++)
      {
        spin = int_01(generator);
        spin = spin*2. - 1.;
        spin_matrix_(i, j) = spin;
      }
    }
  }


  // Initial energy
  E_ = 0.0;

  for(int64_t i =0; i < L_; i++)
  {
    for (int64_t j= 0; j < L_; j++)
    {
      E_ -= (double)spin_matrix_(i,j) * (spin_matrix_(pbc(i,-1),j) + spin_matrix_(i,pbc(j,-1)));
    }
  }


  // Initial Magnetization
  double M_ = 0.0;
  for(int64_t i =0; i < L_; i++)
  {
    for (int64_t j= 0; j < L_; j++)
    {
      M_ += (double)spin_matrix_(i,j);
    }
  }

  // Boltzmann dist
  boltzmann_ = arma::vec(17).fill(0.);
  boltzmann_(0) = exp(8./T_);
  boltzmann_(4) = exp(4./T_);
  boltzmann_(8) = exp(0.);
  boltzmann_(12) = exp(-4./T_);
  boltzmann_(16) = exp(-8./T_);

}


// inline function for periodic boundary conditions
int Lattice::pbc(int64_t i, int64_t add) { 
  return (i+L_+add) % L_;
}



void Lattice::Metropolis()
{
  uniform_real_distribution<double> index(0, 1);
  uniform_int_distribution<int64_t> dist(0, L_-1);
  
  for (int64_t i=0; i<N_; i++)
  {

    int64_t x = dist(generator);
    int64_t y = dist(generator);
    

    int64_t dE =  2 * spin_matrix_(x,y)*
        (spin_matrix_(x, pbc(y, -1))+
        spin_matrix_(x, pbc(y, 1))+
        spin_matrix_(pbc(x, -1), y) +
        spin_matrix_(pbc(x, 1), y));

      
    
    double r = index(generator);
    if (dE <= 0)
    {
      spin_matrix_(x,y) *= -1;
      E_ += (double) dE;
      M_ += (double) 2.0*spin_matrix_(x,y);
    }

    else if(r<boltzmann_(dE+8))
    {
      spin_matrix_(x,y) *= -1;
      E_ += (double) dE;
      M_ += (double) 2.0*spin_matrix_(x,y);
      
    }
  }
}




void Lattice::MCMC()
{
  int64_t burnin = MC_cycles_/10;
  int64_t cyc = MC_cycles_ - burnin;
  double norm = 1./(cyc*N_);
  expectedE_ = 0.0;
  expectedE2_ = 0.0;
  expectedM_ = 0.0;
  expectedM2_ = 0.0;
  epsilon_ = 0.0;

  for (int64_t i=1; i<=burnin; i++)
  {
    Metropolis();
  }

  for (int64_t cycles=1; cycles<=cyc; cycles++)
  {
    Metropolis();
    expectedE_ += E_;
    expectedE2_ += E_*E_;
    expectedM_ += fabs(M_);
    expectedM2_ += M_*M_;
    energy(cycles-1) = E_;
  }
  expectedE_ = expectedE_ * norm;
  expectedE2_ = expectedE2_ * norm;
  expectedM_ = expectedM_ * norm;
  expectedM2_ = expectedM2_ * norm;
  cv_ = (expectedE2_*N_ - expectedE_ * expectedE_ *(N_*N_)) /(T_*T_*N_);
  chi_ = (expectedM2_*N_ - expectedM_ * expectedM_ *(N_*N_)) /(T_*N_);
}



void Lattice::MCMC_burn_in_time_study(string filename)
{
  eps = arma::vec(MC_cycles_).fill(0.);
  mag = arma::vec(MC_cycles_).fill(0.);
  expectedE_ = 0.0;
  expectedM_ = 0.0;

  for (int cycles=1; cycles<=MC_cycles_; cycles++)
  {
    Metropolis();
    expectedE_ += E_;
    expectedM_ += fabs(M_);
    eps(cycles-1) = expectedE_ /(N_*cycles);
    mag(cycles-1) = expectedM_ /(N_*cycles);
    energy(cycles-1) = E_;
  }
   
  
  string fname = filename + ".csv";
  ofstream ofile;
  ofile.open(fname);
  int width = 12;
  int prec = 4;

  ofile << setw(width) << setprecision(prec) << scientific << "Cycles"
  << setw(width) << setprecision(prec) << scientific << "Epsilon/N"
  << setw(width) << setprecision(prec) << scientific << "Magnetization/N"
  << setw(width) << setprecision(prec) << scientific << "Epsilon" << endl;
  for (int i=1; i<=MC_cycles_; i++)
  {
    ofile << setw(width) << setprecision(prec) << scientific << i
    << setw(width) << setprecision(prec) << scientific << eps(i-1)
    << setw(width) << setprecision(prec) << scientific << mag(i-1)
    << setw(width) << setprecision(prec) << scientific << energy(i-1) << endl;
  }

}




void Lattice::write_file_problem4()
{
  
  string filename = "Problem4.txt";

  //Analytical results
  double Z = 2*exp(8/T_) + 2*exp(-8/T_) + 12;
  double epsilon = (4/Z)*(exp(-8/T_) - exp(8/T_));
  double epsilon2 = (8/Z)*(exp(-8/T_) + exp(8/T_));
  double magnetization_abs= abs(2/Z * (exp(8/T_) + 2));
  double magnetization2 = 2/Z * (exp(8/T_) + 1);
  double c_v = L_*L_/(T_*T_) * (epsilon2 - epsilon*epsilon);
  double _chi = L_*L_/T_ * (magnetization2 - magnetization_abs*magnetization_abs);


  ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec = 4;


  ofile << "Number of cycles = " << MC_cycles_ << endl;
  ofile << "Numerical results" << endl;
  ofile << "{epsilon}=" << setw(width) << setprecision(prec) << scientific << expectedE_ << "  "
  << "{|m|}=" << setw(width) << setprecision(prec) << scientific << expectedM_ << "  "
  << "{C_v}=" << setw(width) << setprecision(prec) << scientific << cv_ << "  "
  << "{chi}=" << setw(width) << setprecision(prec) << scientific << chi_ << endl;

  ofile << "Analytical results" << endl;
  ofile << "{epsilon}=" << setw(width) << setprecision(prec) << scientific << epsilon << "  "
  << "{|m|}=" << setw(width) << setprecision(prec) << scientific << magnetization_abs << "  "
  << "{C_v}=" << setw(width) << setprecision(prec) << scientific << c_v << "  "
  << "{chi}=" << setw(width) << setprecision(prec) << scientific << _chi << endl; 
  ofile.close();
}


void Lattice::write_file()
{
  
  string filename = "Problem6.txt";

  ofstream ofile;
  ofile.open(filename);
  int width = 15;
  int prec = 8;

  ofile << "Number of cycles =" << MC_cycles_ << endl;
  ofile << "Numerical results" << endl;
  ofile << "{epsilon}=" << setw(width) << setprecision(prec) << scientific << expectedE_ << " ";
  ofile << "{|m|}=" << setw(width) << setprecision(prec) << scientific << expectedM_ << " ";
  ofile << "{C_v}=" << setw(width) << setprecision(prec) << scientific << cv_ << " ";
  ofile << "{chi}=" << setw(width) << setprecision(prec) << scientific << chi_ << endl;
 
  ofile.close();
}


