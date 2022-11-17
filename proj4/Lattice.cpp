#include "Lattice.hpp"


void Lattice::Initializer(int L, double T, int MC_cycles, bool ordered)
{
  L_ = L;
  T_ = T;
  N_ = L_*L_;
  MC_cycles_ = MC_cycles;
  E_dist_ = new double[N_*MC_cycles_];


  // Initial matrix
  arma::Mat<int> spin_matrix = arma::Mat<int>(L_, L_).fill(1);
  spin_matrix_ = spin_matrix;

  mt19937 generator (123);
  uniform_int_distribution<int> int_01(0.0,1.0);

  if(ordered == false)
  {
    int spin;
    for (int i = 0; i < L_; i ++)
    {
      for (int j = 0; j < L_; j ++)
      {
        spin = int_01(generator);
        spin = spin*2. - 1.;
        this->spin_matrix_(i, j) = spin;
      }
    }
  }


  // Initial energy
  double E_ = 0;
  for(int i =0; i < L_; i++)
  {
    for (int j= 0; j < L_; j++)
    {
      E_ -= spin_matrix_(i,j) * (spin_matrix_(pbc(i,L_,-1),j) + spin_matrix_(i,pbc(j,L_,-1)));
    }
  }



  // Initial Magnetization
  double M_ = 0.;
  for(int i =0; i < L_; i++)
  {
    for (int j= 0; j < L_; j++)
    {
      M_ += spin_matrix_(i,j);
    }
  }

  // Boltzmann dist
  boltzmann_ = new double[17];
  boltzmann_[0] = exp(8/T_);
  boltzmann_[4] = exp(4/T_);
  boltzmann_[8] = exp(0);
  boltzmann_[12] = exp(-4/T_);
  boltzmann_[16] = exp(-8/T_);
 
}


// inline function for periodic boundary conditions
int Lattice::pbc(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}



void Lattice::Metropolis()
{
  uniform_real_distribution<double> index(0, 1);
  uniform_int_distribution<int> dist(0, L_-1);
  
  for (int i=0; i<=N_; i++)
  {

    int x = dist(generator);
    int y = dist(generator);
    

    int dE =  2 * spin_matrix_(x,y)*
        (spin_matrix_(x,pbc(y,L_,-1))+
        spin_matrix_(pbc(x,L_,-1),y) +
        spin_matrix_(x,pbc(y,L_,+1)) +
        spin_matrix_(pbc(x,L_,1),y));
    
    double r = index(generator);
/*
    if (dE <= 0)
    {
      E_ += dE;
      M_ += 2*spin_matrix_(x,y);
      spin_matrix_(x,y) *= -1;     
    }
*/

     if(r<=boltzmann_[dE+8])
    {
      E_ += dE;
      M_ += 2*spin_matrix_(x,y);
      spin_matrix_(x,y) *= -1;
    }
  }
}




void Lattice::MCMC()
{
  double norm = 1./(MC_cycles_*N_);
  //cout << norm << endl;
  expectedE_ = 0.0;
  expectedE2_ = 0.0;
  expectedM_ = 0.0;
  expectedM2_ = 0.0;
  epsilon_ = 0.0;

  for (int cycles=1; cycles<=MC_cycles_; cycles++)
  {
    Metropolis();
    expectedE_ += E_;
    expectedE2_ += E_*E_;
    expectedM_ += abs(M_);
    expectedM2_ += M_*M_;
    epsilon_ += E_;
    E_dist_[cycles] = E_;
    
  }
  //cout << expectedE_ << endl;

  expectedE_ = expectedE_ * norm;
  expectedE2_ = expectedE2_ * norm;
  expectedM_ = expectedM_ * norm;
  expectedM2_ = expectedM2_ * norm;
  epsilon_ /= N_;
  cv_ = (expectedE2_*N_ - expectedE_ * expectedE_ *(N_*N_)) /(T_*T_*N_);
  chi_ = (expectedM2_*N_ - expectedM_ * expectedM_ *(N_*N_)) /(T_*T_*N_);

  
  //cout << expectedE_ << endl;
}



void Lattice::MCMC_burn_in_time(string filename_)
{
  string filename = filename_;
    ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec = 4;
    //ofile << "MC_cycle" << " " << "{epsilon}" << "{|m|}" << " " << "E" <<  endl;

  double norm = 1./(MC_cycles_*N_);
  double *eps = new double[MC_cycles_];
  double *mag = new double[MC_cycles_];
  //cout << norm << endl;
  expectedE_ = 0.0;
  expectedE2_ = 0.0;
  expectedM_ = 0.0;
  expectedM2_ = 0.0;


  for (int cycles=1; cycles<=MC_cycles_; cycles++)
  {
    Metropolis();
    expectedE_ += E_;
    expectedM_ += abs(M_);
    eps[cycles] = expectedE_ * norm;
    mag[cycles] = expectedM_ * norm;
    E_dist_[cycles] = E_;
    //cout << E_dist_[cycles] << endl;

    ofile << setw(width) << setprecision(prec) << scientific << cycles << " " << eps[cycles] << " " << mag[cycles] << " " << expectedE_ << endl;
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
  double cv = L_*L_/(T_*T_) * (epsilon2 - epsilon*epsilon);
  double chi = L_*L_/T_ * (magnetization2 - magnetization_abs*magnetization_abs);


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

  ofile << "Analytical results" << endl;
  ofile << "{epsilon}=" << setw(width) << setprecision(prec) << scientific << epsilon << " ";
  ofile << "{|m|}=" << setw(width) << setprecision(prec) << scientific << magnetization_abs << " ";
  ofile << "{C_v}=" << setw(width) << setprecision(prec) << scientific << cv << " ";
  ofile << "{chi}=" << setw(width) << setprecision(prec) << scientific << chi << endl << " "; 
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


