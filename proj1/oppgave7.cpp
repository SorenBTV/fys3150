#include<iostream>
#include<armadillo>
#include<string>
#include<fstream>
#include<iomanip>

arma::vec f(arma::vec);
void solve(int a, int b, int c, int n, arma::vec x, double h);


arma::vec f(arma::vec x){
  return 100*exp(-10*x);
}
void solve(int a, int b, int c, int n, arma::vec x, double h){
  arma::vec a_vec = arma::vec(n-1).fill(a);
  arma::vec b_vec = arma::vec(n).fill(b);
  arma::vec c_vec = arma::vec(n-1).fill(c);
  arma::vec v = arma::vec(n).fill(0);
  arma::vec g = f(x)*h*h;

  for (int i = 1; i < n; i++){
    b_vec[i] = b_vec[i] - a_vec[i]*c_vec[i-1]/b_vec[i-1];
    g[i] = g[i] - a_vec[i]*g[i-1]/b_vec[i-1];
  }

  v[n-1] = g[n-1]/b_vec[n-1];

  for (int i = n-1; i > 0; i--){
    v[i-1] = (g[i-1] - c_vec[i-1]*v[i])/b_vec[i-1];
  }

  arma::vec resles = arma::vec(n+2).fill(0);
  arma::vec x1 = arma::linspace(0, 1, n+2);

  for (int i = 1; i < (n+2); i++){
    resles[i] = v[i-1];
  }

  if (n==8) {
    std::string filename1 = "x_v10.txt";
    std::ofstream ofile;
    ofile.open (filename1);
    int width = 12;
    int prec = 4;
    for (int i = 0; i <= (n+1); i ++){
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x1(i)
    << std::setw(width) << std::setprecision(prec) << std::scientific << resles(i)
    << std::endl;
    }
    ofile.close();
  } else {
    int k = 0;
  }

  if (n==98) {
    std::string filename2 = "x_v100.txt";
    std::ofstream ofile;
    ofile.open (filename2);
    int width = 12;
    int prec = 4;
    for (int i = 0; i <= (n+1); i ++){
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x1(i)
    << std::setw(width) << std::setprecision(prec) << std::scientific << resles(i)
    << std::endl;
    }
    ofile.close();
  } else {
    int k1 = 0;
  }
  if (n==998) {
    std::string filename3 = "x_v1000.txt";
    std::ofstream ofile;
    ofile.open (filename3);
    int width = 12;
    int prec = 4;
    for (int i = 0; i <= (n+1); i ++){
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x1(i)
      << std::setw(width) << std::setprecision(prec) << std::scientific << resles(i)
      << std::endl;
    }
    ofile.close();
  } else {
    int k2 = 0;
  }
}

int main(){
  int a = -1;
  int b = 2;
  int c = -1;
  int n = 8;
  int n1 = 98;
  int n2 = 998;
  double h = 1./n;
  double h1 = 1./n1;
  double h2 = 1./n2;
  arma::vec x = arma::linspace(0+1./(n+1), 1-1./(n+1), n);
  arma::vec x1 = arma::linspace(0+1./(n1+1), 1-1./(n1+1), n1);
  arma::vec x2 = arma::linspace(0+1./(n2+1), 1-1./(n2+1), n2);
  solve(a, b, c, n , x, h);
  solve(a, b, c, n1 , x1, h1);
  solve(a, b, c, n2 , x2, h2);


/*  n = 100;
  arma::vec res2 = solve(a, b, c, n , x, h);

  resles = arma::vec(n+2).fill(0);

  for (int i = 1; i < (n+1); i++){
    resles[i] = res2[i-1];
  }
  std::string filename2 = "x_v5.txt";

  std::ofstream gfile;
  gfile.open (filename2);
  width = 12;
  prec = 4;
  for (int i = 0; i < (n+2); i ++){
  gfile << std::setw(width) << std::setprecision(prec) << std::scientific << x1(i)
  << std::setw(width) << std::setprecision(prec) << std::scientific << resles(i)
  << std::endl;
  }
  gfile.close();

  n = 1000;
  arma::vec res3 = solve(a, b, c, n , x, h);

  resles = arma::vec(n+2).fill(0);

  for (int i = 1; i < (n+1); i++){
    resles[i] = res3[i-1];
  }
  std::string filename3 = "x_v6.txt";

  std::ofstream rfile;
  rfile.open (filename3);
  width = 12;
  prec = 4;
  for (int i = 0; i < (n+2); i ++){
  rfile << std::setw(width) << std::setprecision(prec) << std::scientific << x1(i)
  << std::setw(width) << std::setprecision(prec) << std::scientific << resles(i)
  << std::endl;
  }
  rfile.close();
*/
  return 0;
}
