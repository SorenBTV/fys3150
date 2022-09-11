#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

arma::vec f(arma::vec x);
arma::vec x, v, x1;
arma::vec solve(arma::vec a, arma::vec b, arma::vec c, arma::vec g, int n);

int main(){

int n = 4;
double n1 = n+1;

arma::vec a = arma::linspace(-1,-1, n);
arma::vec b = arma::linspace(2,2, n);
arma::vec c = arma::linspace(-1,-1, n);
double h = 1/n1;


x = arma::linspace(0.2, 0.8, n);
arma::vec x1 = arma::linspace(0, 1, n+2);
arma::vec g = f(x)*h*h;
arma::vec res = solve(a, b, c, g, n);
v = arma::linspace(0,0,n+2);
for (int i = 0; i < n; i++){
    v[i+1] = res[i];
}
std::cout << v << "\n";

// Creation of .txt file with the outputs from our results.
std::string filename = "x_v.txt";

std::ofstream ofile;
ofile.open(filename);
int width = 12;
int prec = 4;

for (int i = 0; i<=100; i++){
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x1(i)
          << std::setw(width) << std::setprecision(prec) << std::scientific << v(i)
          << std::endl;
    }
    ofile.close();

return 0;
}

// Function f
arma::vec f(arma::vec x){
return 100*exp(-10*x);
}


arma::vec solve(arma::vec a, arma::vec b, arma::vec c, arma::vec g, int n) {
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    g[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        g[i] = (g[i] - a[i]*g[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    g[n] = (g[n] - a[n]*g[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        g[i] -= c[i]*g[i+1];
    }
return g;
}
