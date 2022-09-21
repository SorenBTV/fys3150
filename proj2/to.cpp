#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

// Declaration
arma::vec x, xh, v; 
arma::mat eigvec;
arma::mat matrix(double a, double d, int n, arma::mat);

int main(){
// Defining variables
int n = 6;
int L = 1;
double h = 1./n;
double a = -1/(h*h);
double d = 2/(h*h);

arma::vec x = arma::linspace(0, L, n);
arma::vec xh = arma::vec(n).fill(0.);
arma::vec v = arma::vec(n).fill(0.);

for (int i=0 ; i<=n; i++){
xh[i] = xh[0] + i * h;
}


//Defining matrix and calling function
arma::mat A = arma::mat(n, n);
arma::mat mat = matrix(a, d, n, A);
std::cout << mat << "\n";

arma::vec eigval;
arma::mat eigvec;
arma::eig_sym(eigval, eigvec, mat);

//std::cout << eigval << "\n";

return 0;
}

arma::mat matrix(double a, double d, int n, arma::mat){
//Producing an empty matrix
arma::mat A = arma::mat(n, n);

// Filling empty matrix
A(0, 0) = d;
A(0, 1) = a;
A(n-1,n-1) = d;
A(n-1, n-2) = a;

for (int i = 2; i < n; i++){
    A(i-1, i-2) = a;
    A(i-1, i-1) = d;
    A(i-1, i) = a;
}
return A;
}