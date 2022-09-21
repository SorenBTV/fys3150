#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

#define pi 3.14159265359

// Declaration
arma::vec x, xh, v; 
arma::mat eigvec;
arma::mat matrix(double a, double d, int N, arma::mat);
arma::vec ana_eigval(int N, double a, double d);

int main(){
// Defining variables
int n = 6;
int N = n-1;
int L = 1;
double h = 1./(n);
double a = -1/(h*h);
double d = 2/(h*h);

arma::vec x = arma::linspace(0, L, n);
arma::vec xh = arma::vec(n).fill(0.);
arma::vec v = arma::vec(n).fill(0.);

for (int i=0 ; i<=n; i++){
xh[i] = xh[0] + i * h;
}


//Defining matrix and calling function
arma::mat A = arma::mat(N, N);
arma::mat mat = matrix(a, d, N, A);
//std::cout << mat << "\n";

arma::vec eigval;
arma::mat eigvec;
arma::eig_sym(eigval, eigvec, mat);

std::cout << eigval << "\n";

arma::vec ana_eigv = ana_eigval(n, a, d);
std::cout << ana_eigv << "\n";

return 0;
}

arma::mat matrix(double a, double d, int N, arma::mat){

//Producing an empty matrix
arma::mat A = arma::mat(N, N);

// Filling empty matrix
A(0, 0) = d;
A(0, 1) = a;
A(N-1,N-1) = d;
A(N-1, N-2) = a;

for (int i = 2; i < N; i++){
    A(i-1, i-2) = a;
    A(i-1, i-1) = d;
    A(i-1, i) = a;
}
return A;
}

arma::vec ana_eigval(int N, double a, double d){
    arma::vec eig = arma::vec(N).fill(0);

        for (int i=1 ; i<=N ; i++){
            eig[i] = d + 2*a*std::cos(i*pi/(N+1));
            }
    return eig;
}