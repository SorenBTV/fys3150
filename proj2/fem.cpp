#include "functions.hpp"

int main(){

// Defining variables
int N = 10;
int n = N+1;
int L = 1;
double h = 1./(n);
double a = -1/(h*h);
double d = 2/(h*h);


//Defining matrix and calling function to produce matrix
arma::mat B = arma::mat(N, N);
arma::mat A = matrix(a, d, N, B);
//arma::mat A = arma::mat(N,N).randn();
//A = arma::symmatu(A);


//Producing eigenvalues and eigenvector with the Jacobi solver
arma::vec eigenvalues(N);
arma::mat eigenvectors = arma::mat(N, N, arma::fill::eye);
int maxiter = 10000;
int iterations = 0;
bool converged;
jacobi_eigensolver(A, 1e-8, eigenvalues, eigenvectors, maxiter, iterations, converged);
std::cout << "number of iterations = " << iterations << "\n";

return 0;
}