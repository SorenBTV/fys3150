#include "functions.hpp"

int main(){

// Defining variables
int N = 6;
int n = N+1;
int L = 1;
double h = 1./(n);
double a = -1/(h*h);
double d = 2/(h*h);


//Defining matrix and calling function to produce matrix
arma::mat B = arma::mat(N, N);
arma::mat A = matrix(a, d, N, B);
//std::cout << A << "\n";

//Producing eigenvalues and eigenvector with the Jacobi solver
arma::vec eigenvalues(N);
arma::mat eigenvectors(N,N);
int maxiter = 6000;
int iterations = 0;
bool converged;
jacobi_eigensolver(A, 1e-8, eigenvalues, eigenvectors, maxiter, iterations, converged);
//std::cout << abs(eigenvalues) << "\n";

//Producing analytical eigenvalues and eigenvectors
arma::vec ana_eigval = ana_lambda(N, a, d);
arma::mat ana_eigvec = ana_v(N);
//std::cout << abs(ana_eigval) << "\n";


// Comparing Jacobi and analytical solutions
bool comp = compare(eigenvalues, ana_eigval, eigenvectors, ana_eigvec, N);
if(comp == true){std::cout << "Results agree" << "\n";}
else{std::cout << "Results don't agree" << "\n";}

return 0;
}

