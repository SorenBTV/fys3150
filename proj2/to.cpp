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

//Producing eigenvalues and eigenvector with armadillo
arma::vec eigval;
arma::mat eigvec;
arma::eig_sym(eigval, eigvec, A);
//std::cout << eigvec << "\n";

//Producing analytical eigenvalues and eigenvectors
arma::vec ana_eigval = ana_lambda(N, a, d);
arma::mat ana_eigvec = ana_v(N);
//std::cout << ana_eigvec << "\n";


// Comparing armadillo and analytical solutions
double comp = compare(eigval, ana_eigval, eigvec, ana_eigvec, N);
if(comp == true){std::cout << "Results agree" << "\n";}
else{std::cout << "Results don't agree" << "\n";}

return 0;
}