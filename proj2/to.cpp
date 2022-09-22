#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

#define pi 3.14159265359

// Declaration
arma::mat matrix(double a, double d, int N, arma::mat);
arma::vec ana_lambda(int N, double a, double d);
arma::mat ana_v(int N);
bool compare(arma::vec eigval, arma::vec ana_eigval, arma::mat eigvec, arma::mat ana_eigvec, int N);

int main(){
// Defining variables
int n = 7;
int N = n-1;
int L = 1;
double h = 1./(n);
double a = -1/(h*h);
double d = 2/(h*h);


//Defining matrix and calling function to produce matrix
arma::mat A = arma::mat(N, N);
arma::mat B = matrix(a, d, N, A);
//std::cout << B << "\n";

//Producing eigenvalues and eigenvector with armadillo
arma::vec eigval;
arma::mat eigvec;
arma::eig_sym(eigval, eigvec, B);
//std::cout << eigvec << "\n";

//Producing analytical eigenvalues and eigenvectors
arma::vec ana_eigval = ana_lambda(N, a, d);
arma::mat ana_eigvec = ana_v(N);
//std::cout << ana_eigvec << "\n";


// Comparing armadillo and analytical solutions
double comp = compare(eigval, ana_eigval, eigvec, ana_eigvec, N);
if(comp = true){std::cout << "Results agree" << "\n";}
else{std::cout << "Results don't agree" << "\n";}

return 0;
}


// Functions for making our matrix A, solving for analytical eigenvalues and eigenvectors,
// and for comparing them to armadillos eigenvector and eigenvalues

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

arma::vec ana_lambda(int N, double a, double d){
    arma::vec eig = arma::vec(N).fill(0.);

        for (int i=1 ; i<=N ; i++){
            eig[i-1] = d + 2*a*std::cos(i*pi/(N+1));
            }
    return eig;
}

arma::mat ana_v(int N){
    arma::mat v = arma::mat(N,N).fill(0);

    for(int i=1; i<=N; i++){
        arma::vec v_i = arma::vec(N).fill(0);
        for(int j=1; j<=N; j++){
            v_i[j-1] = std::sin((j*i*pi)/(N+1));
        }
    v.col(i-1) = arma::normalise(v_i);
    }
return v;
}

bool compare(arma::vec eigval, arma::vec ana_eigval, arma::mat eigvec, arma::mat ana_eigvec, int N){
    double tol = 1e-7;

    for (int i=0; i<=N; i++){
        double comp_val = std::abs(eigval(i)) - std::abs(ana_eigval(i)) > tol;

        for (int j=0; j<=N; j++){
            double comp_vec = std::abs(eigvec(i,j)) - std::abs(ana_eigvec(i,j)) > tol;

            if(comp_val, comp_vec = true){
                return true;
            }
            else{ return false;}
        }
    }
    if (true){return true;}
    else{return false;}
}