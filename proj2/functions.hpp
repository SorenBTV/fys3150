#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <math.h>
#define pi 3.14159265359

//Declaration
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);

void max_offdiag_symmetric_test();

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);

arma::mat matrix(double a, double d, int N, arma::mat);

arma::vec ana_lambda(int N, double a, double d);

arma::mat ana_v(int N);

bool compare(arma::vec eigenvalues, arma::vec ana_eigval, arma::mat eigenvectors, arma::mat ana_eigvec, int N);



//Functions
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

bool compare(arma::vec eigenvalues, arma::vec ana_eigval, arma::mat eigenvectors, arma::mat ana_eigvec, int N){
    double tol = 1e-8;

    for (int i=0; i<N; i++){
        bool comp_val = abs(ana_eigval(i)) - abs(eigenvalues(i)) < tol;
        for (int j=0; j<N; j++){
            bool comp_vec = abs(ana_eigvec(i,j)) - abs(eigenvectors(i,j)) < tol;
            if(comp_val == true){continue;}
            else{ return false;}
            if(comp_vec == true){continue;}
            else{ return false;}
        }
    }
    return true;
}

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l){
    int N = A.n_cols;
    double max_value = 0;

    for (int i=0; i<N; i++){

        for (int j = i+1; j<N; j++){

            // Changes max_value to the highest absolute value in the off-diagonal and 
            // sets k and l to the respective row and column of the max value
            if(std::abs(A(i,j)) >= max_value){
                max_value = std::abs(A(i,j));
                k = i;
                l = j;
            }
        }
    }
    return max_value;
    }

void max_offdiag_symmetric_test(){
    //Setting up testing matrix
    arma::mat A = arma::mat(4, 4, arma::fill::eye);
    A(1, 2) = -0.7;
    A(2, 1) = -0.7;
    A(0, 3) = 0.5;
    A(3, 0) = 0.5;
    int k,l;
    double max_value;

    //Calling function and checking the correct answer against a tolerance
    max_value = max_offdiag_symmetric(A, k, l);
    //std::cout << max_value << "\n";
    assert(max_value - 0.7 <= 1e-7);
    assert(k == 1);
    assert(l == 2);
}

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    double eps = 1e-8;
    int N = (int)A.n_cols;
        double t, c, s, tau;

        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if(tau < 0){t = -1/(-tau + sqrt(1 + tau * tau));}
        else{t = 1/(tau + sqrt(1 + tau * tau));}
        //if(tau > 0){t = 1/(tau + sqrt(1 + tau * tau));}

        c = 1/(sqrt(1 + t * t));
        s = c * t;


        double Am = A(k,k);

        A(k,k) = A(k,k) * c * c - 2*A(k,l) * c * s + A(l,l) * s * s;
        A(l,l) = A(l,l) * c * c + 2 * A(k,l) * c * s + Am * s * s;
        A(k,l) = 0;
        A(l,k) = 0;

        for (int i=0; i<N; i++)
        {
            if(i != k && i != l){
                double Am = A(i,k);
                A(i,k) = A(i,k) * c - A(i,l) * s;
                A(k,i) = A(i,k);
                A(i,l) = A(i,l) * c + Am * s;
                A(l,i) = A(i,l);
            }
            double Rm = R(i,k);
            R(i,k) = R(i,k) * c - R(i,l) * s;
            R(i,l) = R(i,l) * c + Rm * s;    
        }  
    R = arma::normalise(R);
}

void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged){

    int k,l;
    double max_value = max_offdiag_symmetric(A,k,l);
    
    while(max_value >= eps){
        jacobi_rotate(A, eigenvectors, k, l);
        iterations += 1;
        max_value = max_offdiag_symmetric(A,k,l);
    }

    //Stops if the number of iterations reaches maxiter
    assert(iterations < maxiter);
    converged = true;
    
    //Defining eigenvalues
    for (int i=0; i<(int)A.n_cols; i++){
        eigenvalues(i) = A(i,i);
    }

    arma::uvec sort = arma::sort_index(eigenvalues);
    arma::vec tempval = eigenvalues;
    arma::mat tempvec = eigenvectors;
    for (int i=0; i<(int)A.n_cols; i++){
        eigenvalues(i) = tempval(sort(i));
        eigenvectors.col(i) = tempvec.col(sort(i));
    }
}