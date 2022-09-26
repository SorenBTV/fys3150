#include "functions.hpp"

int main(int argc, char* argv[]){

// Defining variables
int n = atoi(argv[1]);
int N = n-1;
int L = 1;
double h = 1./(n);
double a = -1/(h*h);
double d = 2/(h*h);

arma::vec xh = arma::linspace(0, 1, n+1);

//Defining matrix and calling function to produce matrix
arma::mat B = arma::mat(N, N);
arma::mat A = matrix(a, d, N, B);


//Producing eigenvalues and eigenvector with the Jacobi solver
arma::vec eigenvalues(N);
arma::mat eigenvectors = arma::mat(N, N, arma::fill::eye);
int maxiter = 100000;
int iterations = 0;
bool converged;
jacobi_eigensolver(A, 1e-8, eigenvalues, eigenvectors, maxiter, iterations, converged);

//Setting start and end-point of eigenvectors to 0
arma::mat eig = arma::mat(n+1, 3).fill(0);
for (int i=0; i<3; i++){
    for (int j=1; j<n; j++){
        eig(j,i) = eigenvectors(i,j-1);
    }
}

 // Creation of .txt file with the outputs from our results.
    std::string filename = "xhat_eigenvectors100.txt";
    
    std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec = 4;

    for (int i = 0; i<=n; i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << xh(i) << std::endl;
    }
    for (int i = 0; i<3; i++){ofile << std::setw(width) << std::setprecision(prec) << std::scientific << eigenvalues(i)*eig.col(i);

    }
    ofile.close();

return 0;
}

