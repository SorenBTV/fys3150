#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

// Declaration
arma::mat matrix(int a, int b, int c,int n, arma::mat);
double solve(int a, int b, int c, int n, arma::mat); 

int main(){
// Defining variables
int n = 4;
int a = -1;
int b = 2;
int c = -1;

//Defining matrix and calling function
arma::mat A = arma::mat(n, n);
arma::mat mat = matrix(a, b, c, n, A);
std::cout << mat << "\n";



return 0;

}

arma::mat matrix(int a, int b, int c, int n, arma::mat){
//Producing an empty matrix
arma::mat A = arma::mat(n, n);

// Filling empty matrix
A(0, 0) = b;
A(0, 1) = c;
A(n-1,n-1) = b;
A(n-1, n-2) = a;

for (int i = 2; i < n; i++){
    A(i-1, i-2) = a;
    A(i-1, i-1) = b;
    A(i-1, i) = c;
}
return A;
}

double solve(int a, int b, int c, int n, arma::mat A){

// backward substitution

}
