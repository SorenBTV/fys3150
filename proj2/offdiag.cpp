#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>

//Declaration
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);
void max_offdiag_symmetric_test();


int main(){

//Calling test function
max_offdiag_symmetric_test();

return 0;
}

//Functions

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l){
    int N = A.n_cols;
    double max_value = 0;

    for (int i=0; i<N; i++){
        for (int j=i+1; j<N; j++){

            if(i==j){max_value = max_value;}

            else{
                if(std::abs(A(i,j)) >= max_value){
                    max_value = std::abs(A(i,j));
                    k=i;
                    l=j;}}
            }}
    return max_value;
    }

void max_offdiag_symmetric_test(){
    arma::mat A = arma::mat(4, 4, arma::fill::eye);
    A(1, 2) = -0.7;
    A(2, 1) = -0.7;
    A(0, 3) = 0.5;
    A(3, 0) = 0.5;
    int k,l;
    double max_value;

    max_value = max_offdiag_symmetric(A, k, l);
    std::cout << max_value << "\n";
    assert(max_value - 0.7 <= 1e-7);
    assert(k == 1);
    assert(l == 2);
}