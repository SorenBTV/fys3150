#include <iostream>
#include <armadillo>

arma::vec u(arma::vec);


int main(){
    
    arma::vec x;
    
    x = arma::linspace(0, 1, 6);
    arma::vec res = u(x);
    std::cout << res <<"\n";
    return 0;
}


arma::vec u(arma::vec x_in){
    return 1-(1-exp(-10)*x_in)-exp(-10*x_in);
}