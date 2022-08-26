#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

arma::vec u(arma::vec);


int main(){
    
    // std::string filename = "x_u.txt";
    
    // std::ofstream ofile;
    // ofile.open(filename);
    
    arma::vec x;
    
    x = arma::linspace(0, 1, 100);
    arma::vec res = u(x);
    std::cout << res <<"\n";
    return 0;
}


arma::vec u(arma::vec x_in){
    return 1-(1-exp(-10)*x_in)-exp(-10*x_in);
}