#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>

arma::vec u(arma::vec);


int main(){
    
    arma::vec x;
    
    x = arma::linspace(0, 1, 101);
    arma::vec res = u(x);
    std::cout << res <<"\n";

    std::string filename = "x_u.txt";
    
    std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec = 4;

    for (int i = 0; i<=100; i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x(i)
          << std::setw(width) << std::setprecision(prec) << std::scientific << res(i)
          << std::endl;
    }
    ofile.close();

    return 0;
}


arma::vec u(arma::vec x_in){
    return 1-(1-exp(-10))*x_in-exp(-10*x_in);
}