#include <iostream>
#include <armadillo>

arma::vec x = arma::vec(5).fill(2.);
//double u(double x); //Declaration

int main(){
   //double n = 5;
   double res =  x(0,5);
   std::cout<<res<<"/n";
   return 0;

}