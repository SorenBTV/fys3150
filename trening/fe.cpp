#include <iostream>
#include <armadillo>
#include <fstream>

double f(double x_1); //Declaration
int main(){
    double x_1 = 0;
    double res = f(x_1);
    std::cout<<res <<"\n";
    return 0;
}


double f(double x_1){
    return 100*exp(-10*x_1);
}



