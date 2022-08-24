#include <iostream>
#include <cmath>
#include <fstream>

double f(double x); //Declaration
main(){
    std::vector<double> x = {0 , 1};
    double res = f(x);
    return 0;
}


double f(double x){
    return 100*exp(-10*x);
}

