#include <iostream>
#include <cmath>
#include <fstream>

double f(double x); //Declaration
main(){
    double x = {0 , 1};
    double res = f(x);
    return 0;
}


double f(double x){
    return 100*e**pow(-10*x);
}

