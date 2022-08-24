#include <iostream>
#include <cmath>

double f(double x, double y); //Declaration

int main() {
    double x = 2;
    double y = 1;
    double res = f(x, y); //Should return 2*2 = 4.
    return 0;
}

double f(double x, double y){
    //Provide definition of function here:
    return x*x+y;
}