#include "../include/Weno432.h"


// Evaluate the WENO polynomial at a point 


double evaluate_weno_gradient_y(Vector<double> coeffs, Vector<double> consts,  Point<2> P) {
    
    double x = P(0); double y = P(1); 
    
    double x0 = consts(0); double y0 = consts(1);
    
    double retval = coeffs(2) + 2*coeffs(4)*(y-y0) + coeffs(5)*(x-x0) + 3.0*coeffs(7)*(y-y0)*(y-y0) + coeffs(8)*(x-x0)*(x-x0) + 
                    coeffs(9)*(x-x0)*2.0*(y-y0);


    return retval; 
}
