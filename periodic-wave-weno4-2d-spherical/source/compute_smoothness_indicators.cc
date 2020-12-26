#include "../include/Weno432.h"

// Find smoothness indicators 

// Second order polynomial 

double compute_second_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS, double h) {

    double u_x   = coeffs(0); double u_y   = coeffs(1); 
    
    double const_term = IS(0)*(h*u_x*u_x + h*u_y*u_y);
    
    return (const_term); 
}

// Third Order polynomial 

double compute_third_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS, double h) {

    double h3 = h*h*h; 
    
    double u_x   = coeffs(0); double u_y   = coeffs(1); double u_xx  = coeffs(2); double u_yy  = coeffs(3);
    double u_xy  = coeffs(4); 
    
    double const_term = IS(0)*(4*h3*u_xx*u_xx + h3*u_xy*u_xy + 4*h3*u_yy*u_yy + h*u_x*u_x + h*u_y*u_y);

    double x_term = IS(1)*(4*h*u_x*u_xx + 2*h*u_xy*u_y);

    double y_term = IS(2)*(2*h*u_x*u_xy + 4*h*u_y*u_yy); 

    double x2_term = IS(3)*(4*h*u_xx*u_xx + h*u_xy*u_xy);

    double y2_term = IS(4)*(h*u_xy*u_xy + 4*h*u_yy*u_yy); 

    double xy_term = IS(5)*(4*h*u_xx*u_xy + 4*h*u_xy*u_yy);
    
    return (const_term + x_term + y_term + x2_term + y2_term + xy_term); 
}

// Fourth order polynomial 

double compute_fourth_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS, double h) {
    
    double h3 = h*h*h; double h5 = h*h*h*h*h;
    
    double u_x   = coeffs(0); double u_y   = coeffs(1); double u_xx  = coeffs(2); double u_yy  = coeffs(3);
    double u_xy  = coeffs(4); double u_xxx = coeffs(5); double u_yyy = coeffs(6); double u_xxy = coeffs(7);
    double u_xyy = coeffs(8);
    
    double const_term = IS(0)*(36*h5*u_xxx*u_xxx + 4*h5*u_xxy*u_xxy + 4*h5*u_xyy*u_xyy + 36*h5*u_yyy*u_yyy + 
                               4*h3*u_xx*u_xx + h3*u_xy*u_xy + 4*h3*u_yy*u_yy + h*u_x*u_x + h*u_y*u_y);

    double x_term = IS(1)*(24*h3*u_xx*u_xxx + 4*h3*u_xxy*u_xy + 8*h3*u_xyy*u_yy + 4*h*u_x*u_xx + 2*h*u_xy*u_y);

    double y_term = IS(2)*(8*h3*u_xx*u_xxy + 4*h3*u_xy*u_xyy + 24*h3*u_yy*u_yyy + 2*h*u_x*u_xy + 4*h*u_y*u_yy); 

    double x2_term = IS(3)*(36*h3*u_xxx*u_xxx + 4*h3*u_xxy*u_xxy + 4*h3*u_xyy*u_xyy + 6*h*u_x*u_xxx + 
                             4*h*u_xx*u_xx + 2*h*u_xxy*u_y + h*u_xy*u_xy);

    double y2_term = IS(4)*(4*h3*u_xxy*u_xxy + 4*h3*u_xyy*u_xyy + 36*h3*u_yyy*u_yyy + 
                            2*h*u_x*u_xyy + h*u_xy*u_xy + 6*h*u_y*u_yyy + 4*h*u_yy*u_yy); 

    double xy_term = IS(5)*(24*h3*u_xxx*u_xxy + 8*h3*u_xxy*u_xyy + 24*h3*u_xyy*u_yyy + 4*h*u_x*u_xxy + 
                             4*h*u_xx*u_xy + 4*h*u_xy*u_yy + 4*h*u_xyy*u_y);

    double x3_term = IS(6)*(12*h*u_xx*u_xxx + 2*h*u_xxy*u_xy); 

    double y3_term = IS(7)*(2*h*u_xy*u_xyy + 12*h*u_yy*u_yyy); 

    double x2y_term = IS(8)*(8*h*u_xx*u_xxy + 6*h*u_xxx*u_xy + 4*h*u_xxy*u_yy + 4*h*u_xy*u_xyy);

    double xy2_term = IS(9)*(4*h*u_xx*u_xyy + 4*h*u_xxy*u_xy + 6*h*u_xy*u_yyy + 8*h*u_xyy*u_yy); 

    double x4_term = IS(10)*(9*h*u_xxx*u_xxx + h*u_xxy*u_xxy);

    double y4_term = IS(11)*(h*u_xyy*u_xyy + 9*h*u_yyy*u_yyy); 

    double x2y2_term = IS(12)*(6*h*u_xxx*u_xyy + 4*h*u_xxy*u_xxy + 6*h*u_xxy*u_yyy + 4*h*u_xyy*u_xyy); 

    double x3y_term = IS(13)*(12*h*u_xxx*u_xxy + 4*h*u_xxy*u_xyy); 

    double xy3_term = IS(14)*(4*h*u_xxy*u_xyy + 12*h*u_xyy*u_yyy);
    
    return (const_term + x_term + y_term + x2_term + y2_term + xy_term + x3_term + 
           y3_term + x2y_term + xy2_term + x4_term + y4_term + x2y2_term + x3y_term + xy3_term); 
    
} 
