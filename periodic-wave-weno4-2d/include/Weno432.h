#ifndef WENO432_H_
#define WENO432_H_

#include "Headers.h"
#include "CLS.h"
#include "LU.h"
#include "cell_properties.h"
#include "gradient_flux.h"

// Complementary Function Definitions

double compute_second_order_smoothness_indicator(Vector<double>, Vector<double>, double);
double compute_third_order_smoothness_indicator(Vector<double>, Vector<double>, double);
double compute_fourth_order_smoothness_indicator(Vector<double>, Vector<double>, double);

double evaluate_weno_gradient_x(Vector<double>, Vector<double>, Point<2>);
double evaluate_weno_gradient_y(Vector<double>, Vector<double>, Point<2>);

// Main Class Declaration

class Weno4_2D
{

    void make_grid();
    void setup_system();
    void compute_cell_properties();
    void compute_weno_polynomial_constants();
    void initialize();
    void precompute_matrices();
    void compute_IS_constants();
    void reconstruct();
    void copy_data();
    void compute_rhs();
    void compute_time_step_based_on_cfl_number(double);
    void L_norm(const double);
    void solve_euler();
    void solve_ssprk33();
    void solve_ssprk54();
    void output_results(unsigned int);

    double initial_condition_u(Point<2>);
    double initial_condition_v(Point<2>);
    double exact_solution(Point<2>, double);

    Triangulation<2> triangulation;

    //solution vector and rhs:
    Vector<double> U;
    Vector<double> V;

    Vector<double> rhs1;
    Vector<double> rhs2;

    //exact solution:
    Vector<double> Uexact;

    //solution error:
    Vector<double> local_difference;

    // Coefficients for WENO polynomials
    std::vector<Vector<double>> coeffs_U;

    // WENO polynomial constants (only depend on mesh)
    std::vector<Vector<double>> WENO_poly_consts;
    std::vector<Vector<double>> IS_constants;

    // Fourth Order Stencil
    std::vector<Constrained_LS> CLS_R4;

    // Third Order Stencils
    std::vector<bool> is_admissible_R3; // Centered third order stencil stencil
    std::vector<Constrained_LS> CLS_R3;

    std::vector<bool> is_admissible_R31; // Third order stencil 2
    std::vector<Constrained_LS> CLS_R31;

    std::vector<bool> is_admissible_R32; // Third order stencil 4
    std::vector<Constrained_LS> CLS_R32;

    std::vector<bool> is_admissible_R33; // Third order stencil 5
    std::vector<Constrained_LS> CLS_R33;

    std::vector<bool> is_admissible_R34; // Third order stencil 7
    std::vector<Constrained_LS> CLS_R34;

    // Second order stencils

    std::vector<LUdcmp> LU_R21;
    std::vector<LUdcmp> LU_R22;
    std::vector<LUdcmp> LU_R23;
    std::vector<LUdcmp> LU_R24;

    unsigned int no_cells_per_block;

    std::vector<cell_properties> Cell;

    double dt;
    double finalTime;
    double cfl;
    double h_min;

    unsigned int cell;

    //wave parameters:
    double kx, ky, omega, C;

    FE_DGQ<2> fv;
    DoFHandler<2> dof_handler;

public:
    Weno4_2D(double, double, unsigned int, double, double, double, double);
    void run();
};

#endif /* WENO43_H_ */
