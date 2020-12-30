#include "../include/Weno432.h"

// Construtor for the WENO4 class

Weno4_2D::Weno4_2D(double ft, double cfl_no, unsigned int n_cell, double c)
    : dt(0.0),
      finalTime(ft),
      cfl(cfl_no),
      cell(n_cell),
      C(c),
      fv(0),
      dof_handler(triangulation)
{
}
