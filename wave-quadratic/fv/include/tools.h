#pragma once

#include "cell_average.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iomanip>

namespace fv
{
  // Linfty error:
  template <int dim, int n_gpts, typename Function, typename Vector>
  inline auto linfty_error(const Vector &u, Function u_exact,
                           const dealii::DoFHandler<dim> &dof_handler)
  {

    Vector error(dof_handler.n_dofs());
    CellAverage<dim, n_gpts> cell_average;
    unsigned int i(0);
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      error[i] = abs(cell_average(u_exact, cell) - u[i]);
      i++;
    }
    return *std::max_element(error.begin(), error.end());
  }

  // L2 error:
  template <int dim, int n_gpts, typename Function, typename Vector>
  inline auto l2_error(const Vector &u, Function u_exact,
                       const dealii::DoFHandler<dim> &dof_handler)
  {
    Vector error(dof_handler.n_dofs());
    CellAverage<dim, n_gpts> cell_average;
    unsigned int i(0);
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      error[i] = pow((cell_average(u_exact, cell) - u[i]), 2.0) * cell->measure();
      i++;
    }
    return sqrt(std::accumulate(error.begin(), error.end(), 0.0));
  }

} // namespace fv
