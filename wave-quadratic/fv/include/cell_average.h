#pragma once

// dealii headers:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

/*Cell Average:
Functor to compute cell average of a given function on a cell:

CellAverage is independent of the function so it is templated:

The Function should have the following signature:
        template<typename Point>
        auto operator()(const Point& p);
*/

/*
   We can use this functor to compute cell average of scalar functions ,
   vector functions and Tensor functions:
   TODO : provide examples:
*/

namespace fv
{
  template <int dim, int n_gpts>
  struct CellAverage
  {

    // default constructor:
    CellAverage()
        : quadrature_formula(n_gpts), fv(0),
          fv_values(fv, quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_JxW_values),
          TotalNgpts(quadrature_formula.size()) {}

    // operator() that computes cell average by taking a function-object and cell
    // iterator: Note : The cell iterator is templated so that it can take both
    // triangulation and dofhandler cell iterator:

    template <typename F, typename CellIter>
    auto operator()(const F &f, const CellIter &cell)
    {
      fv_values.reinit(cell);

      // initialize value and also deduce the type of value:
      auto value = f(fv_values.quadrature_point(0));
      value *= fv_values.JxW(0);

      // loop over quadrature points:
      for (int q = 1; q < TotalNgpts; q++)
      {
        auto tmp = f(fv_values.quadrature_point(q));
        tmp *= fv_values.JxW(q);
        value += tmp;
      }

      // divide it by cell volume
      value /= cell->measure();
      return value;
    }

  private:
    const dealii::QGauss<dim> quadrature_formula;
    const dealii::FE_DGQ<dim> fv;
    dealii::FEValues<dim> fv_values;
    const int TotalNgpts;
  };

  /* project:
Loops over the given range of cells and computes the cell average of a given
function on each cell and writes the output in another container:
The interface is similar to STL transform:*/

  template <int dim, int n_gpts, typename BeginCellItr, typename EndCellItr,
            typename OutputItr, typename Function>
  inline void project(BeginCellItr firstc, EndCellItr lastc, OutputItr first,
                      Function f)
  {

    CellAverage<dim, n_gpts> cell_average;
    for (; firstc != lastc; ++firstc, ++first)
      *first = cell_average(f, firstc);
  }
} // namespace fv
