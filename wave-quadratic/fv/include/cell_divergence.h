#pragma once

// dealii headers:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

/*CellDivergence:
works for both interior and boundary cell:
The divergence functor is independent of Global_Polynomial type and Flux functor
so both are templated: The Global_Polynomial is composed of CellPolynomials eg
GlobalPolynomial = std::vector<CellPolynomial> The Flux functor should have the
following signature:
*/

namespace fv
{
  template <int dim, int n_gpts, typename Flux>
  struct CellDivergence
  {

    // iterator type:
    using CellIter = typename dealii::DoFHandler<dim>::active_cell_iterator;

    // default constructor:
    CellDivergence()
        : n_faces(dealii::GeometryInfo<dim>::faces_per_cell), fv(0),
          face_quadrature_formula(n_gpts),
          fv_face_values(fv, face_quadrature_formula, dealii::update_quadrature_points | dealii::update_normal_vectors | dealii::update_JxW_values),
          n_face_q_points(face_quadrature_formula.size()) {}

    // operator() takes cell iterator and Global polynomial as input:
    template <typename GlobalPolynomial>
    double operator()(const CellIter &cell, const GlobalPolynomial &poly)
    {

      // initialize divergence term:
      double divergence(0.0);

      // initialize owner and neighbor index:
      unsigned int owner = cell->active_cell_index(), neighbor = 0;

      // loop over faces and add all the face flux term:
      // the face flux is computed using the flux functor:
      for (unsigned int f = 0; f < n_faces; ++f)
      {
        // initialize normal vector and quad points:
        fv_face_values.reinit(cell, f);

        // neighbor index:
        neighbor = cell->neighbor_index(f);

        // initialize face flux:
        double face_flux(0.0);

        // check for boundary face:
        if (cell->face(f)->at_boundary())
        {
          // loop over quadrature points on the face and compute flux:
          for (int q = 0; q < n_face_q_points; ++q)
            face_flux += flux(poly[owner], poly[owner], fv_face_values.quadrature_point(q), fv_face_values.normal_vector(q)) * fv_face_values.JxW(q);
        }

        // interior face:
        else
        {
          // loop over quadrature points on the face and compute flux:
          for (int q = 0; q < n_face_q_points; ++q)
            face_flux += flux(poly[owner], poly[neighbor], fv_face_values.quadrature_point(q), fv_face_values.normal_vector(q)) * fv_face_values.JxW(q);
        }

        // add face flux to divergence:
        divergence += face_flux;
      }

      // Divide it by cell volume:
      divergence *= (1.0 / cell->measure());
      return divergence;
    }

  private:
    Flux flux;
    const unsigned int n_faces;
    const dealii::FE_DGQ<dim> fv;
    const dealii::QGauss<dim - 1> face_quadrature_formula;
    dealii::FEFaceValues<dim> fv_face_values;
    const int n_face_q_points;
  };

} // namespace fv
