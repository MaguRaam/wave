#pragma once

#include "cell_divergence.h"

// GradientFlux:
// works for any polynomial but gradient(polynomial) should be defined:
namespace fv
{
  template <int dim>
  struct GradientFlux
  {

    // constructor:
    GradientFlux() = default;

    // face flux:
    template <typename CellPolynomial, typename Point, typename FaceNormal>
    double operator()(const CellPolynomial &owner_poly,
                      const CellPolynomial &neighbor_poly, const Point &point,
                      const FaceNormal &face_normal) const
    {
      // add static_assertions for dimension  check:

      double face_flux(0.0);
      for (unsigned int i = 0; i < FaceNormal::dimension; i++)
        face_flux +=
            0.5 *
            (gradient(owner_poly, point)[i] + gradient(neighbor_poly, point)[i]) *
            face_normal[i];
      return face_flux;
    }
  };

  // CellLaplace
  // This resembles the mathematical notation Divergence.Gradient = Laplacian
  template <int dim, int n_gpts>
  using CellLaplace = CellDivergence<dim, n_gpts, GradientFlux<dim>>;

} // namespace fv
