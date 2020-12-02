#pragma once

// A small Generic Finite-Volume Library written using dealii:
/*
        ~Functors and Functions:
        CellAverage
        CellDivergence
        GradientFlux
        CellLaplace
        CellPolyMatrix
        CellPolyVector
        Reconstruct

        ~Types:
        CellLinearPoly
        GlobalPoly
*/

// TODO Need to add assertions for safety:

// c++ headers:
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <functional>

// dealii headers:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>

using namespace dealii;

template <int dim>
using CellIterator = typename DoFHandler<dim>::active_cell_iterator;

// Computes cell average:
template <int dim, int n_gpts = 4> struct CellAverage {
  // Default constructor:
  CellAverage()
      : quadrature_formula(n_gpts), fv(0),
        fv_values(fv, quadrature_formula,
                  update_quadrature_points | update_JxW_values),
        TotalNgpts(fv_values.get_quadrature().size()) {}

  template <typename Function>
  auto operator()(Function f, const CellIterator<dim> &cell) {
    fv_values.reinit(cell);
    auto value{fv_values.JxW(0) * f(fv_values.quadrature_point(0))};
    for (unsigned int q = 1; q < TotalNgpts; q++)
      value += fv_values.JxW(q) * f(fv_values.quadrature_point(q));
    return value / cell->measure();
  }

private:
  // dealii objects:
  const QGauss<dim> quadrature_formula;
  const FE_DGQ<dim> fv;
  FEValues<dim> fv_values;
  const unsigned int TotalNgpts;
};







// Cell Divergence operator:
// works for both interior and boundary cell:
// works for arbirtary polynomial degrees:
// the flux functor should take two polynomial objects and face normal
// TODO: We use only one gauss point for finding face flux
// extend it for arbirtary number of gauss points:

template <int dim, typename Flux> struct CellDivergence {
  CellDivergence()
      : fv(0), face_quadrature_formula(1),
        fv_face_values(fv, face_quadrature_formula, update_normal_vectors) {}

  template <typename GlobalPolynomial>
  double operator()(const CellIterator<dim> &cell,
                    const GlobalPolynomial &poly) {

    double divergence(0.0);

    for (unsigned int f = 0; f < n_faces; ++f) {
      fv_face_values.reinit(cell, f);
      face_normal = cell->face(f)->measure() * fv_face_values.normal_vector(0);

      // check for boundary face:
      if (cell->face(f)->at_boundary())
        divergence += flux(poly[cell->active_cell_index()], face_normal);
      // interior face:
      else
        divergence += flux(poly[cell->active_cell_index()],
                           poly[cell->neighbor_index(f)], face_normal);
    }
    divergence *= (1.0 / cell->measure());
    return divergence;
  }

private:
  Flux flux;
  const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
  const FE_DGQ<dim> fv;
  const QGauss<dim - 1> face_quadrature_formula;
  FEFaceValues<dim> fv_face_values;
  Tensor<1, dim> face_normal;
};

// Gradient flux at face:
// works for any polynomial degree but gradient(polynomial) shouldbe defined:
struct GradientFlux {

  // Interior face flux:
  template <typename CellPolynomial, typename FaceNormal>
  double operator()(const CellPolynomial &owner_poly,
                    const CellPolynomial &neighbor_poly,
                    const FaceNormal &face_normal) {
    double face_flux(0.0);
    for (unsigned int i = 0; i < FaceNormal::dimension; i++)
      face_flux += 0.5 *
                   (gradient(owner_poly)[i] + gradient(neighbor_poly)[i]) *
                   face_normal[i];
    return face_flux;
  }

  // Boundary face flux:
  template <typename CellPolynomial, typename FaceNormal>
  double operator()(const CellPolynomial &owner_poly,
                    const FaceNormal &face_normal) {
    double face_flux(0.0);
    for (unsigned int i = 0; i < FaceNormal::dimension; i++)
      face_flux += gradient(owner_poly)[i] * face_normal[i];
    return face_flux;
  }
};

// Cell Laplace operator:
template <int dim> using CellLaplace = CellDivergence<dim, GradientFlux>;

// CellPolyMatrix:
// Returns Linear polynomial matrix for a given cell iterator:
// The functor is purely based on the geometry of the cell:
// This works for both interior and boundary cells:

// TODO Works only for Linear polynomials:
template <int dim, int n_gpts = 4> struct CellPolyMatrix {
  CellPolyMatrix() : A{GeometryInfo<dim>::faces_per_cell, dim} {}

  FullMatrix<double> operator()(const CellIterator<dim> &cell) {

    auto xo = cell->barycenter();
    auto f = [&xo](auto x) { return x - xo; };
    Tensor<1, dim, double> row;

    // fill rows of coefficient matrix by looping over faces:
    for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i) {

      // check for boundary face for a given cell:
      if (cell->face(i)->at_boundary())
        row = f(cell->face(i)->center());
      else
        row = integrate(f, cell->neighbor(i));

      // fill columns:
      for (unsigned int j = 0; j < dim; j++)
        A(i, j) = row[j];
    }
    return A;
  }

private:
  FullMatrix<double> A;
  CellAverage<dim, n_gpts> integrate;
};

// TODO Works only for Linear polynomials:
// CellPolyVector:
template <int dim, typename GlobalVector, typename BoundaryMap>
inline void cell_poly_vector(const CellIterator<dim> &cell,
                             const GlobalVector &u, BoundaryMap &u_bcs,
                             Vector<double> &b) {

  // Loop over cell neighbors or faces:
  for (unsigned int i = 0; i < b.size(); ++i) {
    if (cell->face(i)->at_boundary())
      b(i) = u_bcs[cell->face(i)] - u(cell->active_cell_index());
    else
      b(i) = u(cell->neighbor_index(i)) - u(cell->active_cell_index());
  }
}

// TODO Works only for Linear polynomials:
// Reconstruct polynomial:
template <int dim> struct Reconstruct {
  Reconstruct() = default;

  // Create Polynomial matrix for all cells:
  void create_matrix(const DoFHandler<dim> &dof_handler) {
    CellPolyMatrix<dim> cell_matrix;
    for (const auto &cell : dof_handler.active_cell_iterators())
      A.push_back(cell_matrix(cell));
  }

  // Polynomial Reconstruct operator:
  template <typename GlobalVector, typename BoundaryMap,
            typename GlobalPolynomial>
  void operator()(const DoFHandler<dim> &dof_handler, const GlobalVector &u,
                  BoundaryMap &u_bcs, GlobalPolynomial &upoly) {

    Householder<double> householder;
    Vector<double> b(GeometryInfo<dim>::faces_per_cell);
    Vector<double> x(dim);
    unsigned int i = 0;
    for (const auto &cell : dof_handler.active_cell_iterators()) {

      householder.initialize(A[i]);             // initialize A
      cell_poly_vector<dim>(cell, u, u_bcs, b); // initialize b
      // solve Ax = b using Least-Square and compute polynomial coefficeints for
      // each cell
      householder.least_squares(x, b);
      upoly[i].coefficient = x;
      i++;
    }
  }

private:
  std::vector<FullMatrix<double>> A;
};

// Linfty error:
template <int dim, int n_gpts = 4, typename Function, typename Vector>
inline auto linfty_error(const Vector &u, Function u_exact,
                  const DoFHandler<dim> &dof_handler) {
  Vector error(dof_handler.n_dofs());
  CellAverage<dim, n_gpts> cell_average;
  unsigned int i(0);
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    error[i] = abs(cell_average(u_exact, cell) - u[i]);
    i++;
  }
  return *std::max_element(error.begin(), error.end());
}

// L2 error:
template <int dim, int n_gpts = 4, typename Function, typename Vector>
inline auto l2_error(const Vector &u, Function u_exact,
                  const DoFHandler<dim> &dof_handler) {
  Vector error(dof_handler.n_dofs());
  CellAverage<dim, n_gpts> cell_average;
  unsigned int i(0);
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    error[i] = pow((cell_average(u_exact, cell) - u[i]),2.0)*cell->measure();
    i++;
  }
  return sqrt(std::accumulate(error.begin(),error.end(),0.0));
}


// L2 norm of a vector:
template <typename Vector> 
inline auto L2norm(const Vector &v) {  
  return sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0));
}


 

//// Rudimentary Cell Linear Polynomial : TODO throw away soon:
// TODO we rely on compiler generated constructors and assignment operators:
template <int dim> struct CellLinearPoly {
  Vector<double> coefficient{dim};
  friend const Vector<double> &gradient(const CellLinearPoly &poly) {
    return poly.coefficient;
  }
  friend Vector<double> &gradient(CellLinearPoly &poly) {
    return poly.coefficient;
  }
};

// global polynomial defined on a grid:
template <int dim> using GlobalPolynomial = std::vector<CellLinearPoly<dim>>;

// Boundary condition type:
template <int dim>
using BoundaryCondition =
    std::map<typename Triangulation<dim>::face_iterator, double>;
