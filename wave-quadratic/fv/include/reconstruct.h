#pragma once

// dealii headers:
#include <deal.II/lac/householder.h>
#include <stencil.h>

namespace fv
{
  // auxillary functions:

  // return the no of faces sharing the boundary for a given cell-iterator:
  template <typename CellIter>
  inline unsigned int n_faces_on_boundary(const CellIter &cell)
  {
    unsigned int n_faces = 0;
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        n_faces++;
    return n_faces;
  }

  // CellMatrix:
  // Returns polynomial matrix for a given cell iterator and a stencil:
  // The functor is purely based on the geometry of the cell:
  // This works for both interior and boundary cells:
  // The CellMatrix class is independent of what polynomial is chosen so the
  // CellPolynomial is templated See Weno3 for the signature of CellPolynomial
  // type:

  template <int dim, int n_gpts, typename CellPolynomial>
  struct CellMatrix
  {
    // default constructor:
    CellMatrix() = default;

    // construct and return the matrix for each cell:
    template <typename CellIter, typename CellStencil>
    dealii::FullMatrix<double> operator()(const CellIter &cell,
                                          const CellStencil &stencil)
    {
      // basis function:
      CellPolynomial basis;
      basis.reinit(cell);

      // matix size:
      const unsigned int m = stencil.size() + n_faces_on_boundary(cell);
      const unsigned int n = basis.size();

      // initialize cell matrix
      dealii::FullMatrix<double> A(m, n);

      // loop over cells in the stencil, compute the cellaverage of the polynomial
      // basis and fill the matix rows:
      unsigned int i = 0; // row index
      for (const auto &cell_n : stencil)
      {
        // fill the rows by taking the average  of the polynomial basis on
        // neighbor cells:
        auto row = cell_average(basis, cell_n);
        // fill columns:
        for (unsigned int j = 0; j < n; ++j)
          A(i, j) = std::move(row[j]);
        i++;
      }

      // loop over faces, evaluate the polynomial basis at the boundary face
      // centre and fill the matrix:
      for (const auto &face : cell->face_iterators())
      {
        // check for boundary face:
        if (face->at_boundary())
        {
          // evaluate at face quadrature point
          auto row = basis(face->center());
          // fill columns:
          for (unsigned int j = 0; j < n; j++)
            A(i, j) = std::move(row[j]);
          i++;
        }
      }
      return A;
    }

  private:
    CellAverage<dim, n_gpts> cell_average;
  };

  // CellVector:
  template <typename CellIter, typename CellStencil, typename GlobalVector,
            typename BoundaryMap>
  inline dealii::Vector<double>
  cell_vector(const CellIter &cell, const CellStencil &stencil,
              const GlobalVector &u, BoundaryMap &u_bcs)
  {
    // initialize vector:
    const int size = stencil.size() + n_faces_on_boundary(cell);
    dealii::Vector<double> b(size);

    // loop over cells in the stencil or interior cells:
    unsigned int i = 0; // row index
    for (const auto &cell_n : stencil)
    {
      b(i) = u(cell_n->active_cell_index()) - u(cell->active_cell_index());
      i++;
    }

    // loop over boundary faces:
    for (const auto &face : cell->face_iterators())
    {
      // check for boundary face:
      if (face->at_boundary())
      {
        b(i) = u_bcs[face] - u(cell->active_cell_index());
        i++;
      }
    }

    return b;
  }

  // Reconstruct:
  template <int dim, int n_gpts, int n_layers, typename CellPolynomial,
            typename Stencil_ = Stencil<dim, n_layers>>
  struct Reconstruct
  {
    // constructor:
    Reconstruct() = default;

    // Create Polynomial matrix for all cells:
    void create_matrix(const dealii::Triangulation<dim> &tria)
    {
      // create stencil:
      stencil.reinit(tria);

      // build matrix looping over each cell:
      A.reserve(tria.n_active_cells());

      CellMatrix<dim, n_gpts, CellPolynomial> cell_matrix;
      for (const auto &cell : tria.active_cell_iterators())
        A.push_back(cell_matrix(cell, stencil(cell)));
    }

    // Polynomial Reconstruct operator:
    template <typename GlobalVector, typename BoundaryMap,
              typename VectorPolynomial>
    void operator()(const dealii::Triangulation<dim> &tria, const GlobalVector &u,
                    BoundaryMap &u_bcs, VectorPolynomial &upoly)
    {

      // solve Ax = b using Least-Square and compute polynomial coefficeints for
      // each cell
      dealii::Householder<double> householder;
      auto cell = tria.begin_active();
      auto endc = tria.end();
      for (unsigned int i = 0; cell != endc; ++cell, ++i)
      {

        householder.initialize(A[i]); // initialize A
        dealii::Vector<double> b = cell_vector(cell, stencil(cell), u, u_bcs);
        dealii::Vector<double> x(A[i].n());
        householder.least_squares(x, b);
        std::move(x.begin(), x.end(), upoly[i].coefficient().begin());
      }
    }

  private:
    Stencil_ stencil;
    std::vector<dealii::FullMatrix<double>> A;
  };
} // namespace fv
