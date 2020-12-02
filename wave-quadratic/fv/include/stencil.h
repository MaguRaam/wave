#pragma once

// dealii headers:
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// stl:
#include <vector>

#include <cell_neighbors.h>

// Sundar
#include <vertex_neighbors.h>

namespace fv
{

  // Stencil:

  template <int dim, unsigned int n_layers>
  class Stencil
  {
    using CellTriaItr = typename dealii::Triangulation<dim>::active_cell_iterator;

  public:
    Stencil() = default;

    // create stencil:
    void reinit(const dealii::Triangulation<dim> &tria)
    {
      // resize stencil to ncells:
      stencil_.resize(tria.n_active_cells());

      // initialize cell neighbors functor:
      CellNeighbors<dim, n_layers> cell_neighbors(tria);

      // Loop over each cell in the stencil and add cell neighbors:
      auto cell = tria.begin_active();
      auto endc = tria.end();
      for (unsigned int c = 0; cell != endc; ++cell, ++c)
        stencil_[c] = cell_neighbors(cell);
    }

    // acess stencil: given tria cell iterator, the function returns the vector of
    // neighbor tria-cell iterators:
    const std::vector<CellTriaItr> &operator()(const CellTriaItr &cell) const
    {
      return stencil_[cell->active_cell_index()];
    }

  private:
    std::vector<std::vector<CellTriaItr>> stencil_;
  };

} // namespace fv
