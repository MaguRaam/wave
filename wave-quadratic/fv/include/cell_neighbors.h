#pragma once

// dealii headers
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// stl
#include <unordered_set>

// Sundar
#include <vertex_neighbors.h>

namespace fv
{
  // CellNeighbors

  template <int dim, unsigned int n_layers>
  class CellNeighbors
  {
    using CellTriaItr = typename dealii::Triangulation<dim>::active_cell_iterator;

  public:
    explicit CellNeighbors(const dealii::Triangulation<dim> &tria) : V2C(tria)
    {

      // map from global index to cell iter:
      auto cell = tria.begin_active();
      auto endc = tria.end();

      cell_iterator.resize(tria.n_active_cells());
      for (unsigned int c = 0; cell != endc; ++cell, ++c)
        cell_iterator[c] = cell;
    }

    // operator()
    // find n_layers of surrounding neighbor cells for a given owner cell index:
    // returns neighbor cell tria iterators:

    std::vector<CellTriaItr> operator()(const CellTriaItr &cell) const
    {

      // initialize old cell-set with owner cell index:
      std::unordered_set<unsigned int> old_cell_set;
      old_cell_set.insert(cell->active_cell_index());

      // loop over cell layers:
      for (unsigned int l = 0; l < n_layers; l++)
      {

        // create new set of cells from old_cell_set:
        auto new_cell_set(old_cell_set);

        // loop over cells in old cell set:
        for (const auto &c : old_cell_set)
        {

          // get neighbor cells that share the vertex with the current cell;
          auto neighbor_cells = closest_cell_neighbors(c);

          // insert these cells in our new cell_set:
          new_cell_set.insert(neighbor_cells.begin(), neighbor_cells.end());
        }

        old_cell_set = std::move(new_cell_set);
      }

      // remove owner cell index since, we are interested only in neighbors:
      old_cell_set.erase(cell->active_cell_index());

      // return vector of cell iterators of neighbor cells:
      std::vector<CellTriaItr> cell_iterators;
      for (const auto &index : old_cell_set)
        cell_iterators.push_back(cell_iterator[index]);

      return cell_iterators;
    }

  private:
    // given cell index the function returns the cells that shares the vertex with
    // the owner cell:
    std::unordered_set<unsigned int>
    closest_cell_neighbors(const unsigned int cell_index) const
    {

      // initialize the set of neighbor cells:
      std::unordered_set<unsigned int> cells;

      // cell iterator from cell index:
      auto cell_iter = cell_iterator[cell_index];

      // loop over cell vertices:
      for (unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell;
           ++v)
      {

        // global index of vertex v:
        auto vertex = cell_iter->vertex_index(v);

        // no of cells sharing vertex v:
        auto n_cells_sharing_vertex = V2C.n_cells_sharing_vertex(vertex);

        // loop over cells sharing vertex v and collect it in cells set:
        for (unsigned int c = 0; c < n_cells_sharing_vertex; ++c)
          cells.insert(V2C.cell_sharing_vertex(vertex, c));
      }

      return cells;
    }

  private:
    Vert2Cell<dim> V2C;
    std::vector<CellTriaItr> cell_iterator;
  };

} // namespace fv
