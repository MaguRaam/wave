#pragma once
// author : Sunder
// Ref. : Applied CFD Techniques: An Introduction based on Finite Element
// Methods, Rainald Lohner (Wiley 2008)

// Class for getting all the cells surrounding a vertex

// dealii headers:
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

namespace fv
{
  // Forward declaration
  template <unsigned int>
  class Vert2Cell;

  // Class function definitions

  template <unsigned int dim>
  class Vert2Cell
  {
    unsigned int n_vertices;
    unsigned int n_cells;
    int *vert_to_cell1; // These two are linked lists
    int *vert_to_cell2;

  public:
    explicit Vert2Cell();
    Vert2Cell(const dealii::Triangulation<dim> &);
    ~Vert2Cell();

    void initialzie(const dealii::Triangulation<dim> &);

    unsigned int n_cells_sharing_vertex(unsigned int) const;
    unsigned int cell_sharing_vertex(unsigned int, unsigned int) const;
  };

  // Default constructor

  template <unsigned int dim>
  inline Vert2Cell<dim>::Vert2Cell()
  {
    n_cells = 0;
    n_vertices = 0;
    vert_to_cell1 = new int[0];
    vert_to_cell2 = new int[0];
  }

  // Constructor taking the triangulation object as input

  template <unsigned int dim>
  Vert2Cell<dim>::Vert2Cell(const dealii::Triangulation<dim> &tria)
  {

    int ipoi1, ipoin, istor;

    n_cells = tria.n_active_cells();
    n_vertices = tria.n_vertices();

    unsigned int vertices_per_cell = dealii::GeometryInfo<dim>::vertices_per_cell;

    vert_to_cell1 = new int[n_cells * vertices_per_cell];
    vert_to_cell2 = new int[n_vertices + 1];

    // First create a cell to vertex connectivity

    int **cell_to_vert;

    cell_to_vert = new int *[n_cells];

    for (unsigned int c = 0; c < n_cells; ++c)
    {
      cell_to_vert[c] = new int[vertices_per_cell];
    }

    // Get iterators for the cells and fill the cell to vertex connectivity

    auto cell = tria.begin_active();
    auto endc = tria.end();

    for (unsigned int c = 0; cell != endc; ++cell, ++c)
    {
      for (unsigned int v = 0; v < vertices_per_cell; ++v)
      {
        cell_to_vert[c][v] = cell->vertex_index(v);
      }
    }

    // Cell pass 1: Count the number of cells connected to each vertex

    // Initialize

    for (unsigned int v = 0; v < n_vertices + 1; ++v)
    {
      vert_to_cell2[v] = 0;
    }

    for (unsigned int c = 0; c < n_cells; ++c)
    {
      for (unsigned int v = 0; v < vertices_per_cell; ++v)
      {

        ipoi1 = cell_to_vert[c][v] + 1;
        vert_to_cell2[ipoi1] = vert_to_cell2[ipoi1] + 1;
      }
    }

    // Storage/reshuffling pass 1

    for (unsigned int v = 1; v < n_vertices + 1; ++v)
    {
      // Loop over all vertices
      // Update storage counter and store
      vert_to_cell2[v] = vert_to_cell2[v] + vert_to_cell2[v - 1];
    }

    // Cell pass 2: Store the cells in vert_to_cell1

    for (unsigned int c = 0; c < n_cells; ++c)
    {
      for (unsigned int v = 0; v < vertices_per_cell; ++v)
      {
        // Update storage counter, storing in vert_to_cell1
        ipoin = cell_to_vert[c][v];
        istor = vert_to_cell2[ipoin] + 1;
        vert_to_cell2[ipoin] = istor;
        vert_to_cell1[istor - 1] = c;
      }
    }

    // Storage/reshuffling pass 2:
    for (unsigned int v = n_vertices; v >= 1; --v)
    {
      vert_to_cell2[v] = vert_to_cell2[v - 1];
    }

    vert_to_cell2[0] = 0;

    // Subtract 1 from all

    for (unsigned int v = 0; v < n_vertices + 1; ++v)
    {
      vert_to_cell2[v] -= 1;
    }

    // Delete the memory

    for (unsigned int c = 0; c < n_cells; ++c)
    {
      delete[] cell_to_vert[c];
    }
  }

  // Destructor

  template <unsigned int dim>
  Vert2Cell<dim>::~Vert2Cell()
  {
    delete[] vert_to_cell1;
    delete[] vert_to_cell2;
  }

  // Get the number of cells surrounding a vertex

  template <unsigned int dim>
  inline unsigned int
  Vert2Cell<dim>::n_cells_sharing_vertex(unsigned int v) const
  {

    Assert(v < n_vertices, ExcIndexRange(v, 0, n_vertices));

    return (vert_to_cell2[v + 1] - (vert_to_cell2[v] + 1)) + 1;
  }

  // Get the global cell index of the ith cell which shares the vertex v

  template <unsigned int dim>
  inline unsigned int Vert2Cell<dim>::cell_sharing_vertex(unsigned int v,
                                                          unsigned int i) const
  {

    // Assert(v < n_vertices, ExcIndexRange(v,0,n_vertices)); // This assertion is
    // not required as it is already being checked in the next function

    unsigned int total_no_of_cells_sharing_vertex = n_cells_sharing_vertex(v);

    Assert(i < total_no_of_cells_sharing_vertex,
           ExcIndexRange(i, 0, total_no_of_cells_sharing_vertex));

    return vert_to_cell1[vert_to_cell2[v] + 1 + i];
  }

} // namespace fv
