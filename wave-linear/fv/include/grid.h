#pragma once

// dealii headers:
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

namespace fv {

template <int dim>
inline void periodic_hypercube(dealii::Triangulation<dim> &triangulation,
                               const double left, const double right,
                               const int refine) {

  bool colorize = true; // Set boundary ids for the four boundaries
  dealii::GridGenerator::hyper_cube(triangulation, left, right, colorize);
  triangulation.refine_global(refine);
  dealii::GridTools::distort_random(0.15, triangulation, true);

  // Add periodicity to the mesh
  if (dim == 1 || dim == 2 || dim == 3) {

    // in x-direction
    std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::Triangulation<dim>::cell_iterator>>
        periodicity_vector_x;
    // Collect periodic faces
    dealii::GridTools::collect_periodic_faces(triangulation, 0, 1, 0,
                                              periodicity_vector_x);

    // add periodicity to the triangulation
    triangulation.add_periodicity(periodicity_vector_x);
  }

  if (dim == 2 || dim == 3) {

    // in y-direction

    std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::Triangulation<dim>::cell_iterator>>
        periodicity_vector_y;

    // Collect periodic faces

    dealii::GridTools::collect_periodic_faces(triangulation, 2, 3, 1,
                                              periodicity_vector_y);

    // add periodicity to the triangulation

    triangulation.add_periodicity(periodicity_vector_y);
  }

  if (dim == 3) {

    // in y-direction

    std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::Triangulation<dim>::cell_iterator>>
        periodicity_vector_z;

    // Collect periodic faces

    dealii::GridTools::collect_periodic_faces(triangulation, 4, 5, 2,
                                              periodicity_vector_z);

    // add periodicity to the triangulation

    triangulation.add_periodicity(periodicity_vector_z);
  }
}
} // namespace fv
