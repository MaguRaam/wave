

#include "fv_tools.h"


// Compute Laplacian of the given function:
// u(x,0) =sin(πkx)sin(πly)  and x ∈ Ω([−1,1]2)

template <int dim> void laplace(unsigned int refine) {

  // generate grid:
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1., 1., false);
  triangulation.refine_global(refine);
  //GridTools::distort_random(0.1, triangulation, true);
  std::ofstream out("../plot/grid.vtk");
  GridOut grid_out;
  grid_out.write_vtk(triangulation, out);

  // enumerate cells:
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(FE_DGQ<dim>{0});
  unsigned int n_cells = dof_handler.n_dofs();

  // initialize u:
  Vector<double> u(n_cells);
  int k(2);
  auto uf = [k](auto p) { return sin(M_PI * k * p(0)) * sin(M_PI * k * p(1)); };

  CellAverage<dim, 4> cell_average;
  unsigned int i(0);
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    u[i] = cell_average(uf, cell);
    i++;
  }

  // set dirchlet bcs u = 0;
  BoundaryCondition<dim> u_bcs;
  for (auto &face : triangulation.active_face_iterators())
    if (face->at_boundary())
      u_bcs[face] = 0.0;


  // Construct polynomial:
  GlobalPolynomial<dim> u_poly(n_cells);
  Reconstruct<dim> reconstruct;
  reconstruct.create_matrix(dof_handler);
  reconstruct(dof_handler, u, u_bcs, u_poly);

  // Compute Laplacian:
  Vector<double> u_laplace(n_cells);
  CellLaplace<dim> cell_laplace;
  i = 0;
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    u_laplace[i] = cell_laplace(cell, u_poly);
    i++;
  }

  // error:
  auto laplace_function = [k](auto p) {
    return (-2.0 * M_PI * M_PI * k * k) * sin(M_PI * k * p(0)) *
           sin(M_PI * k * p(1));
  };

  Vector<double> uf_laplace(n_cells);
  i = 0;
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    uf_laplace[i] = cell_average(laplace_function, cell);
    i++;
  }

  //compute error for each cell:
  Vector<double> error(n_cells);
  i = 0;
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    error[i] = abs(u_laplace[i] - uf_laplace[i]);
    i++;
  }


  std::cout << "ncells = " << n_cells << "  "
            << "Linfty error = "
            << linfty_error<dim, 4>(u_laplace, laplace_function, dof_handler)
            << "  L2 error = "
            <<l2_error<dim, 4>(u_laplace, laplace_function, dof_handler)<<"\n";

  // plot u:
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(u, "u", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(u_laplace, "u_laplace", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(uf_laplace, "uf_laplace",
                           DataOut<dim>::type_dof_data);
  data_out.add_data_vector(error, "error",
                           DataOut<dim>::type_dof_data);
  data_out.build_patches();
  std::ofstream output(dim == 2 ? "../plot/2d.vtk" : "../plot/3d.vtk");
  data_out.write_vtk(output);
}

int main() {

  for (unsigned int refine = 4; refine < 9; refine++)
    laplace<2>(refine);

  return 0;
}

// for 3d replace:
/*auto uf = [k](auto p) {
  return sin(M_PI * k * p(0)) * sin(M_PI * k * p(1)) * sin(M_PI * k * p(2));};*/

// auto laplace_function = [k](auto p){return (-3.0*M_PI*M_PI*k*k)*sin(M_PI * k
// * p(0)) * sin(M_PI * k * p(1)) * sin(M_PI * k * p(2));};
