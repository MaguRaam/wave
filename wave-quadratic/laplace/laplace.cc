#include "fv.h"

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iomanip>

// Compute Laplacian of the given function:
// u(x,0) =sin(πkx)sin(πly)  and x ∈ Ω([−1,1]2)

//In this program we do the convergence study for the Laplacian opeartor
//varying different polynomials and dimension:

using namespace dealii;
using namespace fv;
using namespace std;

template <int dim>
auto laplace(unsigned int refine)
{
  using Polynomial = polynomial::Quadratic<dim>;

  // generate grid:
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1., 1., false);
  tria.refine_global(refine);
  GridTools::distort_random(0.1, tria, true);

  // enumerate cells:
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_DGQ<dim>{0});
  unsigned int n_cells = dof_handler.n_dofs();

  // initial condition u:
  int k(2);
  auto uf = [k](auto p) { return sin(M_PI * k * p(0)) * sin(M_PI * k * p(1)); };

  // project u on grid:
  Vector<double> u(n_cells);
  project<dim, 2>(tria.begin_active(), tria.end(), u.begin(), uf);

  // set dirchlet bcs u = 0;
  using BoundaryMap = map<typename Triangulation<dim>::face_iterator, double>;
  BoundaryMap u_bcs;
  for (auto &face : tria.active_face_iterators())
    if (face->at_boundary())
      u_bcs[face] = 0.0;

  //create polynomial vector:

  vector<Polynomial> polynomial(n_cells);

  auto cell = tria.begin_active();
  auto endc = tria.end();
  for (unsigned int c = 0; cell != endc; ++cell, ++c)
    polynomial[c].reinit(cell);

  //reconstruct Quadratic polynomial:
  Reconstruct<dim, 2, 1, Polynomial> reconstruct;
  reconstruct.create_matrix(tria);
  reconstruct(tria, u, u_bcs, polynomial);

  // Compute Laplacian:
  Vector<double> laplace(n_cells);
  CellLaplace<dim, 2> cell_laplace;
  unsigned int i = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    laplace[i] = cell_laplace(cell, polynomial);
    i++;
  }

  // exact laplacian:
  auto laplacef = [k](auto p) { return (-2.0 * M_PI * M_PI * k * k) * sin(M_PI * k * p(0)) * sin(M_PI * k * p(1)); };
  Vector<double> laplace_exact(n_cells);
  project<dim, 2>(tria.begin_active(), tria.end(), laplace_exact.begin(), laplacef);

  //L2 and Linfty error:
  auto L2 = l2_error<dim, 2>(laplace, laplacef, dof_handler);
  auto Linfty = linfty_error<dim, 2>(laplace, laplacef, dof_handler);

  //compute error for each cell:
  Vector<double> error(n_cells);
  for (unsigned int i = 0; i < n_cells; ++i)
    error[i] = abs(laplace[i] - laplace_exact[i]);

  // plot u:
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(u, "u", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(laplace, "laplace2", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(laplace_exact, "laplace_exact", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(error, "error", DataOut<dim>::type_dof_data);
  data_out.build_patches();
  ofstream output(dim == 2 ? "../laplace/2d.vtk" : "../laplace/3d.vtk");
  data_out.write_vtk(output);

  return make_tuple(sqrt(n_cells), L2, Linfty);
}

namespace plt = matplotlibcpp;
int main()
{
  vector<double> N, L2, Linfty;

  ofstream File("../laplace/convergence.dat", ios::out);
  File.flags(ios::dec | ios::scientific);
  File.precision(5);

  //convergence table:
  for (int refine = 5; refine < 9; ++refine)
  {
    //compute laplacian:
    auto output = laplace<2>(refine);
    N.push_back(get<0>(output));
    L2.push_back(get<1>(output));
    Linfty.push_back(get<2>(output));

    File << "N = " << get<0>(output) << "  "
         << "Linfty error = " << get<2>(output)
         << "  L2 error = " << get<1>(output) << "\n";
  }
  auto rate = [](auto e1, auto e2) { return (log(e1) - log(e2)) / log(2); };
  File << "Rate of convergence\n\n";

  for (unsigned int i = 1; i < N.size(); ++i)
    File << "N = " << N[i] << "  "
         << "Linfty rate = " << rate(Linfty[i - 1], Linfty[i]) << "  L2 rate = " << rate(L2[i - 1], L2[i]) << "\n";
  File.close();

  //Convergence plot:
  plt::loglog(N, L2, {{"label", "l2"}, {"marker", "o"}});
  plt::loglog(N, Linfty, {{"label", "linfty_error"}, {"marker", "o"}});
  plt::legend();
  plt::xlabel("N");
  plt::ylabel("Error");
  plt::title("Convergence");
  plt::savefig("../laplace/convergence.png");
  plt::show();

  return 0;
}
