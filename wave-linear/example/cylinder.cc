#pragma once
 

#include "fv_tools.h"
#include "ode.h"

template <int dim> class wave {
public:
  wave();
  void make_grid();
  void setup_system();
  void operator()(const Vector<double> &x, Vector<double> &dxdt, double /*t*/);
  void write_data(int timestep_number);

  unsigned int n_cells;
private:
  Triangulation<dim> triangulation;
  DoFHandler<dim> dof_handler;
  Reconstruct<dim> reconstruct;
  CellLaplace<dim> laplace;

  Vector<double> u;
  GlobalPolynomial<dim> u_poly;
  BoundaryCondition<dim> u_bcs;
};

// constructor:
template <int dim> wave<dim>::wave() : dof_handler(triangulation) {}

// create cylinder grid:
template <int dim> void wave<dim>::make_grid() {
  GridGenerator::channel_with_cylinder(triangulation);
  triangulation.refine_global(4);

  //set boundary id = 1 for inlet and id = 0 for others:
  for (auto &face : triangulation.active_face_iterators())
    if (face->at_boundary() && face->center()[0] == 0.0 && face->center()[1] >=0.1 && face->center()[1] <=0.3)
      face->set_boundary_id(1);

  // write mesh file:
  std::ofstream out("../plot/grid.vtk");
  GridOut grid_out;
  grid_out.write_vtk(triangulation, out);
}

template <int dim> void wave<dim>::setup_system() {
  // enumerate finite volume cells:
  dof_handler.distribute_dofs(FE_DGQ<dim>{0});
  n_cells = dof_handler.n_dofs();

  // create polynomial matrix for each cell:
  reconstruct.create_matrix(dof_handler);

  // resize solution and polynomial:
  u.reinit(n_cells);
  u_poly.resize(n_cells);
}

//dx/dt = f(x) <=> du/dt = v , dv/dt = L(u)
template <int dim>
void wave<dim>::operator()(const Vector<double> &x, Vector<double> &dxdt,
                           double t) {
  // du/dt = v:
  using std::copy;
  copy(x.begin() + n_cells, x.end(), dxdt.begin());

  // extract u:
  copy(x.begin(), x.begin() + n_cells, u.begin());

  //set bcs:
  for (auto &face : triangulation.active_face_iterators()) {
    if (face->at_boundary())
      (face->boundary_id() == 1) ? (u_bcs[face] = sin(4 * M_PI * t)) : 0.0;
  }

  // Reconstruct polynomial for u:
  reconstruct(dof_handler, u, u_bcs, u_poly);

  // dv/dt = L(u):
  unsigned int i = n_cells;
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    dxdt[i] = laplace(cell, u_poly);
    i++;
  }
}

// write data:
template <int dim> void wave<dim>::write_data(int timestep_number) {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(u, "U", DataOut<dim>::type_dof_data);
  data_out.build_patches();

  const std::string filename = "../plot/solution-" +
                               Utilities::int_to_string(timestep_number, 3) +
                               ".vtk";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level =
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtk(output);
}


using RK4 = runge_kutta4<Vector<double>>;

int main() {
  wave<2> wave_;
  wave_.make_grid();
  wave_.setup_system();

  //integrate in time:
  RK4 rk4;
  const int steps = 5000;
  const double dt = 0.0001;

  //ics:
  Vector<double> x(2 * wave_.n_cells);
  x = 0.0;

  for (size_t n = 0; n < steps; ++n)
	  {
	      rk4.do_step(wave_, x, n * dt, dt);
	      std::cout << n * dt << '\n';
	      if (n%100 == 0) wave_.write_data(n);
	  }

  return 0;
}
