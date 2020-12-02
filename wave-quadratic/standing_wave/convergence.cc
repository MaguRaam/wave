/*Standing wave in cavity:

  Solve for u(x,t),x ∈ Ω([−1,1]^2),t ∈[0,T]

  ∂u/∂t = v                 x in Ω
  ∂v/∂t = ∇^2u              x in Ω
  
  Boundary condition:
  u(x,t) = 0             x in ∂Ω
  
  Initial condition:
  u(x,0) = sin(πkx)sin(πly)  x in Ω
  v(x,0) = 0             x in Ω

*/
#include "fv.h"

using namespace dealii;
using namespace fv;
using namespace std;
namespace plt = matplotlibcpp;

template <int dim>
using BoundaryMap = map<typename Triangulation<dim>::face_iterator, double>;

//wave:
template <int dim, int n_gpts, int n_layers, typename Polynomial>
class wave
{

  const Triangulation<dim> &tria;
  const DoFHandler<dim> &dof_handler;

  Reconstruct<dim, n_gpts, n_layers, Polynomial> reconstruct;
  CellLaplace<dim, n_gpts> cell_laplace;

  const int n_cells;
  Vector<double> u;
  vector<Polynomial> polynomial;
  BoundaryMap<dim> u_bcs;

public:
  //constructor:
  wave(const Triangulation<dim> &tria,
       const DoFHandler<dim> &dof_handler) : tria(tria),
                                             dof_handler(dof_handler),
                                             n_cells(dof_handler.n_dofs()),
                                             u(n_cells),
                                             polynomial(n_cells)
  {
    // set dirchlet bcs u = 0:
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->at_boundary())
      {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            u_bcs[cell->face(f)] = 0.0;
      }
    }

    // construct polynomial for each cell:
    auto cell = tria.begin_active();
    auto endc = tria.end();
    for (unsigned int c = 0; cell != endc; ++cell, ++c)
      polynomial[c].reinit(cell);

    //Initialize matrix in Reconstructor operator:
    reconstruct.create_matrix(tria);
  }

  //operator():
  void operator()(const Vector<double> &x, Vector<double> &dxdt, double /*t*/)
  {
    // extract u:
    copy(x.begin(), x.begin() + n_cells, u.begin());

    // reconstruct polynomial for u:
    reconstruct(tria, u, u_bcs, polynomial);

    // du/dt = v:
    copy(x.begin() + n_cells, x.end(), dxdt.begin());

    // dv/dt = L(u):
    auto cell = dof_handler.begin_active();
    auto endc = dof_handler.end();
    for (unsigned int i = n_cells; cell != endc; ++cell, ++i)
      dxdt[i] = cell_laplace(cell, polynomial);
  }
};

//write data:
template <int dim>
void write(const DoFHandler<dim> &dof_handler, const Vector<double> &u, int n)
{

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(u, "u", DataOut<dim>::type_dof_data);
  data_out.build_patches();
  const std::string filename = "../standing_wave/plot/solution-" +
                               Utilities::int_to_string(n, 3) +
                               ".vtk";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level =
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtk(output);
}


auto run(int refinement)
{
  const int dim = 2, n_gpts = 3, n_layers = 1;

  using Polynomial = polynomial::Quadratic<dim>;
  //Generate grid:
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1., 1., false);
  tria.refine_global(refinement);
  GridTools::distort_random(0, tria, true);

  //Index cells:
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_DGQ<dim>(0));
  const int n_cells = dof_handler.n_dofs();

  cout << "grid created with " << n_cells << " cells\n";

  //Initial condition:
  const int k = 2;
  Vector<double> u(n_cells);
  project<dim, n_gpts>(dof_handler.begin_active(), dof_handler.end(), u.begin(),[k](const auto &x) { return sin(M_PI * k * x(0)) * sin(M_PI * k * x(1)); });

  std::vector<double> v(n_cells, 0.0);

  //Initialize global state x by concatenating u and v:
  Vector<double> x(2 * n_cells);

  //x = [u0,u1,...,n-1,v0,v1,...,vn-1]
  move(u.begin(), u.end(), x.begin());
  move(v.begin(), v.end(), x.begin() + n_cells);

  //time-step:
  const int steps = 500;
  const double dt = 0.001;

  //solve:

  //Wave
  wave<dim, n_gpts, n_layers, Polynomial> Wave(tria, dof_handler);

  //Ode
  euler<Vector<double>> Euler;

  size_t n = 0;
  for (; n < steps; ++n)
  {

    Euler.do_step(Wave, x, n * dt, dt);
    std::cout << n * dt << '\n';

    //write data:
    if (n % 10 == 0)
    {
      //extract u:
      copy(x.begin(), x.begin() + n_cells, u.begin());
      write(dof_handler, u, n);
    }
  }

  //extract u from x for last time instant:
  copy(x.begin(), x.begin() + n_cells, u.begin());

  //exact solution at last time instant
  double omega = sqrt(2*M_PI*M_PI*k*k);
  double t = n * dt;
  auto u_exact = [k, t, omega](auto p) { return (cos(omega * t)) * sin(M_PI * k * p(0)) * sin(M_PI * k * p(1)); };

  //L2 and Linfty error at last time instant:
  auto L2     = l2_error<dim,n_gpts>(u, u_exact, dof_handler);
  auto Linfty = linfty_error<dim,n_gpts>(u, u_exact, dof_handler); 

  return make_tuple(sqrt(n_cells),L2,Linfty);
}

int main()
{
  vector<double> N, L2, Linfty;

  ofstream data("../standing_wave/convergence.dat", ios::out);
  data.flags(ios::dec | ios::scientific);
  data.precision(5);

  //write covergence:
  for (int refine = 5; refine < 9; ++refine)
  {
    auto output = run(refine);
    N.push_back(get<0>(output));
    L2.push_back(get<1>(output));
    Linfty.push_back(get<2>(output));

    data << "N = " << get<0>(output) << "  "
         << "Linfty error = " << get<2>(output)
         << "  L2 error = " << get<1>(output) << "\n";
  }

  auto rate = [](auto e1, auto e2) { return (log(e1) - log(e2)) / log(2); };
  data << "Rate of convergence\n\n";
  for (unsigned int i = 1; i < N.size(); ++i)   
    data << "N = " << N[i] << "  " << "Linfty rate = " << rate(Linfty[i - 1], Linfty[i]) << "  L2 rate = " << rate(L2[i - 1], L2[i]) << "\n";
  
  data.close();

  //Convergence plot:
  plt::loglog(N, L2, {{"label", "l2"}, {"marker", "o"}});
  plt::loglog(N, Linfty, {{"label", "linfty_error"}, {"marker", "o"}});
  plt::legend();
  plt::xlabel("N");
  plt::ylabel("Error");
  plt::title("Convergence");
  plt::savefig("../standing_wave/convergence.png");

  return 0;
}


 