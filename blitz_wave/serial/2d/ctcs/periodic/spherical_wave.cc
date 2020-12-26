/**
 * @file wave_periodic.cc
 * @author Magu
 * @brief Solves the wave equation in 2d for a periodic boundary condition.
 * @version 0.1
 * @date 2020-12-10
 * @copyright Copyright (c) 2020
 * 
 */

/*Solve for u(x,t),x ∈ Ω([−1,1]^2),t ∈[0,T]

  ∂^2u/∂^2t−c^2∇^2u= 0   x∈Ω

  Periodic Boundary condition:

  Exact solution
  u(x,0) =sin(kx + ky - omegat)  x∈Ω
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <blitz/array.h>

using namespace blitz;

//project functions on grid:
template <typename F>
inline Array<double, 2> project(F f, const Array<double, 1> &x, const Array<double, 1> &y)
{
  int nx(x.numElements()), ny(y.numElements());
  Array<double, 2> u(nx, ny);

  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
      u(i, j) = f(x(i), y(j));
  }

  return u;
}

std::string int_to_string(unsigned int value, const unsigned int digits)
{
  std::string lc_string = std::to_string(value);

  if (lc_string.size() < digits)
  {
    // We have to add the padding zeroes in front of the number
    const unsigned int padding_position = (lc_string[0] == '-')
                                              ? 1
                                              : 0;

    const std::string padding(digits - lc_string.size(), '0');
    lc_string.insert(padding_position, padding);
  }

  return lc_string;
}
//write tecplot
void write_tecplot(int n, double t, const Array<double, 1> &x, const Array<double, 1> &y, const Array<double, 2> &u)
{
  std::ofstream tpl;
  const std::string filename = "../plot/plot_" + int_to_string(n, 3) + ".dat";
  tpl.open(filename);
  tpl.flags(std::ios::dec | std::ios::scientific);
  tpl.precision(6);

  int nx(x.numElements()), ny(y.numElements());
  tpl << "TITLE = \"Wave Equation 2D\" " << std::endl
      << "VARIABLES = \"x\", \"y\", \"u\" " << std::endl;
  tpl << "Zone I = " << ny - 2 << " J = " << nx - 2 << std::endl;
  tpl << "SOLUTIONTIME = " << t << std::endl;

  for (int i = 1; i < nx - 1; i++)
    for (int j = 1; j < ny - 1; j++)
      tpl << x(i) << "\t" << y(j) << "\t" << u(i, j) << std::endl;
}

//2d Wave stencil 2nd order stencil:
BZ_DECLARE_STENCIL4(wave22, P1, P2, P3, CFL)
P3(0, 0) = CFL * CFL * (P2(1, 0) + P2(-1, 0) + P2(0, 1) + P2(0, -1) - 4.0 * P2(0, 0)) + 2.0 * P2(0, 0) - P1(0, 0);
BZ_END_STENCIL

void wave(double tf = 1.0)
{
  auto print = [](const auto &x) { std::cout << x << std::endl; };

  //plane wave parameters:
  const double lambda(0.3);                   //wave length
  const double k(2.0 * M_PI / lambda);        //wave number
  const double T(1.0), omega(2.0 * M_PI / T); //time period and angular frequency:

  //dispersion relation to compute wave speed:
  const double c(omega / k);

  //spherical wave solution:
  double t;
  auto spherical_wave = [&t, omega, k, lambda](double x, double y) {
    double r = sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
    return r < lambda ? (1.0 / r) * cos(omega * t - k * r) : 0.0;
  };

  //finite volume square grid:
  const int N = 64;
  const int Ng = 1;
  const double L = 1.0;
  const double h = L / N;

  //cell coordinates:
  Array<double, 1> x(N + 2 * Ng), y(N + 2 * Ng);
  x(Range(1, N)) = (0.5 + tensor::i) * h;
  x(0) = x(N);
  x(N + 1) = x(1);
  y = x;

  Array<double, 2> CFL(N + 2 * Ng, N + 2 * Ng);

  //time step:
  double cfl = 0.5;
  double dt = cfl * h / c;
  CFL = cfl;

  //initial condition P1:
  int nt = 0;
  t = 0;
  auto P1 = project(spherical_wave, x, y);
  write_tecplot(nt, t, x, y, P1);

  //initial condition P2:
  nt++;
  t += dt;
  auto P2 = project(spherical_wave, x, y);
  write_tecplot(nt, t, x, y, P2);

  //initialize P3;
  Array<double, 2> P3(N + 2 * Ng, N + 2 * Ng);
  P3 = 0.0;

  while (t < tf)
  {
    nt++;
    t += dt;

    //update pressure:
    applyStencil(wave22(), P1, P2, P3, CFL);

    if (nt % 10 == 0)
    {
      write_tecplot(nt, t, x, y, P3);
    }

    //periodic bcs:
    Range all = Range::all();

    //make periodic along j dxn:
    P3(all, 0) = P3(all, N);
    P3(all, N + 1) = P3(all, 1);

    //make periodic along i dxn:
    P3(0, all) = P3(N, all);
    P3(N + 1, all) = P3(1, all);

    print(t);

    cycleArrays(P1, P2, P3);
  }
}

int main()
{
  wave();
  return 0;
}
