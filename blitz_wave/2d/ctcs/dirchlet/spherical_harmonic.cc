/**
 * @file wave.cc
 * @author Magu
 * @brief  Spherical waves in 2d
           using 2nd and 4th order accurate laplacian.
 * @version 0.1
 * @date 2020-12-10
 * @copyright Copyright (c) 2020
 *
 */

/*Solve for u(x,t),x ∈ Ω([−1,1]^2),t ∈[0,T]

  ∂^2u/∂^2t−∇^2u= 0   x∈Ω

  Boundary condition:
  u(x,t) = 0      x∈∂Ω

  Initial condition:
  u(x,0) = ?  x∈Ω
  ∂u/∂t(x,0) = ?      x∈Ω

  using 2nd and 4th order accurate laplacian
*/

#include "tools.h"

//2d Wave stencil 2nd order stencil:
BZ_DECLARE_STENCIL4(wave22, P1, P2, P3, C)
P3(0, 0) = C * (P2(1, 0) + P2(-1, 0) + P2(0, 1) + P2(0, -1) - 4.0 * P2(0, 0)) + 2.0 * P2(0, 0) - P1(0, 0);
BZ_END_STENCIL

//2d Wave stencil 4th order stencil:
BZ_DECLARE_STENCIL4(wave24, P1, P2, P3, C)
P3(0, 0) = C * (-60.0 * P2(0, 0) + 16.0 * (P2(1, 0) + P2(-1, 0) + P2(0, 1) + P2(0, -1)) - (P2(2, 0) + P2(-2, 0) + P2(0, 2) + P2(0, -2))) + 2.0 * P2(0, 0) - P1(0, 0);
BZ_END_STENCIL


//initial pulse:


void wave(double tf = 5.0)
{
  auto print = [](const auto &x) { std::cout << x << std::endl; };

  //spherical wave solution:
  const double lambda(0.3);                            //wave length
  const double k(2.0 * M_PI / lambda);                 //wave number
  const double T(1.0), omega(2.0 * M_PI / T);          //time period and angular frequency:
  const double c = omega/k;                            // wave speed:
  double t;

  auto spherical_wave = [&t,omega,k,lambda](double x, double y){
    double r = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
    return r < lambda ? (1.0/r)*cos(omega*t - k*r) : 0.0;
  };
  

  // finite volume square grid:
  const int N = 300;
  const double L = 1.0;
  const double h = L / N;

  // cell coordinates:
  Array<double, 1> x(N), y(N);
  x = (0.5 + tensor::i) * h;
  y = x;


  // cfl and time step:
  double cfl = 0.5, dt = cfl * h / c;

  //Normalized wave speed:
  Array<double, 2> C(N, N);
  //C = (c*c*dt*dt)/(h*h);                  //2nd order scaling
  C = (c * c * dt * dt) / (12.0 * h * h);   //4th order scaling

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

  //initialize P3:
  Array<double, 2> P3(N, N);
  P3 = 0.0;

  while (t < tf)
  {
    nt++;
    t += dt;

    //update pressure:
    applyStencil(wave24(), P1, P2, P3, C);

	 //auto Pexact = project(spherical_wave, x, y);
    if (nt%10 == 0) write_tecplot(nt, t, x, y, P3);

    print(t);

    cycleArrays(P1, P2, P3);
  }
}

int main()
{
  wave();
  return 0;
}
