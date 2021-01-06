

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <blitz/array.h>

using namespace blitz;


// project functions on grid:
template <typename F>
inline Array<double, 2> project(F f, const Array<double, 1> &x,
                                const Array<double, 1> &y) {
  int nx(x.numElements()), ny(y.numElements());
  Array<double, 2> u(nx, ny);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j)
      u(i, j) = f(x(i), y(j));
  }

  return u;
}

std::string int_to_string(unsigned int value, const unsigned int digits) {
  std::string lc_string = std::to_string(value);

  if (lc_string.size() < digits) {
    // We have to add the padding zeroes in front of the number
    const unsigned int padding_position = (lc_string[0] == '-') ? 1 : 0;

    const std::string padding(digits - lc_string.size(), '0');
    lc_string.insert(padding_position, padding);
  }

  return lc_string;
}
// write tecplot
void write_tecplot(int nt, double t, const Array<double, 1> &x,
                   const Array<double, 1> &y, const Array<double, 2> &u) {
  std::ofstream tpl;
  const std::string filename = "../plot/plot_" + int_to_string(nt, 3) + ".dat";
  tpl.open(filename);
  tpl.flags(std::ios::dec | std::ios::scientific);
  tpl.precision(6);

  int nx(x.numElements()), ny(y.numElements());
  tpl << "TITLE = \"Wave Equation 2D\" " << std::endl
      << "VARIABLES = \"x\", \"y\", \"u\" " << std::endl;
  tpl << "Zone I = " << ny << " J = " << nx << std::endl;
  tpl << "SOLUTIONTIME = " << t << std::endl;

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      tpl << x(i) << "\t" << y(j) << "\t" << u(i, j) << std::endl;
}