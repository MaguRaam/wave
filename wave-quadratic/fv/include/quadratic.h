#pragma once

namespace fv
{
  namespace polynomial
  {
    template <int dim>
    struct Quadratic
    {
    };

    // template specialization for dim = 2;
    template <>
    struct Quadratic<2>
    {
      using CellItr = typename dealii::Triangulation<2>::active_cell_iterator;

      // constructor
      Quadratic() = default;

      void reinit(const CellItr &cell)
      {

        // cell center:
        const auto po = cell->barycenter();

        // compute h:
        h = sqrt(cell->measure());

        // extract coordinates of the cell centre:
        xo = po[0];
        yo = po[1];

        // compute constants c1, c2 and c3:
        CellAverage<2, 4> cell_average;

        c1 = cell_average([this](auto p) {double x = p[0]; return (1.0 / (h * h)) * (x - xo) * (x - xo); }, cell);
        c2 = cell_average([this](auto p) {double y = p[1];return (1.0 / (h * h)) * (y - yo) * (y - yo); }, cell);
        c3 = cell_average([this](auto p) {double x = p[0], y = p[1];return (1.0 / (h * h)) * (x - xo) * (y - yo); }, cell);
      }

      // operator(): returns the basis functions of the polynomial evaluated at
      // point p:
      dealii::Vector<double> operator()(const dealii::Point<2, double> &p) const
      {
        const double x = p[0], y = p[1];

        dealii::Vector<double> basis(5);
        basis[0] = (1.0 / h) * (x - xo);
        basis[1] = (1.0 / h) * (y - yo);
        basis[2] = ((1.0 / (h * h)) * (x - xo) * (x - xo)) - c1;
        basis[3] = ((1.0 / (h * h)) * (y - yo) * (y - yo)) - c2;
        basis[4] = ((1.0 / (h * h)) * (x - xo) * (y - yo)) - c3;

        return basis;
      }

      // return gradient of the polynomial evaluated at cell center:
      friend dealii::Vector<double> gradient(const Quadratic<2> &quadratic,
                                             const dealii::Point<2, double> &p)
      {
        const double h = quadratic.h, xo = quadratic.xo, yo = quadratic.yo;
        const double x = p[0], y = p[1];

        // x and y component of gradient
        const double ux = quadratic.coefficients[0], uy = quadratic.coefficients[1],
                     uxx = quadratic.coefficients[2], uyy = quadratic.coefficients[3],
                     uxy = quadratic.coefficients[4];

        dealii::Vector<double> grad(2);
        grad[0] = ux * (1.0 / h) + uxx * (2.0 / (h * h)) * (x - xo) + uxy * (1.0 / (h * h)) * (y - yo);
        grad[1] = uy * (1.0 / h) + uyy * (2.0 / (h * h)) * (y - yo) + uxy * (1.0 / (h * h)) * (x - xo);

        return grad;
      }

      // returns no of coefficients excluding uo:
      constexpr unsigned int size() const { return 5; }

      // access coefficients of the polynomial
      std::array<double, 5> &coefficient() { return coefficients; }

      // read only access coefficients of the polynomial
      const std::array<double, 5> &coefficient() const { return coefficients; }

    private:
      std::array<double, 5> coefficients; // coefficients excluding uo
      double h;                           // sqrt of cell volume
      double xo, yo;                      // cell centre coordinates
      double c1, c2, c3;                  // polynomial constants
    };

    // template specialization for dim = 3
    template <>
    struct Quadratic<3>
    {
      using CellItr = typename dealii::Triangulation<3>::active_cell_iterator;

      // constructor
      Quadratic() = default;

      void reinit(const CellItr &cell)
      {

        // cell center:
        const auto po = cell->barycenter();

        // compute h:
        h = cbrt(cell->measure());

        // extract coordinates of the cell centre:
        xo = po[0];
        yo = po[1];
        zo = po[2];

        // compute constants c1, c2, c3, c4, c5 and c6:
        CellAverage<3, 4> cell_average;

        c1 = cell_average([this](auto p) {double x = p[0]; return (1.0 / (h * h)) * (x - xo) * (x - xo); }, cell);
        c2 = cell_average([this](auto p) {double y = p[1]; return (1.0 / (h * h)) * (y - yo) * (y - yo); }, cell);
        c3 = cell_average([this](auto p) {double z = p[2]; return (1.0 / (h * h)) * (z - zo) * (z - zo); }, cell);
        c4 = cell_average([this](auto p) {double x = p[0], y = p[1]; return (1.0 / (h * h)) * (x - xo) * (y - yo); }, cell);
        c5 = cell_average([this](auto p) {double y = p[1], z = p[2]; return (1.0 / (h * h)) * (y - yo) * (z - zo); }, cell);
        c6 = cell_average([this](auto p) {double x = p[0], z = p[2]; return (1.0 / (h * h)) * (x - xo) * (z - zo); }, cell);
      }

      // operator(): returns the basis functions of the polynomial evaluated at
      // point p:
      dealii::Vector<double> operator()(const dealii::Point<3, double> &p) const
      {
        const double x = p[0], y = p[1], z = p[2];

        dealii::Vector<double> basis(9);
        basis[0] = (1.0 / h) * (x - xo);
        basis[1] = (1.0 / h) * (y - yo);
        basis[2] = (1.0 / h) * (z - zo);
        basis[3] = ((1.0 / (h * h)) * (x - xo) * (x - xo)) - c1;
        basis[4] = ((1.0 / (h * h)) * (y - yo) * (y - yo)) - c2;
        basis[5] = ((1.0 / (h * h)) * (z - zo) * (z - zo)) - c3;
        basis[6] = ((1.0 / (h * h)) * (x - xo) * (y - yo)) - c4;
        basis[7] = ((1.0 / (h * h)) * (y - yo) * (z - zo)) - c5;
        basis[8] = ((1.0 / (h * h)) * (x - xo) * (z - zo)) - c6;

        return basis;
      }

      // return gradient of the polynomial evaluated at cell center:
      friend dealii::Vector<double> gradient(const Quadratic<3> &quadratic,
                                             const dealii::Point<3, double> &p)
      {
        const double h = quadratic.h, xo = quadratic.xo, yo = quadratic.yo, zo = quadratic.zo;
        const double x = p[0], y = p[1], z = p[2];

        // component of gradient
        double ux = quadratic.coefficients[0], uy = quadratic.coefficients[1],
               uz = quadratic.coefficients[2], uxx = quadratic.coefficients[3],
               uyy = quadratic.coefficients[4], uzz = quadratic.coefficients[5],
               uxy = quadratic.coefficients[6], uyz = quadratic.coefficients[7],
               uzx = quadratic.coefficients[8];

        dealii::Vector<double> grad(3);
        grad[0] = ux * (1.0 / h) + uxx * (2.0 / (h * h)) * (x - xo) + uxy * (1.0 / (h * h)) * (y - yo) + uzx * (1.0 / (h * h)) * (z - zo);
        grad[1] = uy * (1.0 / h) + uyy * (2.0 / (h * h)) * (y - yo) + uxy * (1.0 / (h * h)) * (x - xo) + uyz * (1.0 / (h * h)) * (z - zo);
        grad[2] = uz * (1.0 / h) + uzz * (2.0 / (h * h)) * (z - zo) + uyz * (1.0 / (h * h)) * (y - yo) + uzx * (1.0 / (h * h)) * (x - xo);

        return grad;
      }

      // returns no of coefficients excluding uo:
      constexpr unsigned int size() const { return 9; }

      // access coefficients of the polynomial
      std::array<double, 9> &coefficient() { return coefficients; }

      // read only access coefficients of the polynomial
      const std::array<double, 9> &coefficient() const { return coefficients; }

    private:
      std::array<double, 9> coefficients; // coefficients excluding uo
      double h;                           // sqrt of cell volume
      double xo, yo, zo;                  // cell centre coordinates
      double c1, c2, c3, c4, c5, c6;      // polynomial constants
    };
  } // namespace polynomial

} // namespace fv
