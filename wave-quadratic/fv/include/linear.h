#pragma once

namespace fv
{
  namespace polynomial
  {
    // Linear polynomial object constructed on a cell:
    template <int dim>
    struct Linear
    {
    };

    // template specialization for dim = 2;
    template <>
    struct Linear<2>
    {
      using CellItr = typename dealii::Triangulation<2>::active_cell_iterator;

      Linear() = default;

      // constructor
      void reinit(const CellItr &cell)
      {

        // cell center:
        const auto po = cell->barycenter();

        // compute h:
        h = sqrt(cell->measure());

        // extract coordinates of the cell centre:
        xo = po[0];
        yo = po[1];
      }

      // operator(): returns the basis functions of the polynomial evaluated at
      // point p:
      dealii::Vector<double> operator()(const dealii::Point<2, double> &p) const
      {
        const double x = p[0], y = p[1];

        dealii::Vector<double> basis(2);
        basis[0] = (1.0 / h) * (x - xo);
        basis[1] = (1.0 / h) * (y - yo);
        return basis;
      }

      // return gradient of the polynomial evaluated at cell center:
      friend dealii::Vector<double> gradient(const Linear<2> &linear,
                                             const dealii::Point<2, double> &p)
      {
        // compute h:
        const double h = linear.h;

        // x and y component of gradient
        double ux = linear.coefficients[0], uy = linear.coefficients[1];
        dealii::Vector<double> grad(2);
        grad[0] = ux / h;
        grad[1] = uy / h;

        return grad;
      }

      // returns no of coefficients excluding uo:
      constexpr unsigned int size() const { return 2; }

      // access coefficients of the polynomial
      std::array<double, 2> &coefficient() { return coefficients; }

      // read only access coefficients of the polynomial
      const std::array<double, 2> &coefficient() const { return coefficients; }

    private:
      std::array<double, 2> coefficients; // coefficients excluding uo
      double h;                           // sqrt of cell volume
      double xo, yo;                      // cell centre coordinates
    };

    // template specialization for dim = 3
    template <>
    struct Linear<3>
    {
      using CellItr = typename dealii::Triangulation<3>::active_cell_iterator;

      // constructor
      Linear() = default;

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
      }

      // operator(): returns the basis functions of the polynomial evaluated at
      // point p:
      dealii::Vector<double> operator()(const dealii::Point<3, double> &p) const
      {
        const double x = p[0], y = p[1], z = p[2];

        dealii::Vector<double> basis(3);

        basis[0] = (1.0 / h) * (x - xo);
        basis[1] = (1.0 / h) * (y - yo);
        basis[2] = (1.0 / h) * (z - zo);

        return basis;
      }

      // return gradient of the polynomial evaluated at cell center:
      friend dealii::Vector<double> gradient(const Linear<3> &linear,
                                             const dealii::Point<3, double> &p)
      {
        // compute h:
        const double h = linear.h;

        // component of gradient
        double ux = linear.coefficients[0], uy = linear.coefficients[1],
               uz = linear.coefficients[2];

        dealii::Vector<double> grad(3);
        grad[0] = ux / h;
        grad[1] = uy / h;
        grad[2] = uz / h;

        return grad;
      }

      // returns no of coefficients excluding uo:
      constexpr unsigned int size() const { return 3; }

      // access coefficients of the polynomial
      std::array<double, 3> &coefficient() { return coefficients; }

      // read only access coefficients of the polynomial
      const std::array<double, 3> &coefficient() const { return coefficients; }

    private:
      std::array<double, 3> coefficients; // coefficients excluding uo
      double h;                           // sqrt of cell volume
      double xo, yo, zo;                  // cell centre coordinates
    };
  } // namespace polynomial

} // namespace fv
