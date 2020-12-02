#include "../include/Weno432.h"

// Precompute reconstruction matrices

void Weno4_2D::precompute_matrices() {

  std::cout << "Computing the stencil matrices" << std::endl;

  unsigned int N_gp = 2; // No. of quadrature points
  QGauss<2> quadrature_formula(N_gp);

  FEValues<2> fv_values(fv, quadrature_formula,
                        update_quadrature_points | update_JxW_values);

  Point<2> q_point;
  Point<2> C;
  Point<2> P;

  double V_neighbor;

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();

  FullMatrix<double> A_R4; // Fourth Order Stencil (Least Squares Part)
  FullMatrix<double> C_R4; // Fourth Order Stencil (Constraint Part)

  FullMatrix<double> A_R3; // Least Squares Matrix for r=3 stencil
  FullMatrix<double> C_R3; // Constraint Matrix for r=3 stencil

  // One-sided stencils

  FullMatrix<double> C_R31;
  FullMatrix<double> A_R31;

  FullMatrix<double> C_R32;
  FullMatrix<double> A_R32;

  FullMatrix<double> C_R33;
  FullMatrix<double> A_R33;

  FullMatrix<double> C_R34;
  FullMatrix<double> A_R34;

  FullMatrix<double> A_R21(2, 2); // Matrix for r=2 stencil 1
  FullMatrix<double> A_R22(2, 2); // Matrix for r=2 stencil 2
  FullMatrix<double> A_R23(2, 2); // Matrix for r=2 stencil 3
  FullMatrix<double> A_R24(2, 2); // Matrix for r=2 stencil 4

  unsigned int ROWS, index;

  double x0, y0;

  for (unsigned int c = 0; cell != endc; ++cell, ++c) {

    if (c < no_cells_per_block) {

      // =====================================================================
      // r = 4 stencil
      // =====================================================================

      C_R4.reinit(4, 9);
      C_R4 = 0.0;

      x0 = WENO_poly_consts[c](0);
      y0 = WENO_poly_consts[c](1);

      // Fill C_R4 (Constraint r = 4 matrix)

      // Row 1 (P1 cell)

      fv_values.reinit(cell->neighbor(0));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R4(0, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R4(0, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R4(0, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R4(0, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R4(0, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        C_R4(0, 5) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
             WENO_poly_consts[c](5));
        C_R4(0, 6) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](6));
        C_R4(0, 7) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
             WENO_poly_consts[c](7));
        C_R4(0, 8) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](8));
      }

      // Row 2 (P2 cell)

      fv_values.reinit(cell->neighbor(3));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R4(1, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R4(1, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R4(1, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R4(1, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R4(1, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        C_R4(1, 5) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
             WENO_poly_consts[c](5));
        C_R4(1, 6) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](6));
        C_R4(1, 7) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
             WENO_poly_consts[c](7));
        C_R4(1, 8) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](8));
      }

      // Row 3 (P3 cell)

      fv_values.reinit(cell->neighbor(1));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R4(2, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R4(2, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R4(2, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R4(2, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R4(2, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        C_R4(2, 5) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
             WENO_poly_consts[c](5));
        C_R4(2, 6) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](6));
        C_R4(2, 7) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
             WENO_poly_consts[c](7));
        C_R4(2, 8) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](8));
      }

      // Row 4 (P4 cell)

      fv_values.reinit(cell->neighbor(2));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R4(3, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R4(3, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R4(3, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R4(3, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R4(3, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        C_R4(3, 5) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
             WENO_poly_consts[c](5));
        C_R4(3, 6) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](6));
        C_R4(3, 7) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
             WENO_poly_consts[c](7));
        C_R4(3, 8) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
             WENO_poly_consts[c](8));
      }

      // Least Squares Part

      ROWS = 0;
      index = 0;

      if (cell->neighbor(0)->neighbor_index(0) != -1) {
        ROWS++; // S1 Cell
      }

      if (cell->neighbor(0)->neighbor_index(3) != -1) {
        ROWS++; // S2 cell
      }

      if (cell->neighbor(3)->neighbor_index(3) != -1) {
        ROWS++; // S3 Cell
      }

      if (cell->neighbor(3)->neighbor_index(1) != -1) {
        ROWS++; // S4 Cell
      }

      if (cell->neighbor(1)->neighbor_index(1) != -1) {
        ROWS++; // S5 Cell
      }

      if (cell->neighbor(1)->neighbor_index(2) != -1) {
        ROWS++; // S6 Cell
      }

      if (cell->neighbor(2)->neighbor_index(2) != -1) {
        ROWS++; // S7 Cell
      }

      if (cell->neighbor(2)->neighbor_index(0) != -1) {
        ROWS++; // S8 Cell
      }

      A_R4.reinit(ROWS, 9);
      A_R4 = 0.0;

      // S1 cell

      if (cell->neighbor(0)->neighbor_index(0) != -1) {

        fv_values.reinit(cell->neighbor(0)->neighbor(0));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      // S2 cell

      if (cell->neighbor(0)->neighbor_index(3) != -1) {

        fv_values.reinit(cell->neighbor(0)->neighbor(3));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      // S3 cell

      if (cell->neighbor(3)->neighbor_index(3) != -1) {

        fv_values.reinit(cell->neighbor(3)->neighbor(3));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      // S4 cell

      if (cell->neighbor(3)->neighbor_index(1) != -1) {

        fv_values.reinit(cell->neighbor(3)->neighbor(1));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      // S5 cell

      if (cell->neighbor(1)->neighbor_index(1) != -1) {

        fv_values.reinit(cell->neighbor(1)->neighbor(1));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      // S6 cell

      if (cell->neighbor(1)->neighbor_index(2) != -1) {

        fv_values.reinit(cell->neighbor(1)->neighbor(2));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      // S7 cell

      if (cell->neighbor(2)->neighbor_index(2) != -1) {

        fv_values.reinit(cell->neighbor(2)->neighbor(2));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      // S8 cell

      if (cell->neighbor(2)->neighbor_index(0) != -1) {

        fv_values.reinit(cell->neighbor(2)->neighbor(0));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          A_R4(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(0) - WENO_poly_consts[c](0));
          A_R4(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                            (q_point(1) - WENO_poly_consts[c](1));
          A_R4(index, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          A_R4(index, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          A_R4(index, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
          A_R4(index, 5) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(0) - x0) -
               WENO_poly_consts[c](5));
          A_R4(index, 6) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](6));
          A_R4(index, 7) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) * (q_point(1) - y0) -
               WENO_poly_consts[c](7));
          A_R4(index, 8) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) * (q_point(1) - y0) -
               WENO_poly_consts[c](8));
        }

        index++;
      }

      CLS_R4[c].initialize(A_R4, C_R4);

      // =====================================================================
      // r = 3 stencil (Centered Stencil)
      // =====================================================================

      A_R3.reinit(4, 5);
      A_R3 = 0.0;
      C_R3.reinit(4, 5);
      C_R3 = 0.0;

      // Fill A_R3 (Least squares r = 3 matrix)

      // Row 1 (S2 Cell)

      fv_values.reinit(cell->neighbor(0)->neighbor(3));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        A_R3(0, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        A_R3(0, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        A_R3(0, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        A_R3(0, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        A_R3(0, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      // Row 2 (S4 cell)

      fv_values.reinit(cell->neighbor(3)->neighbor(1));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        A_R3(1, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        A_R3(1, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        A_R3(1, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        A_R3(1, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        A_R3(1, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      // Row 3 (S6 cell)

      fv_values.reinit(cell->neighbor(1)->neighbor(2));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        A_R3(2, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        A_R3(2, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        A_R3(2, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        A_R3(2, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        A_R3(2, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      // Row 4 (S8 cell)

      fv_values.reinit(cell->neighbor(2)->neighbor(0));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        A_R3(3, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        A_R3(3, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        A_R3(3, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        A_R3(3, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        A_R3(3, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      // Fill C_R3 (Constraint r = 3 matrix)

      // Row 1 (P1 cell)

      fv_values.reinit(cell->neighbor(0));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R3(0, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R3(0, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R3(0, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R3(0, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R3(0, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      // Row 2 (P2 cell)

      fv_values.reinit(cell->neighbor(3));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R3(1, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R3(1, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R3(1, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R3(1, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R3(1, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      // Row 3 (P3 cell)

      fv_values.reinit(cell->neighbor(1));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R3(2, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R3(2, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R3(2, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R3(2, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R3(2, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      // Row 4 (P4 cell)

      fv_values.reinit(cell->neighbor(2));

      V_neighbor = 0.0;
      for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
        V_neighbor += fv_values.JxW(i);

      for (unsigned int i = 0; i < N_gp * N_gp; i++) {
        q_point = fv_values.quadrature_point(i);
        C_R3(3, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(0) - WENO_poly_consts[c](0));
        C_R3(3, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                      (q_point(1) - WENO_poly_consts[c](1));
        C_R3(3, 2) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
        C_R3(3, 3) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
        C_R3(3, 4) +=
            (1. / V_neighbor) * fv_values.JxW(i) *
            ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
      }

      CLS_R3[c].initialize(A_R3, C_R3);

      // =====================================================================
      // r = 3 stencil 1
      // =====================================================================

      if (cell->neighbor(0)->neighbor_index(0) != -1 &&
          cell->neighbor(0)->neighbor_index(3) != -1 &&
          cell->neighbor(2)->neighbor_index(0) != -1) {
        is_admissible_R31[c] = true;
      }

      else if (cell->neighbor(0)->neighbor_index(0) != -1 &&
               cell->neighbor(0)->neighbor_index(3) != -1) {
        is_admissible_R31[c] = true;
      }

      else if (cell->neighbor(0)->neighbor_index(0) != -1 &&
               cell->neighbor(2)->neighbor_index(0) != -1) {
        is_admissible_R31[c] = true;
      }

      else {
        is_admissible_R31[c] = false;
      }

      if (is_admissible_R31[c]) {

        // Constraint Matrix

        C_R31.reinit(4, 5);

        C_R31 = 0.0;

        // Fill C_R31 (C_R31onstraint r = 4 matrix)

        // Row 1 (P1 cell)

        fv_values.reinit(cell->neighbor(0));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R31(0, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R31(0, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R31(0, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R31(0, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R31(0, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 2 (P2 cell)

        fv_values.reinit(cell->neighbor(3));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R31(1, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R31(1, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R31(1, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R31(1, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R31(1, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 3 (P4 cell)

        fv_values.reinit(cell->neighbor(2));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R31(2, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R31(2, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R31(2, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R31(2, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R31(2, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 4 (S1 cell)

        fv_values.reinit(cell->neighbor(0)->neighbor(0));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R31(3, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R31(3, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R31(3, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R31(3, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R31(3, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Least Squares Matrix

        index = 0;

        if (cell->neighbor(0)->neighbor_index(3) != -1 &&
            cell->neighbor(2)->neighbor_index(0) != -1) {
          ROWS = 2;
        }

        else {
          ROWS = 1;
        }

        A_R31.reinit(ROWS, 5);

        if (cell->neighbor(0)->neighbor_index(3) != -1) { // S2 cell

          fv_values.reinit(cell->neighbor(0)->neighbor(3));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R31(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R31(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R31(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R31(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R31(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }

          index++;
        }

        if (cell->neighbor(2)->neighbor_index(0) != -1) { // S8 cell

          fv_values.reinit(cell->neighbor(2)->neighbor(0));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R31(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R31(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R31(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R31(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R31(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }
        }

        CLS_R31[c].initialize(A_R31, C_R31);
      }

      // =====================================================================
      // r = 3 stencil 3
      // =====================================================================

      if (cell->neighbor(3)->neighbor_index(3) != -1 &&
          cell->neighbor(0)->neighbor_index(3) != -1 &&
          cell->neighbor(3)->neighbor_index(1) != -1) {
        is_admissible_R32[c] = true;
      }

      else if (cell->neighbor(3)->neighbor_index(3) != -1 &&
               cell->neighbor(0)->neighbor_index(3) != -1) {
        is_admissible_R32[c] = true;
      }

      else if (cell->neighbor(3)->neighbor_index(3) != -1 &&
               cell->neighbor(3)->neighbor_index(1) != -1) {
        is_admissible_R32[c] = true;
      }

      else {
        is_admissible_R32[c] = false;
      }

      if (is_admissible_R32[c]) {

        // Constraint Matrix

        C_R32.reinit(4, 5);

        C_R32 = 0.0;

        // Fill C_R32 (C_R32onstraint r = 4 matrix)

        // Row 1 (P1 cell)

        fv_values.reinit(cell->neighbor(0));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R32(0, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R32(0, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R32(0, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R32(0, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R32(0, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 2 (P2 cell)

        fv_values.reinit(cell->neighbor(3));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R32(1, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R32(1, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R32(1, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R32(1, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R32(1, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 3 (P3 cell)

        fv_values.reinit(cell->neighbor(1));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R32(2, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R32(2, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R32(2, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R32(2, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R32(2, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 4 (S3 cell)

        fv_values.reinit(cell->neighbor(3)->neighbor(3));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R32(3, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R32(3, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R32(3, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R32(3, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R32(3, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Least Squares Matrix

        index = 0;

        if (cell->neighbor(0)->neighbor_index(3) != -1 &&
            cell->neighbor(3)->neighbor_index(1) != -1) {
          ROWS = 2;
        }

        else {
          ROWS = 1;
        }

        A_R32.reinit(ROWS, 5);

        if (cell->neighbor(0)->neighbor_index(3) != -1) { // S2 cell

          fv_values.reinit(cell->neighbor(0)->neighbor(3));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R32(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R32(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R32(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R32(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R32(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }

          index++;
        }

        if (cell->neighbor(3)->neighbor_index(1) != -1) { // S4 cell

          fv_values.reinit(cell->neighbor(3)->neighbor(1));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R32(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R32(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R32(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R32(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R32(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }
        }

        CLS_R32[c].initialize(A_R32, C_R32);
      }

      // =====================================================================
      // r = 3 stencil 3
      // =====================================================================

      if (cell->neighbor(2)->neighbor_index(2) != -1 &&
          cell->neighbor(1)->neighbor_index(2) != -1 &&
          cell->neighbor(2)->neighbor_index(0) != -1) {
        is_admissible_R33[c] = true;
      }

      else if (cell->neighbor(2)->neighbor_index(2) != -1 &&
               cell->neighbor(1)->neighbor_index(2) != -1) {
        is_admissible_R33[c] = true;
      }

      else if (cell->neighbor(2)->neighbor_index(2) != -1 &&
               cell->neighbor(2)->neighbor_index(0) != -1) {
        is_admissible_R33[c] = true;
      }

      else {
        is_admissible_R33[c] = false;
      }

      if (is_admissible_R33[c]) {

        // Constraint Matrix

        C_R33.reinit(4, 5);

        C_R33 = 0.0;

        // Fill C_R33 (C_R33onstraint r = 4 matrix)

        // Row 1 (P1 cell)

        fv_values.reinit(cell->neighbor(0));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R33(0, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R33(0, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R33(0, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R33(0, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R33(0, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 2 (P3 cell)

        fv_values.reinit(cell->neighbor(1));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R33(1, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R33(1, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R33(1, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R33(1, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R33(1, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 3 (P4 cell)

        fv_values.reinit(cell->neighbor(2));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R33(2, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R33(2, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R33(2, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R33(2, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R33(2, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 4 (S7 cell)

        fv_values.reinit(cell->neighbor(2)->neighbor(2));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R33(3, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R33(3, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R33(3, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R33(3, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R33(3, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Least Squares Matrix

        index = 0;

        if (cell->neighbor(1)->neighbor_index(2) != -1 &&
            cell->neighbor(2)->neighbor_index(0) != -1) {
          ROWS = 2;
        }

        else {
          ROWS = 1;
        }

        A_R33.reinit(ROWS, 5);

        if (cell->neighbor(1)->neighbor_index(2) != -1) { // S6 cell

          fv_values.reinit(cell->neighbor(1)->neighbor(2));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R33(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R33(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R33(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R33(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R33(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }

          index++;
        }

        if (cell->neighbor(2)->neighbor_index(0) != -1) { // S8 cell

          fv_values.reinit(cell->neighbor(2)->neighbor(0));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R33(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R33(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R33(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R33(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R33(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }
        }

        CLS_R33[c].initialize(A_R33, C_R33);
      }

      // =====================================================================
      // r = 3 stencil 4
      // =====================================================================

      if (cell->neighbor(1)->neighbor_index(1) != -1 &&
          cell->neighbor(3)->neighbor_index(1) != -1 &&
          cell->neighbor(1)->neighbor_index(2) != -1) {
        is_admissible_R34[c] = true;
      }

      else if (cell->neighbor(1)->neighbor_index(1) != -1 &&
               cell->neighbor(3)->neighbor_index(1) != -1) {
        is_admissible_R34[c] = true;
      }

      else if (cell->neighbor(1)->neighbor_index(1) != -1 &&
               cell->neighbor(1)->neighbor_index(2) != -1) {
        is_admissible_R34[c] = true;
      }

      else {
        is_admissible_R34[c] = false;
      }

      if (is_admissible_R34[c]) {

        // Constraint Matrix

        C_R34.reinit(4, 5);

        C_R34 = 0.0;

        // Fill C_R34 (C_R34onstraint r = 4 matrix)

        // Row 1 (P2 cell)

        fv_values.reinit(cell->neighbor(3));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R34(0, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R34(0, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R34(0, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R34(0, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R34(0, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 2 (P3 cell)

        fv_values.reinit(cell->neighbor(1));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R34(1, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R34(1, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R34(1, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R34(1, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R34(1, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 3 (P4 cell)

        fv_values.reinit(cell->neighbor(2));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R34(2, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R34(2, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R34(2, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R34(2, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R34(2, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Row 4 (S5 cell)

        fv_values.reinit(cell->neighbor(1)->neighbor(1));

        V_neighbor = 0.0;
        for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
          V_neighbor += fv_values.JxW(i);

        for (unsigned int i = 0; i < N_gp * N_gp; i++) {
          q_point = fv_values.quadrature_point(i);
          C_R34(3, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(0) - WENO_poly_consts[c](0));
          C_R34(3, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                         (q_point(1) - WENO_poly_consts[c](1));
          C_R34(3, 2) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(0) - x0) - WENO_poly_consts[c](2));
          C_R34(3, 3) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(1) - y0) * (q_point(1) - y0) - WENO_poly_consts[c](3));
          C_R34(3, 4) +=
              (1. / V_neighbor) * fv_values.JxW(i) *
              ((q_point(0) - x0) * (q_point(1) - y0) - WENO_poly_consts[c](4));
        }

        // Least Squares Matrix

        index = 0;

        if (cell->neighbor(3)->neighbor_index(1) != -1 &&
            cell->neighbor(2)->neighbor_index(1) != -1) {
          ROWS = 2;
        }

        else {
          ROWS = 1;
        }

        A_R34.reinit(ROWS, 5);

        if (cell->neighbor(3)->neighbor_index(1) != -1) { // S4 cell

          fv_values.reinit(cell->neighbor(3)->neighbor(1));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R34(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R34(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R34(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R34(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R34(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }

          index++;
        }

        if (cell->neighbor(2)->neighbor_index(1) != -1) { // S6 cell

          fv_values.reinit(cell->neighbor(2)->neighbor(1));

          V_neighbor = 0.0;
          for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
            V_neighbor += fv_values.JxW(i);

          for (unsigned int i = 0; i < N_gp * N_gp; i++) {
            q_point = fv_values.quadrature_point(i);
            A_R34(index, 0) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(0) - WENO_poly_consts[c](0));
            A_R34(index, 1) += (1. / V_neighbor) * fv_values.JxW(i) *
                               (q_point(1) - WENO_poly_consts[c](1));
            A_R34(index, 2) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(0) - x0) -
                                WENO_poly_consts[c](2));
            A_R34(index, 3) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(1) - y0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](3));
            A_R34(index, 4) += (1. / V_neighbor) * fv_values.JxW(i) *
                               ((q_point(0) - x0) * (q_point(1) - y0) -
                                WENO_poly_consts[c](4));
          }
        }

        CLS_R34[c].initialize(A_R34, C_R34);
      }

      // =====================================================================
      // r = 2 stencil 1 (P1 and P2)
      // =====================================================================

      A_R21(0, 0) = C_R3(0, 0);
      A_R21(0, 1) = C_R3(0, 1); // Row 1 (P1 cell)
      A_R21(1, 0) = C_R3(1, 0);
      A_R21(1, 1) = C_R3(1, 1); // Row 2 (P2 cell)

      LU_R21[c].initialize(A_R21);

      // =====================================================================
      // r = 2 stencil 2 (P2 and P3)
      // =====================================================================

      A_R22(0, 0) = C_R3(1, 0);
      A_R22(0, 1) = C_R3(1, 1); // Row 1 (P2 cell)
      A_R22(1, 0) = C_R3(2, 0);
      A_R22(1, 1) = C_R3(2, 1); // Row 2 (P3 cell)

      LU_R22[c].initialize(A_R22);

      // =====================================================================
      // r = 2 stencil 3 (P3 and P4)
      // =====================================================================

      A_R23(0, 0) = C_R3(2, 0);
      A_R23(0, 1) = C_R3(2, 1); // Row 1 (P3 cell)
      A_R23(1, 0) = C_R3(3, 0);
      A_R23(1, 1) = C_R3(3, 1); // Row 2 (P4 cell)

      LU_R23[c].initialize(A_R23);

      // =====================================================================
      // r = 2 stencil 4 (P4 and P1)
      // =====================================================================

      A_R24(0, 0) = C_R3(3, 0);
      A_R24(0, 1) = C_R3(3, 1); // Row 1 (P4 cell)
      A_R24(1, 0) = C_R3(0, 0);
      A_R24(1, 1) = C_R3(0, 1); // Row 2 (P1 cell)

      LU_R24[c].initialize(A_R24);

    } // End of center block loop

  } // End of cell loop

  std::cout << "Done!" << std::endl;
  std::cout << "===========================" << std::endl;

} // End of function
