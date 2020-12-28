#include "../include/Weno432.h"

// Update the solution using SSPRK(5,4) method

void Weno4_2D::solve_ssprk54()
{

    // Using SSPRK(5,4)

    double time = 0.0;
    unsigned int count = 0;

    Vector<double> U_old(dof_handler.n_dofs());
    Vector<double> V_old(dof_handler.n_dofs());

    Vector<double> U_2(dof_handler.n_dofs());
    Vector<double> V_2(dof_handler.n_dofs());

    Vector<double> U_3(dof_handler.n_dofs());
    Vector<double> V_3(dof_handler.n_dofs());

    Vector<double> rhs_3u(no_cells_per_block);
    Vector<double> rhs_3v(no_cells_per_block);

    //std::cout<<"test1"<<std::endl;
    L_norm(time);
    //std::cout<<"test2"<<std::endl;
    while (time < finalTime)
    {

        compute_time_step_based_on_cfl_number(time);

        int output_count;
        output_count = 0.5 / dt;

        if (count == 0)
        {
            copy_data();
        }

        if (count % output_count == 0)
        {
            output_results(count);
        }

        time += dt;

        std::cout << "time = " << time << ", Final time = " << finalTime << std::endl;
        count++;

        U_old = U;
        V_old = V;

        compute_rhs();

        // SSPRK Stage 1

        compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c)
        {
            if (c < no_cells_per_block)
            {
                U(c) = U(c) + 0.391752226571890 * dt * rhs1(c);
                V(c) = V(c) + 0.391752226571890 * dt * rhs2(c);
            }
        }

        // SSPRK Stage 2

        compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c)
        {
            if (c < no_cells_per_block)
            {
                U(c) = 0.444370493651235 * U_old(c) + 0.555629506348765 * U(c) + 0.368410593050371 * dt * rhs1(c);
                V(c) = 0.444370493651235 * V_old(c) + 0.555629506348765 * V(c) + 0.368410593050371 * dt * rhs2(c);
            }
        }

        U_2 = U;
        V_2 = V;

        // SSPRK Stage 3

        compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c)
        {
            if (c < no_cells_per_block)
            {
                U(c) = 0.620101851488403 * U_old(c) + 0.379898148511597 * U(c) + 0.251891774271694 * dt * rhs1(c);
                V(c) = 0.620101851488403 * V_old(c) + 0.379898148511597 * V(c) + 0.251891774271694 * dt * rhs2(c);
            }
        }

        U_3 = U;
        V_3 = V;

        // SSPRK Stage 4

        compute_rhs();

        rhs_3u = rhs1;
        rhs_3v = rhs2;
        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c)
        {
            if (c < no_cells_per_block)
            {
                U(c) = 0.178079954393132 * U_old(c) + 0.821920045606868 * U(c) + 0.544974750228521 * dt * rhs1(c);
                V(c) = 0.178079954393132 * V_old(c) + 0.821920045606868 * V(c) + 0.544974750228521 * dt * rhs2(c);
            }
        }

        //  SSPRK Stage 5

        compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c)
        {
            if (c < no_cells_per_block)
            {
                U(c) = 0.517231671970585 * U_2(c) + 0.096059710526147 * U_3(c) + 0.386708617503269 * U(c) +
                       0.226007483236906 * dt * rhs1(c) + 0.063692468666290 * dt * rhs_3u(c);
                V(c) = 0.517231671970585 * V_2(c) + 0.096059710526147 * V_3(c) + 0.386708617503269 * V(c) +
                       0.226007483236906 * dt * rhs2(c) + 0.063692468666290 * dt * rhs_3v(c);
            }
        }
        L_norm(time);
    }

    //output_results(count);
}
