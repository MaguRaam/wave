#include "../include/Weno432.h"

// Output the results in data files 

void Weno4_2D::output_results(unsigned int i) {
    
    DataOut<2> data_out1;
    data_out1.attach_dof_handler (dof_handler);
    data_out1.add_data_vector (U, "U", DataOut<2>::type_dof_data);
    data_out1.add_data_vector (local_difference, "error", DataOut<2>::type_dof_data);
    data_out1.add_data_vector (Uexact,"exact_solution",DataOut<2>::type_dof_data);
    data_out1.build_patches ();
    const std::string filename1 = "../plot/U_" + Utilities::int_to_string (i, 5) + ".vtk";
    std::ofstream output1 (filename1.c_str());
    data_out1.write_vtk(output1);
    

} 
