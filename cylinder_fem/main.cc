 
#include <deal.II/grid/grid_out.h>              //graphical output of the grid          
#include <deal.II/grid/grid_in.h>               //import mesh

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/base/utilities.h>



 
  using namespace dealii;



  template <int dim>
  class WaveEquation
  {
  public:
    WaveEquation ();
    void run ();

  private:
    void setup_system ();
    void solve_u ();
    void solve_v ();
    void output_results () const;
 
	 
    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> matrix_u;
    SparseMatrix<double> matrix_v;

    Vector<double>       solution_u, solution_v, wave_speed;
    Vector<double>       old_solution_u, old_solution_v;
    Vector<double>       system_rhs;

    double time_step, time;
    unsigned int timestep_number;
    const double theta;
  };


	//member function definitions:
	#include "functions.h"         
   #include "constructor.h"       
   #include "setup_system.h"      
   #include "solve_U_and_V.h"  
   #include "output.h"
   
   //run calls all the private functions sequentially:  
   #include "run.h"   
   


int main ()
{
   
      using namespace dealii;
       
		 
      WaveEquation<2> wave_equation_solver;
      wave_equation_solver.run();
  
   

  return 0;
}
