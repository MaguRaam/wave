#ifndef HEADERS_H_
#define HEADERS_H_

// C++ STD and STL Headers 
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <ctime>
#include <assert.h> 

// GSL Headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// deal.II Headers 
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/householder.h>
#include <deal.II/numerics/vector_tools.h>

// Namespace declaration 
using namespace dealii; 

#endif /* HEADERS_H_ */
