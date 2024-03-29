#include "../include/Weno432.h"

 



void Weno4_2D::make_grid() {

    std::cout << "===========================" << std::endl;
    std::cout << "Making grid" << std::endl;

    std::vector<unsigned int> repetions(2); // No. of cells in x and y directions 
    repetions[0] = cell;
    repetions[1] = cell; 
    
    bool colorize = true;  // Set boundary ids for the four boundaries 
    
    // Diagonal points of the domain 

    Point<2> P1(-2.0, -2.0);
    Point<2> P2( 2.0,  2.0);
    
    /* Center triangulation */ 
    
    Triangulation<2> triangulation0;
    GridGenerator::subdivided_hyper_rectangle(triangulation0, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation0, true);
    
    
    /* Left center */ 
    P1(0) = -6.0; P1(1) =-2.0;
    P2(0) = -2.0; P2(1) = 2.0;
    Triangulation<2> triangulation1;
    GridGenerator::subdivided_hyper_rectangle(triangulation1, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation1, true);
    //GridTools::shift (Point<2>(-4.0, 0.0), triangulation1);
    
    /* Right center */ 
    P1(0) = 2.0; P1(1) =-2.0;
    P2(0) = 6.0; P2(1) = 2.0;
    Triangulation<2> triangulation2;
    GridGenerator::subdivided_hyper_rectangle(triangulation2, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation2, true);
    //GridTools::shift (Point<2>(4.0, 0.0), triangulation2);
    
    /* Bottom center */ 
    P1(0) = -2.0; P1(1) =-6.0;
    P2(0) = 2.0; P2(1) = -2.0;
    Triangulation<2> triangulation3;
    GridGenerator::subdivided_hyper_rectangle(triangulation3, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation3, true);
    //GridTools::shift (Point<2>(0.0, -4.0), triangulation3);
    
    /* Top center */ 
    P1(0) = -2.0; P1(1) =2.0;
    P2(0) = 2.0; P2(1) = 6.0;
    Triangulation<2> triangulation4;
    GridGenerator::subdivided_hyper_rectangle(triangulation4, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation4, true);
    //GridTools::shift (Point<2>(0.0, 4.0), triangulation4);
    
    /* Top Left */ 
    P1(0) = -6.0; P1(1) = 2.0;
    P2(0) = -2.0; P2(1) = 6.0;
    Triangulation<2> triangulation5;
    GridGenerator::subdivided_hyper_rectangle(triangulation5, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation5, true);
    //GridTools::shift (Point<2>(-4.0, 4.0), triangulation5);
    
    /* Top Right */ 
    P1(0) = 2.0; P1(1) = 2.0;
    P2(0) = 6.0; P2(1) = 6.0;
    Triangulation<2> triangulation6;
    GridGenerator::subdivided_hyper_rectangle(triangulation6, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation6, true);
    //GridTools::shift (Point<2>(4.0, 4.0), triangulation6);
    
    /* Bottom Leftt */
    P1(0) = -6.0; P1(1) =-6.0;
    P2(0) = -2.0; P2(1) = -2.0;
    Triangulation<2> triangulation7;
    GridGenerator::subdivided_hyper_rectangle(triangulation7, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation7, true);
    //GridTools::shift (Point<2>(-4.0, -4.0), triangulation7);
    
    /* Bottom Right */
    P1(0) = 2.0; P1(1) =-6.0;
    P2(0) = 6.0; P2(1) = -2.0;
    Triangulation<2> triangulation8;
    GridGenerator::subdivided_hyper_rectangle(triangulation8, repetions, P1, P2, colorize);
    //GridTools::distort_random (0.15, triangulation8, true);
    //GridTools::shift (Point<2>(4.0, -4.0), triangulation8); 
   

    auto grid_transform = [](const Point<2>& in)
    {
     return Point<2>(in(0) + 0.3*std::sin(2.0*in(1)*M_PI/4.0),in(1)+0.4*std::sin(2.0*in(0)*M_PI/4.0));
    };

    /* Merge all of them */ 
    
    GridGenerator::merge_triangulations (triangulation4, triangulation6, triangulation);
    GridGenerator::merge_triangulations (triangulation5, triangulation, triangulation);
    GridGenerator::merge_triangulations (triangulation2, triangulation, triangulation);
    GridGenerator::merge_triangulations (triangulation1, triangulation, triangulation);
    GridGenerator::merge_triangulations (triangulation8, triangulation, triangulation);
    GridGenerator::merge_triangulations (triangulation3, triangulation, triangulation);
    GridGenerator::merge_triangulations (triangulation7, triangulation, triangulation);
    GridGenerator::merge_triangulations (triangulation0, triangulation, triangulation);

    GridTools::transform(grid_transform, triangulation);

    no_cells_per_block = repetions[0]*repetions[1]; 

	Triangulation<2>::active_cell_iterator  cell = triangulation.begin_active(), 
											endc = triangulation.end();
	for (; cell!=endc; ++cell)
	  	cell->set_all_manifold_ids (0);
 
}


