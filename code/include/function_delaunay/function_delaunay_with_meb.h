#pragma once

#include <mpp_utils/create_graded_matrices_from_simplex_tree.h>
#include <gudhi/Simplex_tree.h>

#include <function_delaunay/Point_with_density.h>
#include <function_delaunay/compute_meb_radii_of_simplex_tree.h>
#include <function_delaunay/Grade_map.h>
#include <function_delaunay/get_simplices_from_triangulation.h>

namespace function_delaunay {


  // Returns the total number of simplices
  template<typename GrMat>
    long function_delaunay_with_meb(std::vector<Point_with_density> &input_points,
				    std::vector<GrMat> &graded_matrices,
				    bool only_return_complex_size=false) {

    if(input_points.size()==0) {
      return 0;
    }
    int d = input_points[0].dimension();

#if FUNCTION_DELAUNAY_TIMERS
    complex_timer.start();
#endif

  
    std::cout << "Start complex" << std::endl;
    
    // The simplices are stored just by boundary vertex indices first
    // The filtration values and the order are computed later!
    
    std::vector<std::vector<int> > simplices;

    // Insert the first simplex by hand
    std::vector<int> initial_simplex;
    for(int i=0;i<=d;i++) {
      initial_simplex.push_back(i);
    }
    simplices.push_back(initial_simplex);
    
    get_simplices_from_triangulation(input_points,simplices);
        
#if FUNCTION_DELAUNAY_TIMERS
    complex_timer.stop();
#endif
    
    std::cout << "Collected " << simplices.size() << " simplices" << std::endl;

    
#if WITH_MEMORY_PROFILE
    std::cout << "Memory after complex: " << mem_info() << std::endl;
#endif

    
#if FUNCTION_DELAUNAY_TIMERS
    face_timer.start();
#endif
  
    
    typedef Gudhi::Simplex_tree<> Simplex_tree;
    
    Simplex_tree simplex_tree;
    
    for(auto simplex : simplices) {
      simplex_tree.insert_simplex_and_subfaces(simplex);
    }
    
#if FUNCTION_DELAUNAY_TIMERS
    face_timer.stop();
#endif
    
    long total_number_of_simplices = simplex_tree.num_simplices();
    
    std::cout << "Simplex tree has " << simplex_tree.num_vertices() << " vertices and " << total_number_of_simplices << " simplices" << std::endl;
    
#if WITH_MEMORY_PROFILE
    std::cout << "Memory after face: " << mem_info() << std::endl;
#endif
    
    if(only_return_complex_size) {
      return total_number_of_simplices;
    }

    
#if FUNCTION_DELAUNAY_TIMERS
    meb_timer.start();
#endif

    compute_meb_radii_of_simplex_tree(simplex_tree,input_points);

    std::cout << "Computed all meb values, now sorting" << std::endl;

#if !NDEBUG || !FUNCTION_DELAUNAY_SMART_MEB_TRAVERSAL

    // This should actually not be necessary, but the numerical issues
    // in meb computation seem to cause problems, in a way that the 
    // meb radius of a simplex can be slightly smaller than for one of
    // its faces
    bool adjustment_needed = simplex_tree.make_filtration_non_decreasing();

    if(adjustment_needed) {
      std::cout << "Some meb value had to be adjusted to fit its cofacet" << std::endl;
    } else {
      std::cout << "No adjustment of meb values was needed" << std::endl;
    }
    
#if FUNCTION_DELAUNAY_SMART_MEB_TRAVERSAL
    assert(!adjustment_needed);
#endif

#endif



#if FUNCTION_DELAUNAY_TIMERS
    meb_timer.stop();
#endif
    
#if WITH_MEMORY_PROFILE
    std::cout << "Memory after bigrade: " << mem_info() << std::endl;
#endif
    
    
    std::cout << "Building graded boundary matrices" << std::endl;
    
#if FUNCTION_DELAUNAY_TIMERS
    graded_matrices_timer.start();
#endif
    
    Grade_map<Simplex_tree> grade_map(simplex_tree,input_points);
    
    // The "false" is for avoiding to compute the grade indices, which are not needed except for mpfree
    mpp_utils::create_graded_matrices_from_simplex_tree(simplex_tree,grade_map,graded_matrices,false);
    
    int simplex_tree_dimension=simplex_tree.dimension();
    
    // The simplex tree is no longer needed at this point
    // Sort of a hack, but it seems to work
    simplex_tree.prune_above_dimension(-1);
  
    for(int i=0;i<=simplex_tree_dimension;i++) {
      std::cout << "Simplices in dimension " << i << ": " << graded_matrices[simplex_tree_dimension-i].get_num_cols() << std::endl;
    }
    
#if FUNCTION_DELAUNAY_TIMERS
    graded_matrices_timer.stop();
#endif

#if WITH_MEMORY_PROFILE
  std::cout << "Memory after boundary: " << mem_info() << std::endl;
#endif
  
  return total_number_of_simplices;
  }

  
}
  
