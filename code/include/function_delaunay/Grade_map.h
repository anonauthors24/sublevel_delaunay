#pragma once

#include <mpp_utils/Graded_matrix.h>
#include <function_delaunay/Point_with_density.h>

namespace function_delaunay {

  // Assign filtration values for simplices in simplex tree
  // Assumes that the meb radius is stored in the simplex tree's filtration value
  // and that the input points are sorted by density
  template<typename SimplexTree>
    class Grade_map {
    
  public: 
    typedef SimplexTree Simplex_tree;
    typedef mpp_utils::Grade_struct<double> Grade;
    
  protected:
    Simplex_tree& simplex_tree;
    std::vector<Point_with_density>& input_points;
    
  public:
    
    Grade_map(SimplexTree& tree, std::vector<Point_with_density>& points) 
      : simplex_tree(tree),input_points(points) {
    }
      
      // Returns a grade, as defined in mpp_utils
      template<typename SimplexHandle> Grade operator() (SimplexHandle& handle) {
	
	// We assume that the meb radius has been computed
	// and is stored as filtration value in the simplex
	double meb_radius = simplex_tree.filtration(handle);
	
	int max_vertex_idx = *std::max_element(simplex_tree.simplex_vertex_range(handle).begin(),
					       simplex_tree.simplex_vertex_range(handle).end());
	double density = input_points[max_vertex_idx].density;
	
	return Grade(meb_radius,density);
	
      }
      
  };
}
