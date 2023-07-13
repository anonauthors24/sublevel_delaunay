#pragma once
//mod from density_delcech Grade_map.h
#include <mpp_utils/Graded_matrix.h>
#include <bifiltered_graph/bifiltered_edge.h>

namespace bifiltered_graph {

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
    std::vector<weightedVertex>& input_points;
    //std::vector<Point_with_density>& input_points;
    
  public:
    
    Grade_map(SimplexTree& tree, std::vector<weightedVertex>& points) 
      : simplex_tree(tree),input_points(points) {
    }
      
      // Returns a grade, as defined in mpp_utils
      template<typename SimplexHandle> Grade operator() (SimplexHandle& handle) {
	
	// We assume that the meb radius has been computed
	// and is stored as filtration value in the simplex
	double graph_radius  = simplex_tree.filtration(handle);
	
	int max_vertex_idx = *std::max_element(simplex_tree.simplex_vertex_range(handle).begin(),
					       simplex_tree.simplex_vertex_range(handle).end());

	double weight = input_points[max_vertex_idx].weight;

        //std::cout << "" <<  max_vertex_idx << " " << graph_radius << " " << weight << std::endl;
	
	return Grade(graph_radius, weight);
	
      }
      
  };
}

