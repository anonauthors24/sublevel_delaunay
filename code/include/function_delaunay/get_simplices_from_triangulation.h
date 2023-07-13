#pragma once

#include <function_delaunay/Delaunay_triangulation_accessors.h>

//#define CGAL_EIGEN3_ENABLED
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation.h>

#define GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(d) \
  typedef CGAL::Dimension_tag<d> Dim_tag; \
  typedef Delaunay_triangulation_accessor_d<Dim_tag> DT_accessor; \
  DT_accessor T(d); \
  _get_simplices_from_triangulation(T,input_points,simplices);

namespace function_delaunay {

  template<typename DelaunayTriangulationAccessor> 
    void _get_simplices_from_triangulation(DelaunayTriangulationAccessor& T ,
					   std::vector<Point_with_density> &input_points,
					   std::vector<std::vector<int> > &simplices) {

    
    typedef DelaunayTriangulationAccessor Delaunay_triangulation_accessor;
    
    typedef typename Delaunay_triangulation_accessor::Vertex_handle Vertex_handle;
    typedef typename Delaunay_triangulation_accessor::Full_cell_handle Full_cell_handle;
    typedef typename Delaunay_triangulation_accessor::Point Point;

    std::vector<Vertex_handle> vertices_by_idx;
    
    int counter=0;

    int no_points_outside_convex_hull=0;

    for(Point_with_density& point : input_points) { 
      //std::cout << "******************************" << std::endl;
      //std::cout << "counter=" << counter << std::endl;
      Point p=T.get_point(point);
      
      Vertex_handle vh;
            
      if(T.is_full_dimensional()) {
	int d = T.current_dimension();
	std::vector<Full_cell_handle> cells;
	Full_cell_handle loc = T.locate(p);
	if(T.is_infinite(loc)) {
	  no_points_outside_convex_hull++;
	}
	T.compute_conflict_zone(p,loc,cells);
	
	for(auto cell:cells) {
	  std::vector<int> cell_idx;
	  for(int i=0;i<=d;i++) {
	    if(! T.is_infinite(cell->vertex(i))) {
	      cell_idx.push_back(T.data_of_vertex(cell->vertex(i)).idx);
	    }
	  }
	  std::sort(cell_idx.begin(),cell_idx.end());
	  cell_idx.push_back(counter);
	  simplices.push_back(cell_idx);
	  
	  /*
	    std::cout << "Cell: ";
	    for(int x: cell_idx) {
	    std::cout <<  x << " ";
	    }
	    std::cout << std::endl;
	  */
	}
	vh = T.insert(p,loc);
      } else {
	vh = T.insert(p);
      }
      T.data_of_vertex(vh).idx = counter;
      //T.data_of_vertex(vh).density=point.density;
      vertices_by_idx.push_back(vh);
      
      counter++;
      
    }
    std::cout << "Number of insertions outside convex hull " << no_points_outside_convex_hull << std::endl;
  }

  void get_simplices_from_triangulation(std::vector<Point_with_density> &input_points,
					std::vector<std::vector<int> > &simplices) {

    if(input_points.size()==0) {
      return;
    }
    int d=input_points[0].dimension();

    if(d==2) {
      // Delaunay_triangulation_2 seems faster than the general one
      std::cout << "In the plane, using CGAL::Delaunay_triangulation_2" << std::endl;
      typedef Delaunay_triangulation_accessor_2 DT_accessor;
      DT_accessor T;
      _get_simplices_from_triangulation(T,input_points,simplices);
    } else if(d==3) {
      // Same for Delaunay triangulation_3
      std::cout << "In space, using CGAL::Delaunay_triangulation_3" << std::endl;
      typedef Delaunay_triangulation_accessor_3 DT_accessor;
      DT_accessor T;
      _get_simplices_from_triangulation(T,input_points,simplices);
    } else if(d==4) {
    // Using fixed dimension tags makes the code slightly faster. We dispatch here by
    // dimension, even if it looks a bit unelegant
      GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(4)
    } else if(d==5) {
      GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(5)      
    } else if(d==6) {
      GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(6)
    } else if(d==7) {
      GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(7)
    } else if(d==8) {
      GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(8)
    } else if(d==9) {
      GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(9)
    } else if(d==10) {
      GET_SIMPLICES_FROM_TRIANGULATION_DISPATCH(10)
    } else {
      // If you really want to do that, use the general version
      typedef CGAL::Dynamic_dimension_tag Dim_tag;
      typedef Delaunay_triangulation_accessor_d<Dim_tag> DT_accessor;
      DT_accessor T(d);
      _get_simplices_from_triangulation(T,input_points,simplices);
    }
  }

  

}
