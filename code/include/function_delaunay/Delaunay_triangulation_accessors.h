#pragma once

//#define CGAL_EIGEN3_ENABLED
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation.h>


namespace function_delaunay {

  struct Vertex_info {
    //double density;
    int idx;
  };

  template<typename DimTag>
    struct Delaunay_triangulation_accessor_d {
      
      typedef DimTag Dim_tag;

      typedef CGAL::Epick_d< Dim_tag >  Kernel;
      typedef CGAL::Triangulation_data_structure<Dim_tag, CGAL::Triangulation_vertex<Kernel,Vertex_info>,CGAL::Triangulation_full_cell<Kernel> > TDS;
      typedef CGAL::Delaunay_triangulation<Kernel,TDS> Triangulation;

      typedef typename Triangulation::Vertex_handle Vertex_handle;
      typedef typename Triangulation::Full_cell_handle Full_cell_handle;
      typedef typename Triangulation::Point Point;
    

      Triangulation T;

      int d;

      Delaunay_triangulation_accessor_d(int d) : T(d), d(d) {

	Vertex_handle inf = T.infinite_vertex();
	inf->data().idx=-1;
	//inf->data().density=-999.999;
      }
      
      Point get_point(Point_with_density& p) {
	return Point(p.x.begin(),p.x.end());
      }
      
      int current_dimension() {
	return T.current_dimension();
      }

      bool is_full_dimensional() {
	return current_dimension()==d;
      }

      Full_cell_handle locate(Point& p) {
	return T.locate(p);
      }

      void compute_conflict_zone(Point p,Full_cell_handle loc,std::vector<Full_cell_handle>& cells) {
	T.compute_conflict_zone(p,loc,std::back_inserter(cells));
      }

      bool is_infinite(Vertex_handle vh) {
	return T.is_infinite(vh);
      }

      bool is_infinite(Full_cell_handle c) {
	return T.is_infinite(c);
      }
      
      Vertex_handle insert(Point& p) {
	return T.insert(p);
      }

      Vertex_handle insert(Point& p, Full_cell_handle c) {
	return T.insert(p,c);
      }

      template<typename InputIterator>
      void insert(InputIterator begin, InputIterator end) {
	T.insert(begin,end);
      }

        
      typename Triangulation::Finite_vertex_iterator vertices_begin() {
	return T.finite_vertices_begin();
      }

      typename Triangulation::Finite_vertex_iterator vertices_end() {
	return T.finite_vertices_end();
      }

      typename Triangulation::Finite_full_cell_iterator finite_full_cells_begin() {
	return T.finite_full_cells_begin();
      }

      typename Triangulation::Finite_full_cell_iterator finite_full_cells_end() {
	return T.finite_full_cells_end();
      }

      // This must be templated because some Vertex-iterators cannot be casted to Vertex_handle
      template<typename VertexHandle>
      Vertex_info& data_of_vertex(VertexHandle vh) {
	return vh->data();
      }

    };

    struct Delaunay_triangulation_accessor_2 {
      
      typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
      typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info,Kernel> Vb;
      typedef CGAL::Triangulation_face_base_2<Kernel> Fb;
      typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
      typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> Triangulation;

      typedef typename Triangulation::Vertex_handle Vertex_handle;
      typedef typename Triangulation::Face_handle Full_cell_handle;
      typedef typename Triangulation::Point Point;
    

      Triangulation T;

      Delaunay_triangulation_accessor_2() {

	Vertex_handle inf = T.infinite_vertex();
	inf->info().idx=-1;
	//inf->info().density=-999.999;
      }
      
      Point get_point(Point_with_density& p) {
	return Point(p.x[0],p.x[1]);
      }
      
      int current_dimension() {
	return T.dimension();
      }

      bool is_full_dimensional() {
	return current_dimension()==2;
      }
      

      Full_cell_handle locate(Point& p) {
	return T.locate(p);
      }

      void compute_conflict_zone(Point p,Full_cell_handle loc,std::vector<Full_cell_handle>& cells) {
	T.get_conflicts(p,std::back_inserter(cells),loc);
      }

      bool is_infinite(Vertex_handle vh) {
	return T.is_infinite(vh);
      }

      bool is_infinite(Full_cell_handle c) {
	return T.is_infinite(c);
      }
      
      Vertex_handle insert(Point& p) {
	return T.insert(p);
      }

      Vertex_handle insert(Point& p, Full_cell_handle c) {
	return T.insert(p,c);
      }

      template<typename InputIterator>
      void insert(InputIterator begin, InputIterator end) {
	T.insert(begin,end);
      }

      typename Triangulation::Finite_vertices_iterator vertices_begin() {
	return T.finite_vertices_begin();
      }

      typename Triangulation::Finite_vertices_iterator vertices_end() {
	return T.finite_vertices_end();
      }

      typename Triangulation::Finite_faces_iterator finite_full_cells_begin() {
	return T.finite_faces_begin();
      }

      typename Triangulation::Finite_faces_iterator finite_full_cells_end() {
	return T.finite_faces_end();
      }

      Vertex_info& data_of_vertex(Vertex_handle vh) {
	return vh->info();
      }

  
    };

    struct Delaunay_triangulation_accessor_3 {
      
      typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
      typedef CGAL::Triangulation_vertex_base_with_info_3<Vertex_info,Kernel> Vb;
      typedef CGAL::Triangulation_cell_base_3<Kernel>               Fb;
      typedef CGAL::Triangulation_data_structure_3<Vb,Fb>      Tds;
      typedef CGAL::Delaunay_triangulation_3<Kernel,Tds,CGAL::Fast_location> Triangulation;

      typedef typename Triangulation::Vertex_handle Vertex_handle;
      typedef typename Triangulation::Cell_handle Full_cell_handle;
      typedef typename Triangulation::Point Point;
      typedef Triangulation::Facet Facet;

      Triangulation T;

      Delaunay_triangulation_accessor_3() {

	Vertex_handle inf = T.infinite_vertex();
	inf->info().idx=-1;
	//inf->info().density=-999.999;
      }
      
      Point get_point(Point_with_density& p) {
	return Point(p.x[0],p.x[1],p.x[2]);
      }
      
      int current_dimension() {
	return T.dimension();
      }

      bool is_full_dimensional() {
	return current_dimension()==3;
      }
      

      Full_cell_handle locate(Point& p) {
	return T.locate(p);
      }

      void compute_conflict_zone(Point p,Full_cell_handle loc,std::vector<Full_cell_handle>& cells) {
	std::vector<Facet> facets;
	T.find_conflicts(p,loc,std::back_inserter(facets),std::back_inserter(cells));
      }

      bool is_infinite(Vertex_handle vh) {
	return T.is_infinite(vh);
      }

      bool is_infinite(Full_cell_handle c) {
	return T.is_infinite(c);
      }

      
      Vertex_handle insert(Point& p) {
	return T.insert(p);
      }

      Vertex_handle insert(Point& p, Full_cell_handle c) {
	return T.insert(p,c);
      }

      template<typename InputIterator>
      void insert(InputIterator begin, InputIterator end) {
	T.insert(begin,end);
      }

      typename Triangulation::Finite_vertices_iterator vertices_begin() {
	return T.finite_vertices_begin();
      }

      typename Triangulation::Finite_vertices_iterator vertices_end() {
	return T.finite_vertices_end();
      }

      typename Triangulation::Finite_cells_iterator finite_full_cells_begin() {
	return T.finite_cells_begin();
      }

      typename Triangulation::Finite_cells_iterator finite_full_cells_end() {
	return T.finite_cells_end();
      }
  
      Vertex_info& data_of_vertex(Vertex_handle vh) {
	return vh->info();
      }


    };

}
