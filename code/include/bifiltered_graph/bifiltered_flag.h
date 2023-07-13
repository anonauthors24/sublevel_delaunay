#include <mpp_utils/create_graded_matrices_from_simplex_tree.h>
#include <gudhi/Simplex_tree.h>

#include <bifiltered_graph/Grade_map.h>
#include <function_delaunay/boost_timers.h>

namespace bifiltered_graph {

    typedef Gudhi::Simplex_tree<> Simplex_tree; 
    void build_flag(std::vector<BifilteredEdge> &edges, Simplex_tree& simplex_tree) {
                       //  std::vector<weightedVertex> &input_points) {

        // Add 1-simplices from edges with radius filtration value
        for (auto& edge : edges) {
            simplex_tree.insert_simplex_and_subfaces(edge.vertices , edge.bifil[0]);
        }
        // Set all vertices to have 0 filtration value
        for (auto simplex : simplex_tree.complex_simplex_range()) {
            if (simplex_tree.dimension(simplex) == 0) {
                simplex_tree.assign_filtration(simplex, 0.0);
             }
        }
        std::cout << "Dimension of complex before expansion " << simplex_tree.dimension() << std::endl ;
        // Expand to flag complex
        simplex_tree.expansion(3);
        std::cout << "Dimension of complex " << simplex_tree.dimension() << std::endl ;

        /*
        for (auto simplex : simplex_tree.complex_simplex_range()) {
            if (simplex_tree.dimension(simplex) == 2) {
                std::cout << "filtration value of 2-simplex " << simplex_tree.filtration(simplex) << std::endl;
             }
        } */
}
    template<typename GrMat>
    void flag_complex_with_bifiltration(std::vector<BifilteredEdge> &edges, 
                                std::vector<weightedVertex> &input_points, 
                                std::vector<GrMat> &graded_matrices){
        typedef Gudhi::Simplex_tree<> Simplex_tree; 
        Simplex_tree simplex_tree;
        complex_timer.start();
        build_flag(edges, simplex_tree);
        complex_timer.stop();
        

        graded_matrices_timer.start();
        Grade_map<Simplex_tree> grade_map(simplex_tree,input_points);
        std::cout << "Dimension of complex: " << simplex_tree.dimension() << std::endl;
        std::cout << "Simplex tree has " << simplex_tree.num_vertices() << " vertices and " << simplex_tree.num_simplices() << " simplices" << std::endl;
        mpp_utils::create_graded_matrices_from_simplex_tree(simplex_tree,grade_map,graded_matrices,false);
        graded_matrices_timer.stop();

        std::cout << "Graded matrices computed from simplex tree" << std::endl;
        int simplex_tree_dimension=simplex_tree.dimension();
    
        // The simplex tree is no longer needed at this point
        // Sort of a hack, but it seems to work
        simplex_tree.prune_above_dimension(-1);
  
        for(int i=0;i<=simplex_tree_dimension;i++) {
            std::cout << "Simplices in dimension " << i << ": " << graded_matrices[simplex_tree_dimension-i].get_num_cols() << std::endl;
        }
    }
}
