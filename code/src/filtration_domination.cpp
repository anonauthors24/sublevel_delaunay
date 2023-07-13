#include <iostream>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

//#include <Python.h>

#include <gudhi/Simplex_tree.h>

#include <function_delaunay/boost_timers.h>
#include <function_delaunay/mem_info.h>


#include <multi_chunk/multi_chunk.h>
#include <mpfree/mpfree.h>
#include <bifiltered_graph/bifiltered_flag.h>

int main(int argc, char* argv[]) {
    initialize_timers();
    overall_timer.start();

    if (argc < 3) {
        std::cout << "Usage: ./main reduced_bifiltered_edges.txt filtered_vertices.txt" << std::endl;
        return 1;
    }
    std::cout << "Computing the minimal presentation of homology dimension 1 induced by a bifiltered graph" << std::endl;
    initial_timer.start();
    std::string filename1 = argv[1];
    std::string filename2 = argv[2];

    std::ifstream inputFile1(filename1);
    if (!inputFile1) {
        std::cout << "Error opening file: " << filename1 << std::endl;
        return 1;
    }

    std::ifstream inputFile2(filename2);
    if (!inputFile2) {
        std::cout << "Error opening file: " << filename2 << std::endl;
        return 1;
    }

    std::vector<bifiltered_graph::BifilteredEdge> edges;
    std::string line;
    while (std::getline(inputFile1, line)) {
        std::stringstream lineStream(line);
        std::vector<int> vertices(2);
        std::vector<double> bifil(2);
        lineStream >> vertices[0] >> vertices[1] >> bifil[0] >> bifil[1];
        bifiltered_graph::BifilteredEdge edge(vertices, bifil);
        edges.push_back(edge);
        vertices.clear();
        bifil.clear();
    }
    
    int d = 0; 
        if (std::getline(inputFile2, line)) {
		std::stringstream lineStream(line);
		std::string value;
		while (lineStream >> value) {
		d = d+1;
		}
	 	d = d-1;
	}
    

    std::vector<bifiltered_graph::weightedVertex> vertices;
    int lineNumber = 1;
    while (std::getline(inputFile2, line)) {
        std::stringstream lineStream(line);
        double lastValue = 0.0;
        std::string value;
        while (lineStream >> value) {
            lastValue = std::stod(value);
        }
        bifiltered_graph::weightedVertex vertex;
        vertex.v = lineNumber;
        vertex.weight = lastValue;
        vertices.push_back(vertex);
        vertex = bifiltered_graph::weightedVertex();
        lineNumber++;
    }
    /*
    std::vector<bifiltered_graph::Point_with_density> vertices;
    int lineNumber = 1;
    while (std::getline(inputFile2, line)) {
        std::stringstream lineStream(line);
        std::vector<double> values;
        double lastValue = 0.0;
        std::string value;
        while (lineStream >> value) {
            lastValue = std::stod(value);
            values.push_back(lastValue);
        }
        bifiltered_graph::Point_with_density vertex(values.begin(), values.end(), lastValue);
        vertices.push_back(vertex);
        lineNumber++;
    }
    */
    
    inputFile1.close();
    inputFile2.close();
    line.clear();

    // Clear and release memory for other variables
    std::string().swap(filename1);
    std::string().swap(filename2);
    std::string().swap(line);
    initial_timer.stop();
    std::cout << "Memory after initial: " << mem_info() << std::endl;
    std::cout << "Number of edges collected " << edges.size() << std::endl;

    typedef mpp_utils::Graded_matrix<> Graded_matrix;

    std::vector<Graded_matrix> graded_matrices;

    bifiltered_graph::flag_complex_with_bifiltration<Graded_matrix>(edges, vertices, graded_matrices);
    
    // no need for edges and vertices
    std::vector<bifiltered_graph::BifilteredEdge>().swap(edges);
    std::vector<bifiltered_graph::weightedVertex>().swap(vertices);


    /*
    multi_chunk_timer.start();
    long number_of_chains_after_multi_chunk=0;
    std::cout << "<<<<<<<<<<<<<<<<<<<< Multi-chunk output start" << std::endl;
    multi_chunk::compress(graded_matrices);
    std::cout << "<<<<<<<<<<<<<<<<<<<< Multi_chunk output stop" << std::endl;
    multi_chunk_timer.stop();
    std::cout << "Memory after multi_chunk: " << mem_info() << std::endl;
    for(int i=0;i<graded_matrices.size();i++) {
      number_of_chains_after_multi_chunk += graded_matrices[i].get_num_cols();
    }
    std::cout << "Number of chains after mult-chunk: " << number_of_chains_after_multi_chunk << std::endl; */
    
// mpfree starts here
    mpfree_timer.start();
    std::vector<Graded_matrix> min_rep;    
    mpfree::verbose=true;
    std::cout << "Mpfree starts" << std::endl;
    int min_pres_dim = 1;
    std::cout << "Homology dimension is " << min_pres_dim << std::endl;
    
    /*
    for (auto& matrix : graded_matrices) {
        std::cout << "Graded Matrix:" << std::endl;
        matrix.print(); // Use the print function of the Graded_matrix class or customize the output here
        std::cout << std::endl; // Add a newline between matrices
    }*/
    int index_min_pres_dim = d-min_pres_dim;
    Graded_matrix& first_mpfree_matrix = graded_matrices[index_min_pres_dim];
    Graded_matrix& second_mpfree_matrix = graded_matrices[index_min_pres_dim + 1];
    


    min_rep.resize(1);

    std::cout << "<<<<<<<<<<<<<<<<<<<< Mpfree output start" << std::endl;

    mpfree::compute_minimal_presentation(first_mpfree_matrix, second_mpfree_matrix,
					 min_rep[0],false);
    std::cout << ">>>>>>>>>>>>>>>>>>>> Mpfree output end" << std::endl;
    mpfree_timer.stop();
    std::cout << "Memory after mpfree: " << mem_info() << std::endl;

    overall_timer.stop();
    print_timers();

    return 0;
}



