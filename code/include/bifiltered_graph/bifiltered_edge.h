namespace bifiltered_graph {

    struct BifilteredEdge {
        std::vector<int> vertices;    // Combined vertices vector
        std::vector<double> bifil;    // Bifil vector
        
    BifilteredEdge(std::vector<int> vertex_list, const std::vector<double>& bifil_values)
        : vertices(vertex_list), bifil(bifil_values) {
    }
};
    
    struct weightedVertex {
        int v;
        double weight;
    };

    struct Sort_by_weight {
        bool operator () (weightedVertex& v, weightedVertex& w) {
            return v.weight < w.weight;
        }
    };
}
