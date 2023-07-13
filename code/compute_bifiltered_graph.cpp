#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

struct Point {
    std::vector<double> coordinates;
    double value;
};

double euclideanDistance(const std::vector<double>& point1, const std::vector<double>& point2) {
    double distance = 0.0;
    size_t dimensions = point1.size();
    for (size_t k = 0; k < dimensions; ++k) {
        double diff = point1[k] - point2[k];
        distance += (diff * diff);
    }
    return std::sqrt(distance);
}

void compute_bifilterd_graph(const std::vector<Point>& points) {
    size_t n = points.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double distance = euclideanDistance(points[i].coordinates, points[j].coordinates);
            double maxFunc = std::max(points[i].value, points[j].value);
            std::cout << i << " " << j << " " << distance << " " << maxFunc << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: ./program <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Error opening file: " << filename << std::endl;
        return 1;
    }

    std::vector<Point> points;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        Point point;
        double coordinate;
        while (iss >> coordinate) {
            if (!iss.eof()) {
                point.coordinates.push_back(coordinate);
            } else {
                point.value = coordinate;
            }
        }
        points.push_back(point);
    }

    file.close();

    compute_bifilterd_graph(points);

    return 0;
}

