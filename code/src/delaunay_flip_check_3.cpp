#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K>  Triangulation;
typedef Triangulation::Point          Point;
typedef Triangulation::Cell_handle          Cell_handle;
typedef Triangulation::Facet Facet;

double random_coordinate(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char** argv)
{

  std::vector<Point> input_points;

  if(argc<2) {
      std::cout << "Not enough arguments" << std::endl;
      std::exit(1);
    }

  if(std::string(argv[1])=="--random") {
    srand(time(NULL));
    if(argc<3) {
      std::cout << "Not enough arguments" << std::endl;
      std::exit(1);
    }
    int n = atoi(argv[2]);
    for(int i=0;i<n;i++) {
      input_points.push_back(Point(random_coordinate(0,1),random_coordinate(0,1),random_coordinate(0,1)));
    }

    std::ofstream ofstr("last_instance.txt");
    ofstr.precision(std::numeric_limits<double>::max_digits10);
    long count=0;
    for(Point p : input_points) {
      ofstr << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " " << CGAL::to_double(p.z()) << " " << count++ << std::endl;
    }
  } else {
    
    std::ifstream in(argv[1]);
    std::istream_iterator<Point> begin(in);
    std::istream_iterator<Point> end;

    std::copy(begin,end,std::back_inserter(input_points));
    std::random_shuffle(input_points.begin(),input_points.end());
  }

  Triangulation T;
  
  long number_of_pentahedra=0;
  
  for(Point p : input_points) {

    std::vector<Cell_handle> cells;
    std::vector<Facet> facets;
    Cell_handle loc = T.locate(p);
    if(T.dimension()==3) {
      T.find_conflicts(p,loc,std::back_inserter(facets),std::back_inserter(cells));
      number_of_pentahedra+=cells.size();
    }
    T.insert(p,loc);
  }

  std::cout << "Number of pentahedra: " << number_of_pentahedra << std::endl;
  
  long n = T.number_of_vertices();
  long m = T.number_of_finite_edges();
  long f = T.number_of_finite_facets();
  long c = T.number_of_finite_cells();

  std::cout << "Number of vertices: " << n << std::endl;

  std::cout << "Number of edges: " << m << std::endl;

  std::cout << "Number of Facets: " << f << std::endl;

  std::cout << "Number of Cells: " << c << std::endl;

  std::cout << "-------------" << std::endl;

  std::cout << "Ratio: " << (double)number_of_pentahedra/(n+m+f+c) << std::endl;

  return 0;
}
