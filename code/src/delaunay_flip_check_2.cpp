#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;
typedef Triangulation::Point          Point;
typedef Triangulation::Face_handle          Face_handle;

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
      input_points.push_back(Point(random_coordinate(0,1),random_coordinate(0,1)));
    }

    std::ofstream ofstr("last_instance.txt");
    ofstr.precision(std::numeric_limits<double>::max_digits10);
    long count=0;
    for(Point p : input_points) {
      ofstr << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " " << count++ << std::endl;
    }
  } else {
    
    std::ifstream in(argv[1]);
    std::istream_iterator<Point> begin(in);
    std::istream_iterator<Point> end;

    std::copy(begin,end,std::back_inserter(input_points));
    std::random_shuffle(input_points.begin(),input_points.end());
  }

  Triangulation T;
  
  long number_of_tetrahedra=0;
  
  for(Point p : input_points) {
    std::vector<Face_handle> faces;
    T.get_conflicts(p,std::back_inserter(faces));
    number_of_tetrahedra+=faces.size();
    T.insert(p);
  }

  std::cout << "Number of tetrahedra: " << number_of_tetrahedra << std::endl;
  
  long n = T.number_of_vertices();
  long m = std::distance(T.finite_edges_begin(),T.finite_edges_end());
  long f = T.number_of_faces();

  std::cout << "Number of vertices: " << n << std::endl;

  std::cout << "Number of edges: " << m << std::endl;

  std::cout << "Number of Faces: " << f << std::endl;

  std::cout << "-------------" << std::endl;

  std::cout << "Ratio: " << (double)number_of_tetrahedra/(n+m+f) << std::endl;

  return 0;
}
