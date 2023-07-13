#include <fstream>

#define CGAL_EIGEN3_ENABLED
#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//#include <CGAL/Min_sphere_of_spheres_d.h>
//#include <CGAL/Min_sphere_of_spheres_d_traits_d.h>

#include <function_delaunay/Compute_meb_radius.h>
#include <function_delaunay/Miniball.hpp>

typedef CGAL::Epeck_d< CGAL::Dynamic_dimension_tag >  K;
typedef K::Point_d Point;

#define DIM 3

struct Compute_meb_radius_with_miniball {

  typedef std::vector<std::vector<double>>::const_iterator PointIterator; 
  typedef std::vector<double>::const_iterator CoordIterator;
  
  typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > 
    MB;

  
  template<typename InputIterator>
  double operator()(InputIterator begin, InputIterator end) {
    int d = begin->size();
    MB mb(d,begin,end);
    return sqrt(mb.squared_radius());
  }

};

/*
struct Compute_meb_radius_with_CGAL {

  // TODO: Fix this mess
 
  typedef CGAL::Min_sphere_of_spheres_d_traits_d<K,K::FT,DIM,CGAL::Tag_true> Min_sphere_traits;
  typedef CGAL::Min_sphere_of_spheres_d<Min_sphere_traits> Min_sphere;
  
  template<typename InputIterator>
  double operator()(InputIterator begin, InputIterator end) {

    typedef typename Min_sphere_traits::Sphere Sphere;
    
    std::vector<Sphere> spheres;
    
    for(auto it = begin;it!=end;it++) {
      spheres.push_back(std::make_pair(*it,0));
    }
    Min_sphere meb(spheres.begin(),spheres.end());

    //std::cout << "Number of support points: " << std::distance(meb.support_begin(),meb.support_end()) << std::endl;
    
    for(auto it = meb.support_begin();it!=meb.support_end();it++) {
      for(int i =0;i<DIM;i++) {
	std::cout << ((*it).first)[i] << " ";
      }
      std::cout << std::endl;
    }
    
    return CGAL::to_double(meb.radius());
  }

};

*/

double random_coordinate(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char** argv) {

  std::cout.precision(std::numeric_limits<double>::max_digits10);
  
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K2;
    
    typedef K2::Point_2 Point_2;
    Point_2 p1(0.17687209517549354, 0.074369245709092008);
    Point_2 p2(0.4430498310565249, 0.21541505782651485);
    Point_2 p3(0.73601323819533604, 0.41892177165435712);
    Point_2 c(-0.97409717988013256,2.5681309853871142);

    std::cout << CGAL::to_double(CGAL::squared_distance(p1,c)) << std::endl;
    std::cout << CGAL::to_double(CGAL::squared_distance(p2,c)) << std::endl;
    std::cout << CGAL::to_double(CGAL::squared_distance(p3,c)) << std::endl;

    
    K2::Triangle_2 tri = K2::Construct_triangle_2()(p1,p2,p3);
    std::cout << "On bounded: " << (K2::Bounded_side_2()(tri,c)!=CGAL::ON_UNBOUNDED_SIDE) << std::endl;
    
  }

  {
    typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kd;
    
    typedef Kd::Point_d Point_d;
    std::vector<Point_d> triangle;
    std::vector<double> p1_coors{0.17687209517549354, 0.074369245709092008};
    triangle.push_back(Point_d(p1_coors.begin(),p1_coors.end()));
    std::vector<double> p2_coors{0.4430498310565249, 0.21541505782651485};
    triangle.push_back(Point_d(p2_coors.begin(),p2_coors.end()));
    std::vector<double> p3_coors{0.73601323819533604, 0.41892177165435712};
    triangle.push_back(Point_d(p3_coors.begin(),p3_coors.end()));
    
    std::vector<double> c_coors{-0.97409717988013256,2.5681309853871142};
    Point_d c(c_coors.begin(),c_coors.end());
  
    std::cout << CGAL::to_double(Kd::Squared_distance_d()(triangle[0],c)) << std::endl;
    std::cout << CGAL::to_double(Kd::Squared_distance_d()(triangle[1],c)) << std::endl;
    std::cout << CGAL::to_double(Kd::Squared_distance_d()(triangle[2],c)) << std::endl;

    std::cout << "On bounded: " << (Kd::Contained_in_simplex_d()(triangle.begin(),triangle.end(),c)!=CGAL::ON_UNBOUNDED_SIDE) << std::endl;
    
  }


  int d;

  std::vector<Point> input_points;
  std::vector<std::vector<double>> input_points_raw;

  if(std::string(argv[1])=="-random") {


    long seed = time(NULL);
    std::cout << "Seed is " << seed << std::endl;
    srand(seed);
    

  
    int number = atoi(argv[2]);
    std::cout << "Creating " << number << " random points" << std::endl;
    
    d = atoi(argv[3]);
    std::cout << "In dimension " << d << std::endl;
    
    for(int i=0;i<number;i++) {
      std::vector<double> p;
      for(int j=0;j<d;j++) {
	p.push_back(random_coordinate(0,1));
      }
      input_points.push_back(Point(p.begin(),p.end()));
      input_points_raw.push_back(p);
    }
    
    // Dump to output file (for debugging)
    std::ofstream ofstr("last_instance.txt");
    ofstr.precision(std::numeric_limits<double>::max_digits10);

    
    for(Point p : input_points) {
      for(int j=0;j<d;j++) {
	ofstr << p[j] << " ";
      }
      ofstr << std::endl;
    }
    ofstr.close();
  } else {
    std::ifstream ifstr(argv[1]);
    
    std::string first_line;
    std::getline(ifstr,first_line);

    std::stringstream sstr(first_line);

    std::vector<double> coors;

    double next;
    
    sstr >> next;
    while(sstr.good()) {
      coors.push_back(next);
      sstr >> next;
    }
    d = coors.size();

    std::cout << "Dimension is " << d << std::endl;

    input_points.push_back(Point(coors.begin(),coors.end()));
    input_points_raw.push_back(coors);

    ifstr >> next;
    while(ifstr.good()) {
      coors.clear();
      for(int i=0;i<d;i++) {
	coors.push_back(next);
	ifstr >> next;
      }
      input_points.push_back(Point(coors.begin(),coors.end()));
      input_points_raw.push_back(coors);
    }

    std::cout << "Read " << input_points.size() << " points" << std::endl;
  
  }
  
  for(Point& p : input_points) {
    for(int i=0;i<d;i++) {
      std::cout << p[i] << " ";
    }
    std::cout << std::endl;
  }

  Compute_meb_radius<K> meb;

  double rad = meb(input_points_raw.begin(),input_points_raw.end());

  std::cout << "MEB radius       = " << rad << std::endl;

  Compute_meb_radius_with_miniball meb_miniball;

  double meb_alt = meb_miniball(input_points_raw.begin(),input_points_raw.end());

  std::cout << "MEB via Miniball = " << meb_alt << std::endl;
  
  return 0;

}
