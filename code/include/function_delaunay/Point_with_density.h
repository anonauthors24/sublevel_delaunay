#pragma once

namespace function_delaunay {

  struct Point_with_density {
    std::vector<double> x;
    double density;
    template<typename Iterator>
    Point_with_density(Iterator begin, Iterator end,double density)
      : density(density) {
      std::copy(begin,end,std::back_inserter(x));
    }
    int dimension() const {
      return x.size();
    }
  };
  
  struct Sort_by_density {
    bool operator() (Point_with_density& a, Point_with_density& b) {
      return a.density<b.density;
    }
  };
}
