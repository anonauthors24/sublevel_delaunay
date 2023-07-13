#include <scc/Scc.h>
#include <fstream>
#include <phat/boundary_matrix.h>
#include <phat/persistence_pairs.h>
#include <phat/compute_persistence_pairs.h>
#include <phat/algorithms/twist_reduction.h>

struct Sort_pairs {

  bool operator() (std::pair<double,double>& a,
		   std::pair<double,double>& b) {
    if(a.first < b.first) {
      return true;
    }
    if(a.first>b.first) {
      return false;
    }
    return a.second < b.second;
  }

};

int main(int argc, char** argv) {
  std::ifstream ifstr(argv[1]);
  double threshold;
  if(argc>2) {
    threshold = atof(argv[2]);
  } else {
    threshold = std::numeric_limits<double>::max();
  }
  std::cout << "threshold is " << threshold << std::endl;
  scc::Scc<> scc(ifstr);

  
  std::cout << "Number of parameters " << scc.number_of_parameters() << std::endl;
  int levels = scc.number_of_levels();
  std::cout << "Number of levels " << levels << std::endl;

  std::vector<std::map<int,int>> index_maps;
  index_maps.resize(levels);

  std::vector<double> filt_values;

  long number_of_cols=0;
  for(int i=levels;i>=1;i--) {
    for(int j=0;j<scc.number_of_generators(i);j++) {

      std::vector<double> grades;
      std::vector<std::pair<int,double> > bd;
      scc.next_column(i,std::back_inserter(grades),std::back_inserter(bd));
      if(grades[1]<=threshold) {
	index_maps[i-1][j]=number_of_cols;
	filt_values.push_back(grades[0]);
	number_of_cols++;
      }
    }
    scc.reset(i);
  }
  std::cout << "Total number of generators " << number_of_cols << std::endl;

  phat::boundary_matrix<> M;
  M.set_num_cols(number_of_cols);
  
  long count=0;
  
  for(int i=levels;i>=1;i--) {
    std::cout << "At level " << i << std::endl;
    for(int j=0;j<scc.number_of_generators(i);j++) {
      std::vector<double> grades;
      std::vector<std::pair<int,double> > bd;
      scc.next_column(i,std::back_inserter(grades),std::back_inserter(bd));
      if(grades[1]<=threshold) {
	phat::column col;
	for(auto p : bd) {
	  assert(i<levels);
	  assert(index_maps[i].find(p.first)!=index_maps[i].end());
	  col.push_back(index_maps[i][p.first]);
	}
	M.set_dim(count,levels-i);
	M.set_col(count,col);
	count++;
      }
    }
  }

  // define the object to hold the resulting persistence pairs
  phat::persistence_pairs pairs;
  
  // choose an algorithm (choice affects performance) and compute the persistence pair
  // (modifies boundary_matrix)
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, M );
    
  std::vector<std::pair<double,double> > ppairs;
  
  for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ ) {
    double birth = filt_values[pairs.get_pair(idx).first];
    double death = filt_values[pairs.get_pair(idx).second];
    if(birth!=death) {
      ppairs.push_back(std::make_pair(birth,death));
    }
  }
  std::cout << "There are " << ppairs.size() << " persistence pairs" << std::endl;
  std::sort(ppairs.begin(),ppairs.end(),Sort_pairs());
  std::ofstream ofstr("Flipper_pers_pairs.txt");

#ifdef HIGH_PRECISION_OUTPUT
  ofstr.precision(std::numeric_limits<double>::max_digits10);
#endif

  for( auto p : ppairs) {
      ofstr <<  p.first << " " << p.second << std::endl;
  }
  ofstr.close();
  
  M.save_ascii("sanity.phat.txt");

  ifstr.close();
  return 0;
}
