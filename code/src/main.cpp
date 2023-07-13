#if FUNCTION_DELAUNAY_TIMERS
#include <function_delaunay/boost_timers.h>
#endif

#if WITH_MEMORY_PROFILE
#include <function_delaunay/mem_info.h>
#endif

#include <function_delaunay/function_delaunay_with_meb.h>
#include <function_delaunay/count_size_of_delaunay_triangulation.h>

#include <multi_chunk/multi_chunk.h>
#include <mpfree/mpfree.h>

#include <cassert>
#include <fstream>
#include <string>

bool is_non_negative_integer(std::string arg) {
  
  for(int i=0;i<arg.length();i++) {
    if(! std::isdigit(arg[i])) {
      return false;
    }
  }
  return true;
}

double random_coordinate(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void print_delaunay_comparison(std::vector<function_delaunay::Point_with_density> &input_points,
			       long total_number_of_simplices,
			       long number_of_chains_after_multi_chunk=-1) {

  std::vector<long> number_of_k_simplices;
 
#if FUNCTION_DELAUNAY_TIMERS
  delaunay_timer.start();
#endif

  long total_del_size = function_delaunay::count_size_of_delaunay_triangulation(input_points,number_of_k_simplices);

#if FUNCTION_DELAUNAY_TIMERS
  delaunay_timer.stop();
#endif

  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Delaunay:" << std::endl;
  std::cout << "Total size of Del: " << total_del_size << std::endl;
  for(int i=0;i<number_of_k_simplices.size();i++) {
    std::cout << "Simplices in dimension " << i << ": " << number_of_k_simplices[i] << std::endl;
  }
  std::cout << "Ratio: " << (double)total_number_of_simplices/total_del_size << std::endl;
  if(number_of_chains_after_multi_chunk>=0){
    std::cout << "Multi-chunk-Ratio: " << (double)number_of_chains_after_multi_chunk/total_del_size << std::endl;
  }
  std::cout << "--------------------------------------------" << std::endl;
}


void print_help_message(char* arg0) {

  std::cout << "Usage: " << arg0 << " [OPTIONS] <INPUT_FILE> <OUTPUT_FILE>\n\n";

  std::cout << "Computes the function-Bowyer-Watson filtration for a set of points. The input file consists of a point set in R^d. Each line represents a point by d+1 coordinates, where the last coordinate is interpreted as the function value. The result is a chain complex of free persistence modules in scc2020 format. \n\n";

  std::cout << "Options:\n";
  std::cout << "--multi-chunk            - post-processes chain complex with multi-chunk optimization\n";
  std::cout << "--minpres p              - only compute a minimal presentation on homology dimension p\n";
  std::cout << "--no-output              - the output data is produced, but not written into a file\n";
  std::cout << "--random n d             - instead of reading the input file, the program generates n points in the d-dimensional unit cube u.a.r with random densities\n";
  std::cout << "--height n d             - instead of reading the input file, the program generates n points in the d-dimensional unit cube u.a.r with the last coordinate as density\n";
  std::cout << "--random-file <filename> - the points generated in a random instance are by default dumped into the file \"last_instance.txt\". Another file name can be specified. Using \"-\", the random input is not stored\n";
  std::cout << "--random-seed s          - uses the specified seed. Otherwise, the current time is used\n";
  std::cout << "--only-complex-size      - Only computes the complex and compares its size with Delaunay triangulation. Ignores most other options.\n";
  std::cout << "--no-delaunay-compare    - Skips the comparison with the Delaunay triangulation in the end.\n";
  std::cout << "-h --help                - prints this message\n\n";
}

int main(int argc, char** argv )
{
#if FUNCTION_DELAUNAY_TIMERS
  initialize_timers();
  overall_timer.start();
#endif


#if FUNCTION_DELAUNAY_TIMERS
  initial_timer.start();
#endif

  // Parse command line arguments
  bool write_output_file=true;
  bool outputfile_read=false;
  std::string outfile;
  
  bool inputfile_read=false;
  std::string infile;

  bool random_instance=false;

  bool random_densities=true;
  int no_random_points=-1;
  int random_dim=-1;

  bool store_random_input=true;
  std::string random_file_name="last_instance.txt";

  int random_seed=time(NULL);

  bool use_multi_chunk=false;

  bool compute_min_pres=false;
  int min_pres_dim=-1;

  bool only_complex_size=false;

  bool compare_with_delaunay=true;

  bool print_help=false;

  int pos=1;
  
  while(pos<argc) {

    std::string arg(argv[pos]);
    
    if(arg=="-h" || arg=="--help") {
      print_help=true;
      break;
    } else if(arg=="--random") {
      random_instance=true;
      random_densities=true;
      inputfile_read=true; // a little hacked...
      if(pos+2<argc && is_non_negative_integer(std::string(argv[pos+1])) && is_non_negative_integer(std::string(argv[pos+2]))) {
	no_random_points=atoi(argv[++pos]);
	random_dim=atoi(argv[++pos]);
      } else {
	std::cout << "--random requires two extra (integer) arguments" << std::endl;
	print_help=true;
	break;
      }
    } else if(arg=="--height") {
      random_instance=true;
      random_densities=false;
      inputfile_read=true; // a little hacked...
      if(pos+2<argc && is_non_negative_integer(std::string(argv[pos+1])) && is_non_negative_integer(std::string(argv[pos+2]))) {
	no_random_points=atoi(argv[++pos]);
	random_dim=atoi(argv[++pos]);
      } else {
	std::cout << "--height requires two extra (integer) arguments" << std::endl;
	print_help=true;
	break;
      }
    } else if(arg=="--random-file") {
      if(pos+1<argc) {
	random_file_name=std::string(argv[++pos]);
	if(random_file_name=="-") {
	  store_random_input=false;
	}
      } else {
	std::cout << "--random-file expects an extra argument for the file name" << std::endl;
	print_help=true;
	break;
      }      
    } else if(arg=="--random-seed") {
      if(pos+1<argc && is_non_negative_integer(std::string(argv[pos+1]))) {
	random_seed=atoi(argv[++pos]);
      } else {
	std::cout << "--random-seed requires an extra (integer) argument" << std::endl;
	print_help=true;
	break;
      }
    } else if(arg=="--no-output") {
      write_output_file=false;
      outputfile_read=true;
    } else if(arg=="--only-complex-size") {
      only_complex_size=true;
    } else if(arg=="--no-delaunay-compare") {
      compare_with_delaunay=false;
    } else if(arg=="--multi-chunk") {
      use_multi_chunk=true;
    } else if(arg=="--minpres") {
      compute_min_pres=true;
      if(pos+1<argc && is_non_negative_integer(std::string(argv[pos+1]))) {
	min_pres_dim=atoi(argv[++pos]);
      } else {
	std::cout << "--minpres requires an extra (integer) argument" << std::endl;
	print_help=true;
	break;
      }
    } else {
      if(arg[0]=='-') {
	std::cout << "Unrecognized option: " << arg << std::endl;
	print_help=true;
	break;
      }
      if(! inputfile_read) {
	infile=arg;
	inputfile_read=true;
      } else if(! outputfile_read) {
	outfile=arg;
	outputfile_read=true;
      } else {
	std::cerr << "Ignoring argument " << arg << std::endl;
      }
    }
    pos++;
  }

  if(!print_help && ! random_instance && !inputfile_read) {
    std::cout << "Input file missing" << std::endl;
    print_help=true;
  }

  if(!print_help && write_output_file && !outputfile_read && !only_complex_size) {
   std::cout << "Output file missing" << std::endl;
    print_help=true;
  } 


  if(print_help) {
    print_help_message(argv[0]);
    std::exit(0);
  }

  typedef function_delaunay::Point_with_density Point_with_density;

  std::vector<Point_with_density> input_points;

  int d;

  if(random_instance) {

    std::cout << "Seed is " << random_seed << std::endl;
    srand(random_seed);
    
    std::cout << "Creating " << no_random_points << " random points" << std::endl;

    std::cout << "In dimension " << random_dim << std::endl;

    for(int i=0;i<no_random_points;i++) {
      std::vector<double> p;
      for(int j=0;j<random_dim;j++) {
      	p.push_back(random_coordinate(0,1));
      }
      double density;
      if(random_densities) {
	density = random_coordinate(0,1);
      } else {
	density = p.back();
      }
      input_points.push_back(Point_with_density(p.begin(),p.end(),density));
    }

    d=random_dim;

    if(store_random_input) {

      // Dump to output file (for debugging)
      std::ofstream ofstr(random_file_name.c_str());
    
#ifdef HIGH_PRECISION_OUTPUT
      ofstr.precision(std::numeric_limits<double>::max_digits10);
#endif
      
      for(Point_with_density p : input_points) {
	for(int j=0;j<random_dim;j++) {
	  ofstr << p.x[j] << " ";
	}
	ofstr << p.density << std::endl;
      }
      ofstr.close();
    }

  } else {
    
    std::ifstream ifstr(infile);
    
    std::string first_line;
    std::getline(ifstr,first_line);

    std::stringstream sstr(first_line);

    std::vector<double> coors;

    double next;
    
    sstr >> next;
    coors.push_back(next);
    while(sstr.good()) {
      sstr >> next;
      coors.push_back(next);
    }
    d = coors.size()-1;

    std::cout << "Dimension is " << d << std::endl;

    double density = coors.back();
    coors.pop_back();
    input_points.push_back(Point_with_density(coors.begin(),coors.end(),density));

    ifstr >> next;
    while(ifstr.good()) {
      coors.clear();
      for(int i=0;i<d;i++) {
	coors.push_back(next);
	ifstr >> next;
      }
      input_points.push_back(Point_with_density(coors.begin(),coors.end(),next));
      ifstr >> next;
    }

    std::cout << "Read " << input_points.size() << " points" << std::endl;
  }
  std::sort(input_points.begin(),input_points.end(),function_delaunay::Sort_by_density());

  /*
  for(auto p : input_points) {
    std::cout << "Coors: ";
    for(auto j : p.x) {
      std::cout << j << " ";
    }
    std::cout << "Density: " << p.density << std::endl;
  }
  */


#if FUNCTION_DELAUNAY_TIMERS
    initial_timer.stop();
#endif
    
#if WITH_MEMORY_PROFILE
    std::cout << "Memory after initial: " << mem_info() << std::endl;
#endif


  typedef mpp_utils::Graded_matrix<> Graded_matrix;

  std::vector<Graded_matrix> graded_matrices;

  long total_number_of_simplices = function_delaunay::function_delaunay_with_meb<Graded_matrix>(input_points,graded_matrices,only_complex_size);

#if FUNCTION_DELAUNAY_TIMERS
  std::cout << "Simplices per second: " << (double)total_number_of_simplices/(double(overall_timer.elapsed().wall)/std::pow(10,9)) << std::endl;
  std::cout << "Time per simplex (in microseconds): " << (double)(double(overall_timer.elapsed().wall)/std::pow(10,3))/total_number_of_simplices << std::endl;
#endif  

  if(only_complex_size) {
    print_delaunay_comparison(input_points,total_number_of_simplices);
#if FUNCTION_DELAUNAY_TIMERS
    overall_timer.stop();
    print_timers();
#endif

    return 0;
  }
  

#if FUNCTION_DELAUNAY_TIMERS
  multi_chunk_timer.start();
#endif


  long number_of_chains_after_multi_chunk=0;

  if(use_multi_chunk) {
    std::cout << "Apply multi_chunk" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<< Multi-chunk output start" << std::endl;
    multi_chunk::compress(graded_matrices);
    std::cout << ">>>>>>>>>>>>>>>>>>>> Multi-chunk output end" << std::endl;

    for(int i=0;i<graded_matrices.size();i++) {
      number_of_chains_after_multi_chunk += graded_matrices[i].get_num_cols();
    }
    std::cout << "Number of chains after mult-chunk: " << number_of_chains_after_multi_chunk << std::endl;

#if WITH_MEMORY_PROFILE
    std::cout << "Memory after multi-chunk: " << mem_info() << std::endl;
#endif

  } else {
    number_of_chains_after_multi_chunk=-1;
  }
  
#if FUNCTION_DELAUNAY_TIMERS
  multi_chunk_timer.stop();
#endif



  // To use mpfree, we have to convert the GrMat into the Graded_matrix type of mpfree.
  // This is slightly annoying and inefficient, mostly because the two types are
  // mostly the same (both derive from mpp_utils::Graded_matrix).

  // Also, needs to be defined even if not using mpfree because the later output
  // command has to know min_rep
  std::vector<Graded_matrix> min_rep;
  
  if(compute_min_pres) {

#if FUNCTION_DELAUNAY_TIMERS
    mpfree_timer.start();
#endif

    mpfree::verbose=true;
    std::cout << "using Mpfree" << std::endl;
    std::cout << "Homology dimension is " << min_pres_dim << std::endl;

    int index_of_min_pres_dim = d-min_pres_dim;

    Graded_matrix& first_mpfree_matrix = graded_matrices[index_of_min_pres_dim];
    Graded_matrix& second_mpfree_matrix = graded_matrices[index_of_min_pres_dim+1];
  
    std::cout << "1st matrix: " << first_mpfree_matrix.num_rows << " x " << first_mpfree_matrix.get_num_cols() << std::endl;
    std::cout << "2nd matrix: " << second_mpfree_matrix.num_rows << " x " << second_mpfree_matrix.get_num_cols() << std::endl;
    
    
    min_rep.resize(1);
    std::cout << "<<<<<<<<<<<<<<<<<<<< Mpfree output start" << std::endl;

    mpfree::compute_minimal_presentation(first_mpfree_matrix,
					 second_mpfree_matrix,
					 min_rep[0],false);
    
    std::cout << ">>>>>>>>>>>>>>>>>>>> Mpfree output end" << std::endl;
#if FUNCTION_DELAUNAY_TIMERS
    mpfree_timer.stop();
#endif
  
#if WITH_MEMORY_PROFILE
    std::cout << "Memory after mpfree: " << mem_info() << std::endl;
#endif

  }



#if FUNCTION_DELAUNAY_TIMERS
  file_output_timer.start();
#endif

  if(write_output_file) {

    std::ofstream ofstr(outfile);
#ifdef HIGH_PRECISION_OUTPUT
    ofstr.precision(std::numeric_limits<double>::max_digits10);
#endif
    std::cout << "Writing to output" << std::endl;
    if(compute_min_pres) {
      mpp_utils::print_in_scc_format(min_rep,ofstr,true);
    } else {
      mpp_utils::print_in_scc_format(graded_matrices,ofstr);
    }
    ofstr.close();

#if WITH_MEMORY_PROFILE
    std::cout << "Memory after file output: " << mem_info() << std::endl;
#endif

  }

#if FUNCTION_DELAUNAY_TIMERS
  file_output_timer.stop();
#endif


  if(compare_with_delaunay) {
    print_delaunay_comparison(input_points,total_number_of_simplices,number_of_chains_after_multi_chunk);
  }

#if WITH_MEMORY_PROFILE
  std::cout << "Memory in the end: " << mem_info() << std::endl;
#endif


#if FUNCTION_DELAUNAY_TIMERS
  overall_timer.stop();
  print_timers();
#endif
  



}
