This folder contains the supplemantary experimental files for the
submission

``Delaunay Bifiltrations of Functions on Point Clouds''

The program to compute the bifiltrations is in the folder ``code''.
Most external libraries are contained, but the program requires
CGAL, GUDHI, EIGEN, and BOOST installed on the machine. The paths to these
libraries can be given within the file code/CMakeLists.txt. After that,
the main program can be compiled via

cd code
cmake .
make main

The program ``main'' can be tested, via

./main --random 1000 3 --no-output

to compute the bifiltration for 1000 random points in 3 dimensions and random
function values. More options are available via

./main --help

The benchmark files are in the ``benchmarks'' folder. Consult the NOTES.txt in there for further information.