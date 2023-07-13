# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( codedd )


if (NOT CMAKE_VERBOSE_MAKEFILE)
   set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (NOT CGAL_DIR)
    set(CGAL_DIR "$ENV{HOME}/CGAL-5.5.1")
endif()

if (NOT EIGEN_INCLUDE_DIR)
    set(EIGEN_INCLUDE_DIR "/usr/include/eigen3")
endif()

if (NOT GUDHI_INCLUDE_DIR)
    set(GUDHI_INCLUDE_DIR "$ENV{HOME}/software/gudhi/include")
endif()

if (NOT SCC_INCLUDE_DIR)
    set(SCC_INCLUDE_DIR "scc_mod/include")
endif()

if (NOT MPP_UTILS_INCLUDE_DIR)
    set(MPP_UTILS_INCLUDE_DIR "mpp_utils_mod/include")
endif()


if (NOT PHAT_INCLUDE_DIR)
    set(PHAT_INCLUDE_DIR "phat/include")
endif()


if (NOT MULTI_CHUNK_INCLUDE_DIR)
    set(MULTI_CHUNK_INCLUDE_DIR "multi_chunk_mod/include")
endif()

if (NOT MPFREE_INCLUDE_DIR)
    set(MPFREE_INCLUDE_DIR "mpfree_mod/include")
endif()


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

find_package(Boost OPTIONAL_COMPONENTS timer)

if(Boost_TIMER_FOUND)
    add_compile_options(-DFUNCTION_DELAUNAY_TIMERS=1)
    link_libraries(Boost::timer)
endif()



if (NOT DEFINED WITH_MEMORY_PROFILE)
   set(WITH_MEMORY_PROFILE ON)
endif()

if (NOT DEFINED FUNCTION_DELAUNAY_ALWAYS_USE_MINIBALL)
   set(FUNCTION_DELAUNAY_ALWAYS_USE_MINIBALL false)
endif()

if (NOT DEFINED FUNCTION_DELAUNAY_SMART_MEB_TRAVERSAL)
   set(FUNCTION_DELAUNAY_SMART_MEB_TRAVERSAL true)
endif()


if(WITH_MEMORY_PROFILE)
  add_compile_options(-DWITH_MEMORY_PROFILE=1)
endif()			    

if(FUNCTION_DELAUNAY_ALWAYS_USE_MINIBALL)
  add_compile_options(-DFUNCTION_DELAUNAY_ALWAYS_USE_MINIBALL=1)
endif()	

if(FUNCTION_DELAUNAY_SMART_MEB_TRAVERSAL)
  add_compile_options(-DFUNCTION_DELAUNAY_SMART_MEB_TRAVERSAL=1)
endif()	


add_compile_options(-DHIGH_PRECISION_OUTPUT=1)

# include for local directory

include_directories(include ${SCC_INCLUDE_DIR} ${PHAT_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${GUDHI_INCLUDE_DIR} ${MULTI_CHUNK_INCLUDE_DIR} ${MPP_UTILS_INCLUDE_DIR} ${MPFREE_INCLUDE_DIR})

# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################


create_single_source_cgal_program( "src/main.cpp" )
create_single_source_cgal_program( "src/persistence_up_to_density.cpp" )
create_single_source_cgal_program( "src/delaunay_flip_check_2.cpp" )
create_single_source_cgal_program( "src/delaunay_flip_check_3.cpp" )
add_executable(filtration_domination "src/filtration_domination.cpp")




