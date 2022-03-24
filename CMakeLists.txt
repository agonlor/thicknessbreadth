if(APPLE)
	set(CMAKE_C_COMPILER "/usr/bin/gcc")
	set(CMAKE_CXX_COMPILER "/usr/local/bin/clang")
endif()
	

PROJECT(tb)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(CheckCXXSourceCompiles)

include_directories (include)

FIND_PACKAGE(OpenMP)
set (CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
CHECK_CXX_SOURCE_COMPILES("
#include <omp.h>
int main() {
#if (_OPENMP >= 200805 || _MSC_VER >= 1500)
  return 0;
#else
  breaks_on_purpose
#endif
}
" OPENMP_VERSION)

if(OPENMP_VERSION)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
message("

OpenMP 3.0 not supported by the compiler (${CMAKE_CXX_COMPILER})!
To get optimal performance for the \"chunk\" algorithm, use a C++ compiler with OpenMP 3.0 support (e.g., GCC >=4.4).
To use a different compiler, pass it to cmake in the variable CMAKE_CXX_COMPILER:
 cmake . -DCMAKE_CXX_COMPILER=g++-4.7

")
endif()

SET(SRC
	src/main.cpp
        src/binobject.h  src/binobject.cpp
	)

ADD_DEFINITIONS("-Wall -std=c++11 -O3" -g)

if(APPLE)
	# For Boost
	INCLUDE_DIRECTORIES(/opt/local/include)
	# For stdlib
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lc++ -lc++abi")
endif()

## Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
##
	
ADD_EXECUTABLE(tb ${SRC})

## Required in DGtal
TARGET_LINK_LIBRARIES(tb ${DGTAL_LIBRARIES})
##
