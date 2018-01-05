cmake_minimum_required(VERSION 2.8)

project(destiny)

# https://github.com/rstudio/rstudio/raw/v1.1.411/cmake/modules/FindLibR.cmake
find_package(LibR REQUIRED)

function(get_r_include varname package)
	execute_process(
		COMMAND ${LIBR_EXECUTABLE} --slave -e "stopifnot(require(${package})); cat(system.file('include', package = '${package}'))"
		OUTPUT_VARIABLE output
	)
	set(${varname} ${output} PARENT_SCOPE)
endfunction(get_r_include)

get_r_include(LIBRCPP_INCLUDE_DIRS      Rcpp)
get_r_include(LIBRCPPEIGEN_INCLUDE_DIRS RcppEigen)

include_directories(BEFORE ${LIBR_INCLUDE_DIRS} ${LIBRCPP_INCLUDE_DIRS} ${LIBRCPPEIGEN_INCLUDE_DIRS})

set(SOURCES
	src/Cover_Tree.h
	src/utils.cpp
	src/utils.h
	src/censoring.cpp
	src/knn.cpp
	src/RcppExports.cpp
)

add_library(destiny SHARED ${SOURCES})
target_link_libraries(destiny ${LIBR_LIBRARIES})

add_executable(destiny-test tests/test_main.cpp)
target_link_libraries(destiny-test destiny)

add_custom_target(destiny-package ALL
	#COMMAND find ${CMAKE_SOURCE_DIR} -name "*.o" -exec rm "{}" "\;"
	#COMMAND find ${CMAKE_SOURCE_DIR} -name "*.so" -exec rm "{}" "\;"
	COMMAND ${LIBR_EXECUTABLE} --slave -e "\"stopifnot(require(roxygen2));roxygenize('${CMAKE_SOURCE_DIR}',roclets=c('rd','collate','namespace'))\""
	COMMAND ${LIBR_EXECUTABLE} CMD INSTALL "${CMAKE_SOURCE_DIR}"
)
