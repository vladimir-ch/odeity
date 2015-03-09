cmake_minimum_required( VERSION 2.6 )

project( odeity )

set( PACKAGE_BUGREPORT "chalupecky@math.kyushu-u.ac.jp" )
set( PACKAGE_NAME "Odeity" )
set( PACKAGE_VERSION "1.0.0" )
set( PACKAGE_STRING "Odeity 1.0.0" )
set( PACKAGE_TARNAME "odeity" )

# Prohibit in-source build
if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
  message( FATAL_ERROR "In-source build prohibited." )
endif()

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

find_package( Matheval REQUIRED ) 
if( Matheval_FOUND )
    include_directories( ${Matheval_INCLUDE_DIRS} )
    link_directories( ${Matheval_LIBRARY_DIRS} )
    set( extraLibs ${extraLibs} ${Matheval_LIBRARIES} )
endif()

find_package( NetCDF REQUIRED ) 
if( NetCDF_FOUND )
    include_directories( ${NetCDF_INCLUDE_DIRS} )
    link_directories( ${NetCDF_LIBRARY_DIRS} )
    # I don't know how to obtain the C++ library via pkg-config
    set( extraLibs ${extraLibs} ${NetCDF_LIBRARIES} netcdf_c++ )
endif()

add_subdirectory( vendor/sundials-2.4.0 )
add_subdirectory( src )

add_subdirectory( examples )

