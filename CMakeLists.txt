# Copyright (c) 2010 Vladimir Chalupecky
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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

