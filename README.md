**This repository has been archived.**

odeity
======

A simple C++ library for solving partial differential equations by the method of lines

-----

## Overview

ODEity is a simple C++ library that implements various time integrators for integrating systems of
ordinary differential equations. The focus of ODEity is on systems that are obtained after
semi-discretization of partial differential equations in space (e.g., by finite differences or
finite elements). It provides its own implementation of several explicit Runge-Kutta methods
with adaptive time step selection. Furthermore, it provides an interface consistent with that of the
explicit solver to higher order integration solvers provided by SUNDIALS library. Data output is
done using the NetCDF library. Some examples are also provided implementing the Allen-Cahn equation
and the Cahn-Hilliard equation in 2D.

## Building

### Prerequisities
  
* gcc and g++ compilers
* installed libraries: libmatheval, netcdf


### Compiling

ODEity uses CMake build system. In order to build ODEity (SUNDIALS will be built at the same time),
type the following commands:
 
```
$ cd odeity
$ mkdir build && cd build
$ cmake ..
$ make
```
