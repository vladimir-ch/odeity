# use pkg-config to get the directories and then use these values
# in the FIND_PATH() and FIND_LIBRARY() calls
find_package( PkgConfig )
if( PKG_CONFIG_FOUND )
    pkg_check_modules( NetCDF netcdf QUIET )
endif()

set( NetCDF_DEFINITIONS ${NetCDF_CFLAGS_OTHER} )

find_path( NetCDF_INCLUDE_DIRS
    NAMES netcdf.h
    HINTS
    ${NetCDF_INCLUDEDIR}
    ${NetCDF_INCLUDE_DIRS}
    )

find_library( NetCDF_LIBRARIES
    NAMES netcdf netcdf_c++
    HINTS
    ${NetCDF_LIBDIR}
    ${NetCDF_LIBRARY_DIRS}
    )

# handle the QUIETLY and REQUIRED arguments and set NetCDF_FOUND to TRUE if 
# all listed variables are TRUE
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( NetCDF DEFAULT_MSG
    NetCDF_LIBRARIES NetCDF_INCLUDE_DIRS)

mark_as_advanced(NetCDF_INCLUDE_DIRS NetCDF_LIBRARIES)
