cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(ecgrid LANGUAGES Fortran)

option(V_POLE "Stagger V on the pole" OFF)

if (V_POLE) 
  add_definitions(-D V_POLE)
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Ofast -fopenmp")
endif ()

find_package(MPI)
include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
  	pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
    	include_directories(${NETCDF_INCLUDEDIR})
    	link_directories(${NETCDF_LIBDIR})
    else ()
    	message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
  	message(FATAL_ERROR "Unable to find Pkg-config library!")
  endif ()
endif ()

find_package(BLAS)
find_package(LAPACK)
if (LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()


add_subdirectory(lib/container)
add_subdirectory(lib/datetime)
add_subdirectory(lib/string)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)
set(source_files
  src/utils/const_mod.F90
  src/utils/namelist_mod.F90
  src/utils/time_mod.F90
  src/utils/allocator_mod.F90
  src/utils/history_mod.F90
  src/utils/parallel_mod.F90
  src/utils/debug_mod.F90
  src/core/mesh_mod.F90
  src/core/state_mod.F90
  src/core/static_mod.F90
  src/core/tend_mod.F90
  src/core/operators_mod.F90
  src/core/pv_mod.F90
  src/core/ke_mod.F90
  src/core/ecgrid_mod.F90
  src/tests/swm/mountain_zonal_flow_test_mod.F90
  src/tests/swm/rossby_haurwitz_wave_test_mod.F90
  src/tests/swm/cross_pole_flow_test_mod.F90
)

add_library(ecgrid ${source_files})
target_link_libraries(ecgrid netcdf netcdff fortran_container fortran_datetime fiona flogger ${MPI_Fortran_LIBRARIES} ${lapackblas_libraries} lapack blas)

add_executable(ecgrid_swm_driver.exe src/drivers/ecgrid_swm_driver.F90)
target_link_libraries(ecgrid_swm_driver.exe ecgrid)
