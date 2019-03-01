cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(tmcore LANGUAGES Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

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
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

set(source_files
  src/CRK3.f90
  src/CRK4.f90
  src/Gill4.f90
  src/L.f90
  src/PC3.f90
  src/diagnostics_tools.f90
  src/integration.f90
  src/mesh_info.f90
  src/module_array.f90
  src/module_para.f90
  src/output_netCDF.f90
  src/shallow_water_waves_test_mod.F90
  src/test_cases.f90
)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/.mods")

add_executable(tmcore_test.exe src/MPAS_SW.f90 ${source_files})
target_link_libraries(tmcore_test.exe netcdf netcdff)
