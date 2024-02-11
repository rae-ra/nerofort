cmake_minimum_required(VERSION 3.14.0)

project(nerofort
	LANGUAGES Fortran
	)


include(${PROJECT_SOURCE_DIR}/cmake/nerofort.cmake)


# Config
#add_subdirectory(config)


# Set compiler options
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-fimplicit-none)
  add_compile_options(-ffree-line-length-132)
  add_compile_options(-Wall)
  add_compile_options(-Wextra)
  add_compile_options(-Wimplicit-procedure)
  add_compile_options(-Wconversion-extra)
  # -pedantic-errors triggers a false positive for optional arguments of elemental functions,
  # see test_optval and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95446
  add_compile_options(-pedantic-errors)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8.0)
    add_compile_options(-std=f2018)
  else()
    add_compile_options(-std=f2008ts)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  add_compile_options(-warn declarations,general,usage,interfaces,unused)
  if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 20.2.1.20200827)
    add_compile_options(-standard-semantics)
  endif()
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 18.0)
    add_compile_options(-stand f15)
  else()
    add_compile_options(-stand f18)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  add_compile_options(-Mdclchk)
endif()

# Add src
add_subdirectory(src
)
# Set include directories
target_include_directories(nerofort PRIVATE src)

# Set the output binary directory
set_target_properties(nerofort PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

find_program(FYPP fypp)
if(NOT FYPP)
  message(FATAL_ERROR "Preprocessor fypp not found!")
endif()

