# Include source files: Common files to all kind of compilations
file(GLOB SRC_FILES *.cpp cpu/*.cpp analytic_integral/*.cpp)

# Include sorce files according to the usage or not of external libraries
if (EXTERNAL)
   message(FATAL_ERROR "Not implemented yet")
else()
   list(APPEND SRC_FILES cioverlap/cioverlap_fake.cpp excited/fake_subs.cpp
                         excited/saverho.cpp libint/g2g_libint_fake.cpp)
endif()

# OpenMP library
find_package(OpenMP)
if (NOT OPENMP_FOUND)
  message(FATAL_ERROR "OpenMP library was not found")
endif()

# Add Warning options
if (WARNING)
  add_compile_options(-Wall -Wextra -Wshadow -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare)
endif()

# Compile g2g library
add_compile_options(-O3)
add_library(g2g SHARED ${SRC_FILES})
target_link_libraries(g2g PUBLIC OpenMP::OpenMP_CXX)

# CPU compilations options
if (CPU)
  target_compile_definitions(g2g PRIVATE "CPU_KERNELS=1")
endif()

# External Libraries
if (EXTERNAL)
   message(STATUS "Using libxc and libint libraries")

   # Setting the LD_LIBRARY_PATH
   string(REPLACE ":" ";" RUNTIME_PATH "$ENV{LD_LIBRARY_PATH}")

   # libxc
   find_library(mylibxc NAME xc HINTS ${RUNTIME_PATH} DOC "LibXC Library")
   if (${mylibxc} STREQUAL mylibxc-NOTFOUND)
      message(FATAL_ERROR "LibXC Library was not found")
   else()
      message(STATUS "LibXC was found as ${mylibxc}")
   endif()

   # libint
   find_library(mylibint NAME int2 HINTS ${RUNTIME_PATH} DOC "LibINT Library")
   if (${mylibint} STREQUAL mylibint-NOTFOUND)
      message(FATAL_ERROR "LibINT Library was not found")
   else()
      message(STATUS "LibINT was found as ${mylibint}")
   endif()

   target_compile_definitions(g2g PRIVATE "USE_LIBXC=1" "LIBXC_CPU=1")
   target_link_libraries(g2g PUBLIC ${mylibxc} ${mylibint})
   target_include_directories(g2g PRIVATE $ENV{LIBINT_HOME}/include $ENV{LIBXC_HOME_CPU}/include)
else()
   message(STATUS "Don't use external libraries")
   target_compile_definitions(g2g PRIVATE "USE_LIBXC=0" "LIBXC_CPU=1")
endif()

# precision
if (DOUBLE)
   message(STATUS "Compiling LIO with double precision")
   target_compile_definitions(g2g PRIVATE "FULL_DOUBLE=1")
else()
   message(STATUS "Compiling LIO with simple precision")
   target_compile_definitions(g2g PRIVATE "FULL_DOUBLE=0")
endif()

# Include headers directories
target_include_directories(g2g PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} cpu analytic_integral cuda)
