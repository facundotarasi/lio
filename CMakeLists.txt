cmake_minimum_required(VERSION 3.10)

project(LIO VERSION 1.0)

# CHECK OPTIONS COMPILATIONS
include(compileOptions.cmake)

# ADD LIBRARIES
add_subdirectory(g2g)
add_subdirectory(lioamber)

# ENABLE FORTRAN
enable_language(Fortran)
add_compile_options(-Wall -O3)

# LIOSOLO
add_executable(liosolo.x liosolo/liosolo.f90)
target_link_libraries(liosolo.x PUBLIC g2g lioamber)
target_include_directories(liosolo.x PUBLIC ${PROJECT_BINARY_DIR}/lioamber)

# LIOMD
add_executable(liomd.x liosolo/liomd.f90)
target_link_libraries(liomd.x PUBLIC g2g lioamber)
target_include_directories(liomd.x PUBLIC ${PROJECT_BINARY_DIR}/lioamber)

# INSTALL FOLDERS
install(TARGETS liosolo.x liomd.x DESTINATION bin)
install(TARGETS g2g lioamber LIBRARY DESTINATION lib)

