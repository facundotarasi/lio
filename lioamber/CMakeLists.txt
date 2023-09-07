enable_language(Fortran)

############################ SOURCE FILES ###############################
file(GLOB SRC_FILES *.f90 liomods/*.f90 liosubs/liosubs.f90 fileio/fileio.f90
                    liosubs_math/liosubs_math.f90 mathsubs/mathsubs.f90
                    typedef_cumat/typedef_cumat.f90 linear_algebra/linear_algebra.f90
                    typedef_operator/typedef_operator.f90 faint_cpu/*.f90
                    typedef_sop/typedef_sop.f90 cdft/cdft_data.f90 cdft/cdft_subs.f90
                    properties/properties_data.f90 properties/properties.f90
                    dftd3/dftd3.f90 ehrensubs/ehrensubs.f90 
                    excitedsubs/excitedsubs.f90 fstshsubs/fstshsubs.f90
                    lennard_jones/lj_switch_data.f90 lennard_jones/lj_switch.f90
                    lennard_jones/ljs_mm_interface.f90 
                    converger_subs/converger_subs.f90 debug_tools/debug_tools.f90)
######################### END SOURCE FILES ###############################

############################ LIBRARIES ###############################
# lapack
find_library(mylapack lapack DOC "Lapack library")
if (${mylapack} STREQUAL mylapack-NOTFOUND)
    message(FATAL_ERROR "no lapack found ${mylapack}")
else()
    message(STATUS "Lapack was found as ${mylapack}")
endif()

# blas
find_library(myblas blas DOC "Blas library")
if (${myblas} STREQUAL myblas-NOTFOUND)
    message(FATAL_ERROR "no blas found ${myblas}")
else()
    message(STATUS "Blas was found as ${myblas}")
endif()
######################### END LIBRARIES ###############################

if (WARNING)
   add_compile_options(-Wall -Wextra)
endif()

add_compile_options(-cpp -O3)
add_library(lioamber SHARED ${SRC_FILES})
target_link_libraries(lioamber PUBLIC g2g ${mylapack} ${myblas})

# Add common flags
target_compile_definitions(lioamber PRIVATE "G2G" "TD_SIMPLE" "CPU_DOUBLE" "pack")

