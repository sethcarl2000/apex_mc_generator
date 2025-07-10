#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16)
project(simc_hrsl
  DESCRIPTION "Test of interface between simc 'mc_hrsl' subroutine with cpp"
  LANGUAGES CXX Fortran C
)

#----------------------------------------------------------------------------
# Include packages we will need for the project
#


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include_cc/*.hh
#  ${PROJECT_SOURCE_DIR}/include_cc/*.c
)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/*.f
  ${PROJECT_SOURCE_DIR}/*.cc
#  ${PROJECT_SOURCE_DIR}/include_cc/*.c
)




#----------------------------------------------------------------------------
# Add the executable, use our local headers, and link it to the Geant4 libraries
#
add_executable(hrsl_exec hrsl_exec.cc ${sources} ${headers})
#target_include_directories(hrsl PRIVATE include)
#target_link_libraries( PRIVATE ${Geant4_LIBRARIES})

