#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(wavy)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)  # => added 20/06/01
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
#include_directories(${PROJECT_SOURCE_DIR}/include 
#                    ${PROJECT_SOURCE_DIR}/shared/include
#                    ${Geant4_INCLUDE_DIR}
#                    ${ROOT_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/shared/include
                    ${Geant4_INCLUDE_DIR}
                    #/gpfs/home/t2k/ogawat/myt2kwork/neut_generator/neut_5.4.0.1/include  
                    /opt/MyNEUT/neut_5.4.0.1/include  
                    ${ROOT_INCLUDE_DIRS})
# >>> neut is added 22/01/07
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_SOURCE_DIR}/shared/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/shared/include/*.hh)

# 1. set library path 
# 2. set library name 
link_directories(
                  #/gpfs/home/t2k/ogawat/myt2kwork/neut_generator/neut_5.4.0.1/lib
                  /opt/MyNEUT/neut_5.4.0.1/lib/
                  /home/fiberstudy/.neutlibs/ # add 23/10/19 
                )

#
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(wavy wavy.cc ${sources} ${headers})

message(${Geant4_LIBRARIES})
message(${ROOT_LIBRARIES})

target_link_libraries(wavy ${Geant4_LIBRARIES}) # => replaced 20/06/01
target_link_libraries(wavy ${ROOT_LIBRARIES} ) # => replaced 20/06/01
#target_link_libraries(wavy ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )
target_link_libraries(wavy 
                      neutctrl neutfsipart neutfsivert neutnucfsistep 
                      neutnucfsivert neutpart neutvect neutvtx ) # >>> neut is added 22/01/07

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build wavy. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(wavy_SCRIPTS
    wavy.out
    wavy.in
    init.in
    vis.mac
    compile.sh
    run.sh
#	 run.mac
    tmp.mac
   )

foreach(_script ${wavy_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS wavy DESTINATION bin)

