################################################################################
##  CMakeLists.txt file for mICAL Code                                         ##
##  @author                                                                   ##
##  @mail                                                                     ##
################################################################################

cmake_minimum_required(VERSION 2.8)

# set the project name
project(anal_ical)


# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++1y -m64")

#including root
find_package(ROOT CONFIG REQUIRED)
find_package(ROOT COMPONENTS Minuit)
# Setup Root include directories and compile definitions
include("${ROOT_USE_FILE}")

#including PostgreSQL
find_package(PostgreSQL REQUIRED)
include_directories (${PostgreSQL_INCLUDE_DIRS})
link_directories(${PostgreSQL_LIBRARY_DIRS})
#find_package(PG REQUIRED)
#include_directories (${PG_INCLUDE_DIRS})
#link_directories(${PG_LIBRARY_DIRS})



#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#-----------------------------------------------------------------------------
include(${Geant4_USE_FILE})

# Setup include directory for this project
include_directories(${PROJECT_SOURCE_DIR}/inc)
#link_directories(${PROJECT_SOURCE_DIR})

add_compile_options(-Wall -Wextra -pedantic)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
# file(GLOB sources1 ${PROJECT_SOURCE_DIR}/src/*.pcm)
file(GLOB headers ${PROJECT_SOURCE_DIR}/inc/*.hh)
file(GLOB headers1 ${PROJECT_SOURCE_DIR}/inc/*.h)
# file(GLOB headers2 ${PROJECT_SOURCE_DIR}/inc/*.cpp)

separate_arguments(ROOT_EXE_LINKER_FLAGS)
add_executable(anal_ical ${headers} ${headers1} ${sources})
# add_executable(anal_ical ${headers} ${headers1} ${headers2} ${sources})
# add_executable(anal_ical ${headers} ${sources})


target_link_libraries(anal_ical PUBLIC ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS} ${Geant4_LIBRARIES} ${PostgreSQL_LIBRARIES} config++ Geom)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(inominiical0_field_SCRIPTS
  )

foreach(_script ${inominiical0_field_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(inotarget DEPENDS anal_ical)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS anal_ical DESTINATION bin)
