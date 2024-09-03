#-----PhysicsClasses-----

# Setup the project
cmake_minimum_required(VERSION 3.20 FATAL_ERROR) #minimum version of CMAKE required
project(PhysicsClasses)

# stablish c++ version
set(CMAKE_CXX_STANDARD 17) # same as ROOT built for Ubuntu
set(CMAKE_CXX_STANDARD_REQUIRED TRUE) # force require this standard
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#ROOT
find_package(ROOT 6.20 CONFIG REQUIRED COMPONENTS MathCore Physics)
include(${ROOT_USE_FILE})

#Add headers
include_directories(${CMAKE_SOURCE_DIR}/inc)

# Load Scripts.cmake
include(Scripts.cmake)
# Set install prefix
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

# Do everything
add_physlibrary(NAME PhysicsClasses LINK)
