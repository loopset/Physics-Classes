# Install script for directory: /home/ivan/Physics-Classes

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/ivan/Physics-Classes/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/ivan/Physics-Classes/install/lib" TYPE SHARED_LIBRARY FILES "/home/ivan/Physics-Classes/libPhysicsClasses.so")
  if(EXISTS "$ENV{DESTDIR}/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so"
         OLD_RPATH "/home/ivan/root/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/ivan/Physics-Classes/CMakeFiles/PhysicsClasses.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ivan/Physics-Classes/install/include/AngComparator.h;/home/ivan/Physics-Classes/install/include/AngDifferentialXS.h;/home/ivan/Physics-Classes/install/include/AngFitter.h;/home/ivan/Physics-Classes/install/include/AngGlobals.h;/home/ivan/Physics-Classes/install/include/AngIntervals.h;/home/ivan/Physics-Classes/install/include/CalibrationRunner.h;/home/ivan/Physics-Classes/install/include/CalibrationSource.h;/home/ivan/Physics-Classes/install/include/FitCompositeTF1.h;/home/ivan/Physics-Classes/install/include/FitData.h;/home/ivan/Physics-Classes/install/include/FitGlobals.h;/home/ivan/Physics-Classes/install/include/FitInterface.h;/home/ivan/Physics-Classes/install/include/FitModel.h;/home/ivan/Physics-Classes/install/include/FitObjective.h;/home/ivan/Physics-Classes/install/include/FitPlotter.h;/home/ivan/Physics-Classes/install/include/FitRunner.h;/home/ivan/Physics-Classes/install/include/FitUtils.h;/home/ivan/Physics-Classes/install/include/Interpolators.h;/home/ivan/Physics-Classes/install/include/ModelPlotter.h;/home/ivan/Physics-Classes/install/include/PhysColors.h;/home/ivan/Physics-Classes/install/include/PhysExperiment.h;/home/ivan/Physics-Classes/install/include/PhysOMP.h;/home/ivan/Physics-Classes/install/include/PhysSF.h;/home/ivan/Physics-Classes/install/include/PhysSM.h;/home/ivan/Physics-Classes/install/include/PhysState.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/ivan/Physics-Classes/install/include" TYPE FILE FILES
    "/home/ivan/Physics-Classes/inc/AngComparator.h"
    "/home/ivan/Physics-Classes/inc/AngDifferentialXS.h"
    "/home/ivan/Physics-Classes/inc/AngFitter.h"
    "/home/ivan/Physics-Classes/inc/AngGlobals.h"
    "/home/ivan/Physics-Classes/inc/AngIntervals.h"
    "/home/ivan/Physics-Classes/inc/CalibrationRunner.h"
    "/home/ivan/Physics-Classes/inc/CalibrationSource.h"
    "/home/ivan/Physics-Classes/inc/FitCompositeTF1.h"
    "/home/ivan/Physics-Classes/inc/FitData.h"
    "/home/ivan/Physics-Classes/inc/FitGlobals.h"
    "/home/ivan/Physics-Classes/inc/FitInterface.h"
    "/home/ivan/Physics-Classes/inc/FitModel.h"
    "/home/ivan/Physics-Classes/inc/FitObjective.h"
    "/home/ivan/Physics-Classes/inc/FitPlotter.h"
    "/home/ivan/Physics-Classes/inc/FitRunner.h"
    "/home/ivan/Physics-Classes/inc/FitUtils.h"
    "/home/ivan/Physics-Classes/inc/Interpolators.h"
    "/home/ivan/Physics-Classes/inc/ModelPlotter.h"
    "/home/ivan/Physics-Classes/inc/PhysColors.h"
    "/home/ivan/Physics-Classes/inc/PhysExperiment.h"
    "/home/ivan/Physics-Classes/inc/PhysOMP.h"
    "/home/ivan/Physics-Classes/inc/PhysSF.h"
    "/home/ivan/Physics-Classes/inc/PhysSM.h"
    "/home/ivan/Physics-Classes/inc/PhysState.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ivan/Physics-Classes/install/lib/G__PhysicsClasses.cxx")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/ivan/Physics-Classes/install/lib" TYPE FILE FILES "/home/ivan/Physics-Classes/G__PhysicsClasses.cxx")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ivan/Physics-Classes/install/lib/libPhysicsClasses_rdict.pcm")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/ivan/Physics-Classes/install/lib" TYPE FILE FILES "/home/ivan/Physics-Classes/libPhysicsClasses_rdict.pcm")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ivan/Physics-Classes/install/lib/libPhysicsClasses.rootmap")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/ivan/Physics-Classes/install/lib" TYPE FILE FILES "/home/ivan/Physics-Classes/libPhysicsClasses.rootmap")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/ivan/Physics-Classes/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
