cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(lidar_mars_registration)

###   VARIABLES   ##############################################################
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_CXX_STANDARD 17) #needed also by the pybind so it can build itself in c++14 mode so we can use py:overload_cast
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -pthread -Werror=return-type -Wno-reorder -O2 -march=native -fopenmp ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -Werror=return-type -Wno-reorder -O3 -march=native -fopenmp ")

option(USE_EASY_PBR "Compile with easy_pbr" OFF)

######   PACKAGES   ############################################################
################################# test from basalt headers ###################
function(report_dependency_location DEP_NAME DEP_TARGET)
  if(${DEP_NAME}_DIR)
    # this is set if dependency is found via find_package
    set(ACTUAL_DEP_CMAKE_DIR ${${DEP_NAME}_DIR})
  else()
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.19.0")
      get_target_property(ACTUAL_DEP_CMAKE_DIR ${DEP_TARGET} SOURCE_DIR)
    else()
      # SOURCE_DIR works for interface targets only since 3.19
      set(ACTUAL_DEP_CMAKE_DIR "")
    endif()
  endif()
  get_target_property(ACTUAL_DEP_INCLUDE_DIR ${DEP_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Using ${DEP_NAME} from: ${ACTUAL_DEP_CMAKE_DIR} (include: ${ACTUAL_DEP_INCLUDE_DIR})")
endfunction()

# clang-specific compile flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-exceptions -fcolor-diagnostics -frelaxed-template-template-args -Wno-error=deprecated-declarations")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    # These are disabled to avoid lot's of warnings in Eigen code with clang 10
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation ")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation -S -mllvm --x86-asm-syntax=intel ")
  endif()
endif()

message(STATUS "Including internal Eigen from submodule")
# Disable Eigen tests
set(BUILD_TESTING_PREV_VALUE ${BUILD_TESTING})  # Note: BUILD_TESTING is always defined because we include CTest above
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
add_subdirectory(deps/eigen EXCLUDE_FROM_ALL)
# Restore previous value of BUILD_TESTING
set(BUILD_TESTING ${BUILD_TESTING_PREV_VALUE} CACHE BOOL "Build the testing tree." FORCE)
#    message(STATUS "Finding external Eigen")
#    find_package(Eigen3 3.4 REQUIRED)
report_dependency_location(Eigen3 Eigen3::Eigen)
include_directories(${EIGEN3_INCLUDE_DIR})
if(NOT EIGEN_ROOT)
  set(EIGEN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen")
endif()

set(EIGEN_INCLUDE_DIR_HINTS ${EIGEN_ROOT})
find_package(Eigen3 3.4.90 EXACT REQUIRED MODULE)
message(STATUS "Found Eigen headers in: ${EIGEN3_INCLUDE_DIR}")
if(NOT EIGEN3_INCLUDE_DIR MATCHES "^${EIGEN_ROOT}")
  message(WARNING "Found Eigen headers are outside of specified EIGEN_ROOT '${EIGEN_ROOT}'")
endif()

################################# end test from basalt headers ###################

#deps
if(USE_EASY_PBR)
  set ( BUILD_SHARED_LIBS_TMP ${BUILD_SHARED_LIBS} )
  set ( BUILD_SHARED_LIBS OFF )
  add_subdirectory(${PROJECT_SOURCE_DIR}/deps/easy_pbr )
  set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_TMP})
endif()

if(NOT USE_EASY_PBR)
  set( BUILD_SHARED_LIBS_TMP ${BUILD_SHARED_LIBS} )
  set( BUILD_SHARED_LIBS OFF )
endif()
set(ABSL_BUILD_TESTING OFF CACHE BOOL "Build absl testing")
add_subdirectory(${PROJECT_SOURCE_DIR}/deps/abseil-cpp )
set(SOPHUS_USE_BASIC_LOGGING ON CACHE BOOL "Sophus dont use fmt")
set(BUILD_SOPHUS_TESTS OFF CACHE BOOL "Build sophus tests")
set(BUILD_SOPHUS_EXAMPLES OFF CACHE BOOL "Build sophus examples")
add_subdirectory(${PROJECT_SOURCE_DIR}/deps/Sophus)
if(NOT USE_EASY_PBR)
  set( BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_TMP} )
endif()

if(USE_EASY_PBR)
  message("Using Easy PBR")
  add_definitions(-DUSE_EASY_PBR)
endif()

find_package(catkin REQUIRED COMPONENTS
   rosbag_storage
   nav_msgs
   sensor_msgs
   std_srvs
   tf2_msgs
   tf2_ros
)
if(NOT USE_EASY_PBR)
   catkin_package()
endif()

find_package(Boost REQUIRED COMPONENTS filesystem )
include_directories(${Boost_INCLUDE_DIR})
#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenCV REQUIRED)
find_package(TBB REQUIRED)

#unset(PYTHONINTERP_FOUND CACHE)
#unset(PYTHON_EXECUTABLE CACHE)
#unset(PYTHON_VERSION_STRING CACHE)
#unset(PYTHON_VERSION_MAJOR CACHE)
#unset(PYTHON_VERSION_MINOR CACHE)
#unset(PYTHON_VERSION_PATCH CACHE)
#unset(PYBIND11_INCLUDE_DIR CACHE)
#unset(PYTHON_INCLUDE_DIRS CACHE)
#unset(PYTHON_LIBRARIES CACHE)
#unset(PYTHON_MODULE_PREFIX CACHE)
#unset(PYTHON_MODULE_EXTENSION CACHE)
#unset(PYTHON_LIBRARY CACHE)

#if(USE_EASY_PBR)
#  find_package(EasyPBR REQUIRED)
#  set(PYBIND11_PYTHON_VERSION 3.6 CACHE STRING "" FORCE)
#  message("PYBIND11_PYTHON_VERSION is" , ${PYBIND11_PYTHON_VERSION})
#  add_subdirectory(${EASYPBR_SRC_PATH}/deps/pybind11  [EXCLUDE_FROM_ALL])
#else()
#  add_subdirectory(${PROJECT_SOURCE_DIR}/deps/pybind11  [EXCLUDE_FROM_ALL])
#endif()

if(USE_EASY_PBR)
  set(EASYPBR_SRC_PATH "${PROJECT_SOURCE_DIR}/deps/easy_pbr")
  find_path(EASYPBR_CORE_INCLUDE_DIR  NAMES easy_pbr/Viewer.h PATHS "${EASYPBR_SRC_PATH}/include" )
  find_path(EASYPBR_EASYGL_INCLUDE_DIR  NAMES Texture2D.h PATHS "${EASYPBR_SRC_PATH}/deps/easy_gl" )
  find_path(EASYPBR_LOGURU_INCLUDE_DIR  NAMES loguru.hpp PATHS "${EASYPBR_SRC_PATH}/deps/loguru" )
  find_path(EASYPBR_CONFIGURU_INCLUDE_DIR  NAMES configuru.hpp PATHS "${EASYPBR_SRC_PATH}/deps/configuru" )
  find_path(EASYPBR_UTILS_INCLUDE_DIR  NAMES Profiler.h PATHS "${EASYPBR_SRC_PATH}/deps/utils/include" )
  find_path(EASYPBR_CONCURRENTQUEUE_INCLUDE_DIR  NAMES concurrentqueue.h PATHS "${EASYPBR_SRC_PATH}/deps/concurrent_queue" )
  find_path(EASYPBR_IMGUI_INCLUDE_DIR  NAMES imgui.h PATHS "${EASYPBR_SRC_PATH}/deps/imgui" )
  find_path(EASYPBR_IMGUIZMO_INCLUDE_DIR  NAMES ImGuizmo.h PATHS "${EASYPBR_SRC_PATH}/deps/imguizmo" )
  find_path(EASYPBR_BETTERENUM_INCLUDE_DIR  NAMES enum.h PATHS "${EASYPBR_SRC_PATH}/deps/better_enums" )

  set( EASYPBR_INCLUDE_DIR ${EASYPBR_SRC_PATH}/extern ${EASYPBR_CORE_INCLUDE_DIR} ${EASYPBR_EASYGL_INCLUDE_DIR} ${EASYPBR_LOGURU_INCLUDE_DIR} ${EASYPBR_CONFIGURU_INCLUDE_DIR} ${EASYPBR_UTILS_INCLUDE_DIR}
                           ${EASYPBR_CONCURRENTQUEUE_INCLUDE_DIR} ${EASYPBR_IMGUI_INCLUDE_DIR} ${EASYPBR_IMGUIZMO_INCLUDE_DIR} ${EASYPBR_BETTERENUM_INCLUDE_DIR} )
  add_definitions( -DDEFAULT_CONFIG="${EASYPBR_SRC_PATH}/config/default_params.cfg" )
endif()

###   INCLUDES   #########################################################
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
                        ${PROJECT_SOURCE_DIR}/deps/abseil-cpp
                        ${PROJECT_SOURCE_DIR}/deps/Sophus
                        ${PROJECT_SOURCE_DIR}/deps/tracy
                        ) # Header folder
include_directories(${PROJECT_INCLUDE_DIR} )
#include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${TBB_INCLUDE_DIRS})

if(USE_EASY_PBR)
  include_directories(${EASYPBR_INCLUDE_DIR})
else()
  include_directories(${PROJECT_SOURCE_DIR}/deps/loguru)
  include_directories(${PROJECT_SOURCE_DIR}/deps/Configuru)
endif()
include_directories(${OPENCV_INCLUDE_DIR})
include_directories(${catkin_INCLUDE_DIRS})

###   SOURCES   #################################################################
set(MY_SRC
    ${PROJECT_SOURCE_DIR}/src/MarsMapWindow.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsAssociator.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsCell.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsMap.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsMapParameters.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsPointTypes.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsSurfel.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsSurfelInfo.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsDenseVoxelGrid.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsSparseVoxelGrid.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsBlockSparseVoxelGrid.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsPermutohedralLattice.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsSemanticClass.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsSplineRegistration.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsSplineRegistrator.cxx
    ${PROJECT_SOURCE_DIR}/src/BagReader.cxx
    ${PROJECT_SOURCE_DIR}/src/BagAdaptor.cxx
    ${PROJECT_SOURCE_DIR}/src/VisUtils.cxx
    ${PROJECT_SOURCE_DIR}/src/MarsSplineRegistrationAdaptor.cxx
    ${PROJECT_SOURCE_DIR}/deps/tracy/TracyClient.cpp
)

if(USE_EASY_PBR)
  set(MY_SRC ${MY_SRC}
    ${PROJECT_SOURCE_DIR}/src/EasyPBRwrapper.cxx
  )
else()
  set(MY_SRC ${MY_SRC}
    ${PROJECT_SOURCE_DIR}/src/loguru.cpp
  )
endif()

set(MY_ABSL_LIBS absl::hash absl::inlined_vector absl::flat_hash_map absl::flat_hash_set absl::node_hash_map absl::node_hash_set )

#if(USE_EASY_PBR)
#  message("Using Easy PBR")
#  add_definitions(-DUSE_EASY_PBR)
#endif()

#pybind
#pybind11_add_module( ${CMAKE_PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/PyBridge.cxx )

###   LIBS   ###############################################
add_library(${CMAKE_PROJECT_NAME}_libs SHARED ${MY_SRC} ${DEPS_SRC})

set(LIBS ${LIBS} -lpthread -ldl ${OpenCV_LIBS} ${Boost_LIBRARIES} ${TBB_LIBRARIES} ${MY_ABSL_LIBS} sophus ${catkin_LIBRARIES} )
if(USE_EASY_PBR)
  set(LIBS ${LIBS} easypbr_cpp )
endif()

target_link_libraries(${CMAKE_PROJECT_NAME}_libs Eigen3::Eigen ${LIBS} )
#target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_PROJECT_NAME}_libs )

#definitions for cmake variables that are necesarry during runtime
#target_compile_definitions(${CMAKE_PROJECT_NAME}_libs PRIVATE TRACY_ENABLE) #point to the cmakelist folder of the easy_pbr
target_compile_definitions(${CMAKE_PROJECT_NAME}_libs PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}") #point to the cmakelist folder of the easy_pbr
target_compile_definitions(${CMAKE_PROJECT_NAME}_libs PRIVATE CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}" ) # points to the CMakeList folder of whichever project included easy_pbr

###   EXECUTABLE   #######################################
if(USE_EASY_PBR)
  add_executable( run_${CMAKE_PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/main.cxx  )
  target_link_libraries(run_${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_PROJECT_NAME}_libs )
  target_compile_definitions(run_${CMAKE_PROJECT_NAME} PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}") #point to the cmakelist folder of the easy_pbr
endif()
add_executable( ${CMAKE_PROJECT_NAME}_node ${PROJECT_SOURCE_DIR}/src/MarsRosNode.cxx  )
target_link_libraries(${CMAKE_PROJECT_NAME}_node PRIVATE ${CMAKE_PROJECT_NAME}_libs )
target_compile_definitions(${CMAKE_PROJECT_NAME}_node PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}") #point to the cmakelist folder of the easy_pbr
