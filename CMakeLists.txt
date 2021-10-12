cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(lidar_mars_registration)

###   VARIABLES   ##############################################################
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17) #needed also by the pybind so it can build itself in c++14 mode so we can use py:overload_cast
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -pthread -Werror=return-type -Wno-reorder -O2 -march=native -fopenmp ")

option(USE_EASY_PBR "Compile with easy_pbr" OFF)

######   PACKAGES   ############################################################
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
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenCV REQUIRED)
find_package(TBB REQUIRED)

#deps
add_subdirectory(${PROJECT_SOURCE_DIR}/deps/abseil-cpp)
add_subdirectory(${PROJECT_SOURCE_DIR}/deps/Sophus)

unset(PYTHONINTERP_FOUND CACHE)
unset(PYTHON_EXECUTABLE CACHE)
unset(PYTHON_VERSION_STRING CACHE)
unset(PYTHON_VERSION_MAJOR CACHE)
unset(PYTHON_VERSION_MINOR CACHE)
unset(PYTHON_VERSION_PATCH CACHE)
unset(PYBIND11_INCLUDE_DIR CACHE)
unset(PYTHON_INCLUDE_DIRS CACHE)
unset(PYTHON_LIBRARIES CACHE)
unset(PYTHON_MODULE_PREFIX CACHE)
unset(PYTHON_MODULE_EXTENSION CACHE)
unset(PYTHON_LIBRARY CACHE)

if(USE_EASY_PBR)
  find_package(EasyPBR REQUIRED)
  set(PYBIND11_PYTHON_VERSION 3.6 CACHE STRING "" FORCE)
  message("PYBIND11_PYTHON_VERSION is" , ${PYBIND11_PYTHON_VERSION})
  add_subdirectory(${EASYPBR_SRC_PATH}/deps/pybind11  [EXCLUDE_FROM_ALL])
else()
  add_subdirectory(${PROJECT_SOURCE_DIR}/deps/pybind11  [EXCLUDE_FROM_ALL])
endif()

###   INCLUDES   #########################################################
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
                        ${PROJECT_SOURCE_DIR}/deps/abseil-cpp
                        ${PROJECT_SOURCE_DIR}/deps/Sophus
                        ${PROJECT_SOURCE_DIR}/deps/tracy
                        ) # Header folder
include_directories(${PROJECT_INCLUDE_DIR} )
include_directories(${Boost_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${TBB_INCLUDE_DIRS})
include_directories(${catkin_INCLUDE_DIRS})

if(USE_EASY_PBR)
  include_directories(${EASYPBR_INCLUDE_DIR})
else()
  include_directories(${PROJECT_SOURCE_DIR}/deps/loguru)
  include_directories(${PROJECT_SOURCE_DIR}/deps/Configuru)
  include_directories(${PROJECT_SOURCE_DIR}/deps/pybind11/include)
endif()

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
    ${PROJECT_SOURCE_DIR}/src/loguru.cpp
    ${PROJECT_SOURCE_DIR}/deps/tracy/TracyClient.cpp
)

if(USE_EASY_PBR)
  set(MY_SRC ${MY_SRC}
    ${PROJECT_SOURCE_DIR}/src/EasyPBRwrapper.cxx
  )
endif()

set(MY_ABSL_LIBS absl::hash absl::inlined_vector absl::flat_hash_map absl::flat_hash_set absl::node_hash_map absl::node_hash_set )

if(USE_EASY_PBR)
  message("Using Easy PBR")
  add_definitions(-DUSE_EASY_PBR)
endif()

#pybind
pybind11_add_module( ${CMAKE_PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/PyBridge.cxx )

###   LIBS   ###############################################
add_library(${CMAKE_PROJECT_NAME}_libs SHARED ${MY_SRC}  ${DEPS_SRC}  )

set(LIBS ${LIBS} -lpthread -ldl ${OpenCV_LIBS} ${Boost_LIBRARIES} ${TBB_LIBRARIES} ${MY_ABSL_LIBS} sophus ${catkin_LIBRARIES} )
if(USE_EASY_PBR)
  set(LIBS ${LIBS} ${EASYPBR_LIBRARY} )
endif()

target_link_libraries(${CMAKE_PROJECT_NAME}_libs ${LIBS} )
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_PROJECT_NAME}_libs )

#definitions for cmake variables that are necesarry during runtime
target_compile_definitions(${CMAKE_PROJECT_NAME}_libs PRIVATE TRACY_ENABLE) #point to the cmakelist folder of the easy_pbr
target_compile_definitions(${CMAKE_PROJECT_NAME}_libs PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}") #point to the cmakelist folder of the easy_pbr
target_compile_definitions(${CMAKE_PROJECT_NAME}_libs PRIVATE CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}" ) # points to the CMakeList folder of whichever project included easy_pbr

###   EXECUTABLE   #######################################
if(USE_EASY_PBR)
  add_executable( run_${CMAKE_PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/main.cxx  )
  target_link_libraries(run_${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_PROJECT_NAME}_libs )
target_compile_definitions(run_${CMAKE_PROJECT_NAME} PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}") #point to the cmakelist folder of the easy_pbr
endif()
add_executable( ${CMAKE_PROJECT_NAME}_node ${PROJECT_SOURCE_DIR}/src/MarsRegistration_Node.cxx  )
target_link_libraries(${CMAKE_PROJECT_NAME}_node PRIVATE ${CMAKE_PROJECT_NAME}_libs )
target_compile_definitions(${CMAKE_PROJECT_NAME}_node PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}") #point to the cmakelist folder of the easy_pbr
