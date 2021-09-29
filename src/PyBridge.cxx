/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration

Copyright (c) 2021, Computer Science Institute VI, University of Bonn.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "PyBridge.h"

//my stuff 
#include "BagAdaptor.h"
#ifdef USE_EASY_PBR
#include "EasyPBRwrapper.h"
#include "easy_pbr/Viewer.h"
#endif
#include "MarsSplineRegistrationAdaptor.h"

// https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html
// PYBIND11_MAKE_OPAQUE(std::vector<int>); //to be able to pass vectors by reference to functions and have things like push back actually work 
// PYBIND11_MAKE_OPAQUE(std::vector<float>, std::allocator<float> >);

namespace py = pybind11;

#ifdef USE_EASY_PBR
using namespace easy_pbr;
#endif
PYBIND11_MODULE(lidar_mars_registration, m) {

#ifdef USE_EASY_PBR
    py::class_<EasyPBRwrapper, std::shared_ptr<EasyPBRwrapper>> (m, "EasyPBRwrapper")
    .def_static("create",  &EasyPBRwrapper::create<const std::string&, const std::shared_ptr<Viewer>& > ) //for templated methods like this one we need to explicitly instantiate one of the arguments
    .def("hasStarted", &EasyPBRwrapper::hasStarted)
    .def("isRunning", &EasyPBRwrapper::isRunning)
    .def("shouldSkipOne", &EasyPBRwrapper::shouldSkipOne)
    .def("shouldPlayOne", &EasyPBRwrapper::shouldPlayOne)
    .def("register_module", py::overload_cast<const std::shared_ptr<BagAdaptor> >(&EasyPBRwrapper::register_module) )
    .def("register_module", py::overload_cast<const std::shared_ptr<MarsSplineRegistrationAdaptor> >(&EasyPBRwrapper::register_module) )
    ;
#endif

    // Quick Hack just for passing this class through python
    py::class_<geometry_msgs::TransformStamped, std::shared_ptr<geometry_msgs::TransformStamped>> (m, "TransformStamped")
    ;
    py::class_<sensor_msgs::Imu, std::shared_ptr<sensor_msgs::Imu>> (m, "Imu")
    ;
    py::class_<nav_msgs::Odometry, std::shared_ptr<nav_msgs::Odometry>> (m, "Odometry")
    ;
    py::class_<Sophus::SE3d, std::shared_ptr<Sophus::SE3d>> (m, "SEd")
    ;


    py::class_<MarsSplineRegistrationAdaptor, std::shared_ptr<MarsSplineRegistrationAdaptor> > ( m, "MarsSplineRegistrationAdaptor")
    .def_static("create", &MarsSplineRegistrationAdaptor::create<const std::string> ) //for templated methods like this one we need to explicitly instantiate one of the arguments
    .def("cloud_msg", &MarsSplineRegistrationAdaptor::cloud_msg )
    .def("imu_msgs", &MarsSplineRegistrationAdaptor::imu_msgs )
    .def("gps_msgs", &MarsSplineRegistrationAdaptor::gps_msgs )
    .def("was_keyframe", &MarsSplineRegistrationAdaptor::was_keyframe )
#ifdef USE_EASY_PBR
    .def("register_module", py::overload_cast<const std::shared_ptr<BagAdaptor> >(&MarsSplineRegistrationAdaptor::register_module) )
#endif
    ;

    // BagAdaptor
    py::class_<BagAdaptor, std::shared_ptr<BagAdaptor> > ( m, "BagAdaptor")
    .def_static("create", &BagAdaptor::create<const std::string> ) //for templated methods like this one we need to explicitly instantiate one of the arguments
    .def("get_next_cloud", &BagAdaptor::get_next_cloud )
    .def("has_next_cloud", &BagAdaptor::has_next_cloud )
    .def("get_next_mav", &BagAdaptor::get_next_mav )
    .def("has_imu_until_next_cloud", &BagAdaptor::has_imu_until_next_cloud )
    .def("has_gps_until_next_cloud", &BagAdaptor::has_gps_until_next_cloud )
    .def("get_next_imu", &BagAdaptor::get_next_imu )
    .def("get_next_gps", &BagAdaptor::get_next_gps )
    .def("restart", &BagAdaptor::restart )
    ;
}
