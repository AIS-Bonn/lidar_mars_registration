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
#include "MarsSplineRegistrationAdaptor.h"

#include <ros/ros.h>
#include <loguru.hpp>

//boost
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

int main(int argc, char **argv)
{
    loguru::init(argc, argv); // "-v INFO" by default used for loguru. Other: FATAL, WARNING or ERROR
    char log_path[PATH_MAX];
    loguru::suggest_log_path(PROJECT_SOURCE_DIR, log_path, sizeof(log_path));
    loguru::add_file(log_path, loguru::FileMode::Truncate, loguru::Verbosity_MAX);

    ros::init(argc, argv, "lidar_mars_registration_node");

    std::string config_file_rel = "./config/live.cfg";
    ros::NodeHandle pnh("~");
    pnh.param<std::string>("config_file_rel", config_file_rel, "./config/live.cfg");

    std::string config_file = fs::canonical(fs::path(PROJECT_SOURCE_DIR) / config_file_rel).string();
    LOG(INFO) << "config file: " << config_file;

    MarsSplineRegistrationAdaptor::Ptr mars_spline_registration_adaptor = MarsSplineRegistrationAdaptor::create(config_file);

    // handle callbacks until shut down
    ros::spin();
    return 0;
}
