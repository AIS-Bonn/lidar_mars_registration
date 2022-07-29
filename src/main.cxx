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
#include <loguru.hpp>
#include "easy_pbr/Viewer.h"
#include "easy_pbr/Scene.h"
#include "EasyPBRwrapper.h"
#include "BagAdaptor.h"
#include "MarsSplineRegistrationAdaptor.h"
//boost
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

using namespace easy_pbr;

int main(int argc, char *argv[]) {

    loguru::init(argc, argv); // "-v INFO" by default used for loguru. Other: FATAL, WARNING or ERROR
    char log_path[PATH_MAX];
    loguru::suggest_log_path(PROJECT_SOURCE_DIR, log_path, sizeof(log_path));
    loguru::add_file(log_path, loguru::FileMode::Truncate, loguru::Verbosity_MAX);

    std::string config_file = fs::canonical(fs::path(PROJECT_SOURCE_DIR) / "./config/urban_loco_hk.cfg").string();

    std::shared_ptr<Viewer> view = Viewer::create(config_file);  //creates it with default params
    std::shared_ptr<EasyPBRwrapper> wrapper = EasyPBRwrapper::create(config_file, view); // not used right now.
    std::shared_ptr<BagAdaptor> bag_player = BagAdaptor::create(config_file);
    std::shared_ptr<MarsSplineRegistrationAdaptor> regis = MarsSplineRegistrationAdaptor::create(config_file);
    wrapper->register_module(bag_player);
    wrapper->register_module(regis);
    regis->register_module(bag_player);

    int cur_cloud_id = 0;
    while (true)
    {
        view->update();
        if ( wrapper->isRunning() && bag_player->has_next_cloud() && ( wrapper->shouldPlayOne() || wrapper->shouldSkipOne() ) )
        {
            ++cur_cloud_id;
            LOG(INFO) <<  "show cloud: " << cur_cloud_id;
            auto cloud = bag_player->get_next_cloud( wrapper->shouldSkipOne() );
            if ( wrapper->shouldSkipOne() || cloud == nullptr )
                continue;
            cloud->m_vis.m_show_points = true;

            while ( bag_player->has_imu_until_next_cloud() )
                regis->imu_msgs( bag_player->get_next_imu() );

            while ( bag_player->has_gps_until_next_cloud() )
                regis->gps_msgs( bag_player->get_next_gps() );

            LOG(INFO) << "register now.";
            regis->cloud_msg( cloud );
            LOG(INFO) << "show now.";
        }
    }
    return 0;
}
