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

#pragma once

#include <vector>
#include <memory>
#include <Eigen/Geometry>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/Twist.h>
#include <std_msgs/UInt8.h>
#include <tf2_ros/buffer.h>
#include "sophus/se3.hpp"

namespace rosbag
{
    class Bag;
}
typedef std::shared_ptr<rosbag::Bag> BagPtr;

namespace bag_reader
{

class BagReader
{
public:
    typedef std::shared_ptr<BagReader> Ptr;
    BagReader ( const std::string & bag_file, const ros::Time & start_time = ros::Time(0), const ros::Time & end_time = ros::Time(0) );

    static Ptr create( const std::string & bag_file, const ros::Time & start_time = ros::Time(0), const ros::Time & end_time = ros::Time(0) )
    {
        return std::make_shared<BagReader>(bag_file, start_time, end_time);
    }

    geometry_msgs::TransformStamped getStaticTransform ( const std::string & target, const std::string & source );

    template<typename T>
    void readMessages( std::vector<typename T::ConstPtr> & data,
                       std::vector<double> & times,
                       const std::string & topic,
                       const int & numMessages = std::numeric_limits<int>::max(),
                       const int & skipFactor = 1,
                       const ros::Time & start_ts = ros::TIME_MIN,
                       const ros::Time & end_ts = ros::TIME_MAX ) const;

    template<typename T>
    static void writeMessages( const std::vector<typename T::ConstPtr> & data,
                               const std::vector<ros::Time>& times,
                               const std::string & bag_file,
                               const std::string & topic, const bool & overwrite ); //const;
    template<typename T>
    static void writeMessages( const std::vector<typename T::ConstPtr> & data,
                               const std::vector<double>& times,
                               const std::string & bag_file,
                               const std::string & topic, const bool & overwrite ) //const
    {
        std::vector<ros::Time> ros_times(times.size());
        for ( size_t idx = 0; idx < times.size(); ++idx )
        {
            ros_times[idx] = ros::Time().fromSec(times[idx]);
        }
        writeMessages<T> ( data, ros_times, bag_file, topic, overwrite);
    }


    void addStaticTransform ( const std::string & target, const std::string & source, const Sophus::SE3d & tf );

    void readClouds( std::vector<sensor_msgs::PointCloud2::ConstPtr> & clouds,
                     const std::string & laser_topic,
                     const int & numClouds = std::numeric_limits<int>::max(),
                     const int & skipFactor = 1,
                     const ros::Time & start_ts = ros::TIME_MIN,
                     const ros::Time & end_ts = ros::TIME_MAX );

    void readCloudBagStamps( std::vector<ros::Time> & bag_stamps, const std::string & laser_topic, const int & numClouds, const int & skipFactor );

    void readCloudsWithPoses( std::vector<sensor_msgs::PointCloud2::ConstPtr> & clouds,
                              std::vector<geometry_msgs::TransformStamped> & poses,
                              const std::string & laser_topic,
                              const std::string & target,
                              const std::string & source,
                              const int & numClouds = std::numeric_limits<int>::max(),
                              const int & skipFactor = 1,
                              std::vector<geometry_msgs::TransformStamped> * ref_poses = nullptr,
                              std::string * ref = nullptr );

    void readPoses( std::vector<ros::Time> & timestamps, std::vector<geometry_msgs::TransformStamped> & poses, const std::string & target, const std::string & source, const bool & use_non_static = false );

    void populateTfBufferFromBag ( const bool & use_additional = false, const std::string & filePath = "");
    void populateTfBufferFromFile ( const std::string & filePath, const std::string & target, const std::string & source, const bool & timeInSec = false, const bool & preClearBuffer = false );

    std::shared_ptr<tf2_ros::Buffer> getTfBuffer ( )
    {
        return tfBufferPtr_;
    }
    ros::Time getMinTimeStamp ( ) const
    {
        return m_min_time_stamp;
    }
    ros::Time getMaxTimeStamp ( ) const
    {
        return m_max_time_stamp;
    }
    bool hasNonStatic() const
    {
        return m_has_non_static;
    }
protected:
    std::shared_ptr<tf2_ros::Buffer> tfBufferPtr_;
    std::string bag_file_;
    bool m_has_non_static = false;
    ros::Time m_max_time_stamp;
    ros::Time m_min_time_stamp;
    ros::Time m_load_start_time;
    ros::Time m_load_end_time;
    BagPtr m_bag = nullptr;
};
}
#include "BagReaderImpl.h"

