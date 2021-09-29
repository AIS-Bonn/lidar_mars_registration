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
#include "BagReader.h"
#include "sophusToTf.h"
#include <tf2_msgs/TFMessage.h>
#include <tf2_ros/buffer.h>
#include <rosbag/view.h>
#include <rosbag/bag.h>
//loguru
//#define LOGURU_REPLACE_GLOG 1
#include <loguru.hpp>

namespace bag_reader
{
BagReader::BagReader ( const std::string & bag_file, const ros::Time & start_time, const ros::Time & end_time)
    : bag_file_(bag_file), m_load_start_time (start_time), m_load_end_time ( end_time )
{
    populateTfBufferFromBag();
}
geometry_msgs::TransformStamped BagReader::getStaticTransform ( const std::string & target, const std::string & source )
{
    if ( ! tfBufferPtr_ ) populateTfBufferFromBag();
    //"arm_wrist_2_link", "laser_scanner_link"
    geometry_msgs::TransformStamped laser_link;
    try{
        laser_link = tfBufferPtr_->lookupTransform( target, source, m_min_time_stamp + ros::Duration((m_max_time_stamp - m_min_time_stamp ).toSec()/2) );
    }
    catch (tf2::TransformException ex){
        ROS_ERROR_STREAM(ex.what());
        ros::Duration(1.0).sleep();
    }
    return laser_link;
}

void BagReader::readCloudBagStamps( std::vector<ros::Time> & bag_stamps, const std::string & laser_topic, const int & numClouds, const int & skipFactor )
{
    if ( ! m_bag ) m_bag = std::make_shared<rosbag::Bag>(bag_file_);
    rosbag::Bag & bag = *m_bag;
    rosbag::View view ( bag, rosbag::TopicQuery(laser_topic) );
    bag_stamps.clear();
    bag_stamps.reserve(view.size());
    int skip = 0;
    const int moduloSkipFactor = std::max<int>(1,skipFactor);
    for ( const rosbag::MessageInstance & m : view )
    {
        if ( int(bag_stamps.size()) >= numClouds ) break;
        if ( skip++ % moduloSkipFactor != 0 ) continue;
        bag_stamps.emplace_back(m.getTime());
    }
    if ( ! m_has_non_static )
    {
        m_max_time_stamp = view.getEndTime();
        m_min_time_stamp = m_min_time_stamp > ros::Time() ? m_min_time_stamp : view.getBeginTime();
    }
}

void BagReader::readClouds( std::vector<sensor_msgs::PointCloud2::ConstPtr> & clouds, const std::string & laser_topic, const int & numClouds, const int & skipFactor, const ros::Time & start_ts, const ros::Time & end_ts )
{
    if ( ! m_bag ) m_bag = std::make_shared<rosbag::Bag>(bag_file_);
    rosbag::Bag & bag = *m_bag;
    rosbag::View view ( bag, rosbag::TopicQuery(laser_topic), start_ts, end_ts );
    clouds.clear();
    clouds.reserve(view.size());
    int skip = 0;
    const int moduloSkipFactor = std::max<int>(1,skipFactor);
    for ( const rosbag::MessageInstance & m : view )
    {
        if ( int(clouds.size()) >= numClouds ) break;
        if ( skip++ % moduloSkipFactor != 0 ) continue;
        sensor_msgs::PointCloud2::ConstPtr cloudPtr = m.instantiate<sensor_msgs::PointCloud2>();
        if ( cloudPtr == nullptr ) continue;
        clouds.emplace_back(cloudPtr);
    }
}

void BagReader::readCloudsWithPoses( std::vector<sensor_msgs::PointCloud2::ConstPtr> & clouds, std::vector<geometry_msgs::TransformStamped> & poses, const std::string & laser_topic, const std::string & target, const std::string & source, const int & numClouds, const int & skipFactor, std::vector<geometry_msgs::TransformStamped> * ref_poses, std::string * ref )
{
    clouds.clear();
    poses.clear();
    if ( ! tfBufferPtr_ ) populateTfBufferFromBag();

    std::vector<sensor_msgs::PointCloud2::ConstPtr> initial_clouds;
    readClouds( initial_clouds, laser_topic, numClouds, skipFactor );

    LOG(INFO) << "found " << initial_clouds.size() << " clouds. non stat: " << m_has_non_static;

    //"arm_wrist_1_link", "arm_wrist_2_link"
    geometry_msgs::TransformStamped tfTransform;
    for ( sensor_msgs::PointCloud2::ConstPtr & cloudPtr : initial_clouds )
    {
        if ( m_has_non_static && (cloudPtr->header.stamp < m_min_time_stamp || cloudPtr->header.stamp > m_max_time_stamp )) continue;
        try{
            tfTransform = tfBufferPtr_->lookupTransform(target, source, cloudPtr->header.stamp);
            if ( ref_poses && ref)
            {
                geometry_msgs::TransformStamped tfTransform2 = tfBufferPtr_->lookupTransform(*ref,target,cloudPtr->header.stamp);
                ref_poses->emplace_back(tfTransform2);
            }
            clouds.emplace_back(cloudPtr);
            poses.emplace_back(tfTransform);
        }
        catch (tf2::TransformException ex){
            ROS_ERROR_STREAM_THROTTLE(1,ex.what());
            continue;
        }
    }
}

void BagReader::readPoses( std::vector<ros::Time> & timestamps, std::vector<geometry_msgs::TransformStamped> & poses, const std::string & target, const std::string & source, const bool & use_non_static )
{
    poses.clear();
    if ( ! tfBufferPtr_ ) { populateTfBufferFromBag(); }

    const std::vector<ros::Time> search_timestamps = timestamps;
    timestamps.clear();

    const bool non_static = (use_non_static && m_has_non_static) || !use_non_static;

    geometry_msgs::TransformStamped tfTransform;
    for ( const ros::Time & stamp : search_timestamps )
    {
        if ( !stamp.isZero() && non_static && (stamp < m_min_time_stamp || stamp > m_max_time_stamp) ) { LOG(INFO) << "nope?" << stamp.toNSec() << " " << m_min_time_stamp.toNSec() << " " << m_max_time_stamp.toNSec() ; continue;}
        try{
            tfTransform = tfBufferPtr_->lookupTransform(target, source, stamp);
            poses.emplace_back(tfTransform);
            timestamps.emplace_back(stamp);
        }
        catch (tf2::TransformException ex){
            LOG(INFO) << "ups... " << stamp.toNSec() << " Err: "<< ex.what();
            continue;
        }
    }
}

void BagReader::populateTfBufferFromBag ( const bool & use_additional, const std::string & filePath )
{
    LOG(INFO) << "populating tf buffer.";
    if ( ! m_bag ) m_bag = std::make_shared<rosbag::Bag>(bag_file_);
    std::shared_ptr<rosbag::Bag> additional_bag = use_additional ? std::make_shared<rosbag::Bag>(filePath) : nullptr;
    rosbag::Bag & bag = use_additional ? *additional_bag : *m_bag;

    std::vector<std::string> topics;
    topics.emplace_back("/tf_static");
    topics.emplace_back("/tf");

    if ( ! use_additional )
    {
        std::shared_ptr<rosbag::View> view = std::make_shared<rosbag::View>(bag); // to get correct begin and end time.
        tfBufferPtr_ = std::make_shared<tf2_ros::Buffer>( view->getEndTime() - view->getBeginTime() );
        m_max_time_stamp = view->getEndTime();
        m_min_time_stamp = view->getBeginTime();
        LOG(INFO) << "View min: " << view->getBeginTime() << " max: " << view->getEndTime();
    }
    std::shared_ptr<rosbag::View> view = std::make_shared<rosbag::View>(bag, rosbag::TopicQuery(topics));

    if ( use_additional )
    {
        view = std::make_shared<rosbag::View>(bag, rosbag::TopicQuery(topics), m_min_time_stamp, m_max_time_stamp);
    }

    if ( ( !m_load_start_time.isZero() || !m_load_end_time.isZero() ) )
    {
        const double begin_stamp = view->getBeginTime().toSec();
        const double end_stamp = view->getEndTime().toSec();
        const double start_time = ( std::max<double>(begin_stamp,std::min<double>(begin_stamp + m_load_start_time.toSec(), end_stamp)) );
        const double end_time = ( std::max<double>(begin_stamp,std::min<double>(begin_stamp + m_load_end_time.toSec(), end_stamp)) );
        view = std::make_shared<rosbag::View>(bag, rosbag::TopicQuery(topics), ros::Time().fromSec(start_time), ros::Time().fromSec(end_time));
    }
    LOG(INFO) << "View min: " << view->getBeginTime() << " max: " << view->getEndTime();
    int static_tf_msgs = 0, non_static_tf_msgs = 0;
    for ( const rosbag::MessageInstance & m : *view )
    {
        tf2_msgs::TFMessage::ConstPtr tfMsg = m.instantiate<tf2_msgs::TFMessage>();
        if ( tfMsg == nullptr ) continue;
        if ( m.getTopic() == topics[0] )
        {
            for(const auto& transform : tfMsg->transforms)
            {
                tfBufferPtr_->setTransform(transform, "bagfile", true);
                ++static_tf_msgs;
            }
        }
        else
        {
            if ( ! use_additional ) m_has_non_static = true;
            for(const auto& transform : tfMsg->transforms)
            {
                //LOG(INFO) << "stamp: " << transform.header.stamp;
                tfBufferPtr_->setTransform(transform, "bagfile");
                ++non_static_tf_msgs;
            }
        }
    }
    LOG(INFO) << "View min: " << view->getBeginTime() << " max: " << view->getEndTime();
    LOG(INFO) << "tf buffer ts -> min: " << getMinTimeStamp() << " max: " << getMaxTimeStamp() << " st: " << static_tf_msgs << " nst: " << non_static_tf_msgs;
    LOG(INFO) << "tf buffer populated: " << tfBufferPtr_->allFramesAsString();
}

void BagReader::populateTfBufferFromFile ( const std::string & filePath, const std::string & target, const std::string & source, const bool & timeInSec, const bool & preClearBuffer )
{
    LOG(INFO) << "populating tf buffer from file: " << filePath;
    std::ifstream file ( filePath );
    if ( ! file.is_open() )
        LOG(FATAL) << "file not opened: " << filePath;

    uint64_t ns = 0;
    double s = 0;
    double tx, ty, tz, qw, qx, qy, qz;

    bool use_inverted = false;
    if ( tfBufferPtr_ )
    {
        std::string parent ="";
        if ( tfBufferPtr_->_getParent (source, ros::Time(0), parent) )
            use_inverted = true;
        LOG(INFO) << "source: " << source << " has parent: " << use_inverted << " p: "<< parent;
    }

    std::vector<geometry_msgs::TransformStamped> transforms;
    while ( file.good() )
    {
        if ( timeInSec )
            file >> s;
        else
            file >> ns;
        file >> tx >> ty >> tz >> qx >> qy >> qz >> qw;


        std::cout << "ns: "<< std::to_string(ns) << " t=[" << tx << "," << ty<<"," << tz<<"] q=["<<qx<<","<<qy<<","<<qz<<","<<qw<<"]\n";
        Sophus::SE3d pose;
        pose.setQuaternion(Eigen::Quaterniond(qw,qx,qy,qz));
        pose.translation() << tx, ty, tz;
        if ( use_inverted ) pose = pose.inverse();

        geometry_msgs::TransformStamped transformStamped;
        transformStamped.header.frame_id = use_inverted ? source : target;
        transformStamped.header.seq = transforms.size();
        transformStamped.child_frame_id = use_inverted ? target : source;

        transformStamped.transform = sophusToTransform(pose);
        transformStamped.header.stamp = timeInSec ? ros::Time().fromSec(s) : ros::Time().fromNSec(ns);

        if ( !preClearBuffer && (transformStamped.header.stamp < m_min_time_stamp || transformStamped.header.stamp > m_max_time_stamp) ) continue;
        transforms.emplace_back(transformStamped);
    }
    if ( transforms.empty() ) return;
    if ( preClearBuffer )
    {
        m_max_time_stamp = transforms.back().header.stamp;
        tfBufferPtr_ = std::make_shared<tf2_ros::Buffer>( transforms.back().header.stamp - transforms.front().header.stamp );
    }
    else
    {
        // TODO: make creation correct
    }
    m_has_non_static = true;

    LOG(INFO) << "MS: " << m_min_time_stamp << " LS: " << m_max_time_stamp;
    LOG(INFO) << "first: " << transforms.front().header.stamp << " last: " << transforms.back().header.stamp;
    LOG(INFO) << "target: " << target << " src: " << source;

    for ( const auto & transformStamped : transforms )
    {
        //std::cout << transformStamped << "\n";
        LOG(INFO) << "t: " << transformStamped.header.stamp.toNSec() << " p: " << transformStamped.transform.translation.x << " " << transformStamped.transform.translation.y << " " << transformStamped.transform.translation.z;
        tfBufferPtr_->setTransform(transformStamped, "file", false);
    }
    LOG(INFO) << "tf buffer populated from file: ";
    LOG(INFO) << tfBufferPtr_->allFramesAsString();
}

void BagReader::addStaticTransform ( const std::string & target, const std::string & source, const Sophus::SE3d & tf )
{
    geometry_msgs::TransformStamped transformStamped;
    transformStamped.header.frame_id = target;
    transformStamped.header.seq = 0;
    transformStamped.child_frame_id = source;
    transformStamped.transform = sophusToTransform(tf);
    transformStamped.header.stamp = ros::Time(0);
    tfBufferPtr_->setTransform(transformStamped, "specific_static", true);
}

} // namespace
