#pragma once

#include <geometry_msgs/Transform.h>
#include <geometry_msgs/TransformStamped.h>
#include "sophus/se3.hpp"

//#include <tf/LinearMath/Transform.h>
//static Sophus::SE3d transformToSophus ( const tf::Transform & transform )
//{
//    Sophus::SE3d pose;
//    const tf::Vector3 origin = transform.getOrigin();
//    const tf::Quaternion orientation = transform.getRotation();
//    pose.translation() << origin.x(), origin.y(), origin.z();
//    pose.setQuaternion( Eigen::Quaterniond(orientation.w(),orientation.x(),orientation.y(),orientation.z()) );
//    return pose;
//}

//static tf::Transform sophusToTransform( const Sophus::SE3d & pose )
//{
//    tf::Transform transform;
//    transform.setIdentity();
//    const Eigen::Vector3d origin = pose.translation();
//    transform.setOrigin(tf::Vector3(origin[0],origin[1],origin[2]));
//    const Eigen::Quaterniond orientation = pose.so3().unit_quaternion();
//    transform.setRotation(tf::Quaternion(orientation.x(),orientation.y(),orientation.z(),orientation.w()));
//    return transform;
//}

static geometry_msgs::Transform sophusToTransform( const Sophus::SE3d & pose )
{
    geometry_msgs::Transform transform;
    transform.translation.x = pose.translation().x();
    transform.translation.y = pose.translation().y();
    transform.translation.z = pose.translation().z();
    transform.rotation.w = pose.unit_quaternion().w();
    transform.rotation.x = pose.unit_quaternion().x();
    transform.rotation.y = pose.unit_quaternion().y();
    transform.rotation.z = pose.unit_quaternion().z();
    return transform;
}

static Sophus::SE3d transformToSophus ( const geometry_msgs::Transform & transform )
{
    const Eigen::Quaterniond q(transform.rotation.w,transform.rotation.x,transform.rotation.y,transform.rotation.z);
    if ( q.norm() < 1e-1 ) return Sophus::SE3d();
    Sophus::SE3d pose;
    pose.translation() << transform.translation.x, transform.translation.y, transform.translation.z;
    pose.setQuaternion( q );
    return pose;
}

static geometry_msgs::TransformStamped toTransformStamped ( const geometry_msgs::Transform & transform, const ros::Time & stamp, const std::string & parentFrame, const std::string & childFrame )
{
    geometry_msgs::TransformStamped stampedTransform;
    stampedTransform.transform = transform;
    stampedTransform.header.stamp = stamp;
    stampedTransform.header.frame_id = parentFrame;
    stampedTransform.child_frame_id = childFrame;
    return stampedTransform;
}
