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
#include <rosbag/view.h>
#include <rosbag/bag.h>

namespace bag_reader
{
//template
template<typename T>
void BagReader::readMessages( std::vector<typename T::ConstPtr> & data,
                              std::vector<double> & times,
                              const std::string & topic,
                              const int & numMessages,
                              const int & skipFactor,
                              const ros::Time & start_ts,
                              const ros::Time & end_ts) const
{
    rosbag::Bag bag ( bag_file_ );
    rosbag::View view ( bag, rosbag::TopicQuery(topic), start_ts, end_ts );
    data.clear();
    times.clear();
    data.reserve(view.size());
    times.reserve(view.size());
    int skip = 0;
    const int moduloSkipFactor = std::max<int>(1,skipFactor);
    for ( const rosbag::MessageInstance & m : view )
    {
        if ( int(data.size()) >= numMessages ) break;
        if ( skip++ % moduloSkipFactor != 0 ) continue;
        typename T::ConstPtr dataPtr = m.instantiate<T>();
        if ( dataPtr == nullptr ) continue;
        times.emplace_back(m.getTime().toSec());
        data.emplace_back(dataPtr);
    }
}

template<typename T>
void BagReader::writeMessages( const std::vector<typename T::ConstPtr> & data,
                               const std::vector<ros::Time>& times,
                               const std::string & bag_file,
                               const std::string & topic,
                               const bool & overwrite ) //const
{
    rosbag::Bag bag ( bag_file, overwrite ? rosbag::bagmode::Write : rosbag::bagmode::Append );
    for ( size_t idx = 0; idx < data.size(); ++idx )
    {
        bag.write(topic,times[idx],data[idx]);
    }
    bag.close();
}

}
