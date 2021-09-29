#!/usr/bin/env python3.6
import sys
import os
import time
from easypbr import *
from geometry_msgs.msg import *
from lidar_mars_registration import *

name="drz"
#name="oxford"
#name="urban_loco"

config_file=os.path.join(os.path.abspath(os.getcwd()),"config/"+name+".cfg")
print("config_file: {}".format(config_file))
view = Viewer.create(config_file)
print("created viewer")
wrapper = EasyPBRwrapper.create(config_file,view)
print("created wrapper")
regis = MarsSplineRegistrationAdaptor.create(config_file)
print("created regis")
bag_player=BagAdaptor.create(config_file)
print("created bag")
wrapper.register_module(bag_player)
wrapper.register_module(regis)
regis.register_module(bag_player)
print("added bag adapt to gui")

num_resets = 0
cur_cloud_id = 0
wasRunning = False

while True:
    view.update()
    if wrapper.isRunning() and bag_player.has_next_cloud() and ( wrapper.shouldPlayOne() or wrapper.shouldSkipOne() ):

        cur_cloud_id = cur_cloud_id + 1
        print("show cloud: {}".format(cur_cloud_id))
        cloud = bag_player.get_next_cloud( wrapper.shouldSkipOne() )
        if wrapper.shouldSkipOne() or cloud is None:
            print("skipping.")
            continue
        cloud.m_vis.m_show_points=True

        while bag_player.has_imu_until_next_cloud():
            regis.imu_msgs( bag_player.get_next_imu() )

        while bag_player.has_gps_until_next_cloud():
            regis.gps_msgs( bag_player.get_next_gps() )

        print("register now.")
        regis.cloud_msg( cloud )
        print("show now.")
        wasRunning = True
