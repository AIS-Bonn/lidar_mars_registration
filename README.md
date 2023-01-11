# Lidar MARS Registration with optional EasyPBR


### [Video-Presentation](http://www.ais.uni-bonn.de/videos/iros2021_quenzel/) | [Paper](http://www.ais.uni-bonn.de/papers/IROS_2021_Quenzel.pdf)

[Real-time Multi-Adaptive-Resolution-Surfel 6D LiDAR Odometry using Continuous-time Trajectory Optimization](http://www.ais.uni-bonn.de/videos/iros2021_quenzel/)<br>
 [Jan Quenzel](https://www.ais.uni-bonn.de/%7Ejquenzel/) <sup>1</sup>,
 [Sven Behnke](https://www.ais.uni-bonn.de/behnke/) <sup>1</sup>,
 <br>
 <sup>1</sup>University of Bonn, Autonomous Intelligent Systems Group

<!--
<p align="middle">
  <img src="imgs/goliath_1_crop.png" width="240"/>
</p>-->

### Clone:
```sh
$ git clone https://github.com/AIS-Bonn/lidar_mars_registration.git --recursive
```
### Install dependencies:
```sh
$ sudo apt-get install python3-catkin-pkg
```
### Build with ROS:
```sh
$ catkin build lidar_mars_registration
```
### Running with ROS:
```sh
$ rosrun lidar_mars_registration lidar_mars_registration_node
```
This runs with parameters specified in ./config/live.cfg If you want to run it with a specific parameter set (multiple defined in config subfolder), then use the following:
```sh
$ rosrun lidar_mars_registration lidar_mars_registration_node _config_file_rel:="./config/urban_loco_ca.cfg"
```
Make sure, that "use_ros" is set to true for the registration_adaptor.
### Build with EasyPBR: 
To build with better visualization, you must have first installed [EasyPBR](https://github.com/JanQuenzel/easy_pbr).<br>
Afterwards this example can be build with 
```sh
$ make
```
### Running with EasyPBR:
After building one can run the executable created in the build folder with 
```sh
$ ./build/temp.linux-x86_64-3.6/run_lidar_mars_registration
```
or
```sh
$ python3 python/registration.py
```
Make sure, that "use_ros" is set to false for the registration_adaptor.
### Citation
```
@inproceedings{quenzel2021mars,
  title={{Real-time Multi-Adaptive-Resolution-Surfel 6D LiDAR Odometry using Continuous-time Trajectory Optimization}},
  author={Quenzel, Jan and Behnke, Sven},
  booktitle={Proceedings of IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
  year={2021}
}
```

### Comparison
The modified A-LOAM and floam used within our paper comparison can be found in the following forks: [floam](https://github.com/JanQuenzel/floam) and [A-LOAM](https://github.com/JanQuenzel/A-LOAM)<br>
For SuMa we first converted the bag files into KITTI binary format and used the original implementation by Behley and Stachniss: [SuMa](https://github.com/jbehley/SuMa)

### License
We make our code available under the BSD 3-clause License.
