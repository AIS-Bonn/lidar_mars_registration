# Lidar MARS Registration with optional EasyPBR
### Clone:
```sh
$ git clone https://git.ais.uni-bonn.de/jquenzel/lidar_mars_registration.git --recursive
```
### Install dependencies:
```sh
$ sudo apt-get install libgoogle-glog-dev libgflags-dev libceres-dev python3-catkin-pkg
```
### Build with ROS:
```sh
$ catkin build lidar_mars_registration
```
### Running with ROS:
```sh
$ rosrun lidar_mars_registration lidar_mars_registration_node
```
### Build with EasyPBR: 
To build the example, you must have first installed EasyPBR. Afterwards this example can be build with 
```sh
$ make
```
### Running with EasyPBR:
After building one can run the executable created in the build folder with 
```sh
$ ./build/temp.linux-x86_64-3.6/run_lidar_mrsctmap_registration
```
or
```sh
$ python3 python/registration.py
```


