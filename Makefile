all:
	echo "Building LidarMarsRegistration"
	python3 -m pip install -v --user --editable ./

clean:
	python3 -m pip uninstall lidar_mars_registration
	rm -rf build *.egg-info build *lidar_mars_registration*.so lib*.so
