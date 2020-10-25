# Push Pull Optimization

## Description

This is an implementation of "A Simple Push-Pull Algorithm for Blue-Noise Sampling" by Ahmed et al. [2016] with several modifications.
It was originally based on [Jianwei Guos implementation (one of the authors of the paper)](https://github.com/jianweiguo/push-pull). However after adding additional features and fixing several bugs in the original implementation it is structured fairly differently by now.

The main differences are:
* Added interface to MATLAB
* Adaptation of the push-pull optimization (PPO) to allow rectangular domains
* Added initialization method to start PPO on an existing point set
* Allows to optimize two tiles with compatible boundary properties (the revised data structures would allow a simple extension to even more tiles)
* Fixed several problems with the original implementation
  * Can handle identical points (in the input point set as well as identical points that occur during the PPO)
  * Points now stay in the tile of their respective replica. In the original code it was possible that a point (and its replicas) would gradually move away from the rest of the point set (only very rarely occurs with inputs that require many iterations until convergence).
  * Allow arbitrary number of neighbors for each point (original code leads to memory access violation for more than 20 neighbors). 
  * ...

## Installation
Prerequisites:
* Compiler with C++17 support
* Matlab installation to build the MEX files
* CGAL library - tested with version 5.0.2 (http://www.cgal.org);
  * The CGAL installer offers to download its two dependencies (GMP and MPFR) automatically. This (as at 15.05.2020) is blocked by the corporate firewall. You therefore have to manually install these libraries.
  The installer will display the URLs from which it tried and failed to download them. You can retrieve them from there manually. Alternatively the libraries that are required by CGAL 5.0.2 are available [here](\\sickcn.net\sick\CDSMI\projekte_mrk\Externe\TobiasBamberger\software\CGAL_dependencies).
* **Optional** (to display the point set):
  * Python installation (tested with 3.6.8) with numpy and matplotlib
    * Please refer to https://github.com/lava/matplotlib-cpp for supported versions of python and these libraries
  * define `PLOT_POINTS` in `pointSet.cpp`

With the prerequisites installed the project can be easily built via CMake (http://www.cmake.org).

The project contains two targets:

* **blueNoiseStandalone:** "normal" executable with several tests. Can also be used as playground to experiment with different parameter settings.
* **blueNoiseTwoTiles:** MATLAB executable (MEX) to run the algorithm via matlab. Cf. `myInterface.cpp` for valid input parameters and provided output.
