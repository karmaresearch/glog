# GLog

This project is a fork of <a href="https://github.com/karmaresearch/vlog">VLog</a>.

[![Build Status](https://travis-ci.org/karmaresearch/vlog.svg?branch=master)](https://travis-ci.org/karmaresearch/vlog)

## Installation 

We used CMake to ease the installation process. To build GLog, the following
commands should suffice:

```
mkdir build
cd build
cmake ..
make
```

External libraries should be automatically downloaded and installed in the same directory. The only library that should be already installed is zlib, which is necessary to read gzip files. This library is usually already present by default.

To enable the web-interface, you need to use the -DWEBINTERFACE=1 option to cmake.

If you want to build the DEBUG version of the program, including the web interface: proceed as follows:

```
mkdir build_debug
cd build_debug
cmake -DWEBINTERFACE=1 -DCMAKE_BUILD_TYPE=Debug ..
make
```

## VLDB 2021 experiments

To facilitate the reproduction of the experiments presented in the paper

```
Tsamoura, Efthymia, David Carral, Enrico Malizia, and Jacopo Urbani. "Materializing knowledge bases via trigger graphs." Proceedings of the VLDB Endowment 14, no. 6 (2021): 943-956
```

we have copied all the datasets, scripts, and other useful resources in the folder:

https://drive.google.com/drive/folders/14ZIn8fWkZ7oCbdGZXBTUI7kiI-1TD_kM?usp=sharing
