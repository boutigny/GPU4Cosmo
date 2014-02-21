# Dependencies

In order to build this software the packages below are required:

* A modern C++ compiler
* ROOT data analysis framework: [http://root.cern.ch](http://root.cern.ch)
* A CUDA-capable GPU and the NVIDIA CUDA software development kit: [https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads)
* Healpix library: [http://healpix.jpl.nasa.gov](http://healpix.jpl.nasa.gov)

If you are using [Scientific Linux 6](https://www.scientificlinux.org), you can install these dependencies using
the commands below:


    $ yum install cuda-repo-rhel6
    $ yum install cuda
    $ yum install root chealpix chealpix-devel healpix-c++ healpix-c++-devel healpix-devel healpix root-graf3d root-graf3d-gl root-physics cfitsio cfitsio-devel root-graf-fitsio

For testing this software we are using:

* ROOT v5.34/11
* CUDA v5.5
* Healpix v3.11
