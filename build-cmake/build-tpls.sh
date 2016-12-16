#!/bin/bash

source module-file

do_parallel_build=-j12

do_build=0

do_build_superlu=$do_build
do_build_xml=$do_build
do_build_boost=$do_build
do_build_yaml=$do_build
do_build_zlib=$do_build
do_build_hdf5=$do_build
do_build_netcdf=$do_build
do_build_pnetcdf=$do_build
do_build_cgns=1 #$do_build
do_build_gtest=$do_build
do_build_opennurbs=$do_build

# do manually
do_build_trilinos=0

do_build_percept=0

################################################################################
#Cmake, Version 3.1.0
################################################################################

if [[ $do_build_cmake -eq 1 ]] && [[ $do_build -eq 1 ]]
then

  cd $percept_build_dir/packages/cmake-3.1.0-rc2
  ./configure --prefix=$percept_build_dir/install
  make clean
  make $do_parallel_build
  make -k install
fi

################################################################################
# gtest
################################################################################

if [ $do_build_gtest -eq 1 ]
then

  cd $percept_build_dir
  mkdir -p packages/googletest/googletest/build
  cp sample-gtest-release.config packages/googletest/googletest/build/gtest-release.config
  cd packages/googletest/googletest/build
  rm -rf cmake_install.cmake CM* install_manifest.txt libgtest* Makefile
  chmod +x ./gtest-release.config
  ./gtest-release.config
  make clean
  make
  make -k install

fi

################################################################################
# OpenNURBS - Percept special version
################################################################################

if [ $do_build_opennurbs -eq 1 ]
then

  cd $percept_build_dir/packages/opennurbs
  make

fi

################################################################################
# SuperLU, Version 2.9.2
################################################################################

#Build

if [ $do_build_superlu -eq 1 ]
then

  cd $percept_build_dir/packages/SuperLU_4.3

  #To find out what the correct platform extension PLAT is:

  # uname -m

  #Edit make.inc as shown below (diffs shown from baselien).

  #PLAT = _x86_64
  #SuperLUroot   = /your_path/install/SuperLU_4.3 i.e., $percept_build_dir/install/SuperLU_4.3
  #BLASLIB       = -L/usr/lib64 -lblas
  #CC           = mpicc
  #FORTRAN            = mpif77

  #Now, make some new directories:

  mkdir -p $percept_build_dir/install/SuperLU_4.3
  mkdir -p $percept_build_dir/install/SuperLU_4.3/lib
  mkdir -p $percept_build_dir/install/SuperLU_4.3/include

  cd $percept_build_dir/packages/SuperLU_4.3
  make  $do_parallel_build SuperLUroot=$percept_build_dir/packages/SuperLU_4.3
  cp SRC/*.h $percept_build_dir/install/SuperLU_4.3/include
  cp lib/*.a $percept_build_dir/install/SuperLU_4.3/lib
fi

################################################################################
#libxml2, Version 2.9.2
################################################################################

if [ $do_build_xml -eq 1 ]
then

  cd $percept_build_dir/packages/libxml2-2.9.2
  CC=mpicc CXX=mpicxx ./configure -without-python --prefix=$percept_build_dir/install
  make clean
  make $do_parallel_build
  make -k install
fi

################################################################################
#boost, Version 1.55.0
################################################################################


if [ $do_build_boost -eq 1 ]
then
  cd $percept_build_dir/packages/boost_1_55_0

  echo "using mpi : `which mpicxx` ;" >> ./tools/build/v2/user-config.jam
  ./bootstrap.sh --prefix=$percept_build_dir/install --with-libraries=signals,regex,filesystem,system,mpi,serialization,thread,program_options,exception,graph,graph_parallel
  ./b2 link=static -j 12 2>&1 | tee boost_build_one
  ./b2 link=static -j 12 install 2>&1 | tee boost_build_intall
fi

################################################################################
#yaml-cpp, Version 0.3.0
################################################################################

if [ $do_build_yaml -eq 1 ]
then

  cd $percept_build_dir/packages/yaml-cpp-0.3.0
  mkdir -p build
  cd build
  cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$percept_build_dir/install ..
  make clean
  make  $do_parallel_build
  make -k install
fi

################################################################################
# zlib, Version 1.2.8
################################################################################

if [ $do_build_zlib -eq 1 ]
then
  cd $percept_build_dir/packages/zlib-1.2.8
  CC=gcc CXX=g++ CFLAGS=-O3 CXXFLAGS=-O3 ./configure --prefix=$percept_build_dir/install/
  make clean
  make $do_parallel_build
  make -k install
fi

################################################################################
#hdf5, Version 1.8.12
################################################################################

if [ $do_build_hdf5 -eq 1 ]
then

  cd $percept_build_dir/packages/hdf5-1.8.12
  ./configure CC=mpicc FC=mpif90 CXX=mpicxx CXXFLAGS="-fPIC -O3" CFLAGS="-fPIC -O3" FCFLAGS="-fPIC -O3" \
      --enable-parallel --enable-shared --enable-static-exec --enable-debug=no --enable-production   \
      --with-zlib=$percept_build_dir/install --prefix=$percept_build_dir/install
  make clean
  make $do_parallel_build
  make -k install
  #make check
fi

################################################################################
#pnetcdf, Version 1.6.1
################################################################################

if [ $do_build_pnetcdf -eq 1 ]
then

  cd $percept_build_dir/packages/parallel-netcdf-1.6.1
  rm -f config.cache
  ./configure --disable-fortran --prefix=$percept_build_dir/install AR_FLAGS='cru' CC=mpicc CPPFLAGS="-DNDEBUG" MPICC=mpicc CFLAGS="-I$percept_build_dir/install/include -fPIC -O3" CXXFLAGS="-fPIC -O3" LDFLAGS=-L$percept_build_dir/install/lib
  make clean
  make $do_parallel_build
  make -k install
  #make check

fi

################################################################################
#netcdf, Version 4.3.3
################################################################################

#Complex Models (expert usage only)

#In netcdf/include/netcdf.h, the following defines need to be changed to support complex models.

#define NC_MAX_DIMS     65536    /* max dimensions per file */
#define NC_MAX_VARS     524288   /* max variables per file */

#For a definiton of Complex Models, please note the following page:


if [ $do_build_netcdf -eq 1 ]
then

  cd $percept_build_dir/packages/netcdf-c-4.3.3.1
  rm -f config.cache
  sed -i -e "s/#define NC_MAX_DIMS.*$/#define NC_MAX_DIMS 65536/g" include/netcdf.h
  sed -i -e "s/#define NC_MAX_VARS.*$/#define NC_MAX_VARS 524288/g" include/netcdf.h
  ./configure --prefix=$percept_build_dir/install CC=mpicc FC=mpif90 CXX=mpicxx CPPFLAGS="-DNDEBUG" CFLAGS="-I$percept_build_dir/install/include -O3" LDFLAGS=-L$percept_build_dir/install/lib \
      --disable-fsync --disable-v2 --disable-dap --disable-doxygen --enable-netcdf-4 --enable-shared --enable-pnetcdf

  make $do_parallel_build
  make -k install

fi

################################################################################
#cgns
################################################################################

if [ $do_build_cgns -eq 1 ]
then

  cd $percept_build_dir
  mkdir -p packages/CGNS/build
  cp sample-cgns-release.config packages/CGNS/build/cgns-release.config
  cd packages/CGNS/build
  rm -rf cmake_install.cmake CMakeCache.txt install_manifest.txt Makefile
  chmod +x cgns-release.config
  ./cgns-release.config
  make clean
  make $do_parallel_build
  make -k install
fi

################################################################################
#Trilinos, head
################################################################################

if [ $do_build_trilinos -eq 1 ]
then

  cd $percept_build_dir/packages/Trilinos
  mkdir -p build

  cp $percept_build_dir/sample-trilinos-release.config trilinos-release.config

  cd $percept_build_dir/packages/Trilinos/build
  chmod +x trilinos-release.config


  ./trilinos-release.config
  make $do_parallel_build
  make -k install

fi

################################################################################
#Percept
################################################################################

if [ $do_build_percept -eq 1 ]
then

  cd $percept_build_dir/packages/Percept/build-cmake

  cp sample-percept-release.config  percept-release.config

  ./percept-relase.config
  make $do_parallel_build

fi

