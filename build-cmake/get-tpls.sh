#!/bin/bash
#
# set variables
source module-file
echo percept_build_dir = $percept_build_dir

do_get=1

do_get_cgns=$do_get
do_get_opennurbs=$do_get
do_get_xml=$do_get
do_get_boost=$do_get
do_get_yaml=$do_get
do_get_zlib=$do_get
do_get_hdf5=$do_get
do_get_netcdf=$do_get
do_get_pnetcdf=$do_get
do_get_gtest=$do_get
do_get_trilinos=$do_get

do_get_percept=0
get_cmake=0

################################################################################
#Cmake
################################################################################

if [ $get_cmake -eq 1 ]
then

  cd $percept_build_dir/packages
  curl -s --retry 3 -o cmake-3.1.0-rc2.tar.gz http://www.cmake.org/files/v3.1/cmake-3.1.0-rc2.tar.gz
  tar zxf cmake-3.1.0-rc2.tar.gz

fi

################################################################################
#CGNS, Version 3.2.1
################################################################################

if [ $do_get_cgns -eq 1 ]
then

  echo Getting CGNS...
  cd $percept_build_dir/packages
  git clone -b master https://github.com/CGNS/CGNS.git

fi

################################################################################
# OpenNURBS - special Percept version
################################################################################

#
if [ $do_get_opennurbs -eq 1 ]
then

  echo Getting OpenNURBS...
  cd $percept_build_dir/packages
  mkdir -p opennurbs
  tar zxf $percept_build_dir/opennurbs-percept.tar.gz

fi


################################################################################
#libxml2, Version 2.9.2
################################################################################

if [ $do_get_xml -eq 1 ]
then

  echo Getting libxml2...
  cd $percept_build_dir/packages
  curl -s --retry 3 -o libxml2-2.9.2.tar.gz http://www.xmlsoft.org/sources/libxml2-2.9.2.tar.gz
  tar zxf libxml2-2.9.2.tar.gz
fi

################################################################################
#boost, Version 1.55.0
################################################################################

if [ $do_get_boost -eq 1 ]
then
  echo Getting boost...
  cd $percept_build_dir/packages
  curl -s --retry 3 -o boost_1_66_0.tar.gz http://iweb.dl.sourceforge.net/project/boost/boost/1.66.0/boost_1_66_0.tar.gz
  tar zxf boost_1_66_0.tar.gz
fi

################################################################################
#yaml-cpp, Version 0.3.0
################################################################################


if [ $do_get_yaml -eq 1 ]
then

  echo Getting yaml-cpp...

  cd $percept_build_dir/packages

  git clone https://github.com/jbeder/yaml-cpp.git
  mv yaml-cpp yaml-cpp-0.5.3
  cd yaml-cpp-0.5.3
  git checkout release-0.5.3

fi


################################################################################
# zlib, Version 1.2.8
################################################################################

if [ $do_get_zlib -eq 1 ]
then
  echo Getting zlib...
  cd $percept_build_dir/packages
  curl -s --retry 3 -o zlib-1.2.8.tar.gz http://zlib.net/fossils/zlib-1.2.8.tar.gz
  tar zxf zlib-1.2.8.tar.gz
fi

################################################################################
#hdf5, Version 1.8.12
################################################################################

if [ $do_get_hdf5 -eq 1 ]
then
  echo Getting hdf5...
  cd $percept_build_dir/packages/
  curl -s --retry 3 -o hdf5-1.8.12.tar.gz https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz
  tar zxf hdf5-1.8.12.tar.gz
fi

################################################################################
#netcdf, Version 4.3.3.1
################################################################################

if [ $do_get_netcdf -eq 1 ]
then
  echo Getting netcdf...
  cd $percept_build_dir/packages/
  curl -s --retry 3 -o netcdf-c-4.3.3.1.tar.gz https://codeload.github.com/Unidata/netcdf-c/tar.gz/v4.3.3.1
  tar zxf netcdf-c-4.3.3.1.tar.gz
fi

################################################################################
#parallel netcdf, Version 1.6.1
################################################################################

if [ $do_get_pnetcdf -eq 1 ]
then
  echo Getting pnetcdf...
  cd $percept_build_dir/packages/
  curl -s --retry 3 -o parallel-netcdf-1.6.1.tar.gz http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.6.1.tar.gz
  tar zxf parallel-netcdf-1.6.1.tar.gz
fi

################################################################################
# googletest (gtest)
################################################################################

if [ $do_get_gtest -eq 1 ]
then
  echo Getting googletest...
  cd $percept_build_dir/packages/
  git clone  https://github.com/google/googletest.git
fi


################################################################################
#Trilinos, head
################################################################################


if [ $do_get_trilinos -eq 1 ]
then
  echo Getting Trilinos...
  cd $percept_build_dir/packages/
  git clone https://github.com/trilinos/Trilinos.git
  cd Trilinos

  git checkout master
  mkdir -p build
  cd build
  cp $percept_build_dir/sample-trilinos-release.config trilinos-release.config

fi


################################################################################
#Percept
################################################################################

if [ $do_get_percept -eq 1 ]
then

  # future
  # git clone https://github.com/percept/percept.git
  ls

fi


