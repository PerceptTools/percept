git clone https://github.com/PerceptTools/percept.git

cd percept/build-cmake

./setup.sh &> LOG_SETUP
./get-tpls.sh &> LOG_GET_TPLS
./build-tpls.sh &> LOG_BUILD_TPLS

# TODO just build Trilinos in build-tpls.sh?
cd packages/Trilinos/build

./trilinos-release.config &> LOG_TRILINOS_CMAKE
make -j24 &> LOG_TRILINOS_BUILD
make install &> LOG_TRILINOS_INSTALL

cd ../../../

./sample-percept-release.config &> LOG_PERCEPT_CMAKE
make -j24 &> LOG_PERCEPT_BUILD
