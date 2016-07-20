Requirements: Boost C++ libraries

UBUNTU: (TODO check pakage names)
sudo apt-get install libboost-dev
sudo apt-get install libboost-programoptions-dev
sudo apt-get install libboost-filesystem-dev

OSX:
0. in ~/.bash_profile add
export WIM2D_PATH=/Users/timill/GITHUB-REPOSITORIES/WIM2d
export BOOST_DIR="/opt/local/boost/"
export OMPI_CC=gcc
export OMPI_CXX=g++
export DYLD_LIBRARY_PATH=$BOOST_DIR/lib:$DYLD_LIBRARY_PATH

1. install gcc (4.7 or later) from macports (not xcode version)
eg sudo port install gcc49
sudo port install openmpi-gcc49

sudo ln -s /opt/local/bin/mpicxx /opt/local/bin/mpic++
sudo port select --set mpi openmpi-gcc49-fortran
sudo port select --set gcc gcc-mp-4.9
sudo port select --set g++ g++-mp-4.9
sudo port select --set c++ c++-mp-4.9

2. Download boost_1_60_0.tar.gz
from
https://sourceforge.net/projects/boost/files/boost/1.60.0/
untar boost_1_60_0.tar.gz
$BOOST0=`readlink -f boost_1_60_0`

3. cd $WIM2D_PATH/CXX/boost_1_60_0
cp bconfigure.sh binstall.sh $BOOST0

4. cd $BOOST0
cp tools/build/example/user-config.jam ~
./bconfigure.sh
./binstall.sh

5. To set up a new experiment in directory $NE_DIR:
cd $NE_DIR
$WIM2D_PATH/CXX/scripts/NewExp.sh
cd Build_cpp
make
cd ..
*Edit wim.cfg
./run_cpp.sh
