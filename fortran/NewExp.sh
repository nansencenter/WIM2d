# run this from where you want to set up Build directory for a given grid/wave config

# make a new area to run fortran WIM from
GS=$WIM2D_PATH/fortran/grid_setup
BUILD=$WIM2D_PATH/fortran/Build

# check infiles for grid
if [ ! -f "infile_grid.txt" ] || [ ! -f "infile_waves.txt" ]
then
   echo "Missing either infiles_grid.txt or infiles_waves.txt"
   echo "- Get templates from $GS/infiles"
   exit
fi

$GS/grid_setup.sh

# get Makefile
ln -s $BUILD/Makefile .
