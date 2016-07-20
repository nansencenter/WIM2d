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

mkdir -p Build_F
mv infile_grid.txt infile_waves.txt Build_F
cd Build_F

# set up grid
$GS/scripts/grid_setup.sh

# get Makefile
ln -s $BUILD/Makefile .

# to enable matlab (mex functions)
cd ..
cp $WIM2D_PATH/fortran/scripts/startup_local.m .
cp $WIM2D_PATH/matlab/main/infiles/infile_matlab.txt .

echo "Compile in Build_F"
echo "Run here with code from $WIM2D_PATH/fortran/scripts"
echo "To use mex functions, edit infile_matlab.txt"
echo "and type run('$WIM2D_PATH/matlab/run_WIM2d.m');"
