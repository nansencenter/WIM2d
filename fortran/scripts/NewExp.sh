# run this from where you want to set up Build directory for a given grid/wave config

# make a new area to run fortran WIM from
GS=$WIM2D_PATH/fortran/grid_setup
BUILD=$WIM2D_PATH/fortran/Build
TOOLS=$WIM2D_PATH/fortran/tools

bdir=Build_F
tdir=scripts

# ========================================================
P=`pwd`
Bdir=$P/$bdir
mkdir -p $Bdir

cd $Bdir
ln -s $BUILD/Makefile Makefile
cp $GS/infiles/grid/infile_grid.txt $GS/infiles/waves/infile_waves.txt .

cd $P
echo " "
echo "*******************************************************"
echo "To set up grid and compile:"
echo "cd $bdir;"
echo "Edit infiles_grid.txt & infiles_waves.txt;"
echo "Type ../$tdir/grid_setup.sh"
echo "Compile with 'make' (pure fortran),"
echo "'make mex1', 'make mex2' or 'make mex3' (mex functions),"
echo "or 'make py' (f2py),"
echo "*******************************************************"
echo " "
# ========================================================

# ========================================================
# link useful scripts here
Tdir=$P/$tdir
mkdir -p $Tdir

# script to set-up grid
ln -s $GS/scripts/grid_setup.sh $Tdir/grid_setup.sh

# plotting tools
for f in $TOOLS/*
do
   g=`basename $f`
   ln -s $f $Tdir/$g
done
# ========================================================


# ========================================================
# for pure fortran
cd $P
cp $WIM2D_PATH/fortran/scripts/infiles/infile_nonIO.txt .
ln -s $WIM2D_PATH/fortran/scripts/run_WIM2d.sh


echo " "
echo "*******************************************************"
echo "To run pure fortran,"
echo "edit options in infile_nonIO.txt"
echo "and run with"
echo "$tdir/run_WIM2D.sh"
echo " "

# to enable matlab (mex functions)
ln -s $WIM2D_PATH/fortran/scripts/startup_local.m .
cp $WIM2D_PATH/matlab/main/infiles/infile_matlab.txt .

echo "To use mex functions, edit infile_matlab.txt"
echo "and (in matlab) type run('$WIM2D_PATH/matlab/run_WIM2d.m');"
echo "for pure matlab set MEX_OPT=0."
echo "*******************************************************"
echo " "
