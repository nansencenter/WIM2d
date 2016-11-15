RUN=$WIM2D_PATH/neXtSIMcoupling/run
cp $RUN/startup_local.m .
cp $RUN/launch_sim.m    . 
# could use/link existing one?

FDIR=$WIM2D_PATH/fortran
$FDIR/scripts/NewExp.sh

cp $FDIR/grid_setup/infiles/grid/infile_grid_xsim.txt Build_F/infile_grid.txt

# cat grid     > infile_dirs.txt
# cat nextwim >> infile_dirs.txt
