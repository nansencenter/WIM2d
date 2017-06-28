1. Set up grid
i.    cd Build
ii.   change grid size by editing infile_grid.txt
iii.  run $WIM2D_PATH/fortran/grid_setup/scripts/grid_setup.sh
iv.   find grid files and a picture of the landmask in ../grid

2. Compile with make
3. cd ..; ./run_test_advect.sh
