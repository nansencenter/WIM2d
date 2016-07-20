1. Compile grid-saving fortran code with:
   cd $WIM2D_PATH/grid_setup; make
   There is also "make py" for a f2py python interface
   - easier to define new grids (more complicated landmasks)
   - see $WIM2D_PATH/grid_setup/README.txt

2. Go to where you want to run the WIM.
   i) Get infiles and edit:
      cp $WIM2D_PATH/fortran/grid_setup/infiles/grid/infile_grid.txt .
      cp $WIM2D_PATH/fortran/grid_setup/infiles/waves/infile_waves.txt .
   ii) Run:
      $WIM2D_PATH/fortran/scripts/NewExp.sh
      Gets makefile and creates header files

3. Compile:
   cd Build_F;
   make exec   :  pure fortran
   make py     :  f2py
   make mex2   :  mex function - simple
   make mex3   :  mex function - spectrum in/out
   make mex4   :  mex function - spectrum in/out, does breaking on mesh
   make All    :  exec, f2py
   make all    :  exec, f2py, mex2, mex3, mex4

6. cd ..
   - make sure bin (has binaries) is in matlab/python paths
   - can now run
