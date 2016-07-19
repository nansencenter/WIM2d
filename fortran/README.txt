1. Compile grid-saving fortran code with:
   cd $WIM2D_PATH/grid_setup; make
   There is also "make py" for a f2py python interface
   - easier to define new grids (more complicated landmasks)
   - see README_2.txt

2. go to where you want to run the WIM:
   mkdir Build; cd Build

   NB if you want to have many different grids,
   make a different "Build" dir for each one, then link to the desired dir
   eg "mkdir Build_1d; ln -s Build_1d Build"

3. copy templates for infile_waves.txt and infile_grid.txt 
   from $WIM2D_PATH/grid_setup/infiles;
   edit

4. Get grid files, header files
   $WIM2D_PATH/NewExp.sh

5. Compile
   make exec   :  pure fortran
   make py     :  f2py
   make mex2   :  mex function - simple
   make mex3   :  mex function - spectrum in/out
   make mex4   :  mex function - spectrum in/out, does breaking on mesh
   make All    :  exec, f2py
   make all    :  exec, f2py, mex2, mex3, mex4

6. cd ..
   - make sure Build (has binaries) is in matlab/python paths
   - can now run
