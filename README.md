WIM2d
=====

2d waves-in-ice module
Stand-alone version(s) for use on a single computer.

NB! All versions still being tested

Recommended compiler: gfortran from gnu 
https://gcc.gnu.org/wiki/GFortranBinaries
is recommended (versions from macports also worked but interfered with the compilation
of some other programs)

Versions:
- 1. matlab
- 2. fortran
- 3. ipython interface for fortran code
     - Mac OSX:
         - Install these python packages with macports (needs xcode with command line tools installed):
            * sudo port install py27-ipython
            * sudo port install py27-numpy
            * sudo port install py27-scipy
            * sudo port install py27-matplotlib

         - These packages probably aren't needed, but can be used with nansat,
           and Sentinel1ice, python tools from NERSC for working with satellite images
               * sudo port install py27-matplotlib-basemap
               * sudo port install py27-gdal
               * sudo port install py27-mahotas
               * sudo port install py27-scikit-image
               * sudo port install py27-scikit-learn
               * sudo port install py27-pil

         - NB Make sure macports (i)python is used (check with "which python"): in /opt/local/bin, do
            * sudo ln -s python2.7  python
            * sudo ln -s ipython2.7 ipython

- 4. matlab interface (mex) for fortran code
      - tested with OSX 10.8, matlab 2013a
      - tested with linux (johansen server), matlab 2012b

SETUP:

1. Go to matlab/main
   * run startup_local.m (set paths)
   * run WIM2d.m
   To test/look at results, go to matlab/main/test
   In matlab:
   * run startup_local.m (set paths)
   * run test_WIM2d.m    (NB edit 1st to set IO_OPT=1)

2. Grid setup:
   - Go to fortran/grid_setup
      * compile with make or make py
      * do "./grid_setup.sh" or "python grid_setup.py" respectively
         - or go into ipython and do "run grid_setup.py"
   - ipython script is easier to change (no need to recompile after changing grid)
      and makes a test plot automatically (plots the land mask)
   - can also look at the land mask in matlab with test/test_grid.m

   Compile/run main code:
   - Go to fortran/Build
      * compile with "make" or "make exe"
   - Go to fortran/run
      * run ./run_WIM2d.sh

   To test/look at results:
   - In matlab:
      * cd fortran/run/test
      * run startup_local.m (set paths)
      * run test_WIM2d_F.m
   - In ipython:
      * cd fortran/run
      * run test_run_WIM2d.py
         (NB change 'testing' option so that RUN_OPT=1
          is used when calling run_WIM2d.do_run)

3. Add the folfders fortran/bin & fortan/py_funs to the variable PYTHONPATH
   eg put this in .bashrc or .bash_profile:
   export WIM2D_PATH=/Users/timill/GITHUB-REPOSITORIES/WIM2d # define location of repo
   export PYTHONPATH=$PYTHONPATH:$WIM2D_PATH/fortran/bin:$WIM2D_PATH/fortran/py_funs

Grid setup:
   - same as for (2)

   Compile/run main code:
   - Go to fortran/Build
      * compile with "make py"
   - Go to fortran/run
     * In ipython: run ./run_WIM2d.py
          - This also has some plotting to look at results

   Can still go to fortran/run/test,
   and in matlab:
   * run startup_local.m (set paths)
   * run test_WIM2d_F.m  (NB edit 1st to set IO_OPT=1)

4. Grid setup:
   - same as for (2)

   Compile/run main code:
   - Go to fortran/Build
   * compile with "make mex2" (or "make mex" - this is more of a test,
     with no inputs/outputs - like running option (2) through matlab)
   - Go to fortran/run
     * In matlab: "startup_local; test_mex_io;"
          - This also has some plotting to look at results.
