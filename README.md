WIM2d
=====

2d waves-in-ice module
Stand-alone version for use on a single computer.

Versions:
- 1. matlab
- 2. fortran
- 3. ipython interface for fortran code

SETUP:
1. Go to matlab/main
   * run startup_local.m (set paths)
   * run WIM2d.m
   To test/look at results, go to matlab/main/test
   In matlab:
   * run startup_local.m (set paths)
   * run test_WIM2d.m

2. Go to fortran/Build
   * compile with "make"
   Go to fortran/run
   * run ./run_WIM2d.sh
   To test/look at results, go to fortran/run/test
   In matlab:
   * run startup_local.m (set paths)
   * run test_WIM2d_F.m
   TODO:
   Enable looking at results in ipython.

3. Go to fortran/Build
   * compile with "make py"
   Go to fortran/run
   In ipython:
   * run ./run_WIM2d.py
     - This has some plotting to look at results,
       but it is not as sophisticated/organised
       as the matlab test code yet.

   Can still go to fortran/run/test, and
   in matlab:
   * run startup_local.m (set paths)
   * run test_WIM2d_F.m
