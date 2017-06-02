WIM2d
=====

2d waves-in-ice module
Stand-alone version(s) for use on a single computer.

1. Matlab version

# ========================================================
A. cd matlab/main
matlab &
# ========================================================


# ========================================================
B. get input file & edit if you want
- can change options like grid size (nx,ny),
  number of frequencies and directions (nw,ndir),
  plotting options, save binary outputs or not
- MEX_OPT=-1 gives a slight speedup compared to MEX_OPT=0
   (uses advection code from a mex file)
- MEX_OPT>0 use mex interfaces to fortran code (more work to compile these)
cp infiles/infile_matlab.save_matlab_figures.txt
infiles/infile_matlab.txt
# ========================================================


# ========================================================
C. run in matlab:
> startup_local
> run_WIM2D

You should get some results in out_m or out_mexadv
- figs &/or binary files, depending on infile options
# ========================================================

2. mex interface to fortran code
- tested with OSX 10.10 & 10.8, matlab 2013a
- tested with linux (johansen server), matlab 2012b
- to compile follow README_fortran.txt (TODO write this!)
