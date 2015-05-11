      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! program to run WIM2d without inpus/outputs
      !! (only inputs: infile.txt; only outputs: saved binary files)
      program main

      use mod_wim2d_run, only: wim2d_run

      implicit none

      call wim2d_run()

      end program main
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
