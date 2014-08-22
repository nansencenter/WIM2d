WIM2d
=====

2d waves-in-ice module
Want to have (a) stand-alone version(s) for use on a laptop/server.

In progress.
So far, just have WENO advection scheme imported from NERSC-HYCOM.

PLANS:
- matlab -> good for initial testing but memory is
  limited (ie might not be able to use many
  directions/frequencies)
- version in c -> Elasto-brittle rheology
- fortran -> testing features to go into NERSC-HYCOM or Wavewatch/WAM
- python? very easy to set up parallisation with that (on
  a server) so could be quite efficient. Also can use
  more memory than matlab.
