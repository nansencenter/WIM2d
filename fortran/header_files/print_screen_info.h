      !!print text to screen during run

      print*,' '
      print*,'******************************************************'
      print*,'CFL number:                       ' ,CFL
      print*,'Maximum wave group velocity (m/s):' ,amax
      print*,'Time step (s):                    ' ,dt
      print*,'Number of time steps:             ' ,nt
      print*,'Time interval (h):                ' ,(nt*dt)/60.0/60.0
      print*,' '
      print*,'Grid dimensions:                  ' ,idm,jdm
      print*,'Spatial resolution (km):          ' ,dx/1.0e3,dy/1.0e3
      print*,'Extent of domain (km):            ' ,x_ext,y_ext
      print*,' '
      print*,'Minimum period (s):               ' ,1.0/freq_vec(nw)
      print*,'Maximum period (s):               ' ,1.0/freq_vec(1)
      print*,'Number of wave frequencies:       ' ,nw
      print*,'Number of wave directions:        ' ,ndir
      print*,'Directional resolution (degrees): ' ,360.0/(1.0*ndir)
      print*,'******************************************************'
      print*,' '
