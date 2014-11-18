      !!print text to log file during run

      print*,'Printing info to log file: ',log_file
      open(unit=3,file=trim(log_file),status = 'replace')

      write(3,'(a)'),' '
      write(3,'(a)'),'*************************************************'
      write(3,'(a)'),'Main parameters:'
      write(3,'(a,i2.2)'),'SOLVER:                           ' ,SOLVER
      write(3,'(a,i2.2)'),'GRID_OPT:                         ' ,GRID_OPT
      write(3,'(a)'),'*************************************************'

      write(3,'(a)'),' '
      write(3,'(a)'),'*************************************************'
      write(3,'(a)'),'WIM parameters:'
      write(3,'(a,f4.2)'),'Brine volume fraction:      ',vbf
      write(3,'(a,e10.3)'),'Youngs modulus (Pa):        ',young
      write(3,'(a,e10.3)'),'Flexural strength (Pa):     ',sigma_c
      write(3,'(a,e10.3)'),'Breaking strain:            ',epsc
      write(3,'(a,f5.2)'),'Damping (Pa.s/m):           ',visc_rp
      write(3,'(a)'),'*************************************************'
      print*,' '

      write(3,'(a)'),' '
      write(3,'(a)'),'*************************************************'
      write(3,'(a)'),'Other Parameters:'
      write(3,'(a,f6.1)'),'Time step (s):                    ' ,dt
      write(3,'(a,f4.3)'),'CFL number:                       ' ,CFL
      write(3,'(a,f5.2)'),'Maximum wave group velocity (m/s):' ,amax
      write(3,'(a,f6.1)'),'Time step (s):                    ' ,dt
      write(3,'(a,i4.4)'),'Number of time steps:             ' ,nt
      write(3,'(a,f5.2)'),'Time interval (h):                '          &
     &    ,(nt*dt)/60.0/60.0
      write(3,'(a)'),'*************************************************'

      write(3,'(a)'),' '
      write(3,'(a)'),'*************************************************'
      write(3,'(a,i4.4,a,i4.4)'),'Grid dimensions:                  '   &
     &                        ,idm,' ',jdm
      write(3,'(a,f4.1,a,f4.1)'),'Spatial resolution (km):          '   &
     &                        ,dx/1.0e3,' ',dy/1.0e3
      write(3,'(a,f6.1,a,f6.1)'),'Extent of domain (km):            '   &
     &                        ,x_ext,' ',y_ext

      write(3,'(a)'),' '
      write(3,'(a,f5.2)'),'Minimum period (s):               '          &
     &                   ,1.0/freq_vec(nw)
      write(3,'(a,f5.2)'),'Maximum period (s):               '          &
     &                   ,1.0/freq_vec(1)
      write(3,'(a,i2.2)'),'Number of wave frequencies:       ' , nw
      write(3,'(a,i3.3)'),'Number of wave directions:        ' , ndir
      write(3,'(a,f5.2)'),'Directional resolution (degrees): '
     &                   ,360.0/(1.0*ndir)
      write(3,'(a)'),'*************************************************'

      write(3,'(a)'),' '
      write(3,'(a)'),'*************************************************'
      write(3,'(a)'),'Diagnostics:'
      write(3,'(a,f6.1)'),'MIZ width (km): ',W_MIZ
      write(3,'(a,f6.1,a,f6.1)'),'Dmax range in MIZ (m): '              &
     &                          ,Dmax_min,' ',Dmax_max
      write(3,'(a,e10.3,a,e10.3)'),'tau_x range (Pa): '                 &
     &                          ,taux_min,' ',taux_max
      write(3,'(a,e10.3,a,e10.3)'),'tau_y range (Pa): '                 &
     &                          ,tauy_min,' ',tauy_max
      write(3,'(a)'),'*************************************************'

      write(3,'(a,f6.1)'),' '
      write(3,'(a)'),'*************************************************'
      write(3,'(a,f7.1)'),'Elapsed time (min):',et/60.0
      write(3,'(a)'),'*************************************************'
      write(3,'(a)'),' '

      !!close file
      close(3)
