      !! MAIN INTEGRATION
      !! - Common to WIM2d_run & WIM2d_run_io subroutines

      !nt = 13!stop at a certain point for testing
      !nt = 10*nt!extend to see true steady state
      reps  = 10!frequency of dumping

      do n=1,nt

         critter  = (DO_CHECK_PROG.eq.1).and.(mod(n-1,reps).eq.0)

         if (critter) then
            !!dump fields to binary
            call check_prog(SOLVER,n-1,0,outdir)
         end if

         print*,'n/nt',n,nt
         call wim_step(SOLVER,ADV_DIM)

      end do!! n - finish time stepping

      !!dump final output
      call check_prog(SOLVER,n-1,1,outdir)

      MIZ_MASK = 0.0
      where ((dfloe.gt.0.0).and.(dfloe.lt.dfloe_pack_init))
         MIZ_MASK = 1.0
      end where

      !!Dmax range in MIZ
      Dmax_max = maxval(MIZ_MASK*dfloe)
      tmp1     = dfloe+1.0e4*(1.0-MIZ_MASK)
      Dmax_min = minval(tmp1)

      W_MIZ = 0.0
      if (.true.) then
         !! NB this definition won't always work
         !! for all configurations
         nMIZ  = sum(MIZ_MASK(:,1))
         W_MIZ = nMIZ*dx/1.0e3
         !!
         print*,' '
         print*,'MIZ width (km)',W_MIZ
      endif
      !!
      print*,'Dmax range in MIZ (m): ',Dmax_min,Dmax_max

      !!range of stresses
      taux_max = maxval(tau_x)
      taux_min = minval(tau_x)
      tauy_max = maxval(tau_y)
      tauy_min = minval(tau_y)
      print*,'tau_x range (Pa)',taux_min,taux_max
      print*,'tau_y range (Pa)',tauy_min,tauy_max
      print*,' '

      et = etime(tictoc)!!finish etime
      print*,' '
      print*,'******************************************************'
      print*,'run of 2d WIM2d finished'
      !print *,'Elapsed time (s):',et
      print *,'Elapsed time (min):',et/60.0                             !&
!    &       ,', user:', tictoc(1),                                     &
!    &       ,', sys:' ,tictoc(2)
      print*,'******************************************************'
      print*,' '
