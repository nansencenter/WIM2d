module mod_year_info
implicit none
   type year_info
      integer iyear           ! year         in YYYYMMDD format
      integer imonth          ! month        in YYYYMMDD format
      integer iday            ! day in month in YYYYMMDD format
      integer ihour           ! hours        in HH:MM:SS time format
      integer iminute         ! minute       in HH:MM:SS time format
      integer isecond         ! seconds      in HH:MM:SS time format

      real*8  dsecond         ! seconds in day (exact)
      integer model_day       ! julian day since refyear/refmonth/refday
      real*8  model_time      ! seconds    since refyear/refmonth/refday

      !characters
      character(len=4) cyear        ! year   in YYYYMMDD format
      character(len=2) cmonth       ! month  in YYYYMMDD format
      character(len=2) cday         ! day    in YYYYMMDD format
      character(len=2) chour        ! hour   in HH:MM:SS time format
      character(len=2) cminute      ! minute in HH:MM:SS time format
      character(len=2) csecond      ! second in HH:MM:SS time format
      character(len=3) month_name   ! 'JAN' etc

      character(len=8) cdate        ! current date (YYYYMMDD format)
      character(len=8) ctime        ! current time (HH:MM:SS format)
      integer days_in_months(12)    ! total Days In Months 
      integer days_in_year          ! total Days In year
   end type year_info

   integer, dimension(12),parameter :: months_standard = &
      (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer, dimension(12),parameter :: months_leapyear = &
      (/31,29,31,30,31,30,31,31,30,31,30,31/)
   integer, dimension(12),parameter :: months_365 = &
      (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer, dimension(12),parameter :: months_366 = &
      (/31,29,31,30,31,30,31,31,30,31,30,31/)

   ! reference time - model_day is relative to this
   integer,parameter :: refyear  = 1900
   integer,parameter :: refmonth = 1
   integer,parameter :: refday   = 1

   public :: year_info,model_time_to_year_info  &
               ,make_year_info                  &
               ,refyear,refmonth,refday

contains


subroutine model_time_to_year_info(tt,model_time)
   ! output type(year_info)   :: tt
   ! input double             :: model_time    - time in seconds since reftime
   implicit none


   type(year_info), intent(out)  :: tt
   real*8, intent(in)            :: model_time

   integer  :: julian_day,iyear,imonth,iday
   real*8   :: dtime_seconds

   julian_day     = floor(model_time/(24.*3600.))     !julian day since reftime
   dtime_seconds  = model_time-julian_day*(24.*3600.) !seconds in the day
   
   call juliantodate(iyear,imonth,iday                      &
                     ,julian_day,refyear,refmonth,refday)

   call make_year_info(tt,iyear,imonth,iday,dtime_seconds)

end subroutine model_time_to_year_info


subroutine make_year_info(tt,iyear,imonth,iday,dtime_seconds)
   implicit none

   integer,intent(in)   :: iyear,imonth,iday
   real*8,intent(in)    :: dtime_seconds
   type(year_info), intent(out)  :: tt

   real*8   :: dsecs


   ! ================================================
   !date
   tt%iyear    = iyear
   tt%imonth   = imonth
   tt%iday     = iday
   write(tt%cyear ,'(i4.4)') tt%iyear
   write(tt%cmonth,'(i2.2)') tt%imonth
   write(tt%cday  ,'(i2.2)') tt%iday
   tt%cdate = tt%cyear//tt%cmonth//tt%cday
   ! ================================================


   ! ================================================
   ! also store model day (relative to refyear/refmonth/refday)
   tt%model_day   = datetojulian(iyear,imonth,iday,refyear,refmonth,refday)

   ! also store model time (seconds relative to refyear/refmonth/refday)
   tt%model_time  = 3600*24*tt%model_day+tt%dsecond
   ! ================================================


   ! ================================================
   !time
   tt%dsecond  = dtime_seconds
   tt%ihour    = floor( dtime_seconds/3600.                           )
   tt%iminute  = floor( (dtime_seconds-3600.*tt%ihour)/60.            )
   dsecs       = floor( (dtime_seconds-3600.*tt%ihour-60.*tt%iminute) )
   tt%isecond  = floor( dsecs )
   write(tt%chour  ,'(i2.2)') tt%ihour
   write(tt%cminute,'(i2.2)') tt%iminute
   write(tt%csecond,'(i2.2)') tt%isecond
   tt%ctime = tt%chour//':'//tt%cminute//':'//tt%csecond
   ! ================================================ 


   ! ================================================ 
   ! month name canbe useful sometimes
   select case (tt%imonth)
   case (1)
      tt%month_name  = 'JAN'
   case (2)
      tt%month_name  = 'FEB'
   case (3)
      tt%month_name  = 'MAR'
   case (4)
      tt%month_name  = 'APR'
   case (5)
      tt%month_name  = 'MAY'
   case (6)
      tt%month_name  = 'JUN'
   case (7)
      tt%month_name  = 'JUL'
   case (8)
      tt%month_name  = 'AUG'
   case (9)
      tt%month_name  = 'SEP'
   case (10)
      tt%month_name  = 'OCT'
   case (11)
      tt%month_name  = 'NOV'
   case (12)
      tt%month_name  = 'DEC'
   end select
   ! ================================================ 


   ! ================================================ 
   !!extra info
   tt%days_in_year   = daysinyear  (iyear)
   tt%days_in_months = monthsinyear(iyear)
   ! ================================================ 


end subroutine make_year_info


integer function datetojulian(year,month,day,ryear,rmonth,rday)
   implicit none
   integer, intent(in) :: year,month,day,ryear,rmonth,rday
   integer :: iyear,sum_days,months(12)

   sum_days=0
   do iyear=ryear,year
      sum_days=sum_days+daysinyear(iyear)
      !print *,sum_days
   enddo

   ! Subtract from start of ref year to reference date
   months=monthsinyear(ryear)
   sum_days=sum_days          &
      - sum(months(1:rmonth)) &
      + months(rmonth) - rday + 1
   !print *,sum_days



   ! Subtract from end date in last year to end of year
   months=monthsinyear(year)
   sum_days=sum_days          &
      - sum(months(month:12)) &
      + day -1
   !print *,sum_days

   datetojulian=sum_days
end  function datetojulian


subroutine juliantodate(jday,year,month,day,ryear,rmonth,rday)
   implicit none
   integer, intent(in) :: jday,ryear,rmonth,rday
   integer, intent(out):: year,month,day
   integer :: iyear,sum_days,months(12),imonth,iday

   sum_days=0

   ! Subtract from start of ref year to reference date
   sum_days=sum_days+daysinyear(ryear)
   !print *,sum_days
   months=monthsinyear(ryear)
   sum_days=sum_days          &
      - sum(months(1:rmonth)) &
      + months(rmonth) - rday + 1
   !print *,sum_days


   ! Add years until beyond julian day
   iyear=ryear+1
   do while (sum_days<jday)
      sum_days=sum_days+daysinyear(iyear)
      iyear=iyear+1
      !print *,sum_days
   enddo
   if (sum_days>jday) then
      iyear=iyear-1
   end if

   imonth=12
   months=monthsinyear(iyear)
   do while (sum_days>jday)
      sum_days=sum_days-months(imonth)
      !print *,sum_days
      imonth=imonth-1
   enddo
   imonth=mod(imonth,12)+1

   iday=1
   do while (sum_days<jday)
      sum_days=sum_days+1
      iday=iday+1
   end do
      
   year=iyear
   month=imonth
   day=iday
end  subroutine juliantodate


integer function daysinyear(iyear)
   implicit none
   integer, intent(in) :: iyear
   !print *,iyear,mod(iyear,4)
   if (mod(iyear,4)==0 ) then
      if (mod(iyear,400)==0) then
         daysinyear=366
      elseif (mod(iyear,100)==0) then
         daysinyear=365
      else
         daysinyear=366
      end if
   else
      daysinyear=365
   end if
end function daysinyear


function monthsinyear(iyear)
   implicit none
   integer :: monthsinyear(12)
   integer, intent(in) :: iyear

   if (daysinyear(iyear)==366) then
      monthsinyear=months_leapyear
   else
      monthsinyear=months_standard
   end if

end function



end module mod_year_info
