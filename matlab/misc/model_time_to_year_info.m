function tt = model_time_to_year_info(model_day_,model_seconds_)

%integer iyear           ! year         in YYYYMMDD format
%integer imonth          ! month        in YYYYMMDD format
%integer iday            ! day in month in YYYYMMDD format
%integer ihour           ! hours        in HH:MM:SS time format
%integer iminute         ! minute       in HH:MM:SS time format
%integer isecond         ! seconds      in HH:MM:SS time format

%real*8  model_seconds   ! seconds in day (exact)
%integer model_day       ! julian day since refyear/refmonth/refday

%!characters
%character(len=4) cyear        ! year   in YYYYMMDD format
%character(len=2) cmonth       ! month  in YYYYMMDD format
%character(len=2) cday         ! day    in YYYYMMDD format
%character(len=2) chour        ! hour   in HH:MM:SS time format
%character(len=2) cminute      ! minute in HH:MM:SS time format
%character(len=2) csecond      ! second in HH:MM:SS time format
%character(len=3) month_name   ! 'JAN' etc

%character(len=8) cdate        ! current date (YYYYMMDD format)
%character(len=8) ctime        ! current time (HH:MM:SS format)
%integer days_in_months(12)    ! total Days In Months 
%integer days_in_year          ! total Days In year

months_standard   = [31,28,31,30,31,30,31,31,30,31,30,31];
months_leapyear   = [31,29,31,30,31,30,31,31,30,31,30,31];

% reference time - model_day is relative to this
refyear  = 1900;
refmonth = 1;
refday   = 1;

% convert model day to date
[iyear,imonth,iday]  = juliantodate(model_day_,refyear,refmonth,refday);

%% output structure
tt = make_year_info(iyear,imonth,iday,model_seconds_,refyear,refmonth,refday);

return;


function tt = make_year_info(iyear,imonth,iday,dtime_seconds,refyear,refmonth,refday)

% ================================================
% store model day (relative to refyear/refmonth/refday)
tt.model_day   = datetojulian(iyear,imonth,iday,refyear,refmonth,refday);

% also store model time (seconds relative to refyear/refmonth/refday)
tt.model_seconds  = dtime_seconds;
% ================================================


% ================================================
%date
tt.iyear    = iyear;
tt.imonth   = imonth;
tt.iday     = iday;
%time
tt.ihour    = floor( dtime_seconds/3600.                         );
tt.iminute  = floor( (dtime_seconds-3600.*tt.ihour)/60.          );
tt.isecond  = round( dtime_seconds-3600.*tt.ihour-60.*tt.iminute );
% ================================================

% ================================================
%%matlab datevec for easy conversion to strings:
date_vector = [tt.iyear,tt.imonth,tt.iday,tt.ihour,tt.iminute,tt.isecond];

% date
tt.cyear    = datestr(date_vector,'yyyy');
tt.cmonth   = datestr(date_vector,'mm');
tt.cday     = datestr(date_vector,'dd');
tt.cdate    = [tt.cyear,tt.cmonth,tt.cday];

% time
tt.chour    = datestr(date_vector,'hh');
tt.cminute  = datestr(date_vector,'MM');
tt.csecond  = datestr(date_vector,'SS');
tt.ctime    = [tt.chour,':',tt.cminute,':',tt.csecond];
% ================================================ 


% ================================================ 
% month name can be useful sometimes
switch (tt.imonth)
case (1)
   tt.month_name  = 'JAN';
case (2)
   tt.month_name  = 'FEB';
case (3)
   tt.month_name  = 'MAR';
case (4)
   tt.month_name  = 'APR';
case (5)
   tt.month_name  = 'MAY';
case (6)
   tt.month_name  = 'JUN';
case (7)
   tt.month_name  = 'JUL';
case (8)
   tt.month_name  = 'AUG';
case (9)
   tt.month_name  = 'SEP';
case (10)
   tt.month_name  = 'OCT';
case (11)
   tt.month_name  = 'NOV';
case (12)
   tt.month_name  = 'DEC';
end
% ================================================ 


% ================================================ 
% extra info
tt.days_in_year   = daysinyear  (iyear);
tt.days_in_months = monthsinyear(iyear);
% ================================================ 

return


function jday = datetojulian(year,month,day,ryear,rmonth,rday)


sum_days = 0;
for iyear=ryear:year
   sum_days=sum_days+daysinyear(iyear);
end

% Subtract from start of ref year to reference date
months   = monthsinyear(ryear);
sum_days = sum_days - sum(months(1:rmonth)) +...
            + months(rmonth) - rday + 1;

% Subtract from end date in last year to end of year
months   = monthsinyear(year);
sum_days = sum_days - sum(months(month:12)) + day -1;

jday  = sum_days;
return;


function [iyear,imonth,iday] =  juliantodate(jday,ryear,rmonth,rday)

sum_days = 0;

% Subtract from start of ref year to reference date
sum_days = sum_days+daysinyear(ryear);
months   = monthsinyear(ryear);
sum_days = sum_days - sum(months(1:rmonth)) +...
            + months(rmonth) - rday + 1;

% Add years until beyond julian day
iyear = ryear+1;
while (sum_days<jday)
   sum_days = sum_days+daysinyear(iyear);
   iyear    =  iyear+1;
end
if (sum_days>jday)
   iyear = iyear-1;
end

imonth   = 12;
months   = monthsinyear(iyear);
while (sum_days>jday)
   sum_days = sum_days-months(imonth);
   imonth   = imonth-1;
end
imonth   = mod(imonth,12)+1;

iday=1;
while (sum_days<jday)
   sum_days = sum_days+1;
   iday     = iday+1;
end

return


function days = daysinyear(iyear)

if (mod(iyear,4)==0 )
   if (mod(iyear,400)==0)
      days  = 366;
   elseif (mod(iyear,100)==0)
      days  = 365;
   else
      days  = 366;
   end
else
   days  = 365;
end
return;


function months = monthsinyear(iyear)

months_standard   = [31,28,31,30,31,30,31,31,30,31,30,31];
months_leapyear   = [31,29,31,30,31,30,31,31,30,31,30,31];

if (daysinyear(iyear)==366)
   months   = months_leapyear;
else
   months   = months_standard;
end
return
