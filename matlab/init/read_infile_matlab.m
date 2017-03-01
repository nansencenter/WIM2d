function params   = read_infile_matlab(infile)

if ~exist('infile','var')
   infile   = 'infile_matlab.txt';
end
infile_version = 14;%%latest infile version

if ~exist(infile)
   %% now need infile to run code
   error([infile,' not present - get example from "matlab/main/infiles" directory'])
else
   fid   = fopen(infile);

   %%check infile version:
   infile_version_   = read_next(fid);
   if infile_version_~=infile_version
      error(['Infile version number is: ',num2str(infile_version_),' - should be: ',num2str(infile_version)]);
   end

   %%read in rest of variables:
   while ~feof(fid)
      [x,name] = read_next(fid);
      if ~isempty(x)
         params.(name)  = x;
      end
   end
   fclose(fid);

   fields   = {'start_year',...
               'start_month',...
               'start_day',...
               'start_hour',...
               'start_minute',...
               'start_second'};
   for j=1:6
      fld            = fields{j};
      date_vector(j) = params.(fld);
      params         = rmfield(params,fld);
   end
   params.year_info_start  = datevec2year_info(date_vector);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,name]  = read_next(fid)
%% read next line in text file

lin   = strtrim(fgets(fid));  %% trim leading white space
lin2  = strsplit(lin);        %%split using spaces
xs    = lin2{1};              %%get 1st thing in line (as string)

if strcmp(xs,'')
   % blank line
   x     = [];
   name  = [];
elseif strcmp(xs,'#')
   % comment
   x     = [];
   name  = [];
else
   % proper variable
   [x,status]  = str2num(xs);
   if status==0
      %% numerical conversion failed
      %% - leave as string
      x  = xs;
   end
   name  = lin2{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
