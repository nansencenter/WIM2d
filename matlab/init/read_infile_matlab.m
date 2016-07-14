function params   = read_infile_matlab(infile,verbosity)

if ~exist('infile','var')
   infile   = 'infile_matlab.txt';
end
if ~exist('verbosity','var')
   verbosity   = 1;
end
infile_version = 8;%%latest infile version

if ~exist(infile)
   %% now need infile to run code
   error([infile,' not present - get example from "matlab/main/infiles" directory'])
else
   if verbosity
      disp('********************************************************')
      disp('reading options from infile:')
      disp(infile)
      disp('********************************************************')
      disp(' ')
   end
   fid   = fopen(infile);

   %%check infile version:
   infile_version_   = read_next(fid);
   if infile_version_~=infile_version
      error(['Infile version number is: ',num2str(infile_version_),' - should be: ',num2str(infile_version)]);
   end

   %%read in rest of variables:
   while ~feof(fid)
      [x,name] = read_next(fid,verbosity);
      if ~isempty(x)
         cmd   = ['params.',name,' = ',num2str(x),';'];
         if verbosity
            disp(cmd);
         end
         eval(cmd);
      end
   end
   fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,name]  = read_next(fid,verbosity)
%% read next line in text file

lin   = strtrim(fgets(fid));  %% trim leading white space
lin2  = strsplit(lin);        %%split using spaces
x     = lin2{1};              %%get 1st thing in line

if strcmp(x,'')
   % blank line
   if verbosity
      disp(' ');
   end
   x     = [];
   name  = [];
elseif strcmp(x,'#')
   % comment
   if verbosity
      disp(lin);
   end
   x     = [];
   name  = [];
else
   % proper variable
   x     = str2num(x);
   name  = lin2{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
