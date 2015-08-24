function params   = read_infile_matlab(infile)

if ~exist('infile','var')
   infile   = 'infile_matlab.txt';
end
infile_version = 7;%%latest infile version

if ~exist(infile)
   %% now need infile to run code
   error([infile,' not present - get example from "matlab/main/infiles" directory'])
else
   disp('********************************************************')
   disp('reading options from infile:')
   disp(infile)
   disp('********************************************************')
   disp(' ')
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
         cmd   = ['params.',name,' = ',num2str(x),';'];
         disp(cmd);
         eval(cmd);
      end
   end
   fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,name]  = read_next(fid)
%% read next line in text file

lin   = strtrim(fgets(fid));  %% trim leading white space
lin2  = strsplit(lin);        %%split using spaces
x     = lin2{1};              %%get 1st thing in line

if strcmp(x,'')
   % blank line
   disp(' ');
   x     = [];
   name  = [];
elseif strcmp(x,'#')
   % comment
   disp(lin);
   x     = [];
   name  = [];
else
   % proper variable
   x     = str2num(x);
   name  = lin2{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=find_num(txt)

x0 = strsplit(txt,'!');
x  = str2num(x0{1});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
