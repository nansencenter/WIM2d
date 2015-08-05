%% fn_bfile_info.m
%% Author: Timothy Williams
%% Date: 20150805, 09:20:58 CEST

function [info,vlist]   = fn_bfile_info(bfile)

bid   = fopen(bfile);

do_vlist = 0;
while ~do_vlist
   lin   = strtrim(fgetl(bid));
   if ~strcmp(lin,'')
      %%ignore blank lines
      slin  = strsplit(lin);
      if strcmp(slin{1},'Record')&strcmp(slin{2},'number')
         do_vlist = 1;
      else
         info.(slin{2}) = str2num(slin{1});
      end
   end
end

Nrec  = info.Nrecs;
for n=1:Nrec
   lin      = strtrim(fgetl(bid));
   slin     = strsplit(lin);
   vlist{n} = slin{2};
end

fclose(bid);
