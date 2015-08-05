function [out,info] = fn_check_prog(outdir,n)

bdir  = [outdir,'/binaries/prog'];
flist = dir([bdir,'/wim_prog*.a']);

if length(flist)==0
   error([bdir,' is empty']);
else
   %% determine the number of digits
   f0 = flist(1).name;
   ss = strsplit(f0,'wim_prog');
   ss = strsplit(ss{2},'.a');
   ndig  = length(ss{1});
   fmt   = sprintf('%%%d.%dd',ndig,ndig);
   cts   = num2str(n,fmt);
end

%% get filename and read it
afile       = [bdir,'/wim_prog',cts,'.a'];
[out,info]  = fn_read_general_binary(afile);
