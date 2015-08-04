function [out,info] = fn_check_final(outdir)

bdir        = [outdir,'/binaries'];
afile       = [bdir,'/wim_out.a'];
[out,info]  = fn_read_general_binary(afile,1);
