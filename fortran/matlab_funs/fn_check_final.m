function out = fn_check_final(outdir)

bdir  = [outdir,'/binaries'];
afile = [bdir,'/wim_out.a'];
out   = fn_read_general_binary(afile);
