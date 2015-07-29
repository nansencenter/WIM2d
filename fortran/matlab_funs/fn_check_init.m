function [grid_prams,ice_fields,wave_fields] = fn_check_init(outdir)

bdir  = [outdir,'/binaries'];
afile = [bdir,'/wim_init.a'];

grid_prams  = fn_get_grid(bdir);
out         = fn_read_general_binary(afile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ice_fields.cice      = out.cice;
ice_fields.hice      = out.hice;
ice_fields.Dmax      = out.Dmax;
ice_fields.ICE_MASK  = 1.0*(ice_fields.cice>0);
%%
wave_fields.Hs          = out.Hs ;
wave_fields.Tp          = out.Tp ;
wave_fields.mwd         = out.mwd;
wave_fields.WAVE_MASK   = 1.0*(wave_fields.Hs>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
