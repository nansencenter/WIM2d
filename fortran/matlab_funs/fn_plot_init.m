function fn_plot_init(grid_prams,outdir)

figure;
fn_fullscreen;
[grid_prams,ice_fields,wave_fields] = fn_check_init(outdir);
fn_plot_ice(grid_prams,ice_fields);

figure;
fn_fullscreen;
fn_plot_waves(grid_prams,wave_fields);
