bin/prog_wim2d_cpp.exec --config-file=wim.cfg

echo ' '
echo 'Make plots of progress with:'
echo '$WIM2D_PATH/fortran/tools/plot_prog.sh 0 out_cpp'

echo ' '
echo 'Make plots of progress & movie of Hs with:'
echo '$WIM2D_PATH/fortran/tools/plot_prog.sh 1 out_cpp'

echo ' '
echo 'Make movies of other variables with:'
echo '$WIM2D_PATH/fortran/tools/prog2mp4.sh'
