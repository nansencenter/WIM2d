# make png figures from progress files
# (and of initial and final, if wanted)
# run from "../run" folder
tools="${WIM2D_PATH}/neXtSIMcoupling/tools"
tools2="${WIM2D_PATH}/fortran/tools"
python $tools/plot_prog.py

if [ $# -eq 1 ]
then
   MKMOV=$1 # 1, make movies; 0, don't
else
   MKMOV=0 # default is not make movie
fi

if [ $MKMOV -eq 1 ]
then
   # make movies of these variables:
   vbl_list="Hs Dmax taux tauy"
   cd test_outputs


   for vbl in $vbl_list
   do
      $tools2/prog2mp4.sh $vbl
   done
fi
