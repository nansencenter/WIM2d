# make png figures from progress files
# (and of initial and final, if wanted)
# run from "../run" folder
tools="${WIM2D_PATH}/fortran/tools"
P=`pwd`
outdir=$P/out_io # root results directory

if [ $# -eq 2 ]
then
   MKMOV=$1 # 1, make movies; 0, don't
   outdir=$2
elif [ $# -eq 1 ]
then
   MKMOV=$1 # 1, make movies; 0, don't
else
   MKMOV=0 # default is not make movie
fi

# make png files
cd $outdir
python $tools/plot_prog.py
cd $P

if [ $MKMOV -eq 1 ]
then
   # make movies of these variables:
   vbl_list="Hs Dmax taux tauy"
   cd $outdir/figs/prog

   for vbl in $vbl_list
   do
      echo $tools/prog2mp4.sh $vbl
      $tools/prog2mp4.sh $vbl
   done
fi
