# make png figures from progress files
# (and of initial and final, if wanted)
# run from "../run" folder
tools="${WIM2D_PATH}/fortran/tools"

if [ $# -eq 0 ]
then
   echo "Usage:"
   echo "plot_prog.sh [1/0: do/don't make movie] [root results directory ie with binaries, log etc]"
   echo "Optional arg's: PLOT_PROG PLOT_INIT PLOT_FINAL"
   exit
fi

MKMOV=$1
outdir=$2
if [ $# -ge 3 ]
then
   PLOT_PROG=$3
else
   PLOT_PROG=1
fi

if [ $# -ge 4 ]
then
   PLOT_INIT=$4
else
   PLOT_INIT=1
fi

if [ $# -ge 5 ]
then
   PLOT_FINAL=$5
else
   PLOT_FINAL=1
fi

# make png files from progress files
# (if they exist)
cd $outdir
outdir=`pwd` # change to full path
bindir=$outdir/binaries/prog

afiles=($bindir/wim_prog*.a)
if [ ${#afiles[@]} -eq 0 ]
then
   echo "No prog files in $bindir"
   exit
fi

echo In `pwd`:
rm -rf figs/prog/*
echo python $tools/plot_prog.py --outdir=$outdir --prog=$PLOT_PROG --init=$PLOT_INIT --final=$PLOT_FINAL
python $tools/plot_prog.py --outdir=$outdir --prog=$PLOT_PROG --init=$PLOT_INIT --final=$PLOT_FINAL

if [ $MKMOV -eq 1 ]
then
   # make movies of these variables:
   # vbl_list="Hs Dmax taux tauy"
   vbl_list="Hs"
   cd $outdir/figs/prog

   for vbl in $vbl_list
   do
      echo $tools/prog2mp4.sh $vbl
      $tools/prog2mp4.sh $vbl
   done
else
   echo To make movie
   echo "$tools/prog2mp4.sh Hs" $outdir/figs/prog
   echo "(can also use eg Dmax,taux,tauy)"
fi
