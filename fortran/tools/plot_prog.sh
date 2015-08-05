# make png figures from progress files
# (and of initial and final, if wanted)
# run from "../run" folder
tools="${WIM2D_PATH}/fortran/tools"
P=`pwd`

if [ $# -eq 0 ]
then
   echo "Usage:"
   echo "plot_prog.sh [1/0: do/don't make movie] [results directory (full path)]"
   exit
fi

MKMOV=$1
outdir=$2

# make png files from progress files
# (if they exist)
cd $outdir
bindir=$outdir/binaries/prog

afiles=($bindir/wim_prog*.a)
if [ ${#afiles[@]} -eq 0 ]
then
   echo "No prog files in $bindir"
   exit
fi

echo In `pwd`:
echo python $tools/plot_prog.py
python $tools/plot_prog.py

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
   echo cd $outdir/figs/prog
   echo "$tools/prog2mp4.sh Hs (or Dmax,taux,tauy)"
fi
