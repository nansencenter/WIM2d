# make png figures from progress files
# (and of initial and final, if wanted)
# run from "../run" folder
tools="${WIM2D_PATH}/fortran/tools"

if [ $# -eq 0 ]
then
   echo "Usage:"
   echo "plot_prog.sh [1/0: do/don't make movie] [root results directory ie with binaries, log etc]"
   echo "Optional arg's: PLOT_PROG PLOT_INIT PLOT_FINAL PLOT_INC vname1 vname2 vname3..."
   exit
fi

MKMOV=$1
outdir=$2
PLOT_PROG=1
PLOT_INIT=1
PLOT_FINAL=1
PLOT_INC=1

n=0
vlist=""
nvbl=0
for var in "$@"
do
   n=$((n+1))
   # echo $n $var
   if [ $n -lt 3 ]
   then
      continue
   elif [ $n -eq 3 ]
   then
      PLOT_PROG=$var
   elif [ $n -eq 4 ]
   then
      PLOT_INIT=$var
   elif [ $n -eq 5 ]
   then
      PLOT_FINAL=$var
   elif [ $n -eq 6 ]
   then
      PLOT_INC=$var
   else
      vlist="$vlist $var"
      vbls[$nvbl]=$var
      nvbl=$((nvbl+1))
   fi
done

# make png files from progress files
# (if they exist)
cd $outdir
outdir=`pwd` # change to full path

echo In `pwd`:
if [ "$vlist" != "" ]
then
   # overwrite
   for vbl in $vlist
   do
      if [ "$PLOT_PROG" -eq 1 ]
      then
         vdir=figs/prog/$vbl
         if [ -d "$vdir" ]
         then
            echo rm -rf $vdir
            rm -rf $vdir
         fi
      fi
      if [ "$PLOT_INIT" -eq 1 ]
      then
         vdir=figs/init/$vbl
         if [ -d "$vdir" ]
         then
            echo rm -rf $vdir
            rm -rf $vdir
         fi
      fi
      if [ "$PLOT_FINAL" -eq 1 ]
      then
         vdir=figs/final/$vbl
         if [ -d "$vdir" ]
         then
            echo rm -rf $vdir
            rm -rf $vdir
         fi
      fi
      if [ "$PLOT_INC" -eq 1 ]
      then
         vdir=figs/incwaves/$vbl
         if [ -d "$vdir" ]
         then
            echo rm -rf $vdir
            rm -rf $vdir
         fi
      fi
   done
fi

echo python $tools/plot_prog.py --outdir=$outdir --prog=$PLOT_PROG --init=$PLOT_INIT --final=$PLOT_FINAL  --incwaves=$PLOT_INC $vlist
python $tools/plot_prog.py --outdir=$outdir --prog=$PLOT_PROG --init=$PLOT_INIT --final=$PLOT_FINAL  --incwaves=$PLOT_INC $vlist

if [ $MKMOV -eq 1 ]
then
   # make movies of these variables:
   # vbl_list="Hs Dmax taux tauy"
   vbl_list="Hs"
   cd $outdir/figs/prog

   for vbl in $vbl_list
   do
      ok=1
      if [ $nvbl -gt 0 ]
      then
         ok=0
         for v in ${vbls[@]}
         do
            if [ $v == $vbl ]
            then
               ok=1
               break
            fi
         done
      fi
      if [ $ok -eq 1 ]
      then
         echo $tools/prog2mp4.sh $vbl
         $tools/prog2mp4.sh $vbl
      fi
   done
else
   echo To make movie
   echo "$tools/prog2mp4.sh Hs" $outdir/figs/prog
   echo "(can also use eg Dmax,tau_x,tau_y)"
fi
