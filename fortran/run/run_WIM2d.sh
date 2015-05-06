# place to put results
outdir="out"
mkdir -p $outdir
mkdir -p $outdir/log
mkdir -p $outdir/binaries
mkdir -p $outdir/binaries/prog
rm -f $outdir/binaries/prog/*

# run model
w2d=$WIM2D_PATH
$w2d/fortran/bin/WIM2d.exe

if [ $# -eq 0 ]
then
   echo To make plots of progress:
   echo $w2d/fortran/tools/plot_prog.sh 0 out
   echo To make plots of progress, and make a movie:
   echo $w2d/fortran/tools/plot_prog.sh 1 out
   echo
   echo Else, enter an argument 0 or 1 with this script ie
   echo "./run_WIM2d.sh 0 (plots, no movie)"
   echo "./run_WIM2d.sh 1 (plots, movie)"
   echo
   exit
else
   MK_MOV=$1
fi

# make plots + movie
echo $w2d/fortran/tools/plot_prog.sh $MK_MOV out
$w2d/fortran/tools/plot_prog.sh $MK_MOV out
