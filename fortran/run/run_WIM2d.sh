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

# make plots + movie
if [ $# -eq 0 ]
then
   MK_MOV=1 # default is make movie
else
   MK_MOV=$1
fi

echo $w2d/fortran/tools/plot_prog.sh $MK_MOV out
$w2d/fortran/tools/plot_prog.sh $MK_MOV out
