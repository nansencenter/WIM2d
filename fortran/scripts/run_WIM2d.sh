# place to put results
w2d=$WIM2D_PATH


# ==================================================================
# input/output dirs; check inputs
notpres=0
ifd="infile_dirs.txt"
if [ ! -f $ifd ]
then
   notpres=1
   echo " "
   echo "$ifd not present"
   echo "using defaults:"
   indir=grid
   outdir=out
   echo "grid in $indir"
   echo "outputs to $outdir"
   echo " "
   echo $indir>$ifd
   echo $outdir>>$ifd
else
   lst=(`cat infile_dirs.txt`)
   indir=${lst[0]}
   outdir=${lst[1]}
fi

if [ ! -f $indir/wim_grid.a ];
then
   echo "$indir/wim_grid.a missing";
   echo "Please run grid_setup.sh or grid_setup.py"
   echo "($w2d/fortran/grid_setup/scripts)"
   echo "to set up grid"
   exit
fi
if [ ! -f $indir/wim_grid.b ];
then
   echo "$indir/wim_grid.b missing";
   echo "Please run grid_setup.sh or grid_setup.py"
   echo "($w2d/fortran/grid_setup/scripts)"
   echo "to set up grid"
   exit
fi
# ==================================================================


# ==================================================================
# set up output dirs
rm -rf   "$outdir"
mkdir -p "$outdir"
mkdir -p "$outdir/binaries"
mkdir -p "$outdir/binaries/prog"
mkdir -p "$outdir/diagnostics"
mkdir -p "$outdir/diagnostics/global"
mkdir -p "$outdir/diagnostics/local"
# ==================================================================


# ==================================================================
# run model
bin/WIM2d.exec
if [ $notpres -eq 1 ]
then
   rm $ifd
fi
# ==================================================================


# ==================================================================
# post-processing
if [ $# -eq 0 ]
then
   echo To make plots of progress:
   echo $w2d/fortran/tools/plot_prog.sh 0 $outdir
   echo To make plots of progress, and make a movie:
   echo $w2d/fortran/tools/plot_prog.sh 1 $outdir
   echo
   echo Else, enter an argument 0 or 1 with this script ie
   echo "$w2d/fortran/scripts/run_WIM2d.sh 0 (plots, no movie)"
   echo "$w2d/fortran/scripts/run_WIM2d.sh 1 (plots, movie)"
   echo
   exit
else
   MK_MOV=$1
fi


# make plots + movie
echo $w2d/fortran/tools/plot_prog.sh $MK_MOV $outdir
$w2d/fortran/tools/plot_prog.sh $MK_MOV $outdir
# ==================================================================
