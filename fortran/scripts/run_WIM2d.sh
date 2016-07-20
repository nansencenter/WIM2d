# place to put results

if [ ! -f infile_dirs.txt ]
then
   echo "infile_dirs.txt not present, copy template from infiles"
   exit
else
   lst=(`cat infile_dirs.txt`)
   indir=${lst[0]}
   outdir=${lst[1]}
   if [ ! -f $indir/wim_grid.a ];
   then
      echo "$indir/wim_grid.a missing";
      echo "Please go to ../grid_setup to set up grid"
      exit
   fi
   if [ ! -f $indir/wim_grid.b ];
   then
      echo "$indir/wim_grid.b missing";
      echo "Please go to ../grid_setup to set up grid"
      exit
   fi
fi

mkdir -p $outdir
mkdir -p $outdir/log
mkdir -p $outdir/binaries
mkdir -p $outdir/binaries/prog
rm -f $outdir/binaries/prog/*
rm -rf $outdir/figs/prog/*

# run model
w2d=$WIM2D_PATH
$w2d/fortran/bin/WIM2d.exec

if [ $# -eq 0 ]
then
   echo To make plots of progress:
   echo $w2d/fortran/tools/plot_prog.sh 0 $outdir
   echo To make plots of progress, and make a movie:
   echo $w2d/fortran/tools/plot_prog.sh 1 $outdir
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
echo $w2d/fortran/tools/plot_prog.sh $MK_MOV $outdir
$w2d/fortran/tools/plot_prog.sh $MK_MOV $outdir
