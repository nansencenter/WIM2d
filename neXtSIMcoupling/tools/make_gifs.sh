if [ $# -lt 2 ]
then
   echo "Usage make_gifs.sh [dir with png files] [mesh_opt]"
   echo "mesh_opt = 1 (save output to gifs_mesh) or 0 (save output to gifs_grid)"
   exit
fi

KERNEL=`uname -s`

if [ $KERNEL == Darwin ]
then
   # install coreutils to get greadlink:
   # sudo port install coreutils
   indir=`greadlink -f $1` # want absolute path
   mesh=$2
else
   indir=`readlink -f $1` # want absolute path
   mesh=$2
fi

if [ $mesh -eq 1 ]
then
   indir1=$indir/figs/simul_out_steps_mesh
   outdir=$indir/figs/gifs_mesh
   mkdir  -p $outdir
   cd $outdir
else
   indir1=$indir/figs/simul_out_steps_grid
   outdir=$indir/figs/gifs_grid
   mkdir  -p $outdir
   cd $outdir
   echo "Saving gif to `pwd`"
fi

cd $indir1
for dd in *
do
   out=$dd.gif
   lst=($indir1/$dd/*.png)
   Nlst=${#lst[@]}
   if [ $Nlst -gt 1 ]
   then
      echo convert -delay 15 -loop 0 $indir1/$dd/'*.png' $out
      convert -delay 15 -loop 0 $indir1/$dd/*.png $out
      mv $out $outdir
   fi
done
