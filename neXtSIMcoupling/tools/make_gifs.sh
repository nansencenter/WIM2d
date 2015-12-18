if [ $# -lt 2 ]
then
   echo "Usage make_gifs.sh [dir with png files] [mesh_opt]"
   echo "mesh_opt = 0 (save output to gifs_mesh) or 1 (save output to gifs_grid)"
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
   outdir=gifs_mesh
   mkdir  -p $outdir
   cd $outdir
else
   outdir=gifs_grid
   mkdir  -p $outdir
   cd $outdir
   echo "Saving gif to `pwd`"
fi

out=`basename $indir`.gif
# echo convert -delay 15 -loop 0 $indir/*.png $out
convert -delay 15 -loop 0 $indir/*.png $out
