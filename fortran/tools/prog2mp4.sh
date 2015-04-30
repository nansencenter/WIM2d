if [ $# -eq 0 ]
then
   vbl='Hs'
else
   vbl=$1
fi

P=`pwd`
outdir=$P/out_io/figs/prog
mkdir tmp

n=-1
for stepdir in $outdir/*
do
   if [ ! -f $stepdir ]
   then
      n=$((n+1))
      cn=`printf "%3.3d" $n`
      f0=$stepdir/$vbl.png
      f=tmp/${vbl}$cn.png
      ln -s $f0 $f
   fi
done

# make movie
fps=5 # frames per second
mov=${vbl}_prog.mp4
cd tmp
echo ffmpeg -framerate $fps -i ${vbl}%03d.png -c:v libx264 -r $fps -pix_fmt yuv420p $mov
ffmpeg -framerate $fps -i ${vbl}%03d.png -c:v libx264 -r $fps -pix_fmt yuv420p $mov

#clean up
outdir2=$outdir/../prog_movies
mkdir -p $outdir2
mv $mov $outdir2
cd ..
rm -r tmp

# play the movie (mac)
echo open -a /Applications/QuickTime\ Player.app "$outdir2/$mov"
open -a /Applications/QuickTime\ Player.app "$outdir2/$mov"
