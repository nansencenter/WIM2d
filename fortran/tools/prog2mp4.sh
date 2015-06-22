P=`pwd` # run from figs/prog/ folder

if [ $# -eq 1 ]
then
   vbl=$1
else
   vbl='Hs'
fi

mkdir tmp

n=-1
for stepdir in $P/*
do
   if [ ! -f $stepdir ]
   then
      n=$((n+1))
      cn=`printf "%3.3d" $n`
      f0=$stepdir/$vbl.png
      echo $f0
      f=tmp/${vbl}$cn.png
      echo $f
      echo ln -s $f0 $f
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
outdir2=$P/../prog_movies
mkdir -p $outdir2
mv $mov $outdir2
cd ..
rm -r tmp

# play the movie (mac)
echo open -a /Applications/QuickTime\\ Player.app "$outdir2/$mov"
open -a /Applications/QuickTime\ Player.app "$outdir2/$mov"
