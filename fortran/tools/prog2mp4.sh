if [ $# -eq 0 ]
then
   echo "Usage : prog2mp4.sh [ variable name eg Hs ]"
   echo "or    : prog2mp4.sh [ variable name eg Hs ] [ path to figs/prog ]"
   echo "NB needs to be run from figs/prog directory"
   echo "(where the png files from plot_prog.py are kept)"
   exit
elif [ $# -eq 1 ]
then
   vbl=$1
elif [ $# -eq 2 ]
then
   vbl=$1
   cd $2
else
   vbl='Hs'
fi

P=`pwd` # run from figs/prog/ folder or pass in path there as 2nd argument
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
fps=7 # frames per second
mov=${vbl}_prog.mp4
cd tmp
echo ffmpeg -framerate $fps -i ${vbl}%03d.png -c:v libx264 -r $fps -pix_fmt yuv420p $mov

#ffmpeg -framerate $fps -i ${vbl}%03d.png -c:v libx264 -r $fps -pix_fmt yuv420p $mov
#ffmpeg -r 1/5 -i ${vbl}%03d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p $mov
#ffmpeg -framerate 25 -i ${vbl}%03d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $mov
#ffmpeg -framerate $fps -i ${vbl}%03d.png -b:v 64k -c:v libx264 -r 24 -pix_fmt yuv420p $mov
#ffmpeg -r 5 -i ${vbl}%03d.png -vf "scale=1920:1080,format=yuv420p" -codec:v libx264 $mov
ffmpeg -framerate $fps -i ${vbl}%03d.png -c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p $mov

#ffmpeg -framerate $fps -i ${vbl}%03d.png -c:v libx264 -r $fps -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" $mov


#clean up
outdir2=$P/../prog_movies
mkdir -p $outdir2
mv $mov $outdir2
cd ..
rm -r tmp

# play the movie (mac)
# echo open -a /Applications/QuickTime\\ Player.app "$outdir2/$mov"
# open -a /Applications/QuickTime\ Player.app "$outdir2/$mov"
echo open "$outdir2/$mov"
open "$outdir2/$mov"
