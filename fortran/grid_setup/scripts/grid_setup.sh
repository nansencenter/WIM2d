# run from Build directory
GS="$WIM2D_PATH/fortran/grid_setup"
Gbin="$GS/bin"
P=`pwd`
bin=`readlink -f $P/../bin`

Gexec=${Gbin}/grid_setup.exec
if [ ! -f $Gexec ]
then
   echo "Compiling grid setup executable"
   cd $GS/Build
   make
   cd $P
fi
echo " " 
echo $Gexec
echo " " 
$Gexec

echo " " 
echo "Now compile in : $P"
echo "Run with $bin in matlab/python path (or just call executables in there)"
echo " " 

kernel=`uname -s`

PLOT=1
if [ $PLOT -eq 1 ]
then
   python $GS/scripts/grid_plot.py
fi
   
gdir=../grid
mkdir -p $gdir
mv wim_grid.a wim_grid.b $gdir
cp *.h $gdir

if [ $PLOT -eq 1 ]
then
   mv out/land_mask.png $gdir
   rm -r out

   pic=$gdir/land_mask.png
   if [ $kernel == "Darwin" ]
   then
      # OSX
      open $pic
   else
      display $pic &
   fi
fi
