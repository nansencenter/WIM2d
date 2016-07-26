GS="$WIM2D_PATH/fortran/grid_setup"
bin="$GS/bin"
P=`pwd`

echo " " 
echo ${bin}/grid_setup.exec
echo " " 
${bin}/grid_setup.exec

echo " " 
echo "Now compile in :  $P"
echo "Run with $bin in path"
echo " " 

kernel=(uname -s)

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
      open $pic
   else
      display $pic &
   fi
fi
