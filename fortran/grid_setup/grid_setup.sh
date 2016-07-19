bin="$WIM2D_PATH/fortran/bin"
GS="$WIM2D_PATH/fortran/grid_setup"
# mkdir -p ../run/inputs

echo " " 
echo ${bin}/grid_setup.exec
echo " " 
${bin}/grid_setup.exec

echo " " 
echo "Now compile in :  ."
echo "Run with '.' in path"
echo " " 

kernel=(uname -s)

if [ 1 -eq 1 ]
then
   python $GS/grid_plot.py

   if [ $kernel == "Darwin" ]
   then
      open out/land_mask.png
   else
      display out/land_mask.png &
   fi
fi
