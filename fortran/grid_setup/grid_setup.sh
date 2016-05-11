bin="../bin"
mkdir -p ../run/inputs

echo " " 
echo ${bin}/grid_setup.exe
echo " " 
${bin}/grid_setup.exec

echo " " 
echo "Now compile in :  ../Build"
echo "Run in         :  ../run"
echo " " 

kernel=(uname -s)

if [ 1 -eq 1 ]
then
   python grid_plot.py

   if [ $kernel == "Darwin" ]
   then
      open test/out/land_mask.png
   else
      display test/out/land_mask.png &
   fi
fi
