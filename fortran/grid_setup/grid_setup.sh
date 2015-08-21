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

if [ 1 -eq 1 ]
then
   python grid_plot.py
   open test/out/land_mask.png
fi
