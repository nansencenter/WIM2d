# place to put results
mkdir -p out
mkdir -p log
rm prog/*

../Build/WIM2d.exe

echo " ******************************************"
echo " Now to check results, do:"
echo " cd test"
echo " run test_WIM2d_F.m in matlab"
echo " ******************************************"
echo " "
