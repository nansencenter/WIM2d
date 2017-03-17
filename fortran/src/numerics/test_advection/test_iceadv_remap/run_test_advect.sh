# place to put test results
mkdir -p out test_out
rm out/* test_out/*

ex2=bin/test_advect.exec      # 2d advection executable
echo running $ex2
$ex2

echo " "
echo "**********************************************************"
echo "Testing options:"
echo "1. run ../test_advect_F.m in matlab from this directory"
echo "2. run ../test_advect_F.py in python from this directory"
echo "**********************************************************"
echo " "
