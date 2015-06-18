# place to put test results
mkdir -p out test_out
rm out/* test_out/*

bin/test_advect.exec

echo " "
echo "**********************************************************"
echo "Testing options:"
echo "1. run test_advect_F.m in matlab"
echo "2. run test_advect_F.py in python"
echo "**********************************************************"
echo " "
