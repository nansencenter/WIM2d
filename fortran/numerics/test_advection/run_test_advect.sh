# place to put test results
mkdir -p out test_out
rm out/* test_out/*

bin/test_advect.exe

echo " "
echo "**********************************************************"
echo "now run test_advect_F.m in matlab to test results"
echo "**********************************************************"
echo " "
