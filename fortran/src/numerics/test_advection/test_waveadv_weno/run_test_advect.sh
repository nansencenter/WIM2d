# place to put test results
mkdir -p out test_out
rm out/* test_out/*

dim=2
if [ $# -eq 1 ]
then
   dim=$1
fi

if [ $dim -eq 2 ]
then
   ex2=bin/test_advect.exec      # 2d advection executable
   echo running $ex2
   $ex2
else
   ex1=bin/test_advect_1d.exec   # 1d advection executable
   echo running $ex1
   $ex1
fi

echo " "
echo "**********************************************************"
echo "Testing options:"
echo "1. run ../test_advect_F.m in matlab from this directory"
echo "2. run ../test_advect_F.py in python from this directory"
echo "**********************************************************"
echo " "
