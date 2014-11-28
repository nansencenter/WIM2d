# place to put results
outdir="out"
mkdir -p $outdir
mkdir -p $outdir/log
mkdir -p $outdir/prog
mkdir -p $outdir/binaries
mkdir -p $outdir/binaries/prog
rm prog/*

../bin/WIM2d.exe

echo " ******************************************"
echo " Now to check results, do:"
echo " cd test"
echo " run test_WIM2d_F.m in matlab"
echo " ******************************************"
echo " "
