if [ $# -eq 0 ]
then
   echo "Usage : rename_files.sh directory"
   echo "or : rename_files.sh directory pattern"
   exit
elif [ $# -eq 1 ]
then
   ddir=$1
fi

cd $ddir
P=`pwd` # full path to dir
mkdir -p tmp

for pattern in mesh field
do
   for ext in bin dat
   do
      n=-1
      for f in $pattern*.$ext
      do
         n=$((n+1))
         f0=$P/`basename $f`
         echo $f0
         g=tmp/${pattern}_$n.$ext
         echo $g
         echo ln -s $f0 $g
         ln -s $f0 $g
      done
   done
done
