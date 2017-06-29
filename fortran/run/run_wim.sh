infile=infile_nonIO.txt
if [ ! -f "$infile" ]
then
   Infile=$WIM2D_PATH/fortran/scripts/infiles/infile_nonIO.txt
   echo Getting $Infile...
   echo Edit it and rerun.
   cp $Infile .
   exit
fi

if [ "WIM2D_PATH" == "" ]
then
   echo Please set \$WIM2D_PATH and rerun
   exit
fi

# run executable
EXEC=../bin/wim2d.exec
if [ ! -f "$EXEC" ]
then
   echo $EXEC not present - compiling
   P=`pwd`
   echo cd ../Build
   cd ../Build
   echo make exec
   make exec || { echo ' ' ; echo 'COMPILATION FAILED' ; echo ' ' ; exit 1 ; }
   cd $P
fi
$EXEC

# plotting program
scr=fortran/tools/plot_prog.sh
echo To make plots run
echo \$WIM2D_PATH/$scr
$WIM2D_PATH/$scr
