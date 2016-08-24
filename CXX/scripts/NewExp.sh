WC=$WIM2D_PATH/CXX

# make directory to compile in, and link to makefile
mkdir -p Build_cpp
cd Build_cpp
ln -s $WC/Build/Makefile Makefile

# get eg config file and run script
cd ..
cp $WC/config_files/wim.cfg .
ln -s $WC/scripts/run_cpp.sh run_cpp.sh

echo "Compile in Build_cpp"
echo "Run here:"
echo "- edit options in wim.cfg and type <./run_cpp.sh>"
