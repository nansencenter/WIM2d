sudo ./b2 -j8 \
	 --layout=system \
	 --debug-configuration \
	 --prefix=$BOOST_PREFIX \
	 toolset=darwin \
	 variant=release \
	 threading=single,multi \
	 link=shared,static \
	 cxxflags="-std=c++11 -fPIC -O3"

sudo ./b2 install --prefix=$BOOST_PREFIX

# linkflags="-stdlib=libc++"
