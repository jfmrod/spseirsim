LDLIBS=`eutils-config --libs` -lavcodec -lavformat -lswscale -lavutil  -lshp -lcairo -lgeotiff -ltiff
LDFLAGS= -Wl,-rpath -Wl,/home/jfmrod/gaia/usr/lib
CXXFLAGS=`eutils-config --cxxflags` -g -O2 -I/usr/include/cairo
CC=g++

all : seirsim spseirsim

run :
	./spseirsim -fintra 0.5 -finter 0.6 -fglobal 0.01 -r0 2.7 -tmstart 130 -fmr 0.3 -tmend 200 -tmax 400

seirsim : seirsim.cpp
spseirsim : spseirsim.o videoenc.o

spseirsim.o : spseirsim.cpp videoenc.h
videoenc.o : videoenc.cpp videoenc.h

