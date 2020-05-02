LDLIBS=`eutils-config --libs` -lavcodec -lavformat -lswscale -lavutil  -lshp -lcairo -lgeotiff -ltiff
LDFLAGS= -Wl,-rpath -Wl,/home/jfmrod/gaia/usr/lib
CXXFLAGS=`eutils-config --cxxflags` -g -O2 -I/usr/include/cairo
CC=g++

all : seirsim spseirsim

run :
	./spseirsim -fintra 0.5 -finter 0.6 -fglobal 0.01 -r0 2.7 -tmstart 130 -fmr 0.3 -tmend 200 -tmax 400

runpt:
	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.001 -ftravel 1.5 -r0 3.5 -foldr 0.1 -tostart 150 -toend 300 -events 50:0.6,60:0.25,125:0.30,150:0.80 -fpop data/portugal.agegroups -fshape data/popdensmaps/gadm36_PRT/gadm36_PRT_0.shp -lonLimit -11.0 -tmax 500 -nthreads 40

runpl:
	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.001 -ftravel 1.5 -r0 3.5 -foldr 0.1 -tostart 150 -toend 300 -events 50:0.6,60:0.25,125:0.30,150:0.80 -fpop data/poland.agegroups -fshape data/popdensmaps/gadm36_POL/gadm36_POL_0.shp -tmax 500 -nthreads 40 -fsinglehh 0.05 -flargehh 0.1

seirsim : seirsim.cpp
spseirsim : spseirsim.o videoenc.o

spseirsim.o : spseirsim.cpp videoenc.h
videoenc.o : videoenc.cpp videoenc.h

