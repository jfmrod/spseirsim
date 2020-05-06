LDLIBS=`eutils-config --libs` -lavcodec -lavformat -lswscale -lavutil  -lshp -lcairo -lgeotiff -ltiff
LDFLAGS= -Wl,-rpath -Wl,/home/jfmrod/gaia/usr/lib
CXXFLAGS=`eutils-config --cxxflags` -g -O2 -I/usr/include/cairo
CC=g++

all : seirsim spseirsim

mobile :
	ffmpeg -i test.mp4 -vcodec libx264 -crf 22 -profile:v baseline -r:v 25 -threads 4 -acodec libvo_aacenc -ar 22050 -ab 63k test.mobile.mp4

run :
	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.001 -ftravel 1.5 -r0 3.5 -foldr 0.01 -tostart 150 -toend 300 -events 50:0.6,60:0.25,125:0.30,150:0.80 -tmax 500 -nthreads 40
#	./spseirsim -fintra 0.5 -finter 0.6 -fglobal 0.01 -r0 2.7 -tmstart 130 -fmr 0.3 -tmend 200 -tmax 400

runpt:
	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.01 -ftravel 1.5 -r0 3.5 -foldr 0.005 -tostart 120 -toend 250 -events 50:0.6,60:0.25,125:0.30,121:0.7,249:0.8,450:0.8 -fpop data/portugal.agegroups -fshape data/popdensmaps/gadm36_PRT/gadm36_PRT_0.shp -lonLimit -11.0 -tmax 500 -nthreads 40 -agethres 55
# scenario for portugal requires pop of >55 instead of >50 to be isolated, during this time there should be some social distancing to ensure ICUs are not overwhelmed
# after the older age group isolation removal, one still needs some social distancing 0.8 to prevent a new wave

#	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.001 -ftravel 1.5 -r0 3.5 -foldr 0.01 -tostart 150 -toend 300 -events 50:0.6,60:0.25,125:0.30,150:0.80 -fpop data/portugal.agegroups -fshape data/popdensmaps/gadm36_PRT/gadm36_PRT_0.shp -lonLimit -11.0 -tmax 500 -nthreads 40
#	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.001 -ftravel 1.5 -r0 3.5 -foldr 0.1 -tostart 150 -toend 300 -events 50:0.6,60:0.25,125:0.30,150:0.80 -fpop data/portugal.agegroups -fshape data/popdensmaps/gadm36_PRT/gadm36_PRT_0.shp -lonLimit -11.0 -tmax 500 -nthreads 40

runpl:
	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.001 -ftravel 1.5 -r0 3.5 -foldr 0.1 -tostart 150 -toend 300 -events 50:0.6,60:0.25,125:0.30,150:0.80 -fpop data/poland.agegroups -fshape data/popdensmaps/gadm36_POL/gadm36_POL_0.shp -tmax 500 -nthreads 40 -fsinglehh 0.05 -flargehh 0.1

runuk:
	./spseirsim -fintra 0.6 -finter 0.7 -fglobal 0.001 -ftravel 1.5 -r0 3.5 -foldr 0.1 -tostart 150 -toend 300 -events 50:0.6,60:0.25,125:0.30,150:0.80 -fpop data/uk.agegroups -fshape data/popdensmaps/gadm36_GBR/gadm36_GBR_0.shp -tmax 500 -nthreads 40 -flargehh 0.08

seirsim : seirsim.cpp
spseirsim : spseirsim.o videoenc.o

spseirsim.o : spseirsim.cpp videoenc.h
videoenc.o : videoenc.cpp videoenc.h

