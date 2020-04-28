LDLIBS=`eutils-config --libs` -lavcodec -lavformat -lswscale -lavutil  -lshp
LDFLAGS= -Wl,-rpath -Wl,/home/jfmrod/gaia/usr/lib
CXXFLAGS=`eutils-config --cxxflags` -g -O2
CC=g++

all : seirsim spseirsim

seirsim : seirsim.cpp
spseirsim : spseirsim.o videoenc.o

spseirsim.o : spseirsim.cpp videoenc.h
videoenc.o : videoenc.cpp videoenc.h

