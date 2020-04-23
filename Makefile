LDLIBS=`eutils-config --libs`
CXXFLAGS=`eutils-config --cxxflags` -g -O2

all : covidsim

covidsim : covidsim.cpp
