# The compiler to be used
CXX = g++-14
LINK = $(CXX)

# Paths
TOOFUSPATH = ~/toofus

# The list of flags passed to the compiler
CXXFLAGS = -Wall -Wextra -pedantic -Wconversion -Wno-unknown-pragmas -O3 -std=c++11 -I $(TOOFUSPATH) -DENABLE_PROFILING
# add -fopenmp to enable OpenMP (but for now it is worst, excepted if the number of cells is very high)

GLUTFLAGS = `pkg-config --cflags glut`
GLUTLINK = `pkg-config --libs glut` -framework OpenGL

SOURCES = Neighbor.cpp Control.cpp Node.cpp Bar.cpp Cell.cpp Lhyphen.cpp
HEADERS = $(SOURCES:%.cpp=%.hpp)
OBJECTS = $(SOURCES:%.cpp=%.o)

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean

all: run see

clean:
	rm -f *~ *.o run see conf2z

run: run.cpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)

see: see.cpp see.hpp run
	$(CXX) $(CXXFLAGS) $(GLUTFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(GLUTLINK)

conf2z: conf2z.cpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)
