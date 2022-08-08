# The compiler to be used
CXX = g++
LINK = $(CXX)

# Paths
TOOFUSPATH = ~/toofus

# The list of flags passed to the compiler
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -I $(TOOFUSPATH)


GLUTFLAGS = `pkg-config --cflags glut`
GLUTLINK = `pkg-config --libs glut` -framework OpenGL


SOURCES = Neighbor.cpp Control.cpp Node.cpp Bar.cpp Cell.cpp Lhyphen.cpp
HEADERS = $(SOURCES:%.cpp=%.hpp)
OBJECTS = $(SOURCES:%.cpp=%.o)


%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean

all: run

clean:
	rm -f *~ *.o run

run: run.cpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)

see: see.cpp see.hpp 
	$(CXX) $(CXXFLAGS) $(GLUTFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(GLUTLINK)