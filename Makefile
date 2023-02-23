# The compiler to be used
CXX = g++-12
LINK = $(CXX)

# Paths
TOOFUSPATH = ~/toofus

# The list of flags passed to the compiler
CXXFLAGS = -Wall -Wextra -Wshadow -pedantic -Wconversion -O3 -std=c++20 -I $(TOOFUSPATH)


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
	rm -f *~ *.o run see

run: run.cpp $(HEADERS) $(SOURCES) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS)

see: see.cpp see.hpp run
	$(CXX) $(CXXFLAGS) $(GLUTFLAGS) -c $<
	$(LINK) $(CXXFLAGS) -o $@ $@.o $(OBJECTS) $(GLUTLINK)
