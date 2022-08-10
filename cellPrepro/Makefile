# The compiler to be used
CXX = g++-11

# The list of flags passed to the compiler
CXXFLAGS = -Wall -Wextra -ansi -O3 -std=c++11 -I ~/toofus

# Notice that you need to install toofus before compiling
# JYD : -I /Users/jyd/Documents/Codes/richefeu-toofus-6c99da3
# VR : -I /usr/local/include/toofus

APP = cellPrepro

.PHONY: all clean

all: $(APP)

clean:
	rm -f *.o $(APP)

$(APP): $(APP).cpp $(APP).hpp
	$(CXX) $(CXXFLAGS) $(APP).cpp -o $(APP)

