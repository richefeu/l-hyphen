# The compiler to be used
CXX = g++

# The list of flags passed to the compiler
CXXFLAGS = -Wall -Wextra -ansi -O3 -Wwrite-strings -Wstrict-prototypes -Wuninitialized \
-Wunreachable-code -Wno-missing-braces  -Wno-missing-field-initializers \
-std=c++0x -I /usr/local/include/toofus
APP = l-hyphen

.PHONY: all clean

all: $(APP)

clean:
	rm -f *.o $(APP)

$(APP): $(APP).cpp $(APP).hpp
	$(CXX) $(CXXFLAGS) $(APP).cpp -o $(APP)

