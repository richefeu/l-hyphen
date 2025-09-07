
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  GNU_GPP := $(shell ls /usr/local/bin/g++-* /opt/homebrew/bin/g++-* 2>/dev/null | grep -Eo 'g\+\+-([0-9]+)' | sort -V | tail -1)
  CXX := $(if $(GNU_GPP),$(GNU_GPP),g++)
  
  CXXFLAGS = -O3 -Wall -Wextra -pedantic -Wno-unknown-pragmas -std=c++17 -I ./toofus
  GLLINK = `pkg-config --libs gl glu glut`
	GLFLAGS = `pkg-config --cflags gl glu glut`	
	# on apple, use brew to install freeglut and mesa-glu
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++17 -I ./toofus
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

# Ensure the toofus directory exists
ifeq ($(wildcard ./toofus),)
$(info Cloning ToOfUs repository)
$(shell git clone https://github.com/richefeu/toofus.git > /dev/null 2>&1)
endif

SOURCES = Neighbor.cpp Control.cpp Node.cpp Bar.cpp Cell.cpp Lhyphen.cpp
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all clean clone_toofus

all: run see

clean:
	@echo "\033[0;32m-> Remove object files\033[0m"
	rm -f *.o
	@echo "\033[0;32m-> Remove compiled applications\033[0m"
	rm -f run see 
	@echo "\033[0;32m-> Remove liblhyphen.a\033[0m"
	rm -f liblhyphen.a

clean+: clean
	@echo "\033[0;32m-> Remove local folder toofus\033[0m"
	rm -rf ./toofus
	
clone_toofus:
	@if [ ! -d "./toofus" ]; then \
		@echo "\033[0;32m-> CLONING ToOfUs (Tools often used)\033[0m"; \
		git clone https://github.com/richefeu/toofus.git; \
	fi

%.o: %.cpp 
	@echo "\033[0;32m-> COMPILING OBJECT" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@

liblhyphen.a: $(OBJECTS)
	@echo "\033[0;32m-> BUILDING LIBRARY" $@ "\033[0m"
	ar rcs $@ $^

run: run.cpp liblhyphen.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o run.o
	$(CXX) -o $@ run.o liblhyphen.a
	
see: see.cpp liblhyphen.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see.o $(GLFLAGS)
	$(CXX) -o $@ see.o liblhyphen.a $(GLLINK)
	
conf2z: conf2z.cpp liblhyphen.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o conf2z.o
	$(CXX) -o $@ conf2z.o liblhyphen.a

