#  Copyright or © or Copr. l-hyphen
#
#  This software is developed for an ACADEMIC USAGE
#
#  This software is governed by the CeCILL-B license under French law and
#  abiding by the rules of distribution of free software.  You can  use,
#  modify and/ or redistribute the software under the terms of the CeCILL-B
#  license as circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and  rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty  and the software's author,  the holder of the
#  economic rights,  and the successive licensors  have only  limited
#  liability.
#
#  In this respect, the user's attention is drawn to the risks associated
#  with loading,  using,  modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean  that it is complicated to manipulate,  and  that  also
#  therefore means  that it is reserved for developers  and  experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or
#  data to be ensured and,  more generally, to use and operate it in the
#  same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL-B license and that you accept its terms.


UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  GNU_GPP := $(shell ls /usr/local/bin/g++-* /opt/homebrew/bin/g++-* 2>/dev/null | grep -Eo 'g\+\+-([0-9]+)' | sort -V | tail -1)
  CXX := $(if $(GNU_GPP),$(GNU_GPP),g++)
  
  CXXFLAGS = -O3 -Wall -Wextra -pedantic -Wno-unknown-pragmas -std=c++17 -I ./toofus
  LDFLAGS = 
  GLLINK = `pkg-config --libs gl glu glut`
  GLFLAGS = `pkg-config --cflags gl glu glut`	
  # on apple, use brew to install freeglut and mesa-glu
	
	GLFWLINK = `pkg-config --libs glfw3`
	GLFWFLAGS = `pkg-config --cflags glfw3`
else
  CXX = g++
  CXXFLAGS = -g -fopenmp -O3 -Wall -std=c++17 -I ./toofus
  LDFLAGS = -lomp
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

# Ensure the toofus directory exists
ifeq ($(wildcard ./toofus),)
$(info Cloning ToOfUs repository)
$(shell git clone https://github.com/richefeu/toofus.git > /dev/null 2>&1)
endif

SOURCES = expressionParser.cpp Neighbor.cpp Control.cpp Node.cpp Bar.cpp Cell.cpp Lhyphen.cpp
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all clean clone_toofus

all: run see2

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
	$(CXX) $(LDFLAGS) -o $@ run.o liblhyphen.a
	
see: see.cpp liblhyphen.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see.o $(GLFLAGS)
	$(CXX) $(LDFLAGS) -o $@ see.o liblhyphen.a $(GLLINK)
	
see2: see2.cpp liblhyphen.a
	@echo "\033[0;32m-> BUILDING APPLICATION" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o see2.o $(GLFWFLAGS) -Wno-missing-field-initializers
	$(CXX) $(LDFLAGS) -o $@ see2.o liblhyphen.a $(GLFWLINK) -framework OpenGL
	

