# nodegen / nodeview : simple compilateur C++17, aucune dépendance
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

all: nodegen nodeview

nodegen: nodegen.cpp common.hpp
	$(CXX) $(CXXFLAGS) -o $@ nodegen.cpp

nodeview: nodeview.cpp common.hpp
	$(CXX) $(CXXFLAGS) -o $@ nodeview.cpp

clean:
	rm -f nodegen nodeview

.PHONY: all clean
