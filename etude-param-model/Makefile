# Outils de l'étude paramétrique : simple compilateur C++17, aucune dépendance
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

all: generate_cases

generate_cases: generate_cases.cpp
	$(CXX) $(CXXFLAGS) -o $@ generate_cases.cpp

clean:
	rm -f generate_cases

.PHONY: all clean
