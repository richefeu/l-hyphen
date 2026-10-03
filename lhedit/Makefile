# lhedit se compile avec un simple compilateur C++17 et la bibliothèque standard,
# de sorte qu'il peut être construit directement sur un cluster :  make && ./lhedit
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

TARGET = lhedit
HDRS   = terminal.hpp text_buffer.hpp utf8.hpp lhyphen_lang.hpp

all: $(TARGET)

$(TARGET): main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ main.cpp

clean:
	rm -f $(TARGET)

.PHONY: all clean
