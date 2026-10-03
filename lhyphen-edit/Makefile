# lhedit se compile avec un simple compilateur C++17 et la bibliothèque standard,
# de sorte qu'il peut être construit directement sur un cluster :  make && ../lhedit
# L'exécutable est copié à la racine de l-hyphen, à côté de run et see2.
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

TARGET = lhedit
HDRS   = terminal.hpp text_buffer.hpp utf8.hpp lhyphen_lang.hpp

all: ../$(TARGET)

$(TARGET): main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ main.cpp

../$(TARGET): $(TARGET)
	cp $(TARGET) ../$(TARGET)

clean:
	rm -f $(TARGET) ../$(TARGET)

.PHONY: all clean
