CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O3 -march=native -Wall -Wextra
TARGET    = baseline
SRC       = main.cpp
HDR       = rng.hpp gaussian.hpp blackscholes.hpp montecarlo.hpp

$(TARGET): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC)

.PHONY: test clean
test: $(TARGET)
	./$(TARGET) selftest

clean:
	rm -f $(TARGET)
