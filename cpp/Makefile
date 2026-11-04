# Plain-make build (no cmake needed). `make` produces build/dhj.
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra
LDFLAGS  ?= -pthread -lz

HDRS := $(wildcard include/dhj/core/*.hpp) $(wildcard include/dhj/dynamics/*.hpp)

.PHONY: all clean
all: build/dhj

build/dhj: src/main.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -rf build
