# Plain-make build (no cmake needed). `make` produces build/dhj.
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra
LDFLAGS  ?= -pthread -lz

HDRS := $(wildcard include/dhj/core/*.hpp) $(wildcard include/dhj/dynamics/*.hpp)

.PHONY: all clean tools
all: build/dhj

build/dhj: src/main.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

tools: build/dhj_compare_checkpoints

build/dhj_compare_checkpoints: tools/compare_checkpoints.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf build
