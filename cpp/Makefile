# Plain-make build (no cmake needed). `make` produces build/dhj.
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra
LDFLAGS  ?= -pthread -lz

HDRS := $(wildcard include/dhj/core/*.hpp) $(wildcard include/dhj/dynamics/*.hpp)

.PHONY: all capd clean
all: build/dhj

build/dhj: src/main.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Experimental: `make capd CAPD_CONFIG=<CAPD>/build/bin/capd-config` builds build/dhj_capd, which
# accepts --validated (CAPD-validated integration error, see include/dhj/core/validated.hpp).
# The bridge gets CAPD's own flags (native intervals, -frounding-math). src/main.cpp keeps the
# flags above, so the numerical integrator is compiled as in build/dhj.
CAPD_CONFIG ?= capd-config

capd: build/dhj_capd

build/capd_bridge.o: src/capd_bridge.cpp include/dhj/core/validated.hpp
	@mkdir -p build
	$(CXX) $(shell $(CAPD_CONFIG) --cflags) -DNDEBUG -DDHJ_WITH_CAPD -Iinclude -c $< -o $@

build/dhj_capd: src/main.cpp $(HDRS) build/capd_bridge.o
	$(CXX) $(CXXFLAGS) -DDHJ_WITH_CAPD $< build/capd_bridge.o -o $@ $(LDFLAGS) $(shell $(CAPD_CONFIG) --libs)

clean:
	rm -rf build
