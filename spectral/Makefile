CXX      ?= g++
CXXFLAGS  = -std=c++17 -O3 -Wall -Wextra

# The Vkl precomputation is parallelised with OpenMP. Apple's clang (the `g++`
# on macOS) does not support -fopenmp; build with real GCC to get the threading,
# e.g.  make CXX=g++-14
# Without OpenMP the code still builds and runs correctly, just single-threaded.
OPENMP   := $(shell $(CXX) -fopenmp -x c++ -E /dev/null >/dev/null 2>&1 && echo -fopenmp)
CXXFLAGS += $(OPENMP)

SRC = BetheSF/SingleFile.cpp BetheSF/SingleFileDerived.cpp BetheSF/SingleFileCorrelations.cpp
DEPS = BetheSF/SingleFile.hpp BetheSF/SingleFileDerived.hpp \
       BetheSF/combinations.hpp BetheSF/integerpartitions.hpp

example: example.cpp $(SRC) $(DEPS)
	$(CXX) -o $@ example.cpp $(SRC) $(CXXFLAGS)
ifeq ($(OPENMP),)
	@echo "note: built WITHOUT OpenMP -- n_threads will have no effect."
endif

clean:
	rm -f example correlation.csv

.PHONY: clean
