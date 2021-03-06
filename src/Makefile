EIGEN_PATH = ${Eigen3_DIR}

CC = g++
CHECK = -Winline -Wall
PAR = -lm -std=c++17
CXXFLAGS = $(CHECK) $(PAR)
RELEASEFLAGS = -O3 -fopenmp -DNDEBUG
WASSMS_LDFLAGS = -L. -lwassms
DECONVMS_LDFLAGS = -L. -ldeconvms

ifeq ($(OS),Windows_NT)
    CXXFLAGS += -D WIN32
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CXXFLAGS += -D LINUX
    endif
    ifeq ($(UNAME_S),Darwin)
        CXXFLAGS += -D OSX
        RELEASEFLAGS += -Wa,-q # Assemble AVX/FMA instructions
    endif
endif

.PHONY: all vectorized release debug clean

all: clean wassms deconvms

wassms: buildLibwassms wassms.o
	$(CC) -o ../$@ wassms.o $(CXXFLAGS) $(WASSMS_LDFLAGS)

deconvms: buildLibdeconvms deconvms.o
	$(CC) -o ../$@ deconvms.o $(CXXFLAGS) $(DECONVMS_LDFLAGS)

release: CXXFLAGS += $(RELEASEFLAGS)
release: clean wassms deconvms

debug: CXXFLAGS += -g -O0
debug: clean wassms deconvms

# link libraries
buildLibwassms: spectrum.o similarity.o io.o
	ar r libwassms.a *.o

buildLibdeconvms: spectrum.o io.o ipm.o ldu.o heuristic.o solver/solver.o
	ar r libdeconvms.a *.o

# Compile source files
wassms.o: wassms.cpp
	$(CC) $(CXXFLAGS) -c wassms.cpp

deconvms.o: deconvms.cpp
	$(CC) $(CXXFLAGS) -c deconvms.cpp -I$(EIGEN_PATH)

spectrum.o: spectrum.cpp spectrum.hpp
	$(CC) $(CXXFLAGS) -c spectrum.cpp

io.o: io.cpp io.hpp
	$(CC) $(CXXFLAGS) -c io.cpp

similarity.o: similarity.cpp similarity.hpp
	$(CC) $(CXXFLAGS) -c similarity.cpp

ipm.o: ipm.cpp ipm.hpp
	$(CC) $(CXXFLAGS) -c ipm.cpp -I$(EIGEN_PATH)

ldu.o: ldu.cpp ldu.hpp
	$(CC) $(CXXFLAGS) -c ldu.cpp -I$(EIGEN_PATH)

heuristic.o: heuristic.cpp heuristic.hpp
	$(CC) $(CXXFLAGS) -c heuristic.cpp -I$(EIGEN_PATH)

solver/solver.o: solver/solver.cpp solver/solver.hpp
	$(CC) $(CXXFLAGS) -c solver/solver.cpp -I$(EIGEN_PATH)

# Clean command
clean:
	clear
	-rm wassms
	-rm deconvms
	-rm libwassms.a
	-rm libdeconvms.a
	-rm *.o
	-rm *.gch
