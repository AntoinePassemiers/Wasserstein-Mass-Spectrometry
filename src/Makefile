EIGEN_PATH = ${Eigen3_DIR}

CC = g++
CHECK = -Winline -Wall
PAR = -lm -std=c++17
CXXFLAGS = $(CHECK) $(PAR)
WASSMS_LDFLAGS = -L. -lwassms
DECONVMS_LDFLAGS = -L. -ldeconvms 

.PHONY: all vectorized release debug clean

all: clean wassms deconvms

wassms: buildLibwassms wassms.o
	$(CC) -o ../$@ wassms.o $(CXXFLAGS) $(WASSMS_LDFLAGS)

deconvms: buildLibdeconvms deconvms.o
	$(CC) -o ../$@ deconvms.o $(CXXFLAGS) $(DECONVMS_LDFLAGS)

release: CXXFLAGS += -O3
release: clean wassms deconvms

debug: CXXFLAGS += -g -O0
debug: clean wassms deconvms

# link libraries
buildLibwassms: spectrum.o wasserstein.o io.o
	ar r libwassms.a *.o

buildLibdeconvms: spectrum.o wasserstein.o io.o ipm.o
	ar r libdeconvms.a *.o

wassms.o: wassms.cpp
	$(CC) $(CXXFLAGS) -c wassms.cpp

deconvms.o: deconvms.cpp
	$(CC) $(CXXFLAGS) -c deconvms.cpp -I$(EIGEN_PATH)

spectrum.o: spectrum.cpp spectrum.hpp
	$(CC) $(CXXFLAGS) -c spectrum.cpp

io.o: io.cpp io.hpp
	$(CC) $(CXXFLAGS) -c io.cpp

wasserstein.o: wasserstein.cpp wasserstein.hpp
	$(CC) $(CXXFLAGS) -c wasserstein.cpp

ipm.o: ipm.cpp ipm.hpp
	$(CC) $(CXXFLAGS) -c ipm.cpp -I$(EIGEN_PATH)

clean:
	clear
	-rm wassms
	-rm deconvms
	-rm libwassms.a
	-rm libdeconvms.a
	-rm *.o
	-rm *.gch
