NOMEFILE = wassms.cpp
CHECK = -Winline -Wall
PAR = -lm -std=c++17
CFLAGS = -L. -lwassms $(CHECK) $(PAR)
CFLAGS += -I eigen

.PHONY: all vectorized release debug clean

all: clean wassms deconvms

wassms: buildLibwassms wassms.cpp
	g++ wassms.cpp -o ../$@ $(CFLAGS)

deconvms: buildLibdeconvms deconvms.cpp
	g++ deconvms.cpp -o ../$@ $(CFLAGS)

release: CFLAGS += -O3
release: clean wassms

debug: CFLAGS += -g -DDEBUG  -O0
debug: clean wassms

# link libraries
buildLibwassms: spectrum.o wasserstein.o io.o
	ar r libwassms.a *.o

buildLibdeconvms: spectrum.o wasserstein.o io.o
	ar r libdeconvms.a *.o

spectrum.o: spectrum.cpp
	g++ $(CFLAGS) $(PAR) -c spectrum.cpp

io.o: io.cpp
	g++ $(CFLAGS) $(PAR) -c io.cpp

wasserstein.o: wasserstein.o
	g++ $(CFLAGS) $(PAR) -c wasserstein.cpp

clean:
	clear
	-rm wassms
	-rm deconvms
	-rm libwassms.a
	-rm libdeconvms.a
	-rm *.o
	-rm *.gch