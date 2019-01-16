NOMEFILE = wassms.cpp
CHECK = -Winline -Wall
PAR = -lm -std=c++17
CFLAGS = -L. -lwassms $(CHECK) $(PAR)

.PHONY: all vectorized release debug clean

all: clean wassms deconvms

wassms: buildLib wassms.cpp
	g++ wassms.cpp -o ../$@ $(CFLAGS)

deconvms: buildLib deconvms.cpp
	g++ deconvms.cpp -o ../$@ $(CFLAGS)

release: CFLAGS += -O3
release: clean wassms

debug: CFLAGS += -g -DDEBUG  -O0
debug: clean wassms


# link libraries
buildLib: spectrum.o wasserstein.o
	ar r libwassms.a *.o

spectrum.o: spectrum.cpp
	g++ $(CFLAGS) $(PAR) -c spectrum.cpp

wasserstein.o: wasserstein.o
	g++ $(CFLAGS) $(PAR) -c wasserstein.cpp

clean:
	clear
	-rm wassms
	-rm deconvms
	-rm libwassms.a
	-rm *.o
	-rm *.gch