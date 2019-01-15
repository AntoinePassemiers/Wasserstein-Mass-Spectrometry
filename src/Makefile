NOMEFILE = wassms.cpp

CHECK = -Winline -Wall

PAR = -lm -std=c++17

CFLAGS = -L. -lwassms $(CHECK) $(PAR)

.PHONY: all vectorized release debug cleanComp clean

all: clean wassms

wassms: buildLib $(NOMEFILE)
	g++ $(NOMEFILE) -o ../$@ $(CFLAGS)

release: CFLAGS += -O3
release: clean wassms

debug: CFLAGS += -g -DDEBUG  -O0
debug:  clean wassms
# clean
# -g      : debug symbols
# -DDEBUG : declare DEBUG macro


# link libraries
buildLib: spectrum.o
	ar r libwassms.a *.o

spectrum.o: spectrum.cpp
	g++ $(CFLAGS) $(PAR) -c spectrum.cpp

cleanComp:
	clear
	-rm wassms
	-rm libwassms.a

clean:
	clear
	-rm wassms
	-rm libwassms.a
	-rm *.o
	-rm *.gch