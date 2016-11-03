CXX	= g++ #Else try, g++-6

#See Makefile in HElib/src for more information on the various compiler options and in case of compiler errors.

#CXXFLAGS = -g -O2 -Wfatal-errors -Wshadow -Wall -I/usr/local/include -I../HElib/src -I..
#CXXFLAGS = -g -O2 -std=c++11 -pthread -DFHE_THREADS -DFHE_DCRT_THREADS -DFHE_BOOT_THREADS -I../HElib/src -I..
CXXFLAGS = -g -O0 -std=c++11 -pthread -I../HElib/src -I..
#CXXFLAGS = -g -O2 -std=c++11 -I../HElib/src -I..

HEADERS = $(wildcard *.h)
OTHER_OBJ = findRing.o flt2int_encode.o

all: helib.a.fxpt helib.a.binary helib.a.ulong

findRing.o: findRing.cpp findRing.h 
	 $(CXX) $(CXXFLAGS) -c -o $@ $<

flt2int_encode.o: flt2int_encode.cpp flt2int_encode.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<


API-Bridge.o.fxpt: API-Bridge.cpp $(OTHER_OBJ) $(HEADERS) 
	$(CXX) $(CXXFLAGS) -DFXPT -c -o $@ $< 

helib.a.fxpt: API-Bridge.o.fxpt
	ar rcs $@ $< $(OTHER_OBJ)

API-Bridge.o.binary: API-Bridge.cpp $(OTHER_OBJ) $(HEADERS)
	$(CXX) $(CXXFLAGS) -DBINARY -c -o $@ $<

helib.a.binary: API-Bridge.o.binary
	ar rcs $@ $< $(OTHER_OBJ)

API-Bridge.o.ulong: API-Bridge.cpp $(OTHER_OBJ) $(HEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

helib.a.ulong: API-Bridge.o.ulong
	ar rcs $@ $< $(OTHER_OBJ)

.PHONY: clean
clean:
	rm -f helib.a* API-Bridge.o*
	rm -f *.o
