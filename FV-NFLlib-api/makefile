# Try various compiler options in case of compiler errors.

#CXX = g++-6
CXX = g++
#CXXFLAGS  = -march=native -mtune=native -DNTT_SSE -std=c++11 -funroll-loops -Wall -g -O3 -I../FV-NFLlib	
# If need arises include paths to "cereal/archives/binary.hpp" and "boost/variant.hpp"
CXXFLAGS = -march=native -mtune=native -DNTT_SSE -std=c++11 -funroll-loops -Wall -g -O3 -I../FV-NFLlib -I/home/vikka/Downloads/cereal-master/include/ -I/home/vikka/Downloads/boost_1_61_0/

HEADERS = $(wildcard *.h) $(wildcard *.hpp)

all: fshe.a.binary fshe.a.ulong fshe.a.mpz

API-Bridge.o.binary: API-Bridge.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -DBINARY -c -o API-Bridge.o.binary API-Bridge.cpp

fshe.a.binary: API-Bridge.o.binary
	ar -q fshe.a.binary API-Bridge.o.binary

API-Bridge.o.ulong: API-Bridge.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c -o API-Bridge.o.ulong API-Bridge.cpp

fshe.a.ulong: API-Bridge.o.ulong
	ar -q fshe.a.ulong API-Bridge.o.ulong

API-Bridge.o.mpz: API-Bridge.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -DMPZ -c -o API-Bridge.o.mpz API-Bridge.cpp

fshe.a.mpz: API-Bridge.o.mpz
	ar -q fshe.a.mpz API-Bridge.o.mpz

.PHONY: clean
clean:
	rm -f fshe.a* API-Bridge.o*
