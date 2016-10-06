CXX	= g++	# Or use g++-6
CXXFLAGS= -Wall -O0 -g -std=c++11
LDFLAGS = -lgmpxx -lgmp 
PARALLEL= yes

ifeq ($(PARALLEL),yes)
  CXXFLAGS+= -fopenmp
  LDFLAGS += -fopenmp
endif

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h) $(wildcard *.hpp)
OBJECTS = $(SOURCES:.cpp=.o) 

LIB_FSHE = FV-NFLlib-api/fshe.a
EXEC_FSHE	= main_fshe
CXXFLAGS_FSHE = -pedantic -march=native -mtune=native -DNTT_SSE -std=c++11 -funroll-loops 
LDFLAGS_FSHE = -lgmpxx -lgmp -lnfllib -lmpfr -lm

LIB_HELIB1 = HElib/src/fhe.a
LIB_HELIB2 = HElib-api/helib.a
EXEC_HELIB	= main_helib
#CXXFLAGS_HELIB = -pthread -DFHE_THREADS  -DFHE_DCRT_THREADS -DFHE_BOOT_THREADS -IHElib/src -IHElib-api
CXXFLAGS_HELIB = -pthread -DFHE_THREADS -IHElib/src -IHElib-api
#CXXFLAGS_HELIB = -IHElib/src -IHElib-api
LDFLAGS_HELIB = -L/usr/local/lib -lntl -lgmpxx -lgmp -lm


all: $(EXEC_FSHE)_mpz $(EXEC_FSHE)_binary $(EXEC_FSHE)_ulong $(EXEC_HELIB)_fxpt $(EXEC_HELIB)_binary $(EXEC_HELIB)_ulong
	rm -rf *.dSYM/

main.o.mpz: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -DMPZ -c -o main.o.mpz main.cpp

main.o.binary: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -DBINARY -c -o main.o.binary main.cpp

main.o.ulong: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -c -o main.o.ulong main.cpp
	
main.o.fxpt: $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -DFXPT -c -o main.o.fxpt main.cpp
	
$(EXEC_FSHE)_mpz: main.o.mpz $(LIB_FSHE).mpz Makefile
	$(CXX) -DMPZ -o $(EXEC_FSHE)_mpz main.o.mpz $(LIB_FSHE).mpz $(LDFLAGS) $(LDFLAGS_FSHE) $(CXXFLAGS) $(CXXFLAGS_FSHE)
$(EXEC_FSHE)_binary: main.o.binary $(LIB_FSHE).binary Makefile
	$(CXX) -DBINARY -o $(EXEC_FSHE)_binary main.o.binary $(LIB_FSHE).binary $(LDFLAGS) $(LDFLAGS_FSHE) $(CXXFLAGS) $(CXXFLAGS_FSHE)
$(EXEC_FSHE)_ulong: main.o.ulong $(LIB_FSHE).ulong Makefile
	$(CXX) -o $(EXEC_FSHE)_ulong main.o.ulong $(LIB_FSHE).ulong $(LDFLAGS) $(LDFLAGS_FSHE) $(CXXFLAGS) $(CXXFLAGS_FSHE)

$(EXEC_HELIB)_fxpt: main.o.fxpt $(LIB_HELIB1) $(LIB_HELIB2).fxpt Makefile
	$(CXX) -o $(EXEC_HELIB)_fxpt main.o.fxpt $(LIB_HELIB2).fxpt $(LIB_HELIB1) $(LDFLAGS) $(LDFLAGS_HELIB)
$(EXEC_HELIB)_binary: main.o.binary $(LIB_HELIB1) $(LIB_HELIB2).binary Makefile
	$(CXX) -o $(EXEC_HELIB)_binary main.o.binary $(LIB_HELIB2).binary $(LIB_HELIB1) $(LDFLAGS) $(LDFLAGS_HELIB)
$(EXEC_HELIB)_ulong: main.o.ulong $(LIB_HELIB1) $(LIB_HELIB2).ulong Makefile
	$(CXX) -o $(EXEC_HELIB)_ulong main.o.ulong $(LIB_HELIB2).ulong $(LIB_HELIB1) $(LDFLAGS) $(LDFLAGS_HELIB)

tests: all
	./$(EXEC_FSHE)_mpz
	./$(EXEC_FSHE)_binary
	./$(EXEC_FSHE)_ulong
	./$(EXEC_HELIB)_fxpt
	./$(EXEC_HELIB)_binary
	./$(EXEC_HELIB)_ulong

.PHONY: clean 
clean:
	rm -rf main.o* $(EXEC_FSHE)* $(EXEC_HELIB)* *.o
