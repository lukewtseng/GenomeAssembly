#
# Edison - NERSC 
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
OPENMP = -openmp #Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -g -std=c++11
LIBS =

TARGETS = mpi

all:	$(TARGETS)

mpi: kmer_mpi.cpp kmer_t.hpp pkmer_t.hpp packing.hpp read_kmers.hpp butil.hpp
	$(MPCC) kmer_mpi.cpp -o $@ $(LIBS) $(MPILIBS) $(CFLAGS)



bin.o: bin.cpp bin.h
	$(CC) -c $(CFLAGS) bin.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

