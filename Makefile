#
# Edison - NERSC 
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
OPENMP = -openmp 
#Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -dynamic -g -std=c++11
LIBS = -tbb

TARGETS = kmer_hash
DATA_DIR = /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets
DATA = test
#DATA = human-chr14-synthetic


all: $(TARGETS)

kmer_hash: kmer_mpi.cpp hashmap_mpi.hpp kmer_t.hpp pkmer_t.hpp packing.hpp read_kmers.hpp butil.hpp
	$(MPCC) kmer_mpi.cpp -o $@ $(LIBS) $(MPILIBS) $(CFLAGS)

bin.o: bin.cpp bin.h
	$(CC) -c $(CFLAGS) bin.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

.ONESHELL:
run:
	module load tbb
	salloc -N 2 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 2 -n 8 ./kmer_hash $(DATA_DIR)/$(DATA).txt test

check:
	cat test*.dat | sort > my_solution.txt
	diff my_solution.txt $(DATA_DIR)/$(DATA)_solution.txt

.ONESHELL:
runsmall:
	module load tbb
	salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 1 ./kmer_hash $(DATA_DIR)/smaller/$(DATA).txt test

checksmall:
	cat test*.dat | sort > my_solution.txt
	diff my_solution.txt $(DATA_DIR)/smaller/$(DATA)_solution.txt

clean:
	rm -f *.o $(TARGETS) *.dat *.txt

