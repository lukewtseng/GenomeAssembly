#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <vector>
#include <list>
#include <set>
#include <unordered_map>
#include <numeric>
#include <cstddef>
#include <mpi.h>
#include <pthread.h>
#include <tbb/task_group.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>

#include "kmer_t.hpp"
#include "hashmap_mpi.hpp"
#include "read_kmers.hpp"

template <typename ...Args>
void print(std::string format, Args... args) {
	fflush(stdout);
	//if (upcxx::rank_me() == 0) {
		printf(format.c_str(), args...);
	//}
	fflush(stdout);
}

void task1(mpi_hashmap &hashmap, tbb::concurrent_vector<kmer_pair>& start_nodes, kmer_pair& kmer) {
		bool success = hashmap.insert(kmer);
		if (!success) {
			throw std::runtime_error("Error: HashMap is full!");
		}
		if (kmer.backwardExt() == 'F') {
			start_nodes.push_back(kmer);
		}
}



int main(int argc, char **argv) {

	int n_proc = 1, rank = 0;
	//MPI_Init( &argc, &argv );
	//MPI_Comm_size( MPI_COMM_WORLD, &n_proc );			    
	//MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	std::string kmer_fname = std::string(argv[1]);
	std::string run_type = "";

	if (argc >= 3) {
		run_type = std::string(argv[2]);
	}
		
	int ks = kmer_size(kmer_fname);

				  
	if (ks != KMER_LEN) {			    
		throw std::runtime_error("Error: " + kmer_fname + " contains " +
				std::to_string(ks) + "-mers, while this binary is compiled for " +
				std::to_string(KMER_LEN) + "-mers.  Modify packing.hpp and recompile.");
	}
	
	size_t n_kmers = line_count(kmer_fname);
	
	size_t hash_table_size = n_kmers * (1.0 / 0.5);
	mpi_hashmap hashmap(hash_table_size);
	
	if (run_type == "verbose") {
		print("Initializing hash table of size %d for %d kmers.\n",
		hash_table_size, n_kmers);
	}

	//std::vector <kmer_pair> kmers = read_kmers(kmer_fname, n_proc, rank);
	std::vector <kmer_pair> kmers = read_kmers(kmer_fname);
	
	if (run_type == "verbose") {
		print("Finished reading kmers.\n");
	}

	//for (auto& k : kmers) 
	//	k.print();

	auto start = std::chrono::high_resolution_clock::now();

	tbb::concurrent_vector <kmer_pair> start_nodes;
	for (auto &kmer : kmers) {
		bool success = hashmap.insert(kmer);
		if (!success) {
			throw std::runtime_error("Error: HashMap is full!");
		}
		if (kmer.backwardExt() == 'F') {
			start_nodes.push_back(kmer);
		}
	}
	print("n_kmers: %zu hashmap size:%zu\n", n_kmers, hashmap.size());

	auto end_insert = std::chrono::high_resolution_clock::now();
	
	double insert_time = std::chrono::duration <double> (end_insert - start).count();
	
	if (run_type != "test") {
		print("Finished inserting in %lf\n", insert_time);
	}


	auto start_read = std::chrono::high_resolution_clock::now();
	std::list <std::list <kmer_pair>> contigs;
	
	for (const auto &start_kmer : start_nodes) {
		std::list <kmer_pair> contig;
		contig.push_back(start_kmer);
		while (contig.back().forwardExt() != 'F') {
			kmer_pair kmer;
			bool success = hashmap.find(contig.back().next_kmer(), kmer);
			if (!success) {
				throw std::runtime_error("Error: k-mer not found in hashmap.");
			}
			contig.push_back(kmer);
		}
		contigs.push_back(contig);
	}
	
	auto end_read = std::chrono::high_resolution_clock::now();
	
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration <double> read = end_read - start_read;
	std::chrono::duration <double> insert = end_insert - start;
	std::chrono::duration <double> total = end - start;


	int numKmers = std::accumulate(contigs.begin(), contigs.end(), 0,
			[] (int sum, const std::list <kmer_pair> &contig) {
				return sum + contig.size(); 
			});

	  if (run_type != "test") {
			print("Assembled in %lf total\n", total.count());
		}
		
		if (run_type == "verbose") {
			printf("Rank %d reconstructed %d contigs with %d nodes from %d start nodes."
					" (%lf read, %lf insert, %lf total)\n", rank, contigs.size(), numKmers, start_nodes.size(), read.count(), insert.count(), total.count());
		}
			  
		if (run_type == "test") {
			std::ofstream fout("test_" + std::to_string(rank) + ".dat");
			for (const auto &contig : contigs) {
				fout << extract_contig(contig) << std::endl;
			}
			fout.close();
		}


	return 0;
}
