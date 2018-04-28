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

#include "kmer_t.hpp"
#include "hashmap_mpi.hpp"
#include "read_kmers.hpp"

/*
template <typename ...Args>
void print(std::string format, Args... args) {
	fflush(stdout);
	printf(format.c_str(), args...);
	fflush(stdout);
}
*/


int main(int argc, char **argv) {

	int n_proc = 1, rank = 0;
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &n_proc );			    
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	//print("World size: %d Rank:%d\n", n_proc, rank);
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
	
	if (rank == 0) {
		print("Total number of kmers: %d\n", n_kmers);
	}

	int bufsize = 2 * n_kmers * sizeof(MYMPI_Msg) / n_proc;
	void *bsend_buf = malloc(bufsize);
	MPI_Buffer_attach(bsend_buf, bufsize);
	
	size_t hash_table_size = n_kmers * (1.0 / 0.5);
	
	MYMPI_Hashmap hashmap(hash_table_size, n_proc, rank);
	
	if (run_type == "verbose" && rank == 0) {
		print("Initializing hash table of size %d for %d kmers.\n",
		hash_table_size, n_kmers);
	}

	std::vector <kmer_pair> kmers = read_kmers(kmer_fname, n_proc, rank);
	
	if (run_type == "verbose" && rank == 0) {
		print("Finished reading kmers.\n");
	}

	auto start = std::chrono::high_resolution_clock::now();

	std::vector <kmer_pair> start_nodes;
	uint64_t outgoing = 0;
	for (auto &kmer : kmers) {
		hashmap.sync_insert();

		bool success = hashmap.insert(kmer, outgoing);
		if (!success) {
			throw std::runtime_error("Error: HashMap is full!");
		}
		if (kmer.backwardExt() == 'F') {
			start_nodes.push_back(kmer);
		}
	}	

	print ("\t rank %d read %zu kmers, local insert %zu remote insert %zu\n", rank, kmers.size(), hashmap.size(), outgoing);
	
	uint64_t rsize = 0;
	while (rsize < n_kmers) {
			if (rank == 0)
				print ("Hashmap size: %zu Reduce size: %zu\n", hashmap.size(), rsize);
			auto size = hashmap.size();
			MPI_Allreduce(&size, &rsize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
			hashmap.sync_insert();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	assert(rsize == n_kmers);
	if (rank == 0)
		print ("Total reduce size %zu\n", rsize);
	print("Rank: %d n_kmers: %zu hashmap size:%zu\n", rank, n_kmers, hashmap.size());
	
	auto end_insert = std::chrono::high_resolution_clock::now();
	
	double insert_time = std::chrono::duration <double> (end_insert - start).count();
	
	if (run_type != "test" && rank == 0) {
		print("Finished inserting in %lf\n", insert_time);
	}


	auto start_read = std::chrono::high_resolution_clock::now();
	std::vector <std::list <kmer_pair>> contigs;
	
	for (const auto &start_kmer : start_nodes) {
		contigs.emplace_back(std::list<kmer_pair> (1, start_kmer));
	}

	int total_done = 0;
	uint64_t done_quest = 0;
	uint64_t quest = start_nodes.size();
	bool ready[quest];
	std::fill_n(ready, quest, true);
	int index;
	
	while (total_done < n_proc) {
		//Check incoming MPI message
		hashmap.update(contigs, total_done, ready);
		
		if (done_quest < quest) {
			for (int i = 0; i < contigs.size(); i ++) {	
				if (ready[i]) {
					if (contigs[i].back().forwardExt() == 'F') {
						ready[i] = false;
						done_quest ++;
					
						if (done_quest == quest) {
							MYMPI_Msg pack;
							MPI_Request request;
							for (int target = 0; target < n_proc; target ++) {
								if (target != rank) {
									MPI_Ibsend(&pack, sizeof(MYMPI_Msg), MPI_BYTE, target, static_cast<int>(Type::done), MPI_COMM_WORLD, &request);
								}
							}
							total_done ++;
							
							break;
						} else {
							continue;
						}
					}

					kmer_pair kmer;
					bool success = hashmap.find(contigs[i].back().next_kmer(), kmer, ready, i);	
					if (success) {
						contigs[i].emplace_back(kmer);
					}
				}
			}
		}
	}

	MPI_Barrier( MPI_COMM_WORLD );
	/*
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
	*/
	
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
	
	
	MPI_Buffer_detach(&bsend_buf, &bufsize);
	free(bsend_buf);
	MPI_Finalize( );

	return 0;
}
