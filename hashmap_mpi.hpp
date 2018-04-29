#pragma once 

#include <cstdio>
#include <unordered_map>
#include <tbb/task_group.h>
#include <tbb/concurrent_unordered_map.h>
//using namespace std;

using namespace tbb;

typedef tbb::concurrent_unordered_multimap<uint64_t,kmer_pair> concurrent_hashmap;

struct mpi_hashmap {

	// Global MPI metadata
	size_t n_proc;
	int rank;
	
	// hashtable metadata 
	size_t total_size;
	size_t local_size;

	// Main data structure
	//std::unordered_multimap<uint64_t, kmer_pair> table;
	concurrent_hashmap table;

	// Constuctor
	mpi_hashmap(size_t size);
	
	// Functions
	bool insert(const kmer_pair &kmer);
	bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

	size_t size();
	void write_slot(uint64_t slot, const kmer_pair &kmer);
	kmer_pair read_slot(uint64_t slot);

	//bool request_slot(uint64_t slot);
	//bool slot_used(uint64_t slot);

};


mpi_hashmap::mpi_hashmap(size_t size) {
	//table.reserve(size);
}

bool mpi_hashmap::insert(const kmer_pair &kmer) {
	uint64_t hash = kmer.hash();
	auto ret = table.insert(std::make_pair(hash, kmer));
	return true;
	/*
  uint64_t probe = 0;
	bool success = false;
	do {
		uint64_t slot = (hash + probe++) % size();
		success = request_slot(slot);
		if (success) {
			write_slot(slot, kmer);
		}
	} while (!success && probe < size());
	return success;
	*/
}

bool mpi_hashmap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
	uint64_t hash = key_kmer.hash();
	
	//typedef std::unordered_multimap<uint64_t, kmer_pair>::iterator MMAPIterator;
	auto result = table.equal_range(hash);

	for (auto it = result.first; it != result.second; it++) {
		if (it->second.kmer == key_kmer) {
			val_kmer = it->second;
			return true;
		}
	}

	return false;
	/*
	uint64_t probe = 0;
	bool success = false;
	do {
		uint64_t slot = (hash + probe++) % size();
		if (slot_used(slot)) {
			val_kmer = read_slot(slot);
			if (val_kmer.kmer == key_kmer) {
				success = true;
			}
		}
	} while (!success && probe < size());
	return success;
	*/
}
size_t mpi_hashmap::size() {
	return table.size();
}
