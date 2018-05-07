#pragma once 

#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <mpi.h>
#include "kmer_t.hpp"

//using namespace std;

template <typename ...Args>
void print(std::string format, Args... args) {
	fflush(stdout);
	printf(format.c_str(), args...);
	fflush(stdout);
}

struct MYMPI_Hashmap {

	// Global MPI metadata
	int n_proc;
	int rank;
	
	// hashtable metadata 
	uint64_t total_size;
	uint64_t local_size;

	// Main data structure
	std::unordered_multimap<uint64_t, kmer_pair> table;

	// Constuctor
	MYMPI_Hashmap(uint64_t size);
	MYMPI_Hashmap(uint64_t size, int npc, int rk);
	
	// Functions
	bool insert(const kmer_pair &kmer, uint64_t& outgoing);
	bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer, bool* ready, uint64_t index);

	uint64_t size();

	// Luke's MPI functionality
	void sync_find(std::vector<std::list<kmer_pair>>& contigs, uint64_t& total_done, bool* ready);
	void sync_insert();
	
	// Unused original functionality from hw3
	void write_slot(uint64_t slot, const kmer_pair &kmer);
	kmer_pair read_slot(uint64_t slot);
	bool request_slot(uint64_t slot);
	bool slot_used(uint64_t slot);

};

// Denote the message type
enum class Type {
	insert,
	request, 
	response,
	done
};

// Standardize message package for MPI
struct MYMPI_Msg {
	uint64_t idx;
	kmer_pair kmer;
	pkmer_t key_kmer;
};


void broadcast_done(int n_proc, int rank);


MYMPI_Hashmap::MYMPI_Hashmap(uint64_t size, int npc, int rk) {
	total_size = size;
	rank = rk;
	n_proc = npc;
	local_size = (size + n_proc - 1) / n_proc;
	table.reserve(local_size);
	if (rank == 0) {
		print ("#### Split %d into world size of %d\n", size, n_proc);
		std::cout << "#### Initialize hashmap with local size of " << local_size << std::endl;
	}	
}

void MYMPI_Hashmap::sync_insert() {
	int flag, byte_count;
	MPI_Status status;
	MYMPI_Msg msg;
	
	do {
		MPI_Iprobe(MPI_ANY_SOURCE, static_cast<int> (Type::insert), MPI_COMM_WORLD, &flag, &status);
		if (!flag) 
			break;
		
		MPI_Get_count(&status, MPI_BYTE, &byte_count);
		MPI_Recv (&msg, byte_count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);	
		
		assert(status.MPI_TAG == static_cast<int>(Type::insert));
		//std::cout << "success pass assertion" << endl;
		auto ret = table.insert(std::make_pair(msg.kmer.hash(), msg.kmer));		
	} while (flag);
}

bool MYMPI_Hashmap::insert(const kmer_pair &kmer, uint64_t& outgoing) {
	uint64_t hash = kmer.hash() % total_size;
	uint64_t local_index = hash % local_size;
	int which_proc = int((hash - local_index) / local_size);
	//std::cout << "hash " << hash << " lindex: " << local_index << " sendto: " << which_proc << std::endl;

	if (which_proc == rank) {
		auto ret = table.insert(std::make_pair(kmer.hash(), kmer));
	} else {
		outgoing ++;
		MPI_Request request;
		MYMPI_Msg msg;
		msg.kmer = kmer;
		MPI_Ibsend(&msg, sizeof(MYMPI_Msg), MPI_BYTE, which_proc, static_cast<int>(Type::insert), MPI_COMM_WORLD, &request);
	}
	return true;
}

void MYMPI_Hashmap::sync_find(std::vector<std::list<kmer_pair>>& contigs, uint64_t& total_done, bool* ready) {
	int flag, byte_count;

	do {
		MPI_Status status;
		MYMPI_Msg msg;
		
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
		if (!flag) 
			break;
		

		MPI_Get_count(&status, MPI_BYTE, &byte_count);
		MPI_Recv (&msg, byte_count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);	

		//assert(msg.act == Type::response || msg.act == Type::done);
		if (status.MPI_TAG == static_cast<int>(Type::done)) {
			//print ("Receive Done msg from %d\n", status.MPI_SOURCE);
			total_done ++;
		
		} else if (status.MPI_TAG == static_cast<int>(Type::response)) {
			//print ("Receive response msg from %d, index: %d\n", status.MPI_SOURCE, msg.idx);
            //TODO
            if(msg.idx < 0 || msg.idx >= local_size){
              print("\t\tWARNING: process %d trying to access index %d\n", rank, msg.idx);
            }

			contigs[msg.idx].emplace_back(msg.kmer);
			//contigs[msg.idx].pop_front();
			
			assert(ready[msg.idx] == false);
			ready[msg.idx] = true;
		
		} else if (status.MPI_TAG == static_cast<int>(Type::request)) {
			//assert(status.MPI_TAG == static_cast<int>(Type::request));
			//print ("Receive request msg from %d\n", status.MPI_SOURCE);
			MYMPI_Msg outgoing_msg;
			outgoing_msg.idx = msg.idx;
			
			auto result = table.equal_range(msg.key_kmer.hash());
			for (auto it = result.first; it != result.second; it++) {
				if (it->second.kmer == msg.key_kmer) {
					outgoing_msg.kmer = std::move(it->second);
					break;
				}
			}
			MPI_Request request;
			MPI_Ibsend(&outgoing_msg, sizeof(MYMPI_Msg), MPI_BYTE, status.MPI_SOURCE, static_cast<int>(Type::response), MPI_COMM_WORLD, &request);
			
		} else {
			std::cout << "\tWarning, receive irrecognizable MPI TAG: " << status.MPI_TAG << std::endl;
			exit(0);
		}
		
	} while (flag);
}

bool MYMPI_Hashmap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer, bool * ready, uint64_t index) {
	uint64_t hash = key_kmer.hash() % total_size;
	uint64_t local_index = hash % local_size;
	int which_proc = int((hash - local_index) / local_size);
    if(index < 0 || index >= local_size){
      print("\t\tWARNING: process %d trying to access index %d\n", rank, index);
    }
	
	//std::cout << "hash " << hash << " lindex: " << local_index << " sendto: " << which_proc << std::endl;
	ready[index] = false;

	if (which_proc == rank) {
		//print ("Search within local process %d\n", which_proc);	
		auto result = table.equal_range(key_kmer.hash());
	
		for (auto it = result.first; it != result.second; it++) {
			if (it->second.kmer == key_kmer) {
				val_kmer = it->second;
				ready[index] = true;
				return true;
			}
		}
		print ("Warning, kmer not found in the local process %d\n", rank);
	}	else {
		MYMPI_Msg msg;
		MPI_Request request;
		msg.key_kmer = key_kmer;
		msg.idx = index;
		MPI_Ibsend(&msg, sizeof(MYMPI_Msg), MPI_BYTE, which_proc, static_cast<int> (Type::request), MPI_COMM_WORLD, &request);
	}

	return false;
}

uint64_t MYMPI_Hashmap::size() {
	return table.size();
}

void broadcast_done(int n_proc, int rank) {
	MYMPI_Msg msg;
	MPI_Request request;
	
	for (int target = 0; target < n_proc; target ++) {
		if (target != rank) {
			MPI_Ibsend(&msg, sizeof(MYMPI_Msg), MPI_BYTE, target, static_cast<int>(Type::done), MPI_COMM_WORLD, &request);
		}
	}
}

