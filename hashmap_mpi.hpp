#pragma once 

#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <mpi.h>
#include <atomic>

#include "kmer_t.hpp"
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_queue.h>

//using namespace std;
struct MYMPI_Msg;
struct PK_Msg;

typedef tbb::concurrent_unordered_multimap<uint64_t,kmer_pair> concurrent_hashmap;

typedef tbb::concurrent_queue<PK_Msg> msg_queue;
typedef tbb::concurrent_queue<MYMPI_Msg> insert_queue;

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
	//std::unordered_multimap<uint64_t, kmer_pair> table;
	concurrent_hashmap table;

	msg_queue messages;
	msg_queue remote_messages;
	insert_queue remote_insert;
	// Constuctor
	MYMPI_Hashmap(uint64_t size);
	MYMPI_Hashmap(uint64_t size, int npc, int rk);
	
	// Functions
	bool insert(const kmer_pair &kmer, uint64_t& outgoing);
	bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer, bool* ready, uint64_t index);

	uint64_t size();

	// Luke's MPI functionality
	void communicate_messages(std::vector<std::list<kmer_pair>>& contigs, std::atomic<int>& total_done, bool* ready);
	void deal_remote_messages(std::vector<std::list<kmer_pair>>& contigs, std::atomic<int>& total_done, bool* ready);
	
	void collect_remote_insert();
	void deal_remote_insert();
	
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

struct PK_Msg {
	MYMPI_Msg data;
	Type tag;
	int target;
};


void broadcast_done(int n_proc, int rank);


MYMPI_Hashmap::MYMPI_Hashmap(uint64_t size, int npc, int rk) {
	total_size = size;
	rank = rk;
	n_proc = npc;
	local_size = (size + n_proc - 1) / n_proc;
	//table.reserve(local_size);
	if (rank == 0) {
		print ("#### Split %d into world size of %d\n", size, n_proc);
		std::cout << "#### Initialize hashmap with local size of " << local_size << std::endl;
	}	
}
void MYMPI_Hashmap::deal_remote_insert() {
	while (1) {
		MYMPI_Msg item;
		if (!remote_insert.try_pop(item))
			break;
		try {
			table.insert(std::make_pair(item.kmer.hash(), item.kmer));
		} catch (std::bad_alloc &ba) {
			std::cerr << "catch bad alloc: " << ba.what() << std::endl;
		}
	}
}

void MYMPI_Hashmap::collect_remote_insert() {
	int flag, byte_count;
	MPI_Status status;
	MYMPI_Msg msg;
	
	do {
		
		while (!messages.empty()) {
			PK_Msg item;
			if (messages.try_pop(item)) {
				MPI_Request request;
				MPI_Ibsend(&item.data, sizeof(MYMPI_Msg), MPI_BYTE, item.target, static_cast<int>(item.tag), MPI_COMM_WORLD, &request);
			}
		}

		MPI_Iprobe(MPI_ANY_SOURCE, static_cast<int> (Type::insert), MPI_COMM_WORLD, &flag, &status);
		if (!flag) 
			break;
	
		MPI_Get_count(&status, MPI_BYTE, &byte_count);
		MPI_Recv (&msg, byte_count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);	
		
		assert(status.MPI_TAG == static_cast<int>(Type::insert));
		//auto ret = table.insert(std::make_pair(msg.kmer.hash(), msg.kmer));		
		remote_insert.emplace(msg);
		deal_remote_insert();

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
		//MPI_Request request;
		PK_Msg pkmsg;
		pkmsg.data.kmer = kmer;
		pkmsg.target = which_proc;
		pkmsg.tag = Type::insert;
		messages.emplace(pkmsg);
		//MYMPI_Msg msg;
		//msg.kmer = kmer;
		//MPI_Ibsend(&msg, sizeof(MYMPI_Msg), MPI_BYTE, which_proc, static_cast<int>(Type::insert), MPI_COMM_WORLD, &request);
	}
	return true;
}
void MYMPI_Hashmap::deal_remote_messages(std::vector<std::list<kmer_pair>>& contigs, std::atomic<int>& total_done, bool* ready) {
	
	while (1) {

	PK_Msg msg;
	if (!remote_messages.try_pop(msg))
		break;

	if (msg.tag == Type::done) {
			//print ("Receive Done msg from %d\n", status.MPI_SOURCE);
			total_done ++;
		
	} else if (msg.tag == Type::response) {
			//print ("Receive response msg from %d, index: %d\n", status.MPI_SOURCE, msg.idx);
			if(msg.data.idx < 0 || msg.data.idx >= local_size) {
				print("\t\tWARNING: process %d trying to access index %d\n", rank, msg.data.idx);
			}
			contigs[msg.data.idx].emplace_back(msg.data.kmer);
			assert(ready[msg.data.idx] == false);
			ready[msg.data.idx] = true;
		
	} else if (msg.tag == Type::request) {
			//print ("Receive request msg from %d\n", status.MPI_SOURCE);
			PK_Msg pkmsg;
			pkmsg.data.idx = msg.data.idx;
			pkmsg.target = msg.target;
			pkmsg.tag = Type::response;
			auto result = table.equal_range(msg.data.key_kmer.hash());
			for (auto it = result.first; it != result.second; it++) {
				if (it->second.kmer == msg.data.key_kmer) {
					pkmsg.data.kmer = it->second;
					break;
				}
			}
			messages.emplace(pkmsg);
	} else {
		std::cout << "\tWarning, receive irrecognizable MPI TAG: " << static_cast<int>(msg.tag) << std::endl;
		exit(0);
	}	
	}
}

void MYMPI_Hashmap::communicate_messages(std::vector<std::list<kmer_pair>>& contigs, std::atomic<int>& total_done, bool* ready) {
	int flag, byte_count;

	do {

		while (!messages.empty()) {
			PK_Msg item;
			if (messages.try_pop(item)) {
				MPI_Request request;
				MPI_Ibsend(&item.data, sizeof(MYMPI_Msg), MPI_BYTE, item.target, static_cast<int>(item.tag), MPI_COMM_WORLD, &request);
			}
		}

		MPI_Status status;
		
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
		if (!flag) 
			continue;
		

		PK_Msg receive;
		MPI_Get_count(&status, MPI_BYTE, &byte_count);
		MPI_Recv (&receive.data, byte_count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);	

		receive.target = status.MPI_SOURCE;
		receive.tag = static_cast<Type>(status.MPI_TAG);
		remote_messages.emplace(receive);
		
		deal_remote_messages(contigs, total_done, ready);
	} while (total_done < n_proc);
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
		PK_Msg pkmsg;
		pkmsg.data.key_kmer = key_kmer;
		pkmsg.data.idx = index;
		pkmsg.target = which_proc;
		pkmsg.tag = Type::request;
		messages.emplace(pkmsg);
	}

	return false;
}

uint64_t MYMPI_Hashmap::size() {
	return table.size();
}

void broadcast_done(int n_proc, int rank, msg_queue& messages) {
	PK_Msg pkmsg;
	//MYMPI_Msg msg;
	MPI_Request request;
	pkmsg.tag = Type::done;
	

	for (int target = 0; target < n_proc; target ++) {
		if (target != rank) {
			//MPI_Ibsend(&msg, sizeof(MYMPI_Msg), MPI_BYTE, target, static_cast<int>(Type::done), MPI_COMM_WORLD, &request);
			pkmsg.target = target;
			messages.emplace(pkmsg);
		}
	}
}

