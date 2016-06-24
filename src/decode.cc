/*
 * decode.cc
 *
 *  Created on: Nov 15, 2015
 *      Author: abdullah
 */

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

const char kAlphabet[] = {'A', 'C', 'G', 'T'};

string BinaryToKmerCount(string& prefix, uint64_t kmer_bin, int rem_kmer_len, uint64_t count_mask) {
	int count = kmer_bin & count_mask;
	string kmer (prefix);
	int shr = 62;
	for (int i = 0; i < rem_kmer_len; ++i) {
		kmer += kAlphabet[(kmer_bin >> shr) & 3];
		shr -= 2;
	}
	return kmer + '\t' + to_string(count);
}

string BinaryToKmer(uint64_t kmer_bin, int kmer_len) {
	string kmer;
	for (int i = (kmer_len - 1) << 1; i >= 0; i -= 2)
		kmer += kAlphabet[(kmer_bin >> i) & 3];
	return kmer;
}


int main(int argc, char** argv) { // ./bin/kcmbt_dump #file
	if (argc != 2) {
		cerr << "./bin/kcmbt_dump number_of_threads_used_in_kcmbt" << endl;
		return EXIT_FAILURE;
	}

	int total_file = atoi(argv[1]);

	ofstream out_file("kmer_list.txt");
	for (int i = 0; i < total_file; ++i) {
		// read prefix binary file
		string pre_name("pre_out_" + to_string(i));
		cout << pre_name << endl;
		FILE* pre_file = fopen(pre_name.c_str() , "rb");

		uint64_t stat;
		fread(&stat, sizeof(uint64_t), 1, pre_file); // first line contains some stats about kmer_len (32 bit) and tree_pow (32 bit)
		int kmer_len = stat >> 32; // k value of kmer
		int tree_pow = stat & 0xFFFFFFFF; // tree power; # of trees = 4^tree_pow


		int total_tree = 1 << (tree_pow << 1);
		int tree_mask = total_tree - 1;
		int rem_kmer_len = kmer_len - tree_pow;
		uint64_t count_mask = (1ULL << (64 - 2 * rem_kmer_len)) - 1;

		//cout << total_tree << " " << tree_mask << " " << rem_kmer_len << " " << count_mask << endl;
		vector<uint64_t> count_arr (total_tree);
		vector<string> prefix_arr (total_tree);
		fread(&count_arr[0], sizeof(uint64_t), total_tree, pre_file);
		for (int i = 0; i < total_tree; ++i)
			for (int j = 2 * tree_pow - 2; j >= 0; j -= 2)
				prefix_arr[i] += kAlphabet[(i >> j) & 3];

		fclose(pre_file);

		cout << "done prefix" << endl;

		// read suffix & count binary file
		string suf_name("suf_out_" + to_string(i));
		FILE* kmer_file = fopen(suf_name.c_str(), "rb");
		uint64_t total_kmer = 0;
		for (int i = 0; i < total_tree; ++i) {
			if (count_arr[i] == 0)
				continue;
			total_kmer += count_arr[i];
			//cout << i << ":\t" << count_arr[i] << "\t" << total_kmer << endl;
			uint64_t* rem_kmer_arr = new uint64_t[count_arr[i]];
			string* kmer_arr = new string[count_arr[i]];

			fread(rem_kmer_arr, sizeof(uint64_t), count_arr[i], kmer_file);
			string out_str;
			for (int j = 0; j < count_arr[i]; ++j)
				out_str += BinaryToKmerCount(prefix_arr[i], rem_kmer_arr[j], rem_kmer_len, count_mask) + "\n";
			out_file << out_str;

			delete[] kmer_arr;
			delete[] rem_kmer_arr;
		}

		fclose(kmer_file);

		cout << "done suffix" << endl;

		string big_name("big_out_" + to_string(i));
		FILE* big_file = fopen(big_name.c_str(), "rb");
		string big_kmer_str;
		uint64_t big_kmer[2];
		while (fread(big_kmer, sizeof(uint64_t), 2, big_file)) {
			big_kmer_str += BinaryToKmer(big_kmer[0], kmer_len);
			big_kmer_str += "\t" + to_string(big_kmer[1]) + "\n";
		}
		out_file << big_kmer_str;
		fclose(big_file);

		cout << "done big" << endl;
	}

	out_file.close();

	return 0;
}




