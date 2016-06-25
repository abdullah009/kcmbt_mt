

#include <iostream>
#include <chrono>
#include <zlib.h>
#include <stdio.h>
#include <ctime>
#include <fstream>
#include <sys/resource.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iterator>
#include <thread>

#include "burst_sort_kmer.h"
#include "kseq.h"
#include "fastq_reader.h"

using namespace std;


// initialize reader for gzip file
KSEQ_INIT(gzFile, gzread)

// using directive
using Clock = std::chrono::high_resolution_clock;
using TimePoint = Clock::time_point;

// constant
const int kMaxBuff = 1 << 10;

// global variables
uint64_t* all_count;
uint64_t* uniq_count; //
uint64_t** all_count_arr;
TrieNode*** root;

string in_file_name, out_file_name ("out");


void HowToUse() {
	cout << "\n\n\n";

	string h("KCMBT (k-mer Counter based on Multiple Burst Trees) 1.0.0");
	cout << h << endl;
	vector<char> c_arr(h.size(), '=');
	copy(c_arr.begin(), c_arr.end(), ostream_iterator<char>(cout));
	cout << endl;

	cout << "Usage: \n\t" << "./bin/kcmbt -k <k-mer length> -i <@file_listing_fastq_files or fastq_file> -t <number_of_threads>" << endl;
	cout << "Example: \n\t" << "./bin/kcmbt -k 28 -i srr.fastq -t 4" << endl;
	cout << "Parameters: \n";
	cout << "\t-k\t\tk-mer length (10 <= k <= 32, default 28)\n";
	cout << "\t-i\t\tinput file in fastq format (start with @ if the file contains a list of fastq files)\n";
	//cout << "\t-o\t\toutput file (a binary file; please run hr_kcmbt to generate human readable output)\n";
	cout << "\t-t\t\tnumber of threads (please use 2^x threads, x = 0, , 2, 3, ..)\n";

	cout << "\n\n\n";
}

void SortKmer(int tree_ind, int pre_ind, FILE* out_file, int th_ind, int total_thread) {
	uint64_t kmer_uniq = 0;
	uint64_t pos_arr[total_thread + 1];
	pos_arr[0] = 0;
	for (int i = 0; i < total_thread; ++i) {
		Traverse(root[i][tree_ind], 0, kmer_uniq, th_ind);
		pos_arr[i + 1] = kmer_uniq;
	}
    //cout << th_ind << "\t" << kmer_uniq << " ";
	int len = pos_arr[total_thread];
	MergeKArr(pos_arr, th_ind, total_thread);
	kmer_uniq = CompactKmer(&kmer_arr[th_ind][0], len, th_ind, tree_ind);
	uniq_count[th_ind] += kmer_uniq;	
	pre_arr[th_ind][pre_ind] = (kmer_uniq | (uint64_t)tree_ind << 32);
    
	uint64_t written = fwrite(kmer_arr[th_ind], sizeof(uint64_t), kmer_uniq, out_file);
    
	//cout << th_ind << " " << tree_ind << " " << kmer_uniq << " written " << written << endl;
}

uint64_t ccc[32];
void SortInsertKmer(int th_ind, int tree_type, int ind, uint64_t* buffer) {
	uint64_t kmer_uniq = 0;
	Traverse(root[th_ind][tree_type * total_tree + ind], 0, kmer_uniq, th_ind);
	ccc[th_ind] += kmer_uniq;
	//delete[] root[th_ind][tree_type * total_tree + ind];
	//root[th_ind][tree_type * total_tree + ind] = nullptr;
	uint64_t k_len =  2 * (kmer_len - tree_pow);
	uint64_t mask_kmer = ((1ULL << k_len) - 1) << (64 - k_len);
	uint64_t tree_mask = total_tree - 1;

	uint64_t pos_arr[total_tree];
	std::memset(pos_arr, 0, total_tree * sizeof(uint64_t));

	for (int i = 0; i < kmer_uniq; ++i) {
		uint64_t k_c = kmer_arr[th_ind][i] & k_count_mask[tree_type];
		uint64_t k_p = kmer_arr[th_ind][i] & k_kmer_mask[tree_type];
		for (int j = 0; j <= tree_type; ++j) {
			uint64_t kmer = ((k_p << 2 * j) & mask_kmer) | k_c;
			int tree_ind = ((ind << 2 * j) | (k_p >> (64 - 2 * j))) & tree_mask;
			++all_count_arr[th_ind][tree_ind];
			buffer[tree_ind * kMaxBuff + (pos_arr[tree_ind]++)] = kmer;
			if (pos_arr[tree_ind] >= kMaxBuff) {
				InsertBatch(root[th_ind][tree_ind], &buffer[tree_ind * kMaxBuff], pos_arr[tree_ind]);
				pos_arr[tree_ind] = 0;
			}
		}
	}
	for (uint64_t i = 0; i < total_tree; ++i) 
		InsertBatch(root[th_ind][i], &buffer[i * kMaxBuff], pos_arr[i]);
}

bool Compare(const BigKmer& k1, const BigKmer& k2) {
	return k1.kmer < k2.kmer;
}

void CountBigKmer(int th_ind) {
	int n = bkmer_arr[th_ind].size(), j = 0;
	//cout << th_ind << " CountBigKmer " << n << endl;
	uint64_t* bkmerw_arr = new uint64_t[2 * n];
	sort(bkmer_arr[th_ind].begin(), bkmer_arr[th_ind].end(), Compare);
	for (int i = 0; i < n; ++i) {
		bkmerw_arr[j++] = bkmer_arr[th_ind][i].kmer;
		bkmerw_arr[j] = bkmer_arr[th_ind][i].count + k_count_mask[0];
		while (i + 1 < n && bkmer_arr[th_ind][i].kmer == bkmer_arr[th_ind][i + 1].kmer) {
			bkmerw_arr[j] += bkmer_arr[th_ind][i + 1].count;
			++i;
		}
		++j;
	}
	string big ("big_" + out_file_name + "_" + to_string(th_ind));
	FILE* big_file = fopen(big.c_str(), "wb");
	fwrite(bkmerw_arr, sizeof(uint64_t), j, big_file);
	fclose(big_file);

#ifdef DEBUG
	cout << n << " big kmer count " << j << endl;
#endif
	uniq_count[th_ind] += (j >> 1);

	delete[] bkmerw_arr;
}




void ComputeKmer(int th_ind, int total_thread, vector<int>& file_ind_arr, vector<string>& in_file_arr, vector<vector<pair<long, long>>>& f_pos_arr) {
	// mapping letters to int
	char codes[256];
	for(int i = 0; i < 256; ++i)
		codes[i] = -1;
	codes['A'] = codes['a'] = 0;
	codes['C'] = codes['c'] = 1;
	codes['G'] = codes['g'] = 2;
	codes['T'] = codes['t'] = 3;
	// ===

	temp_mem[th_ind] = new uint64_t[total_thread * kMaxBucketSize]; // temporary memory for sorting

	uint64_t total_entry = kTotalLayer * total_tree;
	uint64_t* buffer = new uint64_t[total_entry * kMaxBuff];// __attribute__((aligned(64))); // kmer buffers, while full, insert into trees
	int pos_arr[total_entry]; // pos in buffer
	int kmer_count_arr[total_entry]; // count kmers in each tree

	all_count_arr[th_ind] = new uint64_t[total_tree];
	for (int j = 0; j < total_tree; ++j)
		all_count_arr[th_ind][j] = 0;
	
	root[th_ind] = new TrieNode*[total_entry];
	for (int i = 0; i < kTotalLayer; ++i) {
		for (int j = 0; j < total_tree; ++j) {
			uint64_t ind = i * total_tree + j;
			root[th_ind][ind] = new TrieNode(i, th_ind);
			root[th_ind][ind]->InitRoot();
			pos_arr[ind] = 0;
			kmer_count_arr[ind] = 0;
		}
	}

	uint64_t k_mask_arr[kAlphabetSize], r_shift_arr[kAlphabetSize], l_shift_arr[kAlphabetSize];
	for (int i = 0; i < kAlphabetSize; ++i) {
		k_mask_arr[i] = ~0ULL >> ((32 - kmer_len - i) << 1); // mask for kmer+x left aligned
		r_shift_arr[i] = (kmer_len + i - tree_pow) << 1; // right shift for tree kmer+x
		l_shift_arr[i] = (32 - kmer_len - i + tree_pow) << 1; // left shift for tree kmer+x
	}

	uint64_t kmer, kmer_rev, kmer_bit;
	uint64_t len, kmer_count = 0, tree_ind, c_layer, r_pos;
	uint64_t prev_k_less;

	for (int f = 0; f < file_ind_arr.size(); ++f) { // read each file and generate kmers
		gzFile fp;
		kseq_t *seq;
		//cout << th_ind << ":" << in_file_arr[file_ind_arr[f]] << endl;
		fp = gzopen(in_file_arr[file_ind_arr[f]].c_str(), "r");
		seq = kseq_init(fp);

		while (kseq_read(seq) >= 0) { // read read sequences
			prev_k_less = kPrevInit; // 2 means unassigned
			r_pos = 0;
			c_layer = 0;
			len = 0;
			kmer = 0ULL;
			kmer_rev = 0ULL;

			char *s = seq->seq.s;
			for (int i=0; i<seq->seq.l; ++i, ++s) {
				int c = codes[*s];
				if (c < 0) { //'N'
					if (c_layer > 0) {
						c_layer--;
						bool b = prev_k_less == 1? true : false; // true means kmer is minimal, false means rev kmer is small
						tree_ind = b * (kmer >> r_shift_arr[c_layer]) + !b * (kmer_rev >> r_shift_arr[c_layer]);
						tree_ind &= mask;
						kmer_bit = b * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !b * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount);
						buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
						if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { // does not affect that much
							InsertBatch(root[th_ind][(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
							kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
							pos_arr[c_layer * total_tree + tree_ind] = 0;
						}
					}

					prev_k_less = kPrevInit;
					r_pos = 0;
					c_layer = 0;
					len = 0;
					kmer = 0ULL;
					kmer_rev = 0ULL;
					continue;
				}
				++len;
				kmer = (kmer << 2) | c;
				kmer_rev = kmer_rev | ((3ULL - c) << r_pos);
				r_pos += 2;
				if (len >= kmer_len) {
					uint64_t kmer0 = kmer & k_mask;
					uint64_t kmer_rev0 = (kmer_rev >> 2 * c_layer) & k_mask;
					bool b = kmer0 < kmer_rev0;
					bool cond1 = (b && prev_k_less == 1) || (!b && prev_k_less == 0); // same kmer as before
					bool cond2 = prev_k_less > 1 || (c_layer < kTotalLayerMinusOne && cond1 ); // check whether continue with this extended kmer
					c_layer += cond2;

					if (!cond2) {
						c_layer += cond1;
						c_layer--;

						tree_ind = cond1 * (prev_k_less * (kmer >> r_shift_arr[c_layer])   + !prev_k_less * (kmer_rev >> r_shift_arr[c_layer]))
								+ !cond1 * (prev_k_less * ((kmer >> r_shift_arr[c_layer] + 2))  + !prev_k_less * (kmer_rev >> r_shift_arr[c_layer]));
						tree_ind &= mask;
						kmer_bit = cond1 * (prev_k_less * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !prev_k_less * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount))
								+ !cond1 * (prev_k_less * (((kmer << (l_shift_arr[c_layer] - 2)) & k_kmer_mask[c_layer]) | kInitCount) + !prev_k_less * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount));
						buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
						all_count_arr[th_ind][tree_ind] += (c_layer == 0);
						if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { 
							InsertBatch(root[th_ind][(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
							kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
							pos_arr[c_layer * total_tree + tree_ind] = 0;
						}

						kmer = kmer & k_mask;
						kmer_rev = (kmer_rev >> (c_layer + 1) * 2) & k_mask;
						bool cond3 = c_layer == kTotalLayerMinusOne;
						c_layer = 1 - cond3;
						r_pos = (kmer_len - cond3) << 1;
						prev_k_less = b + cond3 * 2;
					}
					else
						prev_k_less = b;

					++kmer_count;
				}
			}

			if (c_layer > 0) {
				c_layer--;
				bool b = prev_k_less == 1? true : false;
				tree_ind = b * (kmer >> r_shift_arr[c_layer]) + !b * (kmer_rev >> r_shift_arr[c_layer]);
				tree_ind &= mask;
				kmer_bit = b * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !b * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount);
				buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
				if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { 
					InsertBatch(root[th_ind][(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
					kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
					pos_arr[c_layer * total_tree + tree_ind] = 0;
				}
			}
		}

		kseq_destroy(seq);
		gzclose(fp);
		//putchar('a' + th_ind * total_thread + f);
		fflush(stdout);
	}

	//cout << file_ind_arr.size() * total_thread << endl;

	// file split	
	char* s = (char*) malloc((1 << 13)* sizeof(char));
	char* seq = s;
	for (int f = file_ind_arr.size() * total_thread; f < in_file_arr.size(); ++f) { // read each file and generate kmers
		FastqReader fr;
		fr.Init(in_file_arr[f].c_str(), f_pos_arr[f][th_ind].first, f_pos_arr[f][th_ind].second);
		
		//cout << th_ind << ":" << in_file_arr[f] << " " << f_pos_arr[f][th_ind].first << " " <<  f_pos_arr[f][th_ind].second << endl;

		s = seq;
		while (fr.ReadSeq(s)) { // read read sequences
			prev_k_less = kPrevInit; // 2 means unassigned
			r_pos = 0;
			c_layer = 0;
			len = 0;
			kmer = 0ULL;
			kmer_rev = 0ULL;
			for (int i = 0; i < fr.len; ++i, ++s) {
				int c = codes[*s];
				if (c < 0) { //'N'
					if (c_layer > 0) {
						c_layer--;
						bool b = prev_k_less == 1? true : false; // true means kmer is minimal, false means rev kmer is small
						tree_ind = b * (kmer >> r_shift_arr[c_layer]) + !b * (kmer_rev >> r_shift_arr[c_layer]);
						tree_ind &= mask;
						kmer_bit = b * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !b * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount);
						buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
						if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { // does not affect that much
							InsertBatch(root[th_ind][(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
							kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
							pos_arr[c_layer * total_tree + tree_ind] = 0;
						}
					}

					prev_k_less = kPrevInit;
					r_pos = 0;
					c_layer = 0;
					len = 0;
					kmer = 0ULL;
					kmer_rev = 0ULL;
					continue;
				}
				++len;
				kmer = (kmer << 2) | c;
				kmer_rev = kmer_rev | ((3ULL - c) << r_pos);
				r_pos += 2;
				if (len >= kmer_len) {
					uint64_t kmer0 = kmer & k_mask;
					uint64_t kmer_rev0 = (kmer_rev >> 2 * c_layer) & k_mask;
					bool b = kmer0 < kmer_rev0;
					bool cond1 = (b && prev_k_less == 1) || (!b && prev_k_less == 0); // same kmer as before
					bool cond2 = prev_k_less > 1 || (c_layer < kTotalLayerMinusOne && cond1 ); // check whether continue with this extended kmer
					c_layer += cond2;

					if (!cond2) {
						c_layer += cond1;
						c_layer--;

						tree_ind = cond1 * (prev_k_less * (kmer >> r_shift_arr[c_layer])   + !prev_k_less * (kmer_rev >> r_shift_arr[c_layer]))
								+ !cond1 * (prev_k_less * ((kmer >> r_shift_arr[c_layer] + 2))  + !prev_k_less * (kmer_rev >> r_shift_arr[c_layer]));
						tree_ind &= mask;
						kmer_bit = cond1 * (prev_k_less * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !prev_k_less * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount))
								+ !cond1 * (prev_k_less * (((kmer << (l_shift_arr[c_layer] - 2)) & k_kmer_mask[c_layer]) | kInitCount) + !prev_k_less * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount));
						buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
						all_count_arr[th_ind][tree_ind] += (c_layer == 0);
						if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { 
							InsertBatch(root[th_ind][(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
							kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
							pos_arr[c_layer * total_tree + tree_ind] = 0;
						}

						kmer = kmer & k_mask;
						kmer_rev = (kmer_rev >> (c_layer + 1) * 2) & k_mask;
						bool cond3 = c_layer == kTotalLayerMinusOne;
						c_layer = 1 - cond3;
						r_pos = (kmer_len - cond3) << 1;
						prev_k_less = b + cond3 * 2;
					}
					else
						prev_k_less = b;

					++kmer_count;
				}
			}

			if (c_layer > 0) {
				c_layer--;
				bool b = prev_k_less == 1? true : false;
				tree_ind = b * (kmer >> r_shift_arr[c_layer]) + !b * (kmer_rev >> r_shift_arr[c_layer]);
				tree_ind &= mask;
				kmer_bit = b * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !b * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount);
				buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
				if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { 
					InsertBatch(root[th_ind][(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
					kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
					pos_arr[c_layer * total_tree + tree_ind] = 0;
				}
			}
			s = seq;
		}

		fr.Destroy();	
		//putchar('a' + th_ind);
		fflush(stdout);
		s = seq;
	}
	free(seq);

	uint64_t max_len = 0;
	for (uint64_t i = 0; i < kTotalLayer; ++i) {
		for (uint64_t j = 0; j < total_tree; ++j) {
			InsertBatch(root[th_ind][i * total_tree + j], &buffer[(i * total_tree + j) * kMaxBuff], pos_arr[i * total_tree + j]);
			kmer_count_arr[i * total_tree + j] += pos_arr[i * total_tree + j];
			if (kmer_count_arr[i * total_tree + j] > max_len)
				max_len = kmer_count_arr[i * total_tree + j];
		}
	}



	ccc[th_ind] = 0;
	all_count[th_ind] = kmer_count;

	uint64_t t_sum = 0;
	for (uint64_t i = 0; i < kTotalLayer; ++i) {
		uint64_t t_s = 0;
		for (uint64_t j = 0; j < total_tree; ++j) {
				t_s += kmer_count_arr[i * total_tree + j];
		}
		//std::cout << i << ": " << t_s << std::endl;
		t_sum += (i + 1) * t_s;
	}

	max_len *= (total_thread * 3);
	if (max_len > 300 * 1024 * 1024)
		max_len = 350 * 1024 * 1024;

	kmer_arr[th_ind] = new uint64_t[max_len];

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	//cout << th_ind << " fi: " << file_ind_arr.size() << "\tmem used: " << kmer_count << "\t" << max_len << " count " << all_count[th_ind] << endl;



	for (uint64_t i = 1; i < kTotalLayer; ++i) {
		for (uint64_t j = 0; j < total_tree; ++j) {
			SortInsertKmer(th_ind, i, j, buffer);
		}
#ifdef DEBUG
		putchar(i + 'a');
		fflush(stdout);
#endif
	}

	//cout << th_ind <<   "  ccc " << ccc[th_ind] << endl;

	delete[] buffer;
}


void TraverseKmer(int th_ind, int total_thread, vector<int>& order) {
	uniq_count[th_ind] = 0;
	//return;
	int tree_per_thread = ceil(1.0 * total_tree / total_thread);
	int s = th_ind * tree_per_thread;
	int e = s + tree_per_thread > total_tree? total_tree : s + tree_per_thread;
	
	pre_arr[th_ind] = new uint64_t[total_tree + 1];
	pre_arr[th_ind][0] = (kmer_len << 32) | tree_pow;
	for (int i = 1; i <= total_tree; ++i)
		pre_arr[th_ind][i] = 0;

	string suf ("suf_" + out_file_name + "_" + to_string(th_ind));
	FILE* out_file = fopen(suf.c_str(), "wb");
	if (!out_file)
		cerr << "ERROR" << endl;
	int tree_c = 0;	
	for (int i = 0; i < order.size(); i++)
		SortKmer(order[i], i + 1, out_file, th_ind, total_thread);
	fclose(out_file);

	string pre ("pre_" + out_file_name + "_" + to_string(th_ind));
	FILE* pre_file = fopen(pre.c_str(), "wb");
	fwrite(pre_arr[th_ind], sizeof(uint64_t), total_tree + 1, pre_file);
	fclose(pre_file);

	delete[] pre_arr[th_ind];
	delete[] kmer_arr[th_ind];
	delete[] temp_mem[th_ind];

	CountBigKmer(th_ind);
}


// main function

int main(int argc, char** argv) {
	if (argc % 2 == 0) {
		HowToUse();
		return EXIT_FAILURE;
	}

	int total_thread = thread::hardware_concurrency(); // total thread supported on the running machine
	for (int i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-k") == 0)
			kmer_len = std::atoi(argv[i + 1]); // command line kmer size. default = 28
		else if (strcmp(argv[i], "-i") == 0)
			in_file_name = argv[i + 1];
		else if (strcmp(argv[i], "-o") == 0)
			out_file_name = argv[i + 1];
		else if (strcmp(argv[i], "-t") == 0)
			total_thread = std::atoi(argv[i + 1]); // command line kmer size. default = 28
	}

	if (in_file_name.empty()) {
		cout << "Please specify input file";
		HowToUse();
		return EXIT_FAILURE;
	}
	

	kTotalMaskBitArr = new uint64_t[kTotalLayer]; // how many bits are left not for kmer
	k_count_mask = new uint64_t[kTotalLayer]; // mask for counter bits right aligned
	k_kmer_mask = new uint64_t[kTotalLayer]; // mask for kmer+x bits left aligned
	for (int i = 0; i < kTotalLayer; ++i) {
		kTotalMaskBitArr[i] = (32 - kmer_len + tree_pow - i) << 1;
		k_count_mask[i] = (1ULL << kTotalMaskBitArr[i]) - 1;
		k_kmer_mask[i] = ~k_count_mask[i];
	}

	mask = total_tree - 1; // tree index mask
	k_mask = (~0ULL) >> ((32 - kmer_len) << 1); // kmer+0 mask right aligned

#ifdef DEBUG
	std::cout << kTotalMaskBitArr[0] << " " << kTotalMaskBitArr[1] << " " << kTotalMaskBitArr[2] << " " << kTotalMaskBitArr[3] << " " <<  k_count_mask[0] << " " <<   k_kmer_mask[0] << " " <<      mask << " " <<  k_mask << std::endl;
#endif

	// time count starts
	TimePoint tp1 = Clock::now();
	std::clock_t cput1 = std::clock();

	// extract input file names
	std::vector<std::string> in_file_arr;
	if (in_file_name[0] == '@') {
		std::ifstream in_file(in_file_name.substr(1, in_file_name.length() - 1));
		std::string file_str;
		while (in_file >> file_str) {
			in_file_arr.push_back(file_str);
			//std::cout << file_str << "\n";
		}
		in_file.close();
	}
	else
		in_file_arr.push_back(in_file_name);

	// ===

	cout << total_thread << " threads" << endl;

	temp_mem = new uint64_t*[total_thread]; // used as temporarry memory for sorting
	root = new TrieNode**[total_thread]; // roots of tries
	kmer_arr = new uint64_t*[total_thread]; // contains kmers
	all_count = new uint64_t[total_thread]; // kmer count for each thread
	all_count_arr = new uint64_t*[total_thread];

	// threads to generate (k+x)-mers, traverse tries, and insert them all in k-mer tries 
	int total_file = in_file_arr.size();
	vector<thread> th_ins_arr;
	int file_per_th = total_file / total_thread;
	int used_file = file_per_th * total_thread;

	vector<vector<int>> file_ind_arr (total_thread);
	for (int i = 0; i < used_file; ++i) { 
		file_ind_arr[i % total_thread].push_back(i);
	}

	vector<vector<pair<long, long>>> f_pos_arr(total_file);
	for (int i = used_file; i < total_file; ++i) {
		f_pos_arr[i].resize(total_thread);
		SplitFile(in_file_arr[i].c_str(), total_thread, f_pos_arr[i]);
	}
	
	for (int i = 0; i < total_thread; ++i) { 
		th_ins_arr.push_back(thread(ComputeKmer, i, total_thread, ref(file_ind_arr[i]), ref(in_file_arr), ref(f_pos_arr))); 
	}

	for (int i = 0; i < total_thread; ++i) 
		th_ins_arr[i].join(); 
	
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	std::cout << "mem used: " << usage.ru_maxrss << std::endl;
	

	// k-mer tree traverse
	TimePoint tp3 = Clock::now();
	std::clock_t cput3 = std::clock();
	std::cout << "insertion time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp3-tp1).count()) / 1000 << "\tcpu " << (cput3 - cput1) / (double)CLOCKS_PER_SEC << "\n";

	// distibution computation of trees 
	struct TempCount {
		uint64_t c;
		int i;
	};
	uint64_t tc_total = 0;
	vector<TempCount> tc_arr (total_tree);
	for (int i = 0; i < total_tree; ++i) {
		tc_arr[i] = {0, i};
		for (int j = 0; j < total_thread; ++j) {
			tc_arr[i].c += all_count_arr[j][i];
		}
		tc_total += tc_arr[i].c;
	}

	sort(tc_arr.begin(), tc_arr.end(), [](const TempCount& a, const TempCount& b){
				return a.c < b.c;
			});

	uint64_t tc_avg = ceil(1.0 * tc_total / total_thread);
	//cout << "total " << tc_total << "\tavg: " << tc_avg << endl;
	uint64_t tc_sum = 0;
	int tc_i = 1;
	vector<vector<int>> order(total_thread);
	for (int i = 0; i < total_tree; ++i) {
	//	cout << i << ":" << tc_sum + tc_arr[i].c << ":" << tc_i * tc_avg << endl;
		if (tc_sum + tc_arr[i].c > tc_i * tc_avg) 
			tc_i++;

		tc_sum += tc_arr[i].c;
		order[tc_i - 1].push_back(tc_arr[i].i);
	}

	uniq_count = new uint64_t[total_thread];
	pre_arr = new uint64_t*[total_thread];
	bkmer_arr.resize(total_thread);

	// threads to traverse k-mer tries 
	vector<thread> th_tra_arr;
	for (int i = 0; i < total_thread; ++i) 
		th_tra_arr.push_back(thread(TraverseKmer, i, total_thread, ref(order[i]))); 

	for (int i = 0; i < total_thread; ++i) 
		th_tra_arr[i].join(); 

	delete[] root;
	delete[] kmer_arr;
	delete[] pre_arr;
	delete[] temp_mem;

	TimePoint tp4 = Clock::now();
	std::clock_t cput4 = std::clock();



#ifdef DEBUG
	std::cout << bucket_fill << " out " << bucket_size << "\t ratio: " << (1.0 * bucket_fill / bucket_size) << "\n";
	std::cout << tc1 << " node: " << num_node << "\t" << num_bucket << std::endl;
	std::cout << "sort: " << time_radm << "\tin_sort: " << tmm1 << "\ttr_sort: " << tmm2 << "\n";
#endif
#ifdef DEBUG
	std::cout << " over " << over_c << "\t" << " :: " << tkc5 << " " <<  tkc4 << " " <<  tkc3 << ":" << tkc2 << ":" << tkc1 << " uniq: " << uniq_count << "\tkmer: " << kmer_count << "\n";
#endif
	
	// aggregate kmer count
	uint64_t total_uniq_count = 0, total_kmer_count = 0;
	for (int i = 0; i < total_thread; ++i) {
		total_uniq_count += uniq_count[i];
		total_kmer_count += all_count[i];
		//cout << i << " all: " << all_count[i]  << "\t uniq: " << uniq_count[i] << endl; 
	}
	delete[] all_count;
	delete[] uniq_count;

	// statistics
	std::cout << "unique kmer: " << total_uniq_count << "\ttotal kmer: " << total_kmer_count << "\n";
	std::cout << "insertion time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp3-tp1).count()) / 1000 << "\tcpu " << (cput3 - cput1) / (double)CLOCKS_PER_SEC << "\n";
//	std::cout << "2nd phase time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp3-tp2).count()) / 1000 << "\tcpu " << (cput3 - cput2) / (double)CLOCKS_PER_SEC << "\n";
	std::cout << "3rd phase time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp4-tp3).count()) / 1000 << "\tcpu " << (cput4 - cput3) / (double)CLOCKS_PER_SEC << "\n";
	std::cout << "total time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp4-tp1).count()) / 1000 << "\tcpu " << (cput4 - cput1) / (double)CLOCKS_PER_SEC << "\n";

#ifdef DEBUG
	for (int i = 0; i < 30; ++i)
		if (depth_c[i])
			std::cout << i << ":\t" << depth_c[i] << std::endl;
#endif
	std::time_t curr_time = std::time(nullptr);
	std::cout << "date_time: " << std::asctime(std::localtime(&curr_time)) << std::endl;
}
