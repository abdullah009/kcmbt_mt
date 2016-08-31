#ifndef COMMON_H_
#define COMMON_H_

#include <cmath>

using namespace std;

#define MALLOC

struct BigKmer {
	uint64_t kmer;
	uint64_t count;
};

vector<vector<BigKmer>> bkmer_arr;

const uint64_t tree_pow = 6;
const uint64_t total_tree = pow(4, tree_pow);
const uint64_t kTotalLayer = 4;

const uint64_t kDepthInc = 1;

const uint64_t kInitCount = 1; // make it 1
uint64_t kmer_len = 28;
const uint64_t kMaxBucketSize = 8 * 256 * 1024;
const uint64_t kRootBucketSize = 2 * 32 * 1024;
const uint64_t kMinBucketSize = 1024;
const uint64_t kAlphabetSize = 4;

const uint64_t kPrevInit = 2;
const uint64_t kTotalLayerMinusOne = kTotalLayer - 1;


uint64_t* kTotalMaskBitArr;
uint64_t mask, k_mask;
uint64_t** kmer_arr;
uint64_t** pre_arr;
uint64_t over_c = 0, an = 0, bn = 0, cn = 0, dn = 0;

uint64_t* k_count_mask;
uint64_t* k_kmer_mask;

uint64_t** temp_mem __attribute((aligned(64)));

#endif
