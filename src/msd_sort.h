

#ifndef MSD_SORT8_H_
#define MSD_SORT8_H_

#include <vector>
#include <cmath>
#include <chrono>
#include <cstring>
#include <algorithm>
#include "common.h"

using namespace std;

double time_radm = 0.0;
chrono::time_point<chrono::high_resolution_clock> tprm1, tprm2;

void DoRadixSort(uint64_t* data_begin, uint64_t* data_end, uint64_t *temp_begin, const uint64_t depth, int type) {
	uint64_t *temp_end = temp_begin + (data_end - data_begin);

	bool odd = 0;
	uint64_t limit = (64 - 2 * depth);
	for (uint64_t i = kTotalMaskBitArr[type]; i < limit; i += 8) {
		uint64_t count_arr[0x100] = {};
		for (uint64_t *p = data_begin; p != data_end; ++p)
			++count_arr[(*p >> i) & 0xFF];
		uint64_t* bucket_arr[0x100];
		uint64_t* q = temp_begin;
		for (uint64_t j = 0; j < 0x100; q += count_arr[j++])
			bucket_arr[j] = q;
		for (uint64_t* p = data_begin; p != data_end; ++p)
			*bucket_arr[(*p >> i) & 0xFF]++ = *p;
		swap(data_begin, temp_begin);
		swap(data_end, temp_end);

		odd ^= 1;
	}

	if (odd) {
		swap(data_begin, temp_begin);
		swap(data_end, temp_end);
		move(temp_begin, temp_end, data_begin);
	}
}

void RadixSortLSD8(uint64_t *begin, uint64_t *end, uint64_t *begin1, uint64_t maxshift, int type) {
	bool odd = 0;
	uint64_t *end1 = begin1 + (end - begin);
	//cout << type << " " << kTotalMaskBitArr[type] << " " << maxshift << endl;
	for (uint64_t shift = kTotalMaskBitArr[type]; shift < maxshift; shift += 8) {
		size_t count[0x100] = {};
		for (uint64_t *p = begin; p != end; p++)
			count[(*p >> shift) & 0xFF]++;
		uint64_t *bucket[0x100], *q = begin1;
		for (int i = 0; i < 0x100; q += count[i++])
			bucket[i] = q;
		for (uint64_t *p = begin; p != end; p++)
			*bucket[(*p >> shift) & 0xFF]++ = *p;
		swap(begin, begin1);
		swap(end, end1);
		odd ^= 1;
	}

	if (odd) {
		swap(begin, begin1);
		swap(end, end1);
	}
	else
		move(begin, end, begin1);
}

void DoMSDSort6(uint64_t *begin, uint64_t *end, uint64_t *begin1, uint64_t depth, int type) {
	uint64_t shift = 58 - 2 * depth;
    uint64_t *end1 = begin1 + (end - begin);
    size_t count[0x40] = {};
    for (uint64_t *p = begin; p != end; p++)
        count[(*p >> shift) & 0x3F]++;
    uint64_t *bucket[0x40], *obucket[0x40], *q = begin1;
    for (int i = 0; i < 0x40; q += count[i++])
        obucket[i] = bucket[i] = q;
    for (uint64_t *p = begin; p != end; p++)
        *bucket[(*p >> shift) & 0x3F]++ = *p;
    for (int i = 0; i < 0x40; ++i)
    	RadixSortLSD8(obucket[i], bucket[i], begin + (obucket[i] - begin1), shift, type);
}

uint64_t MergeAndDivideMSD(uint64_t* bucket, uint64_t start, uint64_t end, uint64_t depth, int type, int th_ind) {
#ifdef TIME1
	tprm1 = chrono::high_resolution_clock::now();
#endif
	//cout << type << " " << kTotalMaskBitArr[type] << " d: " << depth << endl;
	DoMSDSort6(&bucket[start], &bucket[end], temp_mem[th_ind], depth, type);
	//DoRadixSort(&bucket[start], &bucket[end], temp_mem, depth, type);
	//sort(&bucket[start], &bucket[end]);
#ifdef TIME1
	tprm2 = chrono::high_resolution_clock::now();
	time_radm += chrono::duration_cast<chrono::microseconds>(tprm2 - tprm1).count();
#endif

	uint64_t* first = &bucket[start], *last = &bucket[end];
	uint64_t* result = first;
	uint64_t count = *result & k_count_mask[type];
	uint64_t kmer_part = *result & k_kmer_mask[type];
	uint64_t count_mask = k_count_mask[type], kmer_mask = k_kmer_mask[type];

	while (++first != last) {
		uint64_t first_c = (*first & count_mask);
		uint64_t first_k = (*first & kmer_mask);
		if ((kmer_part != first_k) || (count + first_c >= count_mask)) {
			*result = kmer_part | count;
			*(++result) = *first;
			count = first_c;
			kmer_part = first_k;
		}
		else
			count += first_c;
	}
	*result = kmer_part | count;
	return distance(&bucket[0], ++result);
}

uint64_t Merge(uint64_t* bucket, uint64_t start2, uint64_t end2, int type) {
	uint64_t* first1 = &bucket[0], *last1 = &bucket[start2];
	uint64_t* first2 = &bucket[start2], *last2 = &bucket[end2];
	uint64_t* out = new uint64_t[end2];
	uint64_t* result = out;
	uint64_t dist = 0;
	while (true) {
		if (first1 == last1) {
			dist = distance(&out[0], result) + distance(first2, last2);
			copy(first2, last2, result);
			break;
		}
		if (first2 == last2) {
			dist = distance(&out[0], result) + distance(first1, last1);
			copy(first1, last1, result);
			break;
		}
		uint64_t k1 = *first1 & k_kmer_mask[type];
		uint64_t k2 = *first2 & k_kmer_mask[type];
		if (k1 < k2)
			*result++ = *first1++;
		else if (k1 > k2)
			*result++ = *first2++;
		else {
			uint64_t c1 = *first1 & k_count_mask[type];
			uint64_t c2 = *first2 & k_count_mask[type];
			if (c1 + c2 < k_count_mask[type]) {
				*result++ = k1 | (c1 + c2);
				first1++;
				first2++;
			}
			else {
				*result++ = *first1++;
				*result++ = *first2++;
			}
		}
	}
	move(&out[0], &out[dist], &bucket[0]);
	delete[] out;
	out = nullptr;
	return dist;
}



 
struct HeapNode {
	uint64_t kmer;
	int i;
};

class MinHeap{
public:
	MinHeap(HeapNode* arr, int heap_sz) : cont(arr), sz(heap_sz) {
		int i = (sz - 1) >> 1;
		while (i >= 0) {
			MinHeapify(i);
			--i;
		}
	}

	HeapNode GetMin() {
		return cont[0];
	}
	
	void ReplaceMin(HeapNode x) {
		cont[0] = x;
		MinHeapify(0);
	}

private:
	void Swap(HeapNode* x, HeapNode* y) {
		HeapNode t = *x;
		*x = *y;
		*y = t;
	}

	void MinHeapify(int i) {
		int l = (i << 1) + 1;
		int r = l + 1;
		int min_ind = i;
		if (l < sz && cont[l].kmer < cont[i].kmer)
			min_ind = l;
		if (r < sz && cont[r].kmer < cont[min_ind].kmer)
			min_ind = r;
		if (min_ind != i) {
			Swap(&cont[i], &cont[min_ind]);
			MinHeapify(min_ind);
		}
	}

	HeapNode* cont;
	int sz;
};


void MergeKArr(uint64_t* pos_arr, int th_ind, int total_thread) {
	uint64_t* temp_arr = new uint64_t[total_thread + 1];
	for (int i = 0; i <= total_thread; ++i)
		temp_arr[i] = pos_arr[i];

	uint64_t n = pos_arr[total_thread];
	HeapNode* h_arr = new HeapNode[total_thread];
	for (int i = 0; i < total_thread; ++i) {
		h_arr[i].kmer = kmer_arr[th_ind][pos_arr[i]++];
		h_arr[i].i = i;
	}

	MinHeap heap(h_arr, total_thread);
	
	uint64_t* out_arr = new uint64_t[n];
	uint64_t c = 0;
	for (uint64_t i = 0; i < n; ++i) {
		HeapNode r = heap.GetMin();
		out_arr[c++] = r.kmer;

		if (pos_arr[r.i] < temp_arr[r.i + 1]) 
			r.kmer = kmer_arr[th_ind][pos_arr[r.i]++];
		else
			r.kmer = 0xFFFFFFFFFFFFFFFF;
		heap.ReplaceMin(r);	
	}

	delete[] h_arr;
	h_arr = nullptr;
	move(&out_arr[0], &out_arr[n], &kmer_arr[th_ind][0]);
	delete[] out_arr;
	out_arr = nullptr;
	delete[] temp_arr;
	temp_arr = nullptr;
}

uint64_t CompactKmer(uint64_t* bucket, uint64_t len, int th_ind, uint64_t tree_ind) {
	uint64_t count_mask = k_count_mask[0], kmer_mask = k_kmer_mask[0];
	uint64_t* out = new uint64_t[len];
	uint64_t* result = out;

	uint64_t* first1 = &bucket[0], *last1 = &bucket[len];
	uint64_t* temp1 = first1;
	++first1;
	while (first1 != last1) {
		if ((*temp1 & kmer_mask) != (*first1 & kmer_mask)) {
			if ((*temp1 & count_mask) != count_mask)
				*result++ = *temp1;
			temp1 = first1;
			++first1;
		}
		else {
			if (((*temp1 & count_mask) + (*first1 & count_mask)) >= count_mask) { // check & remove this line
				BigKmer bkmer {(tree_ind << ((kmer_len - tree_pow) << 1)) | ((*temp1 & kmer_mask) >> ((32 - kmer_len + tree_pow) << 1)), (*temp1 & count_mask) + (*first1 & count_mask) - count_mask};
				bkmer_arr[th_ind].push_back(bkmer);
				*temp1 = (*temp1 & kmer_mask) | count_mask;
			}
			else
				*temp1 = (*temp1 & kmer_mask) | ((*temp1 & count_mask) + (*first1 & count_mask));
			++first1;
		}
	}
	if ((*temp1 & count_mask) != count_mask)
		*result++ = *temp1;

	uint64_t dist =  distance(&out[0], result);
	move(&out[0], &out[dist], &bucket[0]);
	delete[] out;
	out = nullptr;
	return dist;
}


#endif
