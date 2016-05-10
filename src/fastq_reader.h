
#include <iostream>
#include <cstdio>
#include <limits>
#include <vector>

using namespace std;

const int kBlockSize = 1 << 14;
const int kLineSize = 1 << 14;

class FastqReader {
public:
	void Init(const char* in_file, long s, long e);
	void Destroy();
	bool ReadSeq(char* &s);
	int len;
private:
	FILE* fp;
	char *seq, *sseq;
	long start, end;
	bool done;
};

void FastqReader::Init(const char* in_file, long s, long e) {
	fp = fopen(in_file, "r");
	fseek(fp, s, SEEK_SET);
	start = s;
	end = e;
	done = false;
	seq = (char*) malloc((kBlockSize + 1) * sizeof(char));
	seq[0] = '\0';
	sseq = seq;
}

void FastqReader::Destroy() {
	if (sseq)
		free(sseq);
	fclose(fp);
}


bool FastqReader::ReadSeq(char* &s) {
	int p = 0;
	for (int i = 0; i < 4; ++i) {
		while (*seq && *seq != '\n') {
			if (i == 1) 
				s[p++] = *seq;	
			++seq;
		}
		if (!*seq) {
			if (done) 
				return false;
			seq = sseq;
			long diff = end - ftell(fp);
			int sz = kBlockSize < diff? kBlockSize : diff;
			if (sz == diff)
				done = true;
			fread(seq, sizeof(char), sz, fp);	
			seq[sz] = '\0';
			while (*seq && *seq != '\n') {
				if (i == 1) 
					s[p++] = *seq;	
				++seq;
			}
			if (*seq)
				++seq;
		}
		else 
			++seq;
		if (i == 1) {
			s[p] = '\0';
			len = p;
		}

	}
	return true;
}




void SplitFile(const char* in_file, int total_th, vector<pair<long, long>>& f_pos_arr) {
	FILE* fp = fopen(in_file, "r");
	fseek(fp, 0, SEEK_END);
	long len = ftell(fp);
	
	
	f_pos_arr[0].first = 0;
	for (int i = 1; i < total_th; ++i) {
		long start = i * len / total_th;
		fseek(fp, start, SEEK_SET);
		
		int pos = 0;
		char c_arr[8] = {'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a'};

		char str[kLineSize];
		fgets(str, kLineSize, fp);
		for (int j = 0; j < 8; ++j) {
//		cout << i << ", " << j << ": " << str << endl;
			if (feof(fp)) {
				start = len;
				break;
			}
			int c = fgetc(fp);
			if (c == '@' && pos > 1 && c_arr[pos - 2] == '+') {
				ungetc(c, fp);
				start = ftell(fp);
				break;
			}
			c_arr[pos++] = c;
			fgets(str, kLineSize, fp);
		}
		fgets(str, kLineSize, fp);
		f_pos_arr[i].first = start;
	}
		
	for (int i = 0; i < total_th  - 1; ++i)
		f_pos_arr[i].second = f_pos_arr[i + 1].first;
	f_pos_arr[total_th - 1].second = len;
}

