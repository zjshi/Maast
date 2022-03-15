#if __linux__
#include <linux/version.h>
#if LINUX_VERSION_CODE > KERNEL_VERSION(2,6,22)
#define _MAP_POPULATE_AVAILABLE
#endif
#endif

#ifdef _MAP_POPULATE_AVAILABLE
#define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
#else
#define MMAP_FLAGS MAP_PRIVATE
#endif

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <cstring>

#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

using namespace std;


// this program scans its input (fastq text stream) for forward k mers,

// usage:
//    g++ -O3 --std=c++11 -o vfkmrz_bunion vfkmrz_bunion.cpp
//    ./vfkmrz_bunion -k1 </path/to/kmer_list1> -k2 </path/to/kmer_list2>
//
// standard fastq format only for input, otherwise failure is almost guaranteed. 

// global variable declaration starts here
constexpr auto k = 31;

// set operation mode
// valid values: 0, 1, 2
// 0 is set union operation; 1 is set intersection operation; 2 is set difference([set1-set2]);
constexpr auto s_mod = 0;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 256 * 1024 * 1024;
constexpr auto buffer_size = 256 * 1024 * 1024;

// output file path
constexpr auto out_path = "/dev/stdout";

// get time elapsed since when it all began in milliseconds.
long chrono_time() {
    using namespace chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

// number of bits per single nucleotide base
constexpr int bpb = 2;

size_t get_fsize(const char* filename) {
    struct stat st;
    stat(filename, &st);
    return st.st_size;
}


char* get_ftype(const char* filename) {
    int fn_len = strlen(filename);
    char *ftype = (char *)malloc(5);
    
    for(int i = 0; i < 4; ++i) {
        ftype[i] = filename[fn_len - 4 + i];
    }

    ftype[4] = '\0';
    
    return ftype;
} 


template <class int_type>
int_type bit_encode(const char c) {
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    }

    assert(false);
}


template <class int_type>
char bit_decode(const int_type bit_code) {
    switch (bit_code) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    }
    assert(false);
}

template <class int_type>
void make_code_dict(int_type* code_dict) {
    code_dict['A'] = bit_encode<int_type>('A');
    code_dict['C'] = bit_encode<int_type>('C');
    code_dict['G'] = bit_encode<int_type>('G');
    code_dict['T'] = bit_encode<int_type>('T');
}

template <class int_type>
int_type seq_encode(const char* buf, int len, const int_type* code_dict, const int_type b_mask) {
    int_type seq_code = 0;
    for (int i=0;  i < len;  ++i) {
        const int_type b_code = code_dict[buf[i]];
        seq_code |= ((b_code & b_mask) << (bpb * (len - i - 1)));
    }
    return seq_code;
}

template <class int_type>
void seq_decode(char* buf, const int len, const int_type seq_code, int_type* code_dict, const int_type b_mask) {
    for (int i=0;  i < len-1;  ++i) {
        const int_type b_code = (seq_code >> (bpb * (len - i - 2))) & b_mask;
        buf[i] = bit_decode<int_type>(b_code);
    }

    buf[len-1] = '\0';
}


template <class int_type>
void bit_load(const char* k_path, vector<char>& buffer, vector<tuple<int_type, int_type>>& k_vec, const int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

    char* window = buffer.data();

    uintmax_t n_lines = 0;

    int fd;
    fd = open(k_path, O_RDONLY);

    int cur_pos = 0;
    int snp_pos = 0;

    char seq_buf[k];
	char snp_id[16];

    //auto fh = fstream(out_path, ios::out | ios::binary);

	bool id_switch = false;
    bool has_wildcard = false;

    while (true) {

        const ssize_t bytes_read = read(fd, window, step_size);

        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
            cerr << "unknown fetal error when reading " << k_path << endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0;  i < bytes_read;  ++i) {
            char c = toupper(window[i]);
            if (c == '\n') {
                ++n_lines;

                if (has_wildcard) {
                    has_wildcard = false;
                    continue;    
                }

                auto code = seq_encode<int_type>(seq_buf, k, code_dict, b_mask);

				snp_id[snp_pos] = '\0';
				int_type id_int = stoull(snp_id);
				//int_type id_int = 1;

                k_vec.push_back(tuple<int_type, int_type>(code, id_int));

                cur_pos = 0;
				snp_pos = 0;

				id_switch = false;
            } else if (c == '\t'){
				id_switch = true;
			} else {
                if (c == 'N') {
                    has_wildcard = true;    
                }

				if (id_switch) {
					snp_id[snp_pos++] = c;	
				} else {
					seq_buf[cur_pos++] = c;
				}
            }
        }

        //fh.write(&kmers[0], kmers.size());

        // cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    }

    auto timeit = chrono_time();
}


template <class int_type>
void bit_load(vector<char>& buffer, vector<tuple<int_type, int_type>>& k_vec, const int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

    char* window = buffer.data();

    uintmax_t n_lines = 0;

    int cur_pos = 0;
    int snp_pos = 0;

    char seq_buf[k];
	char snp_id[16];

    //auto fh = fstream(out_path, ios::out | ios::binary);

	bool id_switch = false;
    bool has_wildcard = false;

    while (true) {

        const ssize_t bytes_read = read(STDIN_FILENO, window, step_size);

        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
            cerr << "unknown fetal error when reading from stdin" << endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0;  i < bytes_read;  ++i) {
            char c = toupper(window[i]);
            if (c == '\n') {
                ++n_lines;

                if (has_wildcard) {
                    has_wildcard = false;
                    continue;    
                }

                auto code = seq_encode<int_type>(seq_buf, k, code_dict, b_mask);

				snp_id[snp_pos] = '\0';
				int_type id_int = stoull(snp_id);
				//int_type id_int = 1;

                k_vec.push_back(tuple<int_type, int_type>(code, id_int));

                cur_pos = 0;
				snp_pos = 0;
				id_switch = false;
            } else if (c == '\t'){
				id_switch = true;
			} else {
                if (c == 'N') {
                    has_wildcard = true;    
                }

				if (id_switch) {
					snp_id[snp_pos++] = c;	
				} else {
					seq_buf[cur_pos++] = c;
				}
            }
        }

        //fh.write(&kmers[0], kmers.size());

        // cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    }

    auto timeit = chrono_time();
}


template <class int_type>
void binary_load(const char* k_path, vector<tuple<int_type, int_type>>& k_vec) {
    size_t filesize = get_fsize(k_path);
    //Open file
    int fd = open(k_path, O_RDONLY, 0);
    assert(fd != -1);
    //Execute mmap
    //uint64_t* mmappedData = (uint64_t *) mmap(NULL, filesize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
    int_type* mmappedData = (int_type *) mmap(NULL, filesize, PROT_READ, MMAP_FLAGS, fd, 0);
    assert(mmappedData != MAP_FAILED);
    //Write the mmapped data to stdout (= FD #1)

    // write(1, mmappedData, filesize);

    // char seq_buf[k+1];

    auto l_start = chrono_time();

    for (uint64_t i = 0; i < filesize/8; i=i+2) {
        // seq_decode<uint_fast64_t>(seq_buf, k, mmappedData[i], b_mask);

        auto kmer_int = mmappedData[i];
        auto snp = mmappedData[i+1];

		string snp_str = to_string(snp);
		
		if (snp_str[6] == '2'){
			snp_str[6] = '0';
		} else if (snp_str[6] == '3') {
			snp_str[6] = '1';
		}

        auto k_pair = make_tuple(kmer_int, stoull(snp_str));

        k_vec.push_back(k_pair);
    }

    //Cleanup
    int rc = munmap(mmappedData, filesize);
    assert(rc == 0);
    close(fd);
}


template <class int_type>
bool cmp_tuple(const tuple<int_type, int_type> &a, const tuple<int_type, int_type> &b){
	return get<0>(a) < get<0>(b);
}

template <class int_type>
void multi_btc64() {	
    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);

    vector<tuple<int_type, int_type>> kdb;
    vector<char> buffer(buffer_size);

    bit_load<int_type>(buffer, kdb, code_dict, b_mask);	

    auto timeit = chrono_time();
    sort(kdb.begin(), kdb.end(), cmp_tuple<int_type>);
    // typename vector<int_type>::iterator ip = unique(kdb.begin(), kdb.end());
    // kdb.resize(std::distance(kdb.begin(), ip));
    cerr << "Done!\n" << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the kmer list has " << kdb.size() << " kmers" << endl;

    // char seq_buf[k+1];
	
	vector<int_type> o_buff;

    ofstream fh(out_path, ofstream::out | ofstream::binary);

    for (auto it = kdb.begin(); it != kdb.end(); ++it) {
        // seq_decode(seq_buf, k+1, *it, code_dict, b_mask);    
        // fh << seq_buf << "\n";
        // fh << *it << "\n";

		// cerr << get<0>(*it) << '\t' << get<1>(*it) << '\n';
		o_buff.push_back(get<0>(*it));
		o_buff.push_back(get<1>(*it));
    }

    fh.write((char*)&o_buff[0], o_buff.size() * sizeof(int_type));

    fh.close();
}

template <class int_type>
void multi_btc64(int n_path, char** kpaths) {	
    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);

    vector<tuple<int_type, int_type>> kdb;
    vector<char> buffer(buffer_size);

    for (int i = 1; i < n_path; ++i) {
        cerr << kpaths[i] << endl;

		char* kp_type = get_ftype(kpaths[i]);
		
		if (strcmp(kp_type, ".tsv") == 0) {
			bit_load<int_type>(kpaths[i], buffer, kdb, code_dict, b_mask);	
		} else if (strcmp(kp_type, ".bin") == 0) {
			binary_load<int_type>(kpaths[i], kdb);
		} else {
			assert(false);
		}
    }

    auto timeit = chrono_time();

	sort(kdb.begin(), kdb.end(), cmp_tuple<int_type>);
    // typename vector<int_type>::iterator ip = unique(kdb.begin(), kdb.end());
    // kdb.resize(std::distance(kdb.begin(), ip));
    cerr << "Sorting done! " << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the kmer list has " << kdb.size() << " kmers" << endl;

	char seq_buf[k+1];
    // ofstream fh(out_path, ofstream::out | ofstream::binary);
	vector<int_type> o_buff;

	// move onto checkout when checkout_flag is true
	bool checkout_flag = true;
	vector<tuple<int_type, int_type>> auto_queue;

    cerr << "start to check conflicts" << endl;
    for (auto it = kdb.begin(); it+1 != kdb.end(); ++it) {
		// seq_decode(seq_buf, k+1, get<0>(*it), code_dict, b_mask);    
		// cerr << seq_buf << '\t' << get<1>(*it) << '\n';
		
		if (get<0>(*it) == get<0>(*(it+1))) {
			auto spe1 = stoi(to_string(get<1>(*it)).substr(0, 6));
			auto spe2 = stoi(to_string(get<1>(*(it+1))).substr(0, 6));

			if (spe1 != spe2) {
				checkout_flag = false;			
			}
			
			auto_queue.push_back(*it);	

			continue;
		}

		// check out when code(i) != code(i+1)
		if (!checkout_flag) {
			for(auto iq = auto_queue.begin(); iq != auto_queue.end(); ++iq){
				//cerr << get<0>(*iq) << " - " << get<1>(*iq) << '\n';
			}
			auto_queue.clear();
			checkout_flag = true;
		} else {
			if (auto_queue.size() > 0){
				for(auto iq = auto_queue.begin(); iq != auto_queue.end(); ++iq){
					o_buff.push_back(get<0>(*iq));
					o_buff.push_back(get<1>(*iq));
				}

				auto_queue.clear();
			}
			o_buff.push_back(get<0>(*it));
			o_buff.push_back(get<1>(*it));
		} 
    }

	auto end_ele = kdb.back();
	if (checkout_flag) {
		if (auto_queue.size() > 0){
			for(auto iq = auto_queue.begin(); iq != auto_queue.end(); ++iq){
				o_buff.push_back(get<0>(*iq));
				o_buff.push_back(get<1>(*iq));
			}

			auto_queue.clear();
		}
		o_buff.push_back(get<0>(end_ele));
		o_buff.push_back(get<1>(end_ele));
	}

    cerr << "the kmer list has " << o_buff.size()/2<< " kmers after purging conflicts" << endl;

	vector<tuple<int_type, int_type>>().swap(kdb);

	unordered_map<uint64_t, tuple<int, int>> snp_indx;	

    for (auto it = o_buff.begin(); it != o_buff.end(); it=it+2) {
		auto snp = stoull(to_string(*(it+1)).substr(0, 6) + to_string(*(it+1)).substr(7));
		auto snp_type = stoi(to_string(*(it+1)).substr(6, 1));
		
		if (snp_indx.find(snp) == snp_indx.end()) {
			auto type_pair = make_tuple(0, 0);
			snp_indx.insert({snp, type_pair});
		}

		assert(snp_type == 0 || snp_type == 1);

		if (snp_type == 0) {  
			get<0>(snp_indx[snp]) = 1;
		} else {
			get<1>(snp_indx[snp]) = 1;
		}
	}

	vector<int_type> v_buff;

	for (auto it = o_buff.begin(); it != o_buff.end(); it=it+2) {
		auto snp = stoull(to_string(*(it+1)).substr(0, 6) + to_string(*(it+1)).substr(7));

		assert(snp_indx.find(snp) != snp_indx.end());

		if (get<0>(snp_indx[snp]) + get<1>(snp_indx[snp]) == 2) {
			v_buff.push_back(*it);
			v_buff.push_back(*(it+1));			
		}
	}

    cerr << "the kmer list has " << v_buff.size()/2<< " kmers after purging conflicts" << endl;

    ofstream fh(out_path, ofstream::binary);

    fh.write((char*)&v_buff[0], v_buff.size() * sizeof(int_type));
    fh.close();
}

void display_usage(char *fname){
	cout << "usage: " << fname << " fpath [fpath ...]\n";
}

int main(int argc, char** argv){		
	if (argc == 2 && string(argv[1]) == "-h") {
		display_usage(argv[0]);
	} else if (argc >= 2) {
        multi_btc64<uint64_t>(argc, argv);		
	} else if (argc == 1) {
        multi_btc64<uint64_t>();
    } else {
        cerr << argv[0] << " reads from stdin or takes at least one arguments!" << endl;
		display_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    return 0;
}
