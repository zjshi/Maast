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
#include <thread>

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

void make_comp_map(char* comp_map) {
    comp_map['A'] = 'T';
    comp_map['C'] = 'G';
    comp_map['G'] = 'C';
    comp_map['T'] = 'A';
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
void load_profile(const char* k_path, vector<char>& buffer, vector<int_type>& kv1, vector<int_type>& kv2, vector<int_type>& kv3, vector<int_type>& kv4, vector<tuple<int_type, int_type>>& k_info, int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

    char* window = buffer.data();

    uintmax_t n_lines = 0;

    int fd;
    fd = open(k_path, O_RDONLY);

    int k_cur = 0;
    int snp_cur = 0;
    int pos_cur = 0;

    char kbuf1[k];
    char kbuf2[k];
    char kbuf3[k];
    char kbuf4[k];

	char snp_pos[16];
	char kmer_pos[4];

    //auto fh = fstream(out_path, ios::out | ios::binary);

	int cur_field = 0;
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
            char c = window[i];

            if (c == '\n') {
                ++n_lines;

                if (has_wildcard) {
                    has_wildcard = false;
                    continue;    
                }

                auto code1 = seq_encode<int_type>(kbuf1, k, code_dict, b_mask);
                auto code2 = seq_encode<int_type>(kbuf2, k, code_dict, b_mask);
                auto code3 = seq_encode<int_type>(kbuf3, k, code_dict, b_mask);
                auto code4 = seq_encode<int_type>(kbuf4, k, code_dict, b_mask);

				kv1.push_back(code1);
				kv2.push_back(code2);
				kv3.push_back(code3);
				kv4.push_back(code4);

				snp_pos[snp_cur] = '\0';
				int_type id_int = stoull(snp_pos);
				
				kmer_pos[pos_cur] = '\0';
				int_type k_pos = stoull(kmer_pos);

               	k_info.push_back(tuple<int_type, int_type>(id_int, k_pos));

				k_cur = 0;
                pos_cur = 0;
				snp_cur = 0;

				cur_field = 0;
            } else if (c == '\t'){
				++cur_field;
				k_cur = 0;
			} else {
                if (c == 'N') {
                    has_wildcard = true;    
                }

				if (cur_field == 0) {
					snp_pos[snp_cur++] = c;	
				} else if (cur_field == 1) {
					kmer_pos[pos_cur++] = c;	
				} else if (cur_field == 2) {
					kbuf1[k_cur++] = c;
				} else if (cur_field == 3) {
					kbuf2[k_cur++] = c;
				} else if (cur_field == 4) {
					kbuf3[k_cur++] = c;
				} else if (cur_field == 5) {
					kbuf4[k_cur++] = c;
				} else {
					//do nothing;
				}
            }
        }
    }


	assert(kv1.size() == k_info.size() && kv1.size() == kv2.size());
	assert(kv1.size() == kv3.size() && kv1.size() == kv4.size());

    auto timeit = chrono_time();
    close(fd);
}

template <class int_type>
void fna_load_pool(const char* fna_path, vector<char>& buffer, unordered_map<int_type, int_type>& k_map, const int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

    char comp_map[1 << (sizeof(char) * 8)];
    make_comp_map(comp_map);

    char* window = buffer.data();

    uintmax_t n_lines = 0;

    int fd;
    fd = open(fna_path, O_RDONLY);

    int cur_pos = 0;

	vector<char> seq_buf(10*1000*1000);
	char* bases = seq_buf.data();

	char kmer_buff[k];
	char rckmer_buff[k];

	bool is_base = false;
    bool has_wildcard = false;

    while (true) {

        const ssize_t bytes_read = read(fd, window, step_size);

        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
            cerr << "unknown fetal error when reading " << fna_path << endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0;  i < bytes_read;  ++i) {
			char c = window[i];
            if (c == '\n') {
				if (!is_base) {
					is_base = true;
				}
				continue;
            } else if (c == '>') {
				for (int j = 0; j < cur_pos-k+1; ++j) {
					for (int l = 0; l < k; ++l) {
						kmer_buff[l] = bases[j+l];

						if (kmer_buff[l] == 'N') {
							has_wildcard = true;
							break;						
						}
					}


					if (has_wildcard) {
						has_wildcard = false;
						continue;
					}

					auto kmer_int = seq_encode<int_type>(kmer_buff, k, code_dict, b_mask);
					
					if (k_map.find(kmer_int) == k_map.end()) {
						k_map.insert({kmer_int, 1});
					} else {
						++k_map[kmer_int];
					}

					for (int l = k-1; l >= 0; --l) {
						rckmer_buff[k-1-l] = comp_map[kmer_buff[l]];
					}

					/* not really necessary when rc kmers present
					auto rckmer_int = seq_encode<int_type>(rckmer_buff, k, code_dict, b_mask);

					if (k_map.find(rckmer_int) == k_map.end()) {
						k_map.insert({rckmer_int, 1});
					} else {
						++k_map[rckmer_int];
					}
					*/
				}

				cur_pos = 0;
				is_base = false;

				++n_lines;
			} else {
                if (is_base) {
                    bases[cur_pos++] = toupper(c);
                }
            }
        }

    }

	if (cur_pos >= k) {
		for (int j = 0; j < cur_pos-k+1; ++j) {
			for (int l = 0; l < k; ++l) {
				kmer_buff[l] = bases[j+l];

				if (kmer_buff[l] == 'N') {
					has_wildcard = true;
					break;						
				}
			}

			if (has_wildcard) {
				has_wildcard = false;
				continue;
			}

			auto kmer_int = seq_encode<int_type>(kmer_buff, k, code_dict, b_mask);
			
			if (k_map.find(kmer_int) == k_map.end()) {
				k_map.insert({kmer_int, 1});
			} else {
				++k_map[kmer_int];
			}

			/* not really necessary when rc kmers present
			for (int l = k-1; l >= 0; --l) {
				rckmer_buff[k-1-l] = comp_map[kmer_buff[l]];
			}

			auto rckmer_int = seq_encode<int_type>(rckmer_buff, k, code_dict, b_mask);
			if (k_map.find(rckmer_int) == k_map.end()) {
				k_map.insert({rckmer_int, 1});
			} else {
				++k_map[rckmer_int];
			}
			*/
		}

		cur_pos = 0;
		++n_lines;
	}

	buffer.clear();
	cerr << fna_path << endl;
	cerr << "number of sequences " << n_lines/2 << endl;
	cerr << "number of unique kmers: "<< k_map.size() << endl << endl;
    auto timeit = chrono_time();
    close(fd);
}


template <class int_type>
void bit_load_pool(const char* k_path, vector<char>& buffer, unordered_map<int_type, int_type>& k_map, const int_type* code_dict, const int_type b_mask) {
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

                auto kmer_int = seq_encode<int_type>(seq_buf, k, code_dict, b_mask);

                snp_id[snp_pos] = '\0';
				int_type kcount = stoull(snp_id);
					
				assert(k_map.find(kmer_int) == k_map.end());
				k_map.insert({kmer_int, kcount});

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
    close(fd);
}


template <class int_type>
void bin_load_pool(const char* p_path, unordered_map<int_type, int_type>& k_map) {
    size_t filesize = get_fsize(p_path);
    //Open file
    int fd = open(p_path, O_RDONLY, 0);
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
        auto kcount = mmappedData[i+1];

		assert(k_map.find(kmer_int) == k_map.end());
		k_map.insert({kmer_int, kcount});
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
void set_kvecs(char* db_path, const int kv_n, vector<int_type>* kvecs, vector<tuple<int_type, int_type>>& kinfo) {	
	assert(kv_n == 4);

    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);

    vector<char> buffer(buffer_size);

	load_profile<int_type>(db_path, buffer, kvecs[0], kvecs[1], kvecs[2], kvecs[3], kinfo, code_dict, b_mask);

	cerr << "DB loading OK!" << endl;
}

template <class int_type>
void multi_dbval(int kv_n, vector<int_type>* kvecs, int n_path, vector<string>& kpaths, vector<tuple<int, int, int, int, int>>& prof_vec) {	
	assert(kv_n == 4);

    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);

    vector<char> buffer(buffer_size);

	int64_t prof_size = kvecs[0].size();

	vector<int> lc_vecs[4];

	for (int i = 0; i < 4; ++i) {
		lc_vecs[i].reserve(prof_size);
	}

	prof_vec.resize(prof_size, make_tuple(0,0,0,0,0));

	unordered_map<int_type, int_type> kpool;

    for (int i = 0; i < n_path; ++i) {
		char* kp_type = get_ftype(kpaths[i].c_str());
		

		if (strcmp(kp_type, ".bin") == 0) {
			bin_load_pool<int_type>(kpaths[i].c_str(), kpool);
		} else if (strcmp(kp_type, ".fna") == 0) {
			fna_load_pool<int_type>(kpaths[i].c_str(), buffer, kpool, code_dict, b_mask);	
		} else {
			bit_load_pool<int_type>(kpaths[i].c_str(), buffer, kpool, code_dict, b_mask);	
		}

		// splitted loops 
		for (int j = 0; j < 4; ++j) {
			for(auto it = kvecs[j].begin(); it != kvecs[j].end(); ++it) {
				if(kpool.find(*it) == kpool.end()) {
					lc_vecs[j].push_back(0);
				} else {
					lc_vecs[j].push_back(kpool[*it]);
				}
			}
		}			
		
		assert(lc_vecs[0].size() == lc_vecs[1].size());
		assert(lc_vecs[0].size() == lc_vecs[2].size());
		assert(lc_vecs[0].size() == lc_vecs[3].size());

		const int64_t lc_size = lc_vecs[0].size();

		for (int64_t j = 0; j < lc_size; ++j) {
			auto lc_sum = lc_vecs[0][j] + lc_vecs[1][j] + lc_vecs[2][j] + lc_vecs[3][j];

			auto ref_sum = lc_vecs[0][j] + lc_vecs[2][j];
			auto alt_sum = lc_vecs[1][j] + lc_vecs[3][j];
			
			/*
			if (strcmp(kp_type, ".fna") == 0) {
				lc_sum = lc_sum / 2;
			}
			*/

			if (lc_sum == 0) {
				++get<0>(prof_vec[j]);			
			} else if (lc_sum == 1) {
				++get<1>(prof_vec[j]);			
			} else {
				++get<2>(prof_vec[j]);			
			}

			if (ref_sum > 0) {
				++get<3>(prof_vec[j]);			
			}

			if (alt_sum > 0) {
				++get<4>(prof_vec[j]);			
			}
		}

		for (int j = 0; j < 4; ++j) {
			lc_vecs[j].clear();
		}

		kpool.clear();
    }

    auto timeit = chrono_time();
}

void display_usage(char *fname){
	cout << "usage: " << fname << " -d profile_path -n identifier [-t n_threads] [-o output_path] [-L path to list of input] inpath1 [ inpath2 ...]\n";
}

int main(int argc, char** argv){		
	extern char *optarg;
    extern int optind;

    bool dbflag = false;
    bool inflag = false;
	bool idflag = false;
	bool list_flag = false;


    char* fname = argv[0];
    char* db_path = (char *)"";
	char* list_path = (char *)"";
    char* oname = (char *)"/dev/stdout";
	char* spe_id = (char *)"";

    int n_threads = 1;

    int opt;
    while ((opt = getopt(argc, argv, "d:n:t:L:o:h")) != -1) {
        switch (opt) {
            case 'd':
                dbflag = true;
                db_path = optarg;
                break;
            case 'n':
                idflag = true;
                spe_id = optarg;
                break;
            case 't':
                n_threads = stoi(optarg);
                break;
			case 'L':
				list_flag = true;
				list_path = optarg;
				break;
            case 'o':
                oname = optarg;
                break;
            case 'h': case '?':
                display_usage(fname);
                exit(1);
        }
    }

    cerr << fname << '\t' << db_path << '\t' << n_threads << endl;

    if (!dbflag) {
        cerr << "missing argument: -d <sckmerdb_path: string>\n";
        display_usage(fname);
        exit(1);
    }

	if (!idflag) {
		cerr << "missing argument: -n <species identifier>\n";
        display_usage(fname);
        exit(1);
	}


	int in_pos = optind;

	if (list_flag) {
        cerr << "program reads a list of kmer pools for checking kmer uniqueness: " << list_path << endl;
    } else { 
        if (optind == argc) {
            cerr << "missing argument: input (>1)\n";
            display_usage(fname);
            exit(1);
        }
    }


	vector<uint64_t> kvecs[4]; 
	const int max_size = 100 * 1000 * 1000;

	for (int i = 0; i < 4; ++i) {
        kvecs[i].reserve(max_size);
    }

	vector<tuple<uint64_t, uint64_t>> kinfo;
	
	set_kvecs<uint64_t>(db_path, 4, kvecs, kinfo);

	vector<string> input_array[n_threads];

	auto label = 0;
	if (list_flag) {
        ifstream file(list_path);
        string line;
		int l_count = 0;
        while (getline(file, line)) {
			label = l_count % n_threads;
            string tmp_line = line;
            input_array[label].push_back(tmp_line);
			++l_count;
        }
    }

	if (optind < argc) {
		for(; optind < argc; optind++) {
			auto slabel = (optind - in_pos + label) % n_threads;
			input_array[slabel].push_back(string(argv[optind]));
		}
	}

	vector<tuple<int, int, int, int, int>> prof_vecs[n_threads];
	vector<thread> th_array;
	
	for (int i = 0; i < n_threads; ++i) {
		th_array.push_back(thread(multi_dbval<uint64_t>, 4, kvecs, input_array[i].size(), ref(input_array[i]), ref(prof_vecs[i])));
	}


	for (thread & ith : th_array) {
		ith.join();
	}
	th_array.clear();

	vector<tuple<int, int, int, int, int>> reduced_prof;
	reduced_prof.resize(kvecs[0].size(), make_tuple(0,0,0,0,0));

    uint64_t lsb = 1;
    uint64_t b_mask = (lsb << bpb) - lsb;

    uint64_t code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<uint64_t>(code_dict);

	char sbuf1[k+1];
	char sbuf2[k+1];
	char sbuf3[k+1];
	char sbuf4[k+1];

	ofstream myfile;
	myfile.open(oname);

	for (int j = 0; j < kinfo.size(); ++j) {
		for (int i = 0; i < n_threads; ++i) {
			get<0>(reduced_prof[j]) += get<0>(prof_vecs[i][j]);	
			get<1>(reduced_prof[j]) += get<1>(prof_vecs[i][j]);	
			get<2>(reduced_prof[j]) += get<2>(prof_vecs[i][j]);	
			get<3>(reduced_prof[j]) += get<3>(prof_vecs[i][j]);	
			get<4>(reduced_prof[j]) += get<4>(prof_vecs[i][j]);	
		}

		auto info_pair = kinfo[j];

		seq_decode<uint64_t>(sbuf1, k+1, kvecs[0][j], code_dict, b_mask);
		seq_decode<uint64_t>(sbuf2, k+1, kvecs[1][j], code_dict, b_mask);
		seq_decode<uint64_t>(sbuf3, k+1, kvecs[2][j], code_dict, b_mask);
		seq_decode<uint64_t>(sbuf4, k+1, kvecs[3][j], code_dict, b_mask);

		myfile << get<0>(info_pair) << '\t' << get<1>(info_pair) << '\t' << sbuf1 << '\t' << sbuf2 << '\t' << sbuf3 << '\t' << sbuf4 << '\t' << get<0>(reduced_prof[j]) << '\t' << get<1>(reduced_prof[j]) << '\t'<< get<2>(reduced_prof[j]) << '\t' << spe_id << '\t' << get<3>(reduced_prof[j]) << '\t'<< get<4>(reduced_prof[j]) << '\n';
	}

    return 0;
}
