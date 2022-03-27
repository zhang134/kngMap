#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include "kseq.h"
#include "kvec.h"
#include "edlib.h"


#ifdef __cplusplus
extern "C" {
#endif
using namespace std;
typedef struct {
	int threads; // number of threads
	//int num_mini; // how many minimizer for each read
	float cutoff; // clustering threshold
	short int k;  // kmer length
	short int knum; // kmer num per sequence
	short int w;  // find a minimizer for every $w consecutive k-mers
	unsigned int seq_num;  // total sequence num
	char *input;
	char *reads_file;
	char *output;
	char *pos_file;
	
} option;

typedef struct {
        uint64_t value;  // hash value
	unsigned int id;     // sequence id
        unsigned int pos;
        short int strand;       // strand:
} kmer_t;

typedef struct {
        //uint64_t value;  // hash value
        unsigned int pos_gen;
        unsigned int pos_seq;
	unsigned int id_gen;
	unsigned int id_seq;
        short int strand_gen;
	short int strand_seq;
	uint64_t value;  // hash value
} mm_match_t;

typedef struct {
	int matched_num;
	int id;
} matched_num_t;

typedef struct {
        unsigned int id; // ID number in the sequence file
        unsigned int len; // length of the sequence
        char *header;
	char *seq_is;
	char *seq_reverse;
	bool reverse_calculated;
	//float cutoff; // clustering threshold
	//short int assigned; // assigned: 1, not assigned 0
	//unsigned int is_seed; // is it seed, 1 yes, 0 no  index in a seeds set if it is a seed
	//unsigned int cluster; // cluster id
	unsigned int kmer_num;     // kmer number = len - k + 1
	char *cigar;
	unsigned int aligned_pos; // aligned position in the target genome
	unsigned int target_genome_id;
	//unsigned int seq_pos_start;
	//unsigned int seq_pos_end;
	//unsigned int gen_pos_start;
	//unsigned int gen_pos_end;
	//bool strand;  // plus 1, minus 0
	kmer_t *kmerlist;
} seqs;

// seed kmer, plus the index from 0, 1, ... , to N (seed num)
typedef struct {
        uint64_t value;  // hash value
        unsigned int id;     // sequence id
        unsigned int pos;
        short int strand;       // strand
	//int index;
} kmer_seed_t;

typedef struct {
	unsigned int match_num; // matched kmer num with each seed
	unsigned int id;	// seed id
	short int strand;
} match_t;

int cmp(const void *a,const void *b);
int cmp_matched_num(const void *a,const void *b);
int cmp_id(const void *a,const void *b);
int cmp_len_descending(const void *a,const void *b);
int set_options(int argc, char *argv[], option *opt);
int print_usage (char *arg);
void seq_sort_descending(seqs *seq, long int seqs_num);
//void mm_extract(const char *str, int len, int w, int k, uint32_t rid, kmer_t * arry);
void mm_extract(const char *str, int len, int k, int knum, uint32_t rid, kmer_t * arry);
int kmer_extract_for_one_seq(seqs *seq, int k, int w);
//int seeds_generate(seqs *arry, option opt, unsigned int start_id, kmer_seed_t *seeds_kmer, int seed_num, int *cluster_index);
kmer_seed_t * seeds_generate(seqs *arry, option opt, unsigned int start_id, int seed_num, int *cluster_index);
seqs * seeds_seq_generate(seqs *arry, option opt, unsigned int *start_id, unsigned int end_id, int seed_num, int *actual_seed_num,int *cluster_index);
seqs * seeds_rhat_generate_kmer(seqs *arry, option opt, unsigned int *start_id, unsigned int end_id, int seed_num, int *actual_seed_num, int *cluster_index, unsigned int * seeds_kmer_num);
seqs * seeds_seq_generate_kmer(seqs *arry, option opt, unsigned int *start_id, unsigned int end_id, int seed_num, int *actual_seed_num,int *cluster_index);
kmer_t * rhat_generate(seqs *seeds_seq_set, int seeds_num, unsigned int *point_list, unsigned int seeds_kmer_num, unsigned int kmer_code_num);
int strand_determine(kmer_t *list1, unsigned int n1, kmer_t *list2, unsigned int n2, int *matched_num);
int cluster_one(seqs *arry, option opt, seqs *seeds_seq, unsigned int seed_num, uint64_t id);
int cluster_one_with_rhat(seqs *arry, option opt, seqs *seeds_seq, unsigned int seed_num, uint64_t id, unsigned int *point_list, kmer_t *seeds_kmer_list);
void seq_reverse(seqs *arry);
int cluster_multi_threads(seqs *arry, option opt, unsigned int seqs_total);
int sequence_map_one(option opt, kmer_t *seeds_kmer_list, unsigned int *point_list, seqs *reads, seqs *gen_list);
int sequence_map_one_with_locate(option opt, seqs *reads, seqs *gen_list, short int *array_gen_id, bool *array_plus, unsigned int *array_seq_start, unsigned int *array_seq_end, unsigned int *array_gen_start, unsigned int *array_gen_end);
int sequence_map_write(option opt, seqs *reads);
void printAlignment(const char* query, const char* target, const unsigned char* alignment, const int alignmentLength, const int position, const EdlibAlignMode modeCode);
void bomb_error(const char *message);
