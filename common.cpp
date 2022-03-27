#include <zlib.h>  
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <set>
#include <map>
#include "omp.h"
#include "common.h"
#include "kseq.h"
#include "edlib.h"
//using namespace std;
unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


int cmp(const void *a,const void *b){
	return (*(kmer_t*)a).value <= (*(kmer_t*)b).value ? -1 : 1;
}
int cmp_mm_match_t(const void *a,const void *b){
	return (*(mm_match_t*)a).pos_gen < (*(mm_match_t*)b).pos_gen ? -1 : 1;
}
bool cmp_sort(kmer_t a, kmer_t b){
	return a.value < b.value;
}


int cmp_id(const void *a,const void *b) {
	return (*(kmer_t*)a).id <= (*(kmer_t*)b).id ? -1 : 1;
}
int cmp_matched_num(const void *a,const void *b){
	return (*(matched_num_t*)a).matched_num <= (*(matched_num_t*)b).matched_num ? 1 : -1;
}
int cmp_len_descending(const void *a,const void *b){
	return (*(seqs*)a).len <= (*(seqs*)b).len ? 1 : -1;
}
// sort sequences by descending order
void seq_sort_descending(seqs *seq_list, long int seqs_num){
	qsort(seq_list, seqs_num, sizeof(seq_list[0]), cmp_len_descending);
	//return 0;
}

char txt_option_g[] = "\tinput genome filename in fasta format, required\n";
char txt_option_r[] = "\tinput reads filename in fasta format, required\n";
char txt_option_o[] = "\toutput filename, required\n";
char txt_option_p[] = "\tposition filename, required\n";
char txt_option_c[] =
"\tsequence identity threshold, default 0.9\n \
\tthis is the default cd-hit's \"global sequence identity\" calculated as:\n \
\tnumber of identical amino acids in alignment\n \
\tdivided by the full length of the shorter sequence\n";
char txt_option_t[] = "\tnumber of threads, default 1; with 0, all CPUs will be used\n";
char txt_option_n[] = "\tnumber of reads, required\n";
//char txt_option_N[] = "\tnumber of genomes, optinal\n";
char txt_option_k[] = "\tlength of k-mer, default 20, required\n";

int print_usage (char *arg) {
	cout << "\n\nUsage: "<< arg << " [Options] \n\nOptions\n\n";
	cout << "   -g" << txt_option_g;
	cout << "   -r" << txt_option_r;
	cout << "   -p" << txt_option_p;
	cout << "   -o" << txt_option_o;
	//cout << "   -n" << txt_option_n;
	//cout << "   -k" << txt_option_k;
	cout << endl;
	//cout << "   -c" << txt_option_c;
	//cout << "   -k" << txt_option_k;
	cout << "   -t" << txt_option_t;
	//cout << "   -n" << txt_option_n;
	//cout << "   -N" << txt_option_N;
	cout << "\n";

}


int set_options(int argc, char *argv[], option *opt){
	int i, n;
	if (argc == 1){
		print_usage(argv[0]);
		return 0;
	}

	for (i=1; i+1<argc; i+=2){
		if (strcmp(argv[i], "-g") == 0) {
			opt->input = (char *) malloc(strlen(argv[i + 1]) + 1);
			strcpy(opt->input, argv[i + 1]);
			//opt->input[strlen(argv[i + 1])] = '\0';
			//opt.input = argv[i + 1];
			//seeds_seq_sets[i].seq_is = (char *) malloc(strlen(arry[l].seq_is));
                	//strcpy(seeds_seq_sets[i].seq_is, arry[l].seq_is);
		}
		else if (strcmp(argv[i], "-r") == 0) {
                        opt->reads_file = (char *) malloc(strlen(argv[i + 1]));
                        strcpy(opt->reads_file, argv[i + 1]);
                }
		else if (strcmp(argv[i], "-p") == 0) {
			opt->pos_file = (char *) malloc(strlen(argv[i + 1]));
                        strcpy(opt->pos_file, argv[i + 1]);
		}
		else if (strcmp(argv[i], "-o") == 0) {
			opt->output = (char *) malloc(strlen(argv[i + 1]));
			strcpy(opt->output, argv[i + 1]);
		}
		else if (strcmp(argv[i], "-t") == 0) opt->threads = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "-k") == 0) opt->k = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "-c") == 0) opt->output = argv[i + 1];
		else if (strcmp(argv[i], "-n") == 0) opt->seq_num = atoi(argv[i + 1]);
		else print_usage(argv[0]);
	}
	return 1;
}

// extract the miniminzers for each sequence
// str	sequence is
// len	sequence length
// w	find a minimizer for every $w consecutive k-mers
// k	kmer length
// knum	kmer num for each sequence
// rid	sequence id
// arry	kmer list
void mm_extract(const char *str, int len, int k, int knum, uint32_t rid, kmer_t * arry) {
	int i, c, z, l = 0, lianxu = 0, nn = 0, k_calculated = 0;
	int w = floor((len - k + 1)/knum);
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	uint64_t indexx, kmer_min = UINT64_MAX;
	kmer_t kmer_curr, kmer_min_t;
	for (i = 0; i < len; ++i) {
		++l;
		c = seq_nt4_table[(uint8_t)str[i]];
		if (c < 4) { // not an ambiguous base
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		}
		//if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand, plus 0, minus 1
		if (l >= k){
			//kmer_curr.value = kmer[z];
			//kmer_curr.strand = z;
			//kmer_curr.id = rid;
			if (kmer[z] < kmer_min){
				kmer_min_t.value = kmer[z];
				kmer_min_t.strand = z;
				kmer_min_t.id = rid;
				kmer_min_t.pos = i;
				kmer_min = kmer[z];

			}
			++lianxu;
		}
		else {
			continue;
		}
		
		if (lianxu == w) { // special case for the first window - because identical k-mers are not stored yet
			indexx = rid * knum + k_calculated;
			arry[indexx].value  = kmer_min_t.value;
			arry[indexx].strand =  kmer_min_t.strand;
			arry[indexx].id = kmer_min_t.id;
			arry[indexx].pos = kmer_min_t.pos;
			//printf("sequence id: %5d, kmer index: %8d\n", rid, rid * knum + k_calculated);
			kmer_min = UINT64_MAX;
			++k_calculated;
			lianxu = 0;
			if (k_calculated  == 50)
				break;
		}

	}
}


int kmer_extract_for_one_seq(seqs *seq, int k, int w){
	unsigned int kmer_num = seq->len - k + 1;
	seq->kmer_num = kmer_num;
	seq->kmerlist = (kmer_t*)malloc(kmer_num * sizeof(kmer_t));
	unsigned int i;
	int c, z, l = 0, lianxu = 0, nn = 0, k_calculated = 0;
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
        uint64_t indexx, kmer_min = UINT64_MAX;
        kmer_t kmer_curr, kmer_min_t;
	for (i = 0; i < seq->len; ++i){
                ++l;
                c = seq_nt4_table[(uint8_t)seq->seq_is[i]];
                if (c < 4) { // not an ambiguous base
                        kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
                        kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
                }
                //if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
                z = kmer[0] < kmer[1]? 0 : 1; // strand, plus 0, minus 1
                if (l >= k){
                        //kmer_curr.value = kmer[z];
                        //kmer_curr.strand = z;
                        //kmer_curr.id = rid;
                        if (kmer[z] < kmer_min){
                                kmer_min_t.value = kmer[z];
                                kmer_min_t.strand = z;
                                //kmer_min_t.id = rid;
                                kmer_min_t.pos = i;
                                kmer_min = kmer[z];

                        }
                        ++lianxu;
                }
                else {
                        continue;
                }
		if (lianxu == w) { // special case for the first window - because identical k-mers are not stored yet
                        indexx = k_calculated;
                        seq->kmerlist[indexx].value  = kmer_min_t.value;
                        seq->kmerlist[indexx].strand =  kmer_min_t.strand;
                        seq->kmerlist[indexx].id = seq->id;//kmer_min_t.id;
                        seq->kmerlist[indexx].pos = kmer_min_t.pos;
                        //printf("sequence id: %5d, kmer index: %8d\n", rid, rid * knum + k_calculated);
                        kmer_min = UINT64_MAX;
                        ++k_calculated;
                        lianxu = 0;
                }
	}
	return 0;
}

// determain plus direction or minus direction and matched kmer num for seq2 to seq1
// list1: the first kmer list
// n1:    the number of the kmer in list1
// matched_num: the identical kmer num
// return 0 plus direction, 1: minus direction for seq2
int strand_determine(kmer_t *list1, unsigned int n1, kmer_t *list2, unsigned int n2, int *matched_num){
	unsigned int n = n1 + n2, indexx = 0, i, curr_id, pre_id, plus_num = 0, minus_num = 0, same_num = 0;
	short int curr_strand, pre_strand;
	uint64_t pre, curr;
	kmer_t *combine = (kmer_t*)malloc(n * sizeof(kmer_t));
	for (i = 0; i < n1; ++i){
		combine[i].value = list1[i].value;
		combine[i].strand = list1[i].strand;
		combine[i].id = 0;
	}
	for (; i < n; ++i){
		combine[i].value = list2[i - n1].value;
		combine[i].strand = list2[i - n1].strand;
		combine[i].id = 1;
	}
	qsort(combine, n, sizeof(combine[0]), cmp);
	pre = combine[0].value;
	pre_id = combine[0].id;
	pre_strand = combine[0].strand;
	for (i = 1; i < n; ++i){
		curr = combine[i].value;
		curr_id = combine[i].id;
		curr_strand = combine[i].strand;
		if (curr == pre){
			if (curr_id != pre_id){
				++same_num;
				if (curr_strand == pre_strand){
					//switch (curr_strand){
					++plus_num;
					//	case 1: ++minus_num;
					//}
				}
				else
					++minus_num;
			}
			else
				continue;
		}
		else{
			pre = curr;
			pre_id = curr_id;
			pre_strand = curr_strand;
		}
	}
	free(combine);
	*matched_num = same_num; //plus_num + minus_num;
	return plus_num >= minus_num ? 0 : 1; 
}





// Gnerate seeds sets for the first time
//

// seed_num: how many seeds you wante to generate
// cluster_index: the cluster index that has been created, which is a global variate, start from 0
//////kmer_seed_t * seeds_generate(seqs *arry, option opt, unsigned int start_id, int seed_num, int *cluster_index){
//////	unsigned int seed_index = 0, i, j, k, l, m;
//////	unsigned int seed_kmer_all_index = 0; // the sum of kmer for all seeds
//////	int strand, matched_num = 0;
//////	unsigned int *seed_seq_id = (unsigned int*)malloc(seed_num * sizeof(unsigned int));
//////	float sim;
//////	EdlibAlignResult result;//, result1;
//////	for (i = start_id; ; ++i){
//////		if (arry[i].assigned == 0){
//////			if (seed_index ==0){
//////				seed_seq_id[seed_index] = i;
//////				//for (j = 0; j < arry[i].kmer_num; ++j){
//////				//	seeds_kmer[seed_kmer_all_index].value = arry[start_id].kmerlist->value;
//////				//	++seed_kmer_all_index;
//////				//}
//////				++seed_index;	
//////				++(*cluster_index);
//////				arry[i].cluster = *cluster_index;
//////				arry[i].assigned = 1;
//////				seed_kmer_all_index += arry[i].kmer_num;
//////			}
//////			else {
//////				for (k = 0; k < seed_index; ++k){
//////					l = seed_seq_id[k];
//////					strand = strand_determine(arry[l].kmerlist, arry[l].kmer_num, arry[i].kmerlist, arry[i].kmer_num, &matched_num);
//////					//result1 = edlibAlign(arry[i].seq_is, arry[i].len, arry[l].seq_is, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
//////					if (strand == 0)// plus direction
//////					result = edlibAlign(arry[i].seq_is, arry[i].len, arry[l].seq_is, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
//////					else{ // minus direction
//////						if (arry[l].reverse_calculated == 1) {
//////							result = edlibAlign(arry[i].seq_is, arry[i].len, arry[l].seq_reverse, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
//////						}
//////						else{
//////							arry[i].seq_reverse = (char *) malloc(sizeof(char) * arry[i].len + 1);
//////							for (m = 0; m < arry[i].len; ++m){
//////								switch(arry[i].seq_is[arry[i].len - 1 - m]){
//////									case 'A': 
//////										arry[i].seq_reverse[m] = 'T';
//////										break;
//////									case 'C': 
//////										arry[i].seq_reverse[m] = 'G';
//////										break;
//////									case 'G': 
//////										arry[i].seq_reverse[m] = 'C';
//////										break;
//////									case 'T': 
//////										arry[i].seq_reverse[m] = 'A';
//////										break;
//////									default:
//////										arry[i].seq_reverse[m] = arry[i].seq_is[arry[i].len - 1 - m];
//////
//////								}
//////							}
//////							result = edlibAlign(arry[i].seq_reverse, arry[i].len, arry[l].seq_is, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
//////							arry[i].reverse_calculated = 1;
//////
//////						}
//////					}
//////					
//////					if (result.status == EDLIB_STATUS_OK) {
//////	                        		//printf("%d", result.editDistance);
//////						sim = 1.0 - result.editDistance / (float) result.alignmentLength;
//////					}
//////
//////					if (sim >= opt.cutoff){
//////						arry[i].assigned = 1;
//////						arry[i].cluster = arry[l].cluster;
//////						break;
//////					}
//////					//else { // a new seed is formed
//////					//	++seed_index;
//////					//	seed_seq_id[seed_index] = i;
//////                                	//	++(*cluster_index);
//////                               		//	arry[i].cluster = *cluster_index;
//////					//}
//////
//////				}
//////				if (arry[i].assigned == 0) { // a new seed is formed
//////					//++seed_index;
//////					seed_seq_id[seed_index] = i;
//////					++seed_index;
//////					++(*cluster_index);
//////					arry[i].cluster = *cluster_index;
//////					arry[i].assigned = 1;
//////					seed_kmer_all_index += arry[i].kmer_num;
//////				}
//////
//////			}
//////			if (seed_index >= seed_num)
//////				break;
//////			
//////
//////		}
//////	}
//////	edlibFreeAlignResult(result);
//////	//edlibFreeAlignResult(result1);
//////	kmer_seed_t *seeds_kmer = (kmer_seed_t*)malloc(seed_kmer_all_index * sizeof(kmer_seed_t));
//////	k = 0;
//////	for (i = 0; i < seed_index; ++i){
//////		l = seed_seq_id[i];
//////		for (j = 0; j < arry[l].kmer_num; ++j){
//////			seeds_kmer[k].value  = arry[l].kmerlist[j].value;
//////			seeds_kmer[k].strand = arry[l].kmerlist[j].strand;
//////			seeds_kmer[k].pos = arry[l].kmerlist[j].pos;
//////			seeds_kmer[k].id = l;
//////			seeds_kmer[k].index = i;
//////			++k;
//////		}
//////	}
//////	free(seed_seq_id);
//////	return seeds_kmer;
//////}
//////}

void seq_reverse(seqs *arry){
	arry->seq_reverse = (char *) malloc(sizeof(char) * arry->len);
        for (int m = 0; m < arry->len; ++m){
        	switch(arry->seq_is[arry->len - 1 - m]){
                	case 'A':
                        	arry->seq_reverse[m] = 'T';
                                break;
                        case 'C':
				arry->seq_reverse[m] = 'G';
                                break;
                        case 'G':
				arry->seq_reverse[m] = 'C';
                                break;
                        case 'T':
				arry->seq_reverse[m] = 'A';
                                break;
                        default:
                                arry->seq_reverse[m] = arry->seq_is[arry->len - 1 - m];
		}
	}
	arry->reverse_calculated = 1;
}

/*

// Gnerate seeds sequence sets just based on edlibAlign
// arry:      sequence list
// opt:       option paremeter
// start_id:  start sequence id
// end_id:    end sequence id
// seed_num: how many seeds you wante to generate
// actual_seed_num: the actual seed num that generated
// cluster_index: the cluster index that has been created, which is a global variate, start from 0
// return seeds sequence, then free sequence in the sequence list to save memory
seqs * seeds_seq_generate(seqs *arry, option opt, unsigned int *start_id, unsigned int end_id, int seed_num, int *actual_seed_num, int *cluster_index){
	unsigned int seed_index = 0, i, j, k, l, m;
	unsigned int seed_kmer_all_index = 0; // the sum of kmer for all seeds
	int strand, matched_num = 0;
	unsigned int *seed_seq_id = (unsigned int*)malloc(seed_num * sizeof(unsigned int));
	float sim;
	EdlibAlignResult result;//, result1;
	for (i = *start_id; i < end_id; ++i){
		if (arry[i].assigned == 0){
			if (seed_index ==0){
				seed_seq_id[seed_index] = i;
				//for (j = 0; j < arry[i].kmer_num; ++j){
				//	seeds_kmer[seed_kmer_all_index].value = arry[start_id].kmerlist->value;
				//	++seed_kmer_all_index;
				//}
				++seed_index;	
				++(*cluster_index);
				arry[i].cluster = *cluster_index;
				arry[i].assigned = 1;
				arry[i].cutoff = 0;
				arry[i].is_seed = 1;
				//seed_kmer_all_index += arry[i].kmer_num;
			}
			else {
				for (k = 0; k < seed_index; ++k){
					l = seed_seq_id[k];
					strand = strand_determine(arry[l].kmerlist, arry[l].kmer_num, arry[i].kmerlist, arry[i].kmer_num, &matched_num);
					//result1 = edlibAlign(arry[i].seq_is, arry[i].len, arry[l].seq_is, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
					if (strand == 0)// plus direction
						result = edlibAlign(arry[i].seq_is, arry[i].len, arry[l].seq_is, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
					else{ // minus direction
						if (arry[i].reverse_calculated == 1) {
							result = edlibAlign(arry[i].seq_reverse, arry[i].len, arry[l].seq_is, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
						}
						else{
							seq_reverse(&arry[i]);
							result = edlibAlign(arry[i].seq_reverse, arry[i].len, arry[l].seq_is, arry[l].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
						}
					}
					
					if (result.status == EDLIB_STATUS_OK) {
	                        		//printf("%d", result.editDistance);
						sim = 1.0 - result.editDistance / (float) result.alignmentLength;
					}

					if (sim >= opt.cutoff){
						arry[i].cutoff   = sim;
						arry[i].assigned = 1;
						arry[i].cluster  = arry[l].cluster;
						free(arry[i].seq_is);
						if (arry[i].reverse_calculated == 1)
							free(arry[i].seq_reverse);
						break;
					}
					//else { // a new seed is formed
					//	++seed_index;
					//	seed_seq_id[seed_index] = i;
                                	//	++(*cluster_index);
                               		//	arry[i].cluster = *cluster_index;
					//}

				}
				if (arry[i].assigned == 0) { // a new seed is formed
					//++seed_index;
					seed_seq_id[seed_index] = i;
					++seed_index;
					++(*cluster_index);
					arry[i].cluster = *cluster_index;
					arry[i].assigned = 1;
					arry[i].is_seed = 1;
					//seed_kmer_all_index += arry[i].kmer_num;
				}

			}
			if (seed_index >= seed_num){
				*actual_seed_num = seed_index;
				*start_id = i;
				break;
			}
			

		}
	}
	if (i == end_id){
		*actual_seed_num = seed_index;
		*start_id = i - 1;
		free(seed_seq_id);
		//free(seeds_seq_sets);
		return NULL;
	}
	edlibFreeAlignResult(result);
	//edlibFreeAlignResult(result1);
	seqs *seeds_seq_sets = (seqs *)malloc(seed_index * sizeof(seqs));
	k = 0;
	for (i = 0; i < seed_index; ++i){
		l = seed_seq_id[i];
		seeds_seq_sets[i].id       = arry[l].id;
		seeds_seq_sets[i].len      = arry[l].len;
		seeds_seq_sets[i].cutoff   = arry[l].cutoff;
		seeds_seq_sets[i].cluster  = arry[l].cluster;
		seeds_seq_sets[i].assigned = arry[l].assigned;
		seeds_seq_sets[i].is_seed  = arry[l].is_seed;
		seeds_seq_sets[i].header = (char *) malloc(strlen(arry[l].header));
		strcpy(seeds_seq_sets[i].header, arry[l].header);
		seeds_seq_sets[i].seq_is = (char *) malloc(strlen(arry[l].seq_is));
		strcpy(seeds_seq_sets[i].seq_is, arry[l].seq_is);
		seeds_seq_sets[i].kmer_num = arry[l].kmer_num;
		seeds_seq_sets[i].kmerlist = (kmer_t *)malloc(arry[l].kmer_num * sizeof(seqs));
		for (j = 0; j < arry[l].kmer_num; ++j){
			seeds_seq_sets[i].kmerlist[j].value  = arry[l].kmerlist[j].value;
			seeds_seq_sets[i].kmerlist[j].id     = arry[l].kmerlist[j].id;
			seeds_seq_sets[i].kmerlist[j].pos    = arry[l].kmerlist[j].pos;
			seeds_seq_sets[i].kmerlist[j].strand = arry[l].kmerlist[j].strand;
		}
		free(arry[l].header);
		free(arry[l].seq_is);
		free(arry[l].kmerlist);
	}
	free(seed_seq_id);
	return seeds_seq_sets;
}

*/

/*

// Gnerate seeds sequence sets also based on kmer num
// arry:      sequence list
// opt:       option paremeter
// start_id:  start sequence id
// end_id:    end sequence id
// seed_num: how many seeds you wante to generate
// actual_seed_num: the actual seed num that generated
// cluster_index: the cluster index that has been created, which is a global variate, start from 0
// return seeds sequence, then free sequence in the sequence list to save memory
seqs * seeds_seq_generate_kmer(seqs *arry, option opt, unsigned int *start_id, unsigned int end_id, int seed_num, int *actual_seed_num, int *cluster_index){
	unsigned int seed_index = 0, i, j, k, l, m, matched_num_id;
	unsigned int seed_kmer_all_index = 0; // the sum of kmer for all seeds
	int strand, matched_num = 0, matched_num_max = -1, matched_num_strand;
	unsigned int *seed_seq_id = (unsigned int*)malloc(seed_num * sizeof(unsigned int));
	float sim;
	EdlibAlignResult result;//, result1;
	for (i = *start_id; i < end_id; ++i){
		if (arry[i].assigned == 0){
			if (seed_index ==0){
				seed_seq_id[seed_index] = i;
				//for (j = 0; j < arry[i].kmer_num; ++j){
				//	seeds_kmer[seed_kmer_all_index].value = arry[start_id].kmerlist->value;
				//	++seed_kmer_all_index;
				//}
				++seed_index;	
				++(*cluster_index);
				arry[i].cluster = *cluster_index;
				arry[i].assigned = 1;
				arry[i].cutoff = 0;
				arry[i].is_seed = 1;
				//seed_kmer_all_index += arry[i].kmer_num;
				continue;
			}
			else {
				for (k = 0; k < seed_index; ++k){
					l = seed_seq_id[k];
					strand = strand_determine(arry[l].kmerlist, arry[l].kmer_num, arry[i].kmerlist, arry[i].kmer_num, &matched_num);
					if (matched_num > matched_num_max){
						matched_num_max = matched_num;
						matched_num_strand = strand;
						matched_num_id = l;
					}

				}
			}
			if (matched_num_strand == 0)// plus direction
				result = edlibAlign(arry[i].seq_is, arry[i].len, arry[matched_num_id].seq_is, arry[matched_num_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                        else if (arry[i].reverse_calculated == 1){ // minus direction
				//if (arry[i].reverse_calculated == 1) {
				result = edlibAlign(arry[i].seq_reverse, arry[i].len, arry[matched_num_id].seq_is, arry[matched_num_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
			}
			else{
				seq_reverse(&arry[i]);
                                result = edlibAlign(arry[i].seq_reverse, arry[i].len, arry[matched_num_id].seq_is, arry[matched_num_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
			}
			if (result.status == EDLIB_STATUS_OK) {
				//printf("%d", result.editDistance);
                                sim = 1.0 - result.editDistance / (float) result.alignmentLength;
			}
                        if (sim >= opt.cutoff){
                                arry[i].cutoff   = sim;
                                arry[i].assigned = 1;
                                arry[i].cluster  = arry[matched_num_id].cluster;
                                free(arry[i].seq_is);
                                if (arry[i].reverse_calculated == 1){
                                       	free(arry[i].seq_reverse);
                                }
			}
			else { // a new seed is formed
				//++seed_index;
				seed_seq_id[seed_index] = i;
				++seed_index;
				++(*cluster_index);
				arry[i].cluster = *cluster_index;
				arry[i].assigned = 1;
				arry[i].is_seed = 1;
				//seed_kmer_all_index += arry[i].kmer_num;
			}
			matched_num_max = -1;
			if (seed_index >= seed_num){
				*actual_seed_num = seed_index;
				*start_id = i;
				break;
			}
		}
	}
	edlibFreeAlignResult(result);
	if (i == end_id){
		*actual_seed_num = seed_index;
		*start_id = i - 1;
		free(seed_seq_id);
		//free(seeds_seq_sets);
		return NULL;
	}
	*actual_seed_num = seed_index;
	//edlibFreeAlignResult(result);
	//edlibFreeAlignResult(result1);
	seqs *seeds_seq_sets = (seqs *)malloc(seed_index * sizeof(seqs));
	k = 0;
	for (i = 0; i < seed_index; ++i){
		l = seed_seq_id[i];
		seeds_seq_sets[i].id       = arry[l].id;
		seeds_seq_sets[i].len      = arry[l].len;
		seeds_seq_sets[i].cutoff   = arry[l].cutoff;
		seeds_seq_sets[i].cluster  = arry[l].cluster;
		seeds_seq_sets[i].assigned = arry[l].assigned;
		seeds_seq_sets[i].is_seed  = arry[l].is_seed;
		seeds_seq_sets[i].header = (char *) malloc(strlen(arry[l].header));
		strcpy(seeds_seq_sets[i].header, arry[l].header);
		seeds_seq_sets[i].seq_is = (char *) malloc(strlen(arry[l].seq_is));
		strcpy(seeds_seq_sets[i].seq_is, arry[l].seq_is);
		seeds_seq_sets[i].kmer_num = arry[l].kmer_num;
		seeds_seq_sets[i].kmerlist = (kmer_t *)malloc(arry[l].kmer_num * sizeof(kmer_t));
		for (j = 0; j < arry[l].kmer_num; ++j){
			seeds_seq_sets[i].kmerlist[j].value  = arry[l].kmerlist[j].value;
			seeds_seq_sets[i].kmerlist[j].id     = arry[l].kmerlist[j].id;
			seeds_seq_sets[i].kmerlist[j].pos    = arry[l].kmerlist[j].pos;
			seeds_seq_sets[i].kmerlist[j].strand = arry[l].kmerlist[j].strand;
		}
		//printf("%d\n", i);
		//free(arry[l].header);
		//arry[l].header = NULL;
		free(arry[l].seq_is);
		arry[l].seq_is = NULL;
		free(arry[l].kmerlist);
		arry[l].kmerlist = NULL;
	}
	free(seed_seq_id);
	return seeds_seq_sets;
}

*/

/*
// Gnerate seeds sequence sets also based on kmer num, build rht at the same time, but this is wrong, can not build at the same time, just can return the seqs.
// arry:      sequence list
// opt:       option paremeter
// start_id:  start sequence id
// end_id:    end sequence id
// seed_num: how many seeds you wante to generate
// actual_seed_num: the actual seed num that generated
// cluster_index: the cluster index that has been created, which is a global variate, start from 0
// point_list:    a table having 4*exp(k) cells. Each cell records a pointer corresponding to the list of the windows for a certain k-mer
// seeds_kmer_list: sorted k-mers. This and point_list are refered to the rHAT mapper's idea
// return seeds sequence, then free sequence in the sequence list to save memory
seqs * seeds_rhat_generate_kmer(seqs *arry, option opt, unsigned int *start_id, unsigned int end_id, int seed_num, int *actual_seed_num, int *cluster_index, unsigned int *seeds_kmer_num){
	unsigned int seed_index = 0, i, j, k, l, m, x, matched_num_id;
        *seeds_kmer_num = 0;
	unsigned int seed_kmer_all_index = 0; // the sum of kmer for all seeds
	int strand, matched_num = 0, matched_num_max = -1, matched_num_strand;
	unsigned int *seed_seq_id = (unsigned int*)malloc(seed_num * sizeof(unsigned int));
	float sim;
	EdlibAlignResult result;//, result1;
	for (i = *start_id; i < end_id; ++i){
		if (arry[i].assigned == 0){
			if (seed_index ==0){
				seed_seq_id[seed_index] = i;
				//for (j = 0; j < arry[i].kmer_num; ++j){
				//	seeds_kmer[seed_kmer_all_index].value = arry[start_id].kmerlist->value;
				//	++seed_kmer_all_index;
				//}
				++seed_index;	
				++(*cluster_index);
				arry[i].cluster = *cluster_index;
				arry[i].assigned = 1;
				arry[i].cutoff = 0;
				arry[i].is_seed = 1;
				*seeds_kmer_num += arry[i].kmer_num;
				//seed_kmer_all_index += arry[i].kmer_num;
				continue;
			}
			else {
				for (k = 0; k < seed_index; ++k){
					l = seed_seq_id[k];
					strand = strand_determine(arry[l].kmerlist, arry[l].kmer_num, arry[i].kmerlist, arry[i].kmer_num, &matched_num);
					if (matched_num > matched_num_max){
						matched_num_max = matched_num;
						matched_num_strand = strand;
						matched_num_id = l;
					}

				}
			}
			if (matched_num_strand == 0)// plus direction
				result = edlibAlign(arry[i].seq_is, arry[i].len, arry[matched_num_id].seq_is, arry[matched_num_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                        else if (arry[i].reverse_calculated == 1){ // minus direction
				//if (arry[i].reverse_calculated == 1) {
				result = edlibAlign(arry[i].seq_reverse, arry[i].len, arry[matched_num_id].seq_is, arry[matched_num_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
			}
			else{
				seq_reverse(&arry[i]);
                                result = edlibAlign(arry[i].seq_reverse, arry[i].len, arry[matched_num_id].seq_is, arry[matched_num_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
			}
			if (result.status == EDLIB_STATUS_OK) {
				//printf("%d", result.editDistance);
                                sim = 1.0 - result.editDistance / (float) result.alignmentLength;
			}
                        if (sim >= opt.cutoff){
                                arry[i].cutoff   = sim;
                                arry[i].assigned = 1;
                                arry[i].cluster  = arry[matched_num_id].cluster;
                                free(arry[i].seq_is);
                                if (arry[i].reverse_calculated == 1){
                                       	free(arry[i].seq_reverse);
                                }
			}
			else { // a new seed is formed
				//++seed_index;
				seed_seq_id[seed_index] = i;
				++seed_index;
				++(*cluster_index);
				arry[i].cluster = *cluster_index;
				arry[i].assigned = 1;
				arry[i].is_seed = 1;
				(*seeds_kmer_num) += arry[i].kmer_num;
				//seed_kmer_all_index += arry[i].kmer_num;
			}
			matched_num_max = -1;
			if (seed_index >= seed_num){
				*actual_seed_num = seed_index;
				*start_id = i;
				break;
			}
		}
	}
	edlibFreeAlignResult(result);
	if (i == end_id){
		*actual_seed_num = seed_index;
		*start_id = i - 1;
		free(seed_seq_id);
		//free(seeds_seq_sets);
		return NULL;
	}
	*actual_seed_num = seed_index;
	//edlibFreeAlignResult(result);
	//edlibFreeAlignResult(result1);seeds_kmer_num += arry[i].kmer_num;
	//unsigned int kmer_code_num = 1ULL<< (2 * opt.k);
	//point_list = (unsigned int *)malloc((kmer_code_num + 1) * sizeof(unsigned int));
	//point_list[kmer_code_num] = seeds_kmer_num;
	//memset(point_list, 0,10*sizeof(char));
	//seeds_kmer_list = (kmer_t *)malloc(seeds_kmer_num * sizeof(kmer_t));
	//point_list = (kmer_t *)malloc((seeds_kmer_num + 1) * sizeof(kmer_t));
	seqs *seeds_seq_sets = (seqs *)malloc(seed_index * sizeof(seqs));
	k = 0;
	//x = 0;
	for (i = 0; i < seed_index; ++i){
		l = seed_seq_id[i];
		seeds_seq_sets[i].id       = arry[l].id;
		seeds_seq_sets[i].len      = arry[l].len;
		seeds_seq_sets[i].cutoff   = arry[l].cutoff;
		seeds_seq_sets[i].cluster  = arry[l].cluster;
		seeds_seq_sets[i].assigned = arry[l].assigned;
		seeds_seq_sets[i].is_seed  = arry[l].is_seed;
		seeds_seq_sets[i].header = (char *) malloc(strlen(arry[l].header));
		strcpy(seeds_seq_sets[i].header, arry[l].header);
		seeds_seq_sets[i].seq_is = (char *) malloc(strlen(arry[l].seq_is));
		strcpy(seeds_seq_sets[i].seq_is, arry[l].seq_is);
		seeds_seq_sets[i].kmer_num = arry[l].kmer_num;
		seeds_seq_sets[i].kmerlist = (kmer_t *)malloc(arry[l].kmer_num * sizeof(kmer_t));
		for (j = 0; j < arry[l].kmer_num; ++j){
			seeds_seq_sets[i].kmerlist[j].value  = arry[l].kmerlist[j].value;
			seeds_seq_sets[i].kmerlist[j].id     = arry[l].kmerlist[j].id;
			seeds_seq_sets[i].kmerlist[j].pos    = arry[l].kmerlist[j].pos;
			seeds_seq_sets[i].kmerlist[j].strand = arry[l].kmerlist[j].strand;
			//seeds_kmer_list[x].value  = arry[l].kmerlist[j].value;
			//seeds_kmer_list[x].id     = arry[l].kmerlist[j].id;
			//seeds_kmer_list[x].pos    = arry[l].kmerlist[j].pos;
			//seeds_kmer_list[x].strand = arry[l].kmerlist[j].strand;
			//++x;
		}
		//printf("%d\n", i);
		//free(arry[l].header);
		//arry[l].header = NULL;
		free(arry[l].seq_is);
		arry[l].seq_is = NULL;
		free(arry[l].kmerlist);
		arry[l].kmerlist = NULL;
	}
	free(seed_seq_id);
	//qsort(seeds_kmer_list, seeds_kmer_num, sizeof(seeds_kmer_list[0]), cmp);
	//uint64_t pre,  curr, kk, pre_value, cur_value;
	//unsigned int pre_id;
	//short int pre_strand;
	//pre = seeds_kmer_list[0].value;
	//point_list[pre] = 0;
	//point_list[kmer_code_num] = seeds_kmer_num;
	//for (kk = 0; kk < pre; ++kk){
	//	point_list[kk] = point_list[pre];
	//}
	//pre_value = pre;
        //pre_id = seeds_kmer_list[0].id;
        //pre_strand = seeds_kmer_list[0].strand;
        //for (i = 1; i < seeds_kmer_num; ++i){
        //        curr = seeds_kmer_list[i].value;
                //curr_id = seeds_kmer_list[i].id;
                //curr_strand = seeds_kmer_list[i].strand;
        //        if (curr == pre){
                        //if (curr_id != pre_id){
                        //        ++same_num;
                        //        if (curr_strand == pre_strand){
                        //                //switch (curr_strand){
                        //                ++plus_num;
                        //                //      case 1: ++minus_num;
                        //                //}
                        //        }
                        //        else
                        //                ++minus_num;
                        //}
                        //else
         //                       continue;
         //       }
         //       else{
         //               pre = curr;
	//		point_list[pre] = i;
	//		for (kk = pre_value + 1; kk < curr; ++kk){
	//			point_list[kk] = point_list[pre];
	//		}
	//		pre_value = curr;
                        //pre_id = curr_id;
                        //pre_strand = curr_strand;
        //        }
       // }
        //free(combine);
	//for (kk = pre_value + 1; kk < kmer_code_num; ++kk){
	//	point_list[kk] = point_list[kmer_code_num];
	//}
	return seeds_seq_sets;
}


*/
kmer_t * rhat_generate(seqs *seeds_seq_sets, int seeds_num, unsigned int *point_list, unsigned int seeds_kmer_num, unsigned int kmer_code_num){
	unsigned int i, j, x = 0;
	kmer_t *seeds_kmer_list;
	seeds_kmer_list = (kmer_t *)malloc(seeds_kmer_num * sizeof(kmer_t));
	//unsigned int j;
	for (i = 0; i < seeds_num; ++i){
                //j = seeds_seq_sets[i].kmer_num;
                for (j = 0; j < seeds_seq_sets[i].kmer_num; ++j){
                        //seeds_seq_sets[i].kmerlist[j].value  = arry[l].kmerlist[j].value;
                        //seeds_seq_sets[i].kmerlist[j].id     = arry[l].kmerlist[j].id;
                        //seeds_seq_sets[i].kmerlist[j].pos    = arry[l].kmerlist[j].pos;
                        //seeds_seq_sets[i].kmerlist[j].strand = arry[l].kmerlist[j].strand;
                        seeds_kmer_list[x].value  = seeds_seq_sets[i].kmerlist[j].value;
                        seeds_kmer_list[x].id     = i;// seed index      seeds_seq_sets[i].kmerlist[j].id;
                        seeds_kmer_list[x].pos    = seeds_seq_sets[i].kmerlist[j].pos;
                        seeds_kmer_list[x].strand = seeds_seq_sets[i].kmerlist[j].strand;
                        ++x;
                }
                //printf("%d\n", i);
                //free(arry[l].header);
                //arry[l].header = NULL;
                //free(arry[l].seq_is);
                //arry[l].seq_is = NULL;
	}
	time_t now = time(0);
        char* dt = ctime(&now);
        cout << "[INFO] Sort started at: " << dt << endl;
	sort(seeds_kmer_list, seeds_kmer_list + seeds_kmer_num, cmp_sort);
	now = time(0);
	dt = ctime(&now);
	cout << "[INFO] Sort ended at: " << dt << endl;
	//qsort(seeds_kmer_list, seeds_kmer_num, sizeof(seeds_kmer_list[0]), cmp);
	//now = time(0);
        //dt = ctime(&now);
        //cout << "[INFO] Qsort ended at: " << dt << endl;
	uint64_t pre,  curr, kk, pre_value, cur_value;
	//unsigned int pre_id;
	//short int pre_strand;
	pre = seeds_kmer_list[0].value;
	point_list[pre] = 0;
	//point_list[kmer_code_num] = seeds_kmer_num;
	for (kk = 0; kk < pre; ++kk){
		point_list[kk] = point_list[pre];
	}
	pre_value = pre;
        //pre_id = seeds_kmer_list[0].id;
        //pre_strand = seeds_kmer_list[0].strand;
        for (i = 1; i < seeds_kmer_num; ++i){
                curr = seeds_kmer_list[i].value;
		//cout << i << "\t" << curr << endl;
		//printf( "\r i = %11i, hash value = %ld ", i, curr);
		//printf( "\r%4.1f%%", p );
		//fflush( stdout );
                //curr_id = seeds_kmer_list[i].id;
                //curr_strand = seeds_kmer_list[i].strand;
                if (curr == pre){
			continue;
                }
                else{
                        pre = curr;
			point_list[pre] = i;
			for (kk = pre_value + 1; kk < curr; ++kk){
				point_list[kk] = point_list[pre];
			}
			pre_value = curr;
                        //pre_id = curr_id;
                        //pre_strand = curr_strand;
                }
        }
        //free(combine);
	for (kk = pre_value + 1; kk < kmer_code_num; ++kk){
		point_list[kk] = point_list[kmer_code_num];
	}
	return seeds_kmer_list;
}





/*

// Just cluster one sequence based the existing seeds
// seq:		sequence that will be processed
// opt:		option paremeters
// seeds_seq:	seeds sequences list
// seed_num:	number of seed sequence
// id:		sequence id for cluster
// Just based on the top number of identical kmers, not based on the k-mer distance like ESPRIT

int cluster_one(seqs *arry, option opt, seqs *seeds_seq, unsigned int seed_num, uint64_t id){
	
	//unsigned int *point_list;
        //kmer_t *seeds_kmer_list;
	//unsigned int kmer_code_num = 1ULL<< (2 * opt.k);
        //point_list = (unsigned int *)malloc((kmer_code_num + 1) * sizeof(unsigned int));
        //point_list[kmer_code_num] = seeds_kmer_num;
        //memset(point_list, 0,10*sizeof(char));
        //seeds_kmer_list = (kmer_t *)malloc(seeds_kmer_num * sizeof(kmer_t));

	if (arry->assigned == 1)
		return 1;
	int m, matched_max = -1, matched_max_seed_id, strand, matched_max_strand, matched_num; // restore the max matched kmer number to all seeds
	float sim;// *sim_list, *sim_list_rev;
	//int *matched_num_list, *strand_list;
	EdlibAlignResult result;
	//sim_list = (float *)malloc(seed_num * sizeof(float));
	//sim_list_rev = (float *)malloc(seed_num * sizeof(float));
	//strand_list = (int *)malloc(seed_num * sizeof(int));
	//matched_num_list = (int *)malloc(seed_num * sizeof(int));
	for (int i = 0; i < seed_num; ++i){
		//result = edlibAlign(arry->seq_is, arry->len, seeds_seq[i].seq_is, seeds_seq[i].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		//if (result.status == EDLIB_STATUS_OK) {
                	//printf("%d", result.editDistance);
                //	sim = 1.0 - result.editDistance / (float) result.alignmentLength;
        	//}
		//sim_list[i] = sim;
		//seq_reverse(arry);
		//result = edlibAlign(arry->seq_reverse, arry->len, arry[i].seq_is, arry[i].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		//if (result.status == EDLIB_STATUS_OK) {
                        //printf("%d", result.editDistance);
                //        sim = 1.0 - result.editDistance / (float) result.alignmentLength;
                //}
		//sim_list_rev[i] = sim;
		
		strand = strand_determine(arry->kmerlist, arry->kmer_num, seeds_seq[i].kmerlist, seeds_seq[i].kmer_num, &matched_num);
		//matched_num_list[i] = matched_num;
		//strand_list[i] = strand;
		if (matched_num > matched_max){
			matched_max = matched_num;
			matched_max_strand = strand;
			matched_max_seed_id = i;

		}
	}
	if (matched_max_strand == 0)
		result = edlibAlign(arry->seq_is, arry->len, seeds_seq[matched_max_seed_id].seq_is, seeds_seq[matched_max_seed_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	else{ // minus direction
		if (arry->reverse_calculated == 1) {
			result = edlibAlign(arry->seq_reverse, arry->len, seeds_seq[matched_max_seed_id].seq_is, seeds_seq[matched_max_seed_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		}
		else{
			seq_reverse(arry);
			result = edlibAlign(arry->seq_reverse, arry->len, seeds_seq[matched_max_seed_id].seq_is, seeds_seq[matched_max_seed_id].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		}
	}
	if (result.status == EDLIB_STATUS_OK) {
		//printf("%d", result.editDistance);
		sim = 1.0 - result.editDistance / (float) result.alignmentLength;
	}
	if (sim >= opt.cutoff){
		arry->cutoff   = sim;
		arry->assigned = 1;
		arry->cluster  = seeds_seq[matched_max_seed_id].cluster;
		free(arry->seq_is);
		arry->seq_is = NULL;
		if (arry->reverse_calculated == 1){
			free(arry->seq_reverse);
			arry->seq_reverse = NULL;
		}
		free(arry->kmerlist);
		arry->kmerlist = NULL;
		//break;
	}
	//unsigned int seed_id, indexx = 0, flag = 0;
	//match_t *match_list = (match_t*)malloc(seed_num * sizeof(match_t));
	//kmer_t  *kmer_list_id = (kmer_t*)malloc(knumm * sizeof(kmer_t));
        //mm_extract(array[i].seq_is, array[i].len, opt.k, opt.knum, i, kmer_list_id);
	//qsort(kmer_list_id, knumm, sizeof(kmer_list_id[0]), cmp);
	//if (kmer_list_id[0] > seed_kmer[seed_kmer_num - 1] || kmer_list_id[0] < seed_kmer[0)
	//	return 0; // not assigned
	//for (uint64_t i = 0; i < knumm; ++i){
	//	for (uint64_t j = 0; j < seed_kmer_num; ++j){
	//	if (kmer_list_id[i].value == seed_kmer[j].value && flag = 0){
	//			seed_id = seed_kmer[j].id;
	//			match_list[indexx].match_num++;
	//			match_list[indexx].
	//		}
	//}
	//int firstseed =
}

*/

/*

// Just cluster one sequence based the existing seeds and the regional hash table
// seq:		sequence that will be processed
// opt:		option paremeters
// seeds_seq:	seeds sequences list
// seed_num:	number of seed sequence
// id:		sequence id for cluster
// point_list:  pointer list
// seeds_kmer_list: 
// Just based on the top number of identical kmers, not based on the k-mer distance like ESPRIT

int cluster_one_with_rhat(seqs *arry, option opt, seqs *seeds_seq, unsigned int seed_num, uint64_t id, unsigned int *point_list, kmer_t *seeds_kmer_list){
	

	if (arry->assigned == 1)
		return 1;
	int i, m, kmer_num, matched_max = -1, matched_max_seed_id, strand, matched_max_strand, matched_num; // restore the max matched kmer number to all seeds
	matched_num_t *matched_list; // store the matched k-mer num with each seed
	int *plus_num_list, *minus_num_list, seed_idd;
	uint64_t kmer_is;
	unsigned int start, endd, k, seed_id;
	matched_list   = (matched_num_t *)malloc(seed_num * sizeof(matched_num_t));
	for (i = 0; i < seed_num; ++i){
		matched_list[i].matched_num = 0;
		matched_list[i].id = i;
	}
	plus_num_list  = (int *)malloc(seed_num * sizeof(int));
	minus_num_list = (int *)malloc(seed_num * sizeof(int));
	memset(plus_num_list, 0, seed_num * sizeof(int));
	memset(minus_num_list, 0, seed_num * sizeof(int));
	float sim;// *sim_list, *sim_list_rev;
	uint64_t pre, curr;
	EdlibAlignResult result;
	kmer_num = arry->kmer_num;
	qsort(arry->kmerlist, kmer_num, sizeof(arry->kmerlist[0]), cmp);
	pre = arry->kmerlist[0].value;
	kmer_is = pre;
	start = point_list[kmer_is];
	endd  = point_list[kmer_is + 1] - 1;
	for (k = start; k <= endd; ++k){
		seed_id = seeds_kmer_list[k].id;
                if (arry->kmerlist[i].strand == seeds_kmer_list[k].strand)
                	plus_num_list[seed_id] += 1;
                else
			minus_num_list[seed_id] += 1;

                matched_list[seed_id].matched_num += 1;
	}

	for (i = 1; i < kmer_num; ++i){
		kmer_is = arry->kmerlist[i].value;
		curr = kmer_is;
		if (curr == pre)
			continue;
		start = point_list[kmer_is];
		endd  = point_list[kmer_is + 1] - 1;
		for (k = start; k <= endd; ++k){
			seed_id = seeds_kmer_list[k].id;
			//printf ("seed_id = %d, kmer_is = %ld, id = %ld\n", seed_id, kmer_is, id);
			//if (id == 103 && kmer_is == 976)
			//	int aaaa = 1;
			if (arry->kmerlist[i].strand == seeds_kmer_list[k].strand)
				plus_num_list[seed_id] += 1;
			else
				minus_num_list[seed_id] += 1;

			matched_list[seed_id].matched_num += 1;
		}
		pre = kmer_is;
	}
	qsort(matched_list, seed_num, sizeof(matched_list[0]), cmp_matched_num);
	for (i = 0; i < 2; ++i){
		seed_idd = matched_list[i].id;
		if (plus_num_list[seed_idd] > minus_num_list[seed_idd])
			result = edlibAlign(arry->seq_is, arry->len, seeds_seq[seed_idd].seq_is, seeds_seq[seed_idd].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		else{
			if (arry->reverse_calculated == 1) {
                        result = edlibAlign(arry->seq_reverse, arry->len, seeds_seq[seed_idd].seq_is, seeds_seq[seed_idd].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                	}
                	else{
                        	seq_reverse(arry);
                        	result = edlibAlign(arry->seq_reverse, arry->len, seeds_seq[seed_idd].seq_is, seeds_seq[seed_idd].len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                	}

		}
		if (result.status == EDLIB_STATUS_OK) {
                	//printf("%d", result.editDistance);
                	sim = 1.0 - result.editDistance / (float) result.alignmentLength;
        	}
		if (sim >= opt.cutoff){
                	arry->cutoff   = sim;
                	arry->assigned = 1;
                	arry->cluster  = seeds_seq[seed_idd].cluster;
                	free(arry->seq_is);
                	arry->seq_is = NULL;
                	if (arry->reverse_calculated == 1){
                        	free(arry->seq_reverse);
                        	arry->seq_reverse = NULL;
                	}
                	free(arry->kmerlist);
                	arry->kmerlist = NULL;
                	break;
        	}

	}
	free(matched_list); // store the matched k-mer num with each seed
        free(plus_num_list);
	free(minus_num_list);
}

int cluster_multi_threads(seqs *array, option opt, unsigned int seqs_total){
	unsigned int i, j, end_id;
	unsigned int start_id = 0;
	int seed_num = 50, actual_seed_num = 0, cluster_index = -1;
	unsigned int seeds_kmer_num;
	unsigned int *point_list;
	kmer_t *seeds_kmer_list;
	unsigned int kmer_code_num = 1ULL<< (2 * opt.k);
	float p1, p0 = 0;
	clock_t t;
	int tid;
	end_id = seqs_total;
	point_list = (unsigned int *)malloc((kmer_code_num + 1) * sizeof(unsigned int));
	for (i = 0; i < seqs_total; ++i){
		// first generate seeds set
		seqs *seeds_seqs;
        	//point_list = (unsigned int *)malloc((kmer_code_num + 1) * sizeof(unsigned int));
        	//point_list[kmer_code_num] = seeds_kmer_num;
        	//memset(point_list, 0,10*sizeof(char));
        	//seeds_kmer_list = (kmer_t *)malloc(seeds_kmer_num * sizeof(kmer_t));
		//kmer_t *seeds_kmer_list;
		t = clock();
		//seeds_seqs = seeds_seq_generate_kmer(array, opt, &start_id, end_id, seed_num, &actual_seed_num, &cluster_index);
		seeds_seqs = seeds_rhat_generate_kmer(array, opt, &start_id, end_id, seed_num, &actual_seed_num, &cluster_index, &seeds_kmer_num);
		printf("\n[INFO] Generating seeds set time: %.3f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
                printf( "\r%d  sequences finished, new seeds number: %9d\n", start_id + 1, actual_seed_num);
		if (start_id == seqs_total - 1)
			return 0;
		point_list[kmer_code_num] = seeds_kmer_num;
		seeds_kmer_list = rhat_generate(seeds_seqs, actual_seed_num, point_list, seeds_kmer_num, kmer_code_num);
		//point_list[kmer_code_num] = seeds_kmer_num;
		//printf("\n[INFO] Generating seeds set time: %.3f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		//printf( "\r%d  sequences finished, new seeds number: %9d\n", start_id + 1, actual_seed_num);
		// then cluster each sequence based on the seeds set
		#pragma omp parallel for
		for (j = start_id + 1; j < seqs_total; ++j){
			//cluster_one(&array[j], opt, seeds_seqs, actual_seed_num, j);
			cluster_one_with_rhat(&array[j], opt, seeds_seqs, actual_seed_num, j, point_list, seeds_kmer_list);
			//p1 = (100.0 * j) / seqs_total;
			tid = omp_get_thread_num();
			if (tid == 1){
				p1 = (100.0 * j) / seqs_total;
				if (p1 > p0+1E-2){
					printf( "\r%5.2f%%  sequences finished", p1);
                			fflush( stdout );
					p0 = p1;
				}
			}
			//printf("j = %d\n", j);
		}
		start_id += 1;
		i = start_id;
		free(seeds_seqs);
	
	// then generate seeds set again
	}
	return 0;
}
*/

int sequence_map_one_with_locate(option opt, seqs *reads, seqs *gen_list, short int *array_gen_id, bool *array_plus, unsigned int *array_seq_start, unsigned int *array_seq_end, unsigned int *array_gen_start, unsigned int *array_gen_end){
	unsigned int i, jj, j, k, c = 0, tolerant_length = reads->len, id_gen;
	unsigned int seq_pos_end, gen_pos_end, gen_pos_start, seq_pos_start, len_seq_for_align, len_gen_for_align;
	int ii;
	EdlibAlignResult result;
        string numm, temp_str, cigar_reverse;
        char* cigar;
        string cigar_str, algined = "";
        id_gen = array_gen_id[reads->id];
        reads->target_genome_id = id_gen;
        gen_pos_end = array_gen_end[reads->id];
        gen_pos_start = array_gen_start[reads->id];
        len_gen_for_align = gen_pos_end - gen_pos_start + 1;
	if (array_plus[reads->id]){// plus direction
		seq_pos_end = array_seq_end[reads->id];
		seq_pos_start = array_seq_start[reads->id];
		len_seq_for_align = seq_pos_end - seq_pos_start + 1;
		if (len_seq_for_align <= len_gen_for_align){
			result = edlibAlign(reads->seq_is + seq_pos_start, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_start, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
			cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			cigar_str = cigar;
			algined += cigar_str;
               		//printf("%s\n", cigar);
               		//printAlignment(reads->seq_is + seq_pos_start, gen_list[id_gen].seq_is + gen_pos_start, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_HW);
               		//////if (*result.endLocations != len_gen_for_align - 1){
                       	//////	algined += to_string(len_gen_for_align - 1 - *result.endLocations);
                       	//////	algined += "D";
               		//////}
		}
		else{
			result = edlibAlign(gen_list[id_gen].seq_is + gen_pos_start, len_gen_for_align, reads->seq_is + seq_pos_start, len_seq_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
			cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			//printf("%s\n", cigar);
                        //printAlignment(gen_list[id_gen].seq_is + gen_pos_start, reads->seq_is + seq_pos_start, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_HW);
			// convert cigar to query status, not genome status.
			cigar_str = cigar;
			for(ii = 0; ii < cigar_str.length(); ++ii){
				if (cigar_str[ii] == 'I'){
					cigar_str[ii] = 'D';
				}
				else if (cigar_str[ii] == 'D'){
					cigar_str[ii] = 'I';
				}
			}
			algined += cigar_str;
		}
		// align the two ends region
		// first the right end
		if (seq_pos_end + 1 < reads->len) {
			len_seq_for_align = reads->len - seq_pos_end - 1;
			len_gen_for_align = 2 * len_seq_for_align;
			if (gen_pos_end + 1 + len_gen_for_align >= gen_list[id_gen].len){
				len_gen_for_align = gen_list[id_gen].len - 1 - (gen_pos_end + 1);
			}
			result = edlibAlign(reads->seq_is + seq_pos_end + 1, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_end + 1, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			cigar_str = cigar;
                        algined += cigar_str;
			//printf("%s\n", cigar);
                        //printAlignment(reads->seq_is + seq_pos_end + 1, gen_list[id_gen].seq_is + gen_pos_end + 1, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);

		}
		// then the left end region
		if (seq_pos_start != 0){
			len_seq_for_align = seq_pos_start;
			len_gen_for_align = 2 * len_seq_for_align;
			char * start_seq_char, * start_gen_char;
		       	start_seq_char = (char *)malloc(len_seq_for_align * sizeof(char));
			start_gen_char = (char *)malloc(len_gen_for_align * sizeof(char));
			for(i = 0; i < len_seq_for_align; ++i){
				start_seq_char[i] = reads->seq_is[len_seq_for_align - i - 1];
				//start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
			}
			if (gen_pos_start  < len_gen_for_align){
				len_gen_for_align = gen_pos_start;
				//for (i = 0; i < len_gen_for_align; ++i){
				//	start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
				//}
			}
			for(i = 0; i < len_gen_for_align; ++i){
				start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
			}
			result = edlibAlign(start_seq_char, len_seq_for_align, start_gen_char, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
			reads->aligned_pos = gen_pos_start - *result.endLocations;
			cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			cigar_str = cigar;
			//printf("%s\n", cigar);
                        //printAlignment(start_seq_char, start_gen_char, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);
			numm = "";
			temp_str = "";
			cigar_reverse = "";
			for(ii = cigar_str.length() - 1; ii >= 0; --ii){
				if (cigar_str[ii] == 'I'){

					cigar_reverse += numm + temp_str;
					temp_str = "I";
					numm = "";
				}
				else if (cigar_str[ii] == 'D'){
					cigar_reverse += numm + temp_str;
					temp_str = "D";
					numm = "";
				}
				else if (cigar_str[ii] == 'M'){
					cigar_reverse += numm + temp_str;
					temp_str = "M";
					numm = "";
				}
				else{
					numm = cigar_str[ii] + numm;
					//cigar_reverse = cigar_str[i] + cigar_reverse;
				}
				
			}

			cigar_reverse += numm + temp_str;
			algined = cigar_reverse + algined;
			reads->cigar =  (char *) malloc(sizeof(char) * (algined.length()+ 1));
                        strcpy(reads->cigar, algined.c_str());
			//reads->aligned_pos = gen_pos_start - *result.endLocations - 1;
			free(start_seq_char);
			free(start_gen_char);
		}
		else{
			reads->cigar =  (char *) malloc(sizeof(char) * (algined.length()+ 1));
			strcpy(reads->cigar, algined.c_str());
			reads->aligned_pos = gen_pos_start + 1;
			
		}
	}
	else{ // minus direction
		//////FILE *fr;
	        //////fr=fopen(opt.output,"w");
        	//////if (fr == NULL){
                //////	printf("File not open!/n");
        	//////}
        	//////else{
                //////	for (ii = start_i; ii <= end_i; ++ii){
                //////        	fprintf(fr, "%d\t%d\n", mm_list[ii].pos_gen, mm_list[ii].pos_seq);
                //////	}

        	//////}
		//////fclose(fr);

		reads->seq_reverse = (char *) malloc(sizeof(char) * (reads->len + 1));
		for (int m = 0; m < reads->len; ++m){
			switch(reads->seq_is[reads->len - 1 - m]){
				case 'A':
					reads->seq_reverse[m] = 'T';
					break;
				case 'C':
					reads->seq_reverse[m] = 'G';
                                        break;
                                case 'G':
					reads->seq_reverse[m] = 'C';
                                        break;
                                case 'T':
					reads->seq_reverse[m] = 'A';
                                        break;
                                default:
		 			reads->seq_reverse[m] = reads->seq_is[reads->len - 1 - m];
			}
                }
		reads->seq_reverse[reads->len] = '\0';
		reads->reverse_calculated = 1;
		seq_pos_start = array_seq_start[reads->id];
                seq_pos_end = array_seq_end[reads->id];
                len_seq_for_align = seq_pos_end - seq_pos_start + 1;
		if (len_seq_for_align <= len_gen_for_align){
                        result = edlibAlign(reads->seq_reverse + seq_pos_start, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_start, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        //////if (*result.startLocations != 0){
                        //////  algined += to_string(*result.startLocations);
                        //////  algined += "D";
                        //////}
                        cigar_str = cigar;
                        algined += cigar_str;
                        //printf("%s\n", cigar);
                        //printAlignment(reads->seq_reverse + seq_pos_start, gen_list[id_gen].seq_is + gen_pos_start, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_HW);
                        //////if (*result.endLocations != len_gen_for_align - 1){
                        //////  algined += to_string(len_gen_for_align - 1 - *result.endLocations);
                        //////  algined += "D";
                        //////}
                }
		else{
                        result = edlibAlign(gen_list[id_gen].seq_is + gen_pos_start, len_gen_for_align, reads->seq_reverse + seq_pos_start, len_seq_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        //printf("%s\n", cigar);
                        //printAlignment(gen_list[id_gen].seq_is + gen_pos_start, reads->seq_is + seq_pos_start, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_HW);
                        cigar_str = cigar;
			for(ii = 0; ii < cigar_str.length(); ++ii){
                                if (cigar_str[ii] == 'I'){
                                        cigar_str[ii] = 'D';
                                }
                                else if (cigar_str[ii] == 'D'){
                                        cigar_str[ii] = 'I';
                                }
                        }
                        algined += cigar_str;
                }

		// align the two ends region
                // first the right end
                if (seq_pos_end + 1 < reads->len) {
                        len_seq_for_align = reads->len - seq_pos_end - 1;
                        len_gen_for_align = 2 * len_seq_for_align;
                        if (gen_pos_end + 1 == gen_list[id_gen].len){
				cigar_str = to_string(len_seq_for_align) + "I";
			}
			else {
				if (gen_pos_end + 1 + len_gen_for_align >= gen_list[id_gen].len){
                                	len_gen_for_align = gen_list[id_gen].len - 1 - (gen_pos_end + 1);
                        	}
                        	result = edlibAlign(reads->seq_reverse + seq_pos_end + 1, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_end + 1, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
                        	cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        	cigar_str = cigar;
			}
                        algined += cigar_str;
                        //printf("%s\n", cigar);
                        //printAlignment(reads->seq_reverse + seq_pos_end + 1, gen_list[id_gen].seq_is + gen_pos_end + 1, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);

                }

		// then the left end region
                if (seq_pos_start != 0){
                        len_seq_for_align = seq_pos_start;
                        len_gen_for_align = 2 * len_seq_for_align;
                        char * start_seq_char, * start_gen_char;
                        start_seq_char = (char *)malloc(len_seq_for_align * sizeof(char));
                        start_gen_char = (char *)malloc(len_gen_for_align * sizeof(char));
                        for(i = 0; i < len_seq_for_align; ++i){
                                start_seq_char[i] = reads->seq_reverse[len_seq_for_align - i - 1];
                                //start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
                        }
                        if (gen_pos_start  < len_gen_for_align){
                                len_gen_for_align = gen_pos_start;
                                //for (i = 0; i < len_gen_for_align; ++i){
                                //      start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
                                //}
                        }
                        for(i = 0; i < len_gen_for_align; ++i){
                                start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
                        }
			result = edlibAlign(start_seq_char, len_seq_for_align, start_gen_char, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        cigar_str = cigar;
                        //printf("%s\n", cigar);
                        //printAlignment(start_seq_char, start_gen_char, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);
			numm = "";
			temp_str = "";
                        cigar_reverse = "";
                        for(ii = cigar_str.length() - 1; ii >= 0; --ii){
                                if (cigar_str[ii] == 'I'){

                                        cigar_reverse += numm + temp_str;
                                        temp_str = "I";
                                        numm = "";
                                }
                                else if (cigar_str[ii] == 'D'){
                                        cigar_reverse += numm + temp_str;
                                        temp_str = "D";
                                        numm = "";
                                }
                                else if (cigar_str[ii] == 'M'){
                                        cigar_reverse += numm + temp_str;
                                        temp_str = "M";
                                        numm = "";
                                }
                                else{
                                        numm = cigar_str[ii] + numm;
                                        //cigar_reverse = cigar_str[i] + cigar_reverse;
                                }

                        }
                        cigar_reverse += numm + temp_str;
                        algined = cigar_reverse + algined;
			reads->cigar =  (char *) malloc(sizeof(char) * (algined.length()+ 1));
                        strcpy(reads->cigar, algined.c_str());
                        free(start_seq_char);
			free(start_gen_char);
			reads->aligned_pos = gen_pos_start - *result.endLocations;
		}
		else{
			reads->cigar =  (char *) malloc(sizeof(char) * (algined.length()+ 1));
			strcpy(reads->cigar, algined.c_str());
                        reads->aligned_pos = gen_pos_start + 1;
                }		
		//free(reads->seq_reverse);
		
	}
	
	return 1;
}

int sequence_map_one(option opt, kmer_t *seeds_kmer_list, unsigned int *point_list, seqs *reads, seqs *gen_list){
	unsigned int i, jj, j, kmer_match_num = 0, k, c = 0, tolerant_length = reads->len;
	int ii;
	uint64_t hashh;//, hash_now, hash_pre;
	//float peak_cover = 0.0, covered;
	//set<unsigned int> set_gen, set_seq;
	//set<uint64_t> set_hash;
	//map<uint64_t, unsigned int> map_gen, map_seq;
	unsigned int peak_i, peak_ii, gen_position_i, gen_position_ii, seq_position_i, seq_position_ii, start_ii, end_i, end_i1, start_i, start_i1, peak_max = 0, peak_max2 = 0, chain_num1 = 0, chain_num2 = 0, match_selected_num, id_gen, diff_pos_gen, diff_pos_seq;
	j = kmer_extract_for_one_seq(reads, opt.k, opt.w);
	for (i = 0; i < reads->kmer_num; ++i){
		hashh = reads->kmerlist[i].value;
		kmer_match_num += point_list[hashh + 1] - point_list[hashh];
	}
	mm_match_t *mm_list = (mm_match_t*)malloc(kmer_match_num * sizeof(mm_match_t));
	for (i = 0; i < reads->kmer_num; ++i){
		hashh = reads->kmerlist[i].value;
		k = point_list[hashh + 1] - point_list[hashh];
		for (j = point_list[hashh]; j < point_list[hashh + 1]; ++j){
			mm_list[c].value = hashh;
			mm_list[c].pos_gen = seeds_kmer_list[j].pos;
			mm_list[c].pos_seq = reads->kmerlist[i].pos;
			mm_list[c].strand_gen = seeds_kmer_list[j].strand;
			mm_list[c].strand_seq = reads->kmerlist[i].strand;
			mm_list[c].id_gen = seeds_kmer_list[j].id;
			++c;
		}
	}
	qsort(mm_list, kmer_match_num, sizeof(mm_list[0]), cmp_mm_match_t);
	for (i = 0; i < kmer_match_num; ++i){
		peak_i = 0;
		peak_ii = 0;
		gen_position_i = mm_list[i].pos_gen;
		seq_position_i = mm_list[i].pos_seq;
		//pre_pos = gen_position_i;
		for (ii = i - 1; ii > 0; ii--){
			gen_position_ii = mm_list[ii].pos_gen;
			diff_pos_gen = gen_position_i - gen_position_ii;
			if (diff_pos_gen > tolerant_length){
				break;
                        }
			else{
				seq_position_ii = mm_list[ii].pos_seq;
				diff_pos_seq = abs((int)(seq_position_i - seq_position_ii));
				if (abs((int)(diff_pos_seq - diff_pos_gen)) <= 200){
					peak_ii++;
				}
                                peak_i++;
                                start_ii = ii;
                        }
		}
		if (peak_i > peak_max){
                	start_i = start_ii;
                        end_i   = i;
                        peak_max = peak_i;
                }
		if (peak_ii > peak_max2) {
			start_i1 = start_ii;
			end_i1   = i;
			peak_max2 = peak_ii; //density = (set_gen.size() + set_seq.size()) * 1.0 / reads->len ;
		}
	}
	// next step: align the middle region contains seeds
	//////match_selected_num = end_i - start_i + 1;
	int strand, plus_num = 0, minus_num = 0, match_n;
	//////unsigned int *chain1 = new unsigned int [match_selected_num]; // chain1 store the tem chain;
       	//////unsigned int *chain2 = new unsigned int [match_selected_num]; // chain2 store the final chain;
	unsigned int seq_pos_end, gen_pos_end, gen_pos_start, seq_pos_start, len_seq_for_align, len_gen_for_align;
	for (jj = start_i; jj < end_i; ++jj){  // determine the plus or minus direction
	//////	chain_num1 = 0;
	//////	chain1[chain_num1] = ii;
		if (mm_list[jj].strand_gen == mm_list[jj].strand_seq) {
			++plus_num;
		}
		else{
			++minus_num;
		}
	//////	for (jj = ii + 1; jj <= end_i; ++jj){
	//////		if (mm_list[jj].pos_seq >= mm_list[ii].pos_seq){
	//////			chain_num1 += 1;
	//////			chain1[chain_num1] = jj;
	//////			if (mm_list[jj].strand_gen == mm_list[jj].strand_seq){
	//////				++plus_num;
	//////			}
	//////			else{
	//////				++minus_num;
	//////			}
	//////		}
	//////	}
	//////	if (chain_num1 > chain_num2){
	//////		for (jj = 0; jj <= chain_num1; ++jj){
	//////			chain2[jj] = chain1[jj];
	//////		}
	//////		chain_num2 = chain_num1;
		if (plus_num >= minus_num){
			strand = 0;
		}
		else {
			strand = -1;
		}
	//////	}
	//////	if (chain_num1 + 1 == match_selected_num){
	//////		break;
	//////	}
	}
	// next step: get the alignment using edlib for the region from first seed to last seed
	EdlibAlignResult result;
	string numm, temp_str, cigar_reverse;
	match_n = opt.k;
	char* cigar;
	string cigar_str, algined = "";
	id_gen = mm_list[start_i1].id_gen;
	reads->target_genome_id = id_gen;
	gen_pos_end = mm_list[end_i1].pos_gen;
	//seq_pos_end = mm_list[end_i].pos_seq;
	gen_pos_start = mm_list[start_i1].pos_gen + 1 - opt.k;
        //seq_pos_start = mm_list[start_i].pos_seq + 1 - opt.k;
	//len_seq_for_align = seq_pos_end - seq_pos_start + 1;
	len_gen_for_align = gen_pos_end - gen_pos_start + 1;
	if (strand == 0){// plus direction
		seq_pos_end = mm_list[end_i1].pos_seq;
		seq_pos_start = mm_list[start_i1].pos_seq + 1 - opt.k;
		len_seq_for_align = seq_pos_end - seq_pos_start + 1;
		if (len_seq_for_align <= len_gen_for_align){
			result = edlibAlign(reads->seq_is + seq_pos_start, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_start, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
			cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			//////if (*result.startLocations != 0){
                       	//////	algined += to_string(*result.startLocations);
                       	//////	algined += "D";
               		//////}
			cigar_str = cigar;
			algined += cigar_str;
               		printf("%s\n", cigar);
               		printAlignment(reads->seq_is + seq_pos_start, gen_list[id_gen].seq_is + gen_pos_start, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_HW);
               		//////if (*result.endLocations != len_gen_for_align - 1){
                       	//////	algined += to_string(len_gen_for_align - 1 - *result.endLocations);
                       	//////	algined += "D";
               		//////}
		}
		else{
			result = edlibAlign(gen_list[id_gen].seq_is + gen_pos_start, len_gen_for_align, reads->seq_is + seq_pos_start, len_seq_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
			cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			printf("%s\n", cigar);
                        printAlignment(gen_list[id_gen].seq_is + gen_pos_start, reads->seq_is + seq_pos_start, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_HW);
			cigar_str = cigar;
		}
		// align the two ends region
		// first the right end
		if (seq_pos_end + 1 < reads->len) {
			len_seq_for_align = reads->len - seq_pos_end - 1;
			len_gen_for_align = 2 * len_seq_for_align;
			if (gen_pos_end + 1 + len_gen_for_align >= gen_list[id_gen].len){
				len_gen_for_align = gen_list[id_gen].len - 1 - (gen_pos_end + 1);
			}
			result = edlibAlign(reads->seq_is + seq_pos_end + 1, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_end + 1, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			cigar_str = cigar;
                        algined += cigar_str;
			printf("%s\n", cigar);
                        printAlignment(reads->seq_is + seq_pos_end + 1, gen_list[id_gen].seq_is + gen_pos_end + 1, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);

		}
		// then the left end region
		if (seq_pos_start != 0){
			len_seq_for_align = seq_pos_start + 1 - opt.k;
			len_gen_for_align = 2 * len_seq_for_align;
			char * start_seq_char, * start_gen_char;
		       	start_seq_char = (char *)malloc(len_seq_for_align * sizeof(char));
			start_gen_char = (char *)malloc(len_gen_for_align * sizeof(char));
			for(i = 0; i < len_seq_for_align; ++i){
				start_seq_char[i] = reads->seq_is[len_seq_for_align - i - 1];
				//start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
			}
			if (gen_pos_start  < len_gen_for_align){
				len_gen_for_align = gen_pos_start;
				//for (i = 0; i < len_gen_for_align; ++i){
				//	start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
				//}
			}
			for(i = 0; i < len_gen_for_align; ++i){
				start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
			}
			result = edlibAlign(start_seq_char, len_seq_for_align, start_gen_char, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
			reads->aligned_pos = gen_pos_start - *result.endLocations;
			cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			cigar_str = cigar;
			printf("%s\n", cigar);
                        printAlignment(start_seq_char, start_gen_char, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);
			numm = "";
			temp_str = "";
			cigar_reverse = "";
			for(ii = cigar_str.length() - 1; ii >= 0; --ii){
				if (cigar_str[ii] == 'I'){

					cigar_reverse += numm + temp_str;
					temp_str = "I";
					numm = "";
				}
				else if (cigar_str[ii] == 'D'){
					cigar_reverse += numm + temp_str;
					temp_str = "D";
					numm = "";
				}
				else if (cigar_str[ii] == 'M'){
					cigar_reverse += numm + temp_str;
					temp_str = "M";
					numm = "";
				}
				else{
					numm = cigar_str[ii] + numm;
					//cigar_reverse = cigar_str[i] + cigar_reverse;
				}
				
			}

			cigar_reverse += numm + temp_str;
			algined = cigar_reverse + algined;
			//reads->aligned_pos = gen_pos_start - *result.endLocations - 1;
			free(start_seq_char);
			free(start_gen_char);
		}
		else{
			reads->aligned_pos = gen_pos_start + 1;
			
		}
	}
	else{ // minus direction
		//////FILE *fr;
	        //////fr=fopen(opt.output,"w");
        	//////if (fr == NULL){
                //////	printf("File not open!/n");
        	//////}
        	//////else{
                //////	for (ii = start_i; ii <= end_i; ++ii){
                //////        	fprintf(fr, "%d\t%d\n", mm_list[ii].pos_gen, mm_list[ii].pos_seq);
                //////	}

        	//////}
		//////fclose(fr);

		reads->seq_reverse = (char *) malloc(sizeof(char) * reads->len + 1);
		for (int m = 0; m < reads->len; ++m){
			switch(reads->seq_is[reads->len - 1 - m]){
				case 'A':
					reads->seq_reverse[m] = 'T';
					break;
				case 'C':
					reads->seq_reverse[m] = 'G';
                                        break;
                                case 'G':
					reads->seq_reverse[m] = 'C';
                                        break;
                                case 'T':
					reads->seq_reverse[m] = 'A';
                                        break;
                                default:
		 			reads->seq_reverse[m] = reads->seq_is[reads->len - 1 - m];
			}
                }
		reads->reverse_calculated = 1;
		seq_pos_start = reads->len - mm_list[start_i1].pos_seq - 1;
                seq_pos_end = reads->len - 1 - (mm_list[end_i1].pos_seq + 1 - opt.k);
                len_seq_for_align = seq_pos_end - seq_pos_start + 1;
		if (len_seq_for_align <= len_gen_for_align){
                        result = edlibAlign(reads->seq_reverse + seq_pos_start, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_start, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        //////if (*result.startLocations != 0){
                        //////  algined += to_string(*result.startLocations);
                        //////  algined += "D";
                        //////}
                        cigar_str = cigar;
                        algined += cigar_str;
                        //printf("%s\n", cigar);
                        //printAlignment(reads->seq_reverse + seq_pos_start, gen_list[id_gen].seq_is + gen_pos_start, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_HW);
                        //////if (*result.endLocations != len_gen_for_align - 1){
                        //////  algined += to_string(len_gen_for_align - 1 - *result.endLocations);
                        //////  algined += "D";
                        //////}
                }
		// align the two ends region
                // first the right end
                if (seq_pos_end + 1 < reads->len) {
                        len_seq_for_align = reads->len - seq_pos_end - 1;
                        len_gen_for_align = 2 * len_seq_for_align;
                        if (gen_pos_end + 1 + len_gen_for_align >= gen_list[id_gen].len){
                                len_gen_for_align = gen_list[id_gen].len - 1 - (gen_pos_end + 1);
                        }
                        result = edlibAlign(reads->seq_reverse + seq_pos_end + 1, len_seq_for_align, gen_list[id_gen].seq_is + gen_pos_end + 1, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        cigar_str = cigar;
                        algined += cigar_str;
                        //printf("%s\n", cigar);
                        //printAlignment(reads->seq_reverse + seq_pos_end + 1, gen_list[id_gen].seq_is + gen_pos_end + 1, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);

                }

		// then the left end region
                if (seq_pos_start != 0){
                        len_seq_for_align = seq_pos_start;
                        len_gen_for_align = 2 * len_seq_for_align;
                        char * start_seq_char, * start_gen_char;
                        start_seq_char = (char *)malloc(len_seq_for_align * sizeof(char));
                        start_gen_char = (char *)malloc(len_gen_for_align * sizeof(char));
                        for(i = 0; i < len_seq_for_align; ++i){
                                start_seq_char[i] = reads->seq_reverse[len_seq_for_align - i - 1];
                                //start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
                        }
                        if (gen_pos_start  < len_gen_for_align){
                                len_gen_for_align = gen_pos_start;
                                //for (i = 0; i < len_gen_for_align; ++i){
                                //      start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
                                //}
                        }
                        for(i = 0; i < len_gen_for_align; ++i){
                                start_gen_char[i] = gen_list[id_gen].seq_is[gen_pos_start - i - 1];
                        }
			result = edlibAlign(start_seq_char, len_seq_for_align, start_gen_char, len_gen_for_align, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0));
                        cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        cigar_str = cigar;
                        //printf("%s\n", cigar);
                        //printAlignment(start_seq_char, start_gen_char, result.alignment, result.alignmentLength, *result.endLocations, EDLIB_MODE_SHW);
			numm = "";
			temp_str = "";
                        cigar_reverse = "";
                        for(ii = cigar_str.length() - 1; ii >= 0; --ii){
                                if (cigar_str[ii] == 'I'){

                                        cigar_reverse += numm + temp_str;
                                        temp_str = "I";
                                        numm = "";
                                }
                                else if (cigar_str[ii] == 'D'){
                                        cigar_reverse += numm + temp_str;
                                        temp_str = "D";
                                        numm = "";
                                }
                                else if (cigar_str[ii] == 'M'){
                                        cigar_reverse += numm + temp_str;
                                        temp_str = "M";
                                        numm = "";
                                }
                                else{
                                        numm = cigar_str[ii] + numm;
                                        //cigar_reverse = cigar_str[i] + cigar_reverse;
                                }

                        }
                        cigar_reverse += numm + temp_str;
                        algined = cigar_reverse + algined;
                        free(start_seq_char);
			free(start_gen_char);
			reads->aligned_pos = gen_pos_start - *result.endLocations;
		}
		else{
                        reads->aligned_pos = gen_pos_start + 1;
                }		
		//free(reads->seq_reverse);
		
	}
	//////delete[] chain1;
	//////delete[] chain2;
	//////FILE *fr;
	/////fr=fopen(opt.output,"w");
	//////if (fr == NULL){
	//////	printf("File not open!/n");
	//////}
	//////else{
	//////	for (ii = start_i; ii <= end_i; ++ii){
        //////                fprintf(fr, "%d\t%d\n", mm_list[ii].pos_gen, mm_list[ii].pos_seq);
        //////        }

	//////}
	reads->cigar = new char[algined.size() + 1];
	memcpy(reads->cigar, algined.c_str(), algined.size());
	//reads->cigar = aligned;
	edlibFreeAlignResult(result);
	free(cigar);
	//free(reads->seq_is);
	free(reads->kmerlist);
	//fclose(fr);
	return 1;
}
int sequence_map_write(option opt, seqs *reads){
	ofstream outfile2;
	const char* filename = opt.output;
	outfile2.open(filename, ios::out);
	if (!outfile2.is_open()){
                bomb_error("Open output file failure");
	}
	outfile2 << "Query:        " << reads->header;
	outfile2 << "strand:        " << reads->reverse_calculated;
	outfile2 << "genome_id:        " << reads->target_genome_id;
	outfile2 << "cigar: " << reads->cigar;
	outfile2.close();
	return 1;
}
void printAlignment(const char* query, const char* target, const unsigned char* alignment, const int alignmentLength, const int position, const EdlibAlignMode modeCode) {
        int tIdx = -1;
        int qIdx = -1;
        if (modeCode == EDLIB_MODE_HW) {
                tIdx = position;
                for (int i = 0; i < alignmentLength; i++) {
                        if (alignment[i] != EDLIB_EDOP_INSERT)
                                tIdx--;
                }
        }
	for (int start = 0; start < alignmentLength; start += 50) {
                // target
                printf("T: ");
                int startTIdx;
                for (int j = start; j < start + 50 && j < alignmentLength; j++) {
                        if (alignment[j] == EDLIB_EDOP_INSERT)
                                printf("-");
                        else
                                printf("%c", target[++tIdx]);
                        if (j == start)
                                startTIdx = tIdx;
                }
                printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);

                // match / mismatch
                printf("   ");
                for (int j = start; j < start + 50 && j < alignmentLength; j++) {
                        printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
                }
                printf("\n");

                // query
                printf("Q: ");
                int startQIdx = qIdx;
		for (int j = start; j < start + 50 && j < alignmentLength; j++) {
                        if (alignment[j] == EDLIB_EDOP_DELETE)
                                printf("-");
                        else
                                printf("%c", query[++qIdx]);
                        if (j == start)
                                startQIdx = qIdx;
                }
                printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
        }
}
}

void bomb_error(const char *message)
{
        fprintf(stderr, "\nFatal Error:\n%s\nProgram halted !!\n\n", message);
        //temp_files.Clear();
        exit(1);
} // END void bomb_error
