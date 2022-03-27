#include <zlib.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include "common.h"
#include "kseq.h"
#include "kvec.h"
// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)  
	
int main(int argc, char *argv[]){  
	int j;
	unsigned int N;
	gzFile fp, fp2;  
    	kseq_t *seq, *seq2;
    	unsigned int seqs_num = 0, total_num = 0;
    	long bases = 0;  
	int l;
    	unsigned int i, len, chunk = 50000, start_id = 0;
	clock_t t;
	option opt;
	l = set_options(argc, argv, &opt);
	if (l == 0){
		return 0;
	}
	ofstream outfile2;
        const char* filename = opt.output;
        outfile2.open(filename, ios::out);
        if (!outfile2.is_open()){
                bomb_error("Open output file failure");
        }

	//t = clock();
	time_t now = time(0);
	char* dt = ctime(&now);
	cout << "\n\n[INFO] Program started at: " << dt << endl;
	
	// ========================================================================
	// read position file
	string line;
	const char* posfilename = opt.pos_file;
	short int *array_gen_id = new short int[opt.seq_num];
	//bool *array_plus = new bool[opt.seq_num];
	bool array_plus[opt.seq_num] = { 0 };
	unsigned int *array_seq_start = new unsigned int [opt.seq_num];
	unsigned int *array_seq_end   = new unsigned int [opt.seq_num];
	unsigned int *array_gen_start = new unsigned int [opt.seq_num];
	unsigned int *array_gen_end   = new unsigned int [opt.seq_num];
	ifstream libfile(posfilename);
	unsigned int index, genome_id, seq_start, seq_end, gen_start, gen_end;
	short int plus;
	if (!libfile.is_open())
		bomb_error("Open position file failure, exit!");
	while (getline(libfile, line)) {
		istringstream iss(line);
		iss >> index >> genome_id >> plus >> seq_start >> seq_end >> gen_start >> gen_end;
		array_gen_id[index]    = genome_id;
		array_plus[index]      = plus;
		array_seq_start[index] = seq_start;
		array_seq_end[index]   = seq_end;
		array_gen_start[index] = gen_start;
		array_gen_end[index]   = gen_end;
		total_num++;
		printf("\rReading the %10d -th position...", total_num);
                fflush( stdout );

	}
	libfile.close();
	now = time(0);
        dt = ctime(&now);
	cout << "\n";
        cout << "[INFO] Position file reading finished at:         " << dt;
	cout << "[INFO] Reading the genome file...         " << endl;
	//opt.seq_num = total_num;
	if (opt.seq_num > chunk){
		N = chunk;
	}
	else {
		N = opt.seq_num;
	}
	seqs  *array = (seqs*)malloc(200 * sizeof(seqs));
    	fp = gzopen(opt.input, "r"); // STEP 2: open the file handler  
    	seq = kseq_init(fp); // STEP 3: initialize seq 
    	while ((l = kseq_read(seq)) >= 0){ // STEP 4: read sequence  
		//printf("\rThe %10ld -th seq len: %d", seqs_num, seq->seq.l);
		//fflush( stdout );
		bases += seq->seq.l;
        	//seqs_num += 1;
		array[seqs_num].header = (char *) malloc(strlen(seq->name.s) + 1);
		strcpy(array[seqs_num].header, seq->name.s);
		array[seqs_num].id  = seqs_num;
		array[seqs_num].len = seq->seq.l;
		//array[seqs_num].is_seed = 0;
		array[seqs_num].reverse_calculated = 0;
		//array[seqs_num].assigned = 0;
		if (seq->seq.l > 0){
			array[seqs_num].seq_is = (char *) malloc(sizeof(char) * seq->seq.l + 1); 
			strcpy(array[seqs_num].seq_is, seq->seq.s);
		}
		seqs_num += 1;
    	}  
    	printf("[INFO] Total genomes: %d\n", seqs_num);
    	printf("[INFO] Total genome bases: %ld\n", bases);
	now = time(0);
        dt = ctime(&now);
        cout << "[INFO] Genome reading finished at:         " << dt;	
	outfile2 << "@HD        VN:1.5  SO:unsorted\n";
        for (i = 0; i < seqs_num; ++i) {
                outfile2 << "@SQ        SN:" << array[i].header << "\tLN:" << array[i].len << "\n";
        }
	outfile2 << "@PG        ID:mapper	PN:mapper	VN:0.0.1	CL:" << *argv << "\n";
	//printf("[INFO] Reading genome time: %.3f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);	
    	kseq_destroy(seq); // STEP 5: destroy seq  
    	gzclose(fp); // STEP 6: close the file handler  
	//now = time(0);
        //dt = ctime(&now);
        //cout << "[INFO] Genome k-mer extracted finished at: " << dt << endl;
	//printf("[INFO] Extracting k-mer time: %.3f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	//now = time(0);
        //dt = ctime(&now);
	//cout << "[INFO] hash table point finished at: " << dt << endl;	
	// read the sequence file, or you can read one chunk, such as 50000 sequences, for large data
        bases = 0;
	seqs_num = 0;
	total_num = 0;
	fp2 = gzopen(opt.reads_file, "r"); // STEP 2: open the file handler  
        seq2 = kseq_init(fp2); // STEP 3: initialize seq
	seqs  *reads = (seqs*)malloc(chunk * sizeof(seqs));	
	while ((l = kseq_read(seq2)) >= 0){ // STEP 4: read sequence  
                //printf("\rThe %10ld -th seq len: %d", seqs_num, seq->seq.l);
                //fflush( stdout );
                bases += seq2->seq.l;
                //seqs_num += 1;
                reads[seqs_num].header = (char *) malloc(strlen(seq2->name.s) + 1);
                strcpy(reads[seqs_num].header, seq2->name.s);
                reads[seqs_num].id  = total_num;
                reads[seqs_num].len = seq2->seq.l;
                reads[seqs_num].reverse_calculated = 0;
		reads[seqs_num].target_genome_id   = array_gen_id[seqs_num];
                //reads[seqs_num].strand             = array_plus[seqs_num];
                //reads[seqs_num].seq_pos_start      = array_seq_start[seqs_num];
                //reads[seqs_num].seq_pos_end        = array_seq_end[seqs_num];
                //reads[seqs_num].gen_pos_start      = array_gen_start[seqs_num];
                //reads[seqs_num].gen_pos_end        = array_gen_end[seqs_num];
                //reads[seqs_num].assigned = 0;
                if (seq2->seq.l > 0){
                        reads[seqs_num].seq_is = (char *) malloc(sizeof(char) * seq2->seq.l + 1);
                        strcpy(reads[seqs_num].seq_is, seq2->seq.s);
                }
                seqs_num += 1;
		total_num += 1;
		printf("\rThe %10d -th seq len: %d", total_num, seq2->seq.l);
                fflush( stdout );
		if (seqs_num >= chunk){
			cout << "\n" << endl;
			#pragma omp parallel for
			for (i = 0; i < seqs_num; ++i){
				//if (reads[i].id < 14859)
				//	continue;
				sequence_map_one_with_locate(opt, &reads[i], array, array_gen_id, array_plus, array_seq_start, array_seq_end, array_gen_start, array_gen_end);
				printf("\rThe %10d -th seq mapped", reads[i].id);
				fflush( stdout );

			}
			printf("\n");
			// write to output file
			for (i = 0; i < seqs_num; ++i){
				//if (i == 39){
				//	int ff = 0;
				//}
				outfile2 << reads[i].header;
				if (reads[i].reverse_calculated){
                                	outfile2 << "\t16";
                        	}
                        	else{
                                	outfile2 << "\t0";
                        	}
        			//outfile2 << "\t" << reads[i].reverse_calculated ? "16" : "0";
        			outfile2 << "\t" << array[reads[i].target_genome_id].header;
        			outfile2 << "\t" << reads[i].aligned_pos;
				outfile2 << "\t60";
				outfile2 << "\t" << reads[i].cigar;
				outfile2 << "\t*";
			        outfile2 << "\t0";
				outfile2 << "\t0";
				if (reads[i].reverse_calculated){
                                	outfile2 << "\t" << reads[i].seq_reverse;
                                	free(reads[i].seq_reverse);
                                	reads[i].reverse_calculated = 0;
                        	}
                        	else{
                                	outfile2 << "\t" << reads[i].seq_is;
                        	}

				//outfile2 << "\t" << reads[i].seq_is;
				outfile2 << "\t*" << endl;
                                free(reads[i].header);
				free(reads[i].seq_is);
				free(reads[i].cigar);
				if (reads[i].reverse_calculated){
					free(reads[i].seq_reverse);
				}
				else{
					reads[i].reverse_calculated = 0;
				}


      			}
			seqs_num = 0;
		}
		
	}

	if (seqs_num > 0){
		printf("\n");
		#pragma omp parallel for
		for (i = 0; i < seqs_num; ++i){
                	//sequence_map_one(opt, seeds_kmer_list, point_list, &reads[i], array);
			sequence_map_one_with_locate(opt, &reads[i], array, array_gen_id, array_plus, array_seq_start, array_seq_end, array_gen_start, array_gen_end);
			printf("\rThe %10d -th seq mapped", reads[i].id);
                        fflush( stdout );
                }
		// write to output file
               	for (i = 0; i < seqs_num; ++i){
                	outfile2 << reads[i].header;
			if (reads[i].reverse_calculated){
				outfile2 << "\t16";
			}
			else{
				outfile2 << "\t0";
			}
                        //outfile2 << "\t" << reads[i].reverse_calculated ? "16" : "0";
                        outfile2 << "\t" << array[reads[i].target_genome_id].header;
                        outfile2 << "\t" << reads[i].aligned_pos;
                        outfile2 << "\t60";
                        outfile2 << "\t" << reads[i].cigar;
                        outfile2 << "\t*";
                        outfile2 << "\t0";
                        outfile2 << "\t0";
			if (reads[i].reverse_calculated){
				outfile2 << "\t" << reads[i].seq_reverse;
				free(reads[i].seq_reverse);
				reads[i].reverse_calculated = 0;
			}
			else{
				outfile2 << "\t" << reads[i].seq_is;
			}
                        outfile2 << "\t*" << endl;
                        free(reads[i].header);
                        free(reads[i].seq_is);
                        free(reads[i].cigar);
                        free(reads[i].kmerlist);
                        //if (reads[i].reverse_calculated){
			//	free(reads[i].seq_reverse);
			//}
			//else{
			//	reads[i].reverse_calculated = 0;
			//}
		}

	}
	outfile2.close();
	//free(array_plus);
        delete [] array_seq_start;
        delete [] array_seq_end;
     	delete [] array_gen_start;
	delete [] array_gen_end;
	printf("[INFO] Total sequences: %d\n", total_num);
        printf("[INFO] Total sequence bases: %ld\n", bases);
        now = time(0);
        dt = ctime(&now);
        cout << "[INFO] Sequence file reading finished at:         " << dt;

	
	now = time(0);
        dt = ctime(&now);
        cout << "[INFO] Program completed at: " << dt << endl;
	return 0;
}
}
