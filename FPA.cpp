// Updated on 11/17/21

/* 

Program FPA.cpp to identify private alleles from high-throughput sequencing data of diploid individuals 
from multiple populations. The allele frequencies necessary for the analysis are estimated beforehand by 
GFE in the p mode.  The private allele is found only when a site is polymorphic in the total population.  
Alleles are examined only from populations with ML estimates and effective number of sampled chromosomes 
equal to or greater than twenty.  Statistical significance of the polymorphism in a deme is examined when 
deciding to add a new allele.  The input file contains error-rate estimates.  The total coverage is simply 
the sum of population coverage over the populations.  The probability of finding the private allele is 
reported in the output.  The frequency of the private allele in the focal population and that in the total 
population are reported in the output.  The latter is calculated as the arithmetic mean across populations.

Input: Combined GFE output files in the p mode with per-site annotations.     

*/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
using namespace std;

int main(int argc, char *argv[])
{
	// Default values of the options
	const char* in_file_name = {"In_FPA.txt"};
	const char* out_file_name = {"Out_FPA.txt"};
	double min_Nc = 20.0;
	double cv = 5.991;
	int print_help = 0;
	
	int argz = 1; // argument counter

	// Read specified settings
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-in") == 0) {
			in_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-out") == 0) {
			out_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-min_Nc") == 0) {
			sscanf(argv[++argz], "%lf", &min_Nc);
		} else if (strcmp(argv[argz], "-cv") == 0) {
			sscanf(argv[++argz], "%lf", &cv);
		} else {
			fprintf(stderr, "unknown option %s\n", argv[argz]);
			print_help = 1;
			break;
		}
		argz++;
	}
	if (print_help) { // print error/usage message ?
		fprintf(stderr, "USAGE: %s {<options>}\n", argv[0]);
		fprintf(stderr, "	options:\n");
		fprintf(stderr, "	-h: print the usage message\n");
		fprintf(stderr, "	-in <s>: specify the input file name\n");
		fprintf(stderr, "       -out <s>: specify the output file name\n");
		fprintf(stderr, "	-min_Nc <f>: specify the minimum effective number of sampled chromosomes required in a deme\n");
		fprintf(stderr, "       -cv <f>: specify the chi-square critical value for the polymorphism test\n");
		exit(1);
	}

	string line; // String buffer
	
	ifstream inputFile(in_file_name); // Try to open the input file
	if ( !inputFile.is_open() ) { // Exit on failure
		fprintf(stderr, "Cannot open %s for reading.\n", in_file_name);
		exit(1);
	}
	
	// Read the header
	string h_scaf, h_site, h_ref_nuc;
	vector <string> pop_info; // Stores population-specific labels.
	getline(inputFile, line);
	istringstream ss(line);
	ss >> h_scaf >> h_site >> h_ref_nuc;
	string str; // Temporarily stores population-specific labels
	pop_info.clear();
	while (true) {
		ss >> str;
		pop_info.push_back(str);
		if ( ss.eof() ) {
			break;
		}
	}
	int num_pops = (int)pop_info.size()/9;
	printf("%d populations to be analyzed\n", num_pops);

	FILE *outstream;

	// Open the output file
	outstream = fopen(out_file_name, "w");
	if (outstream == NULL ) { // Exit on failure
		fprintf(stderr, "Cannot open %s for writing.\n", out_file_name); 
		exit(1);
	}
	
	// Print out the field names
	fprintf(outstream, "scaffold\tsite\tref_nuc\ttot_cov\tne_pops\tnum_alleles\tprivate_allele\tid_pop\tfocal_frequency\ttotal_frequency\tlog_prob_pa\tMAF\n");
	// printf("scaffold\tsite\tref_nuc\ttot_cov\tne_pops\tnum_alleles\tprivate_allele\tid_pop\tfocal_frequency\ttotal_frequency\tlog_prob_pa\tMAF\n");
	
	// Read the main data
	string scaffold, ref_nuc, n1[num_pops+1], n2[num_pops+1], s_Nc[num_pops+1], s_best_p[num_pops+1], s_best_q[num_pops+1], best_error, s_best_H[num_pops+1], s_pol_llstat[num_pops+1];
	double Nc[num_pops+1], pol_llstat[num_pops+1]; 
	int site, pop_cov[num_pops+1], num_alleles;
	int tot_cov;		// total coverage (sum of the coverage across the populations)
	int ne_pops;            // total number of populations with data
	double sum_Nc;		// Sum of the effective number of sampled chromosomes over populations with data
	vector <string> alleles;   // store allele identities
	int ag;		// allele counter
	int pg;		// population counter
	vector <int> id_pop_a;	// id of the population with an allele
	vector <double> freq_a;		// frequency of the allele in the population
	int num_pops_a;			// number of populations that have an allele
	double sum_freq_a;		// sum of the frequencies of the allele over populations with data
	double mean_freq_a;           // mean of the frequencies of the allele over populations with data
	double maf_total;		// minor allele frequency in the total population
	double Nc_focal;		// effective number of sampled chromosomes in the focal population
	double Nc_other;		// effective number of sampled chromosomes in the other populations
	vector <string> private_allele;	// private allele
	vector <int> id_pop_pa;		// id of the population with the private allele
	vector <double> focal_paf;	// private-allele frequency in the focal population
	vector <double> total_paf;	// private-allele frequency in the total population
	double t_prob_pa;		// temporarily stores the probability of the private allele
	vector <double> log_prob_pa;	// logarithm of the probability of the private allele
	int num_pa;			// number of private alleles

	while ( getline(inputFile, line) ) {
		istringstream ss(line);
		alleles.clear();
		tot_cov = 0;
		ne_pops = 0;
		sum_Nc = 0.0;
		ss >> scaffold >> site >> ref_nuc;
		for (pg = 1; pg <= num_pops; pg++) {
			ss >> n1[pg] >> n2[pg] >> pop_cov[pg] >> s_Nc[pg] >> s_best_p[pg] >> s_best_q[pg] >> best_error >> s_best_H[pg] >> s_pol_llstat[pg];
			tot_cov = tot_cov + pop_cov[pg];
			// printf("site: %d\tpop: %d\tn1: %s\tn2: %s\n", site, pg, n1[pg].c_str(), n2[pg].c_str());
			// fprintf(outstream, "site: %d\tpop: %d\tn1: %s\tn2: %s\n", site, pg, n1[pg].c_str(), n2[pg].c_str());
			if (n1[pg] != "NA") {
				Nc[pg] = atof(s_Nc[pg].c_str());
				if (Nc[pg] >= min_Nc) {
					ne_pops = ne_pops + 1;
					sum_Nc = sum_Nc + Nc[pg];
                                	if ( find(alleles.begin(), alleles.end(), n1[pg]) == alleles.end() ) {
						// printf("%s\n", n1[pg].c_str());
                                        	alleles.push_back(n1[pg]);
                                	}
                        		if (n2[pg] != "NA") {
						pol_llstat[pg] = atof(s_pol_llstat[pg].c_str());
						if (pol_llstat[pg] > cv) {
                                			if ( find(alleles.begin(), alleles.end(), n2[pg]) == alleles.end() ) {
								// printf("%s\n", n2[pg].c_str());
                                        			alleles.push_back(n2[pg]);
							}
                                		}
					}
				}
                        }
		}
		// Count the number of alleles segregating in the population sample
		num_alleles = alleles.size();
		/*
		printf("site: %d\tnum_alleles: %d\n", site, num_alleles);
		fprintf(outstream, "site: %d\tnum_alleles: %d\n", site, num_alleles);
		for (ag=0; ag<num_alleles; ag++) {
			printf("%s\n", alleles.at(ag).c_str());
		}
		*/
		// clear the vectors for private alleles
		private_allele.clear();
		id_pop_pa.clear();
		focal_paf.clear();
		total_paf.clear();
		log_prob_pa.clear();
		for (ag = 0; ag < num_alleles; ag++) {		// Examine each of the alleles
			sum_freq_a = 0.0;
			// Clear the vectors on alleles
			id_pop_a.clear();
			freq_a.clear();
			for (pg = 1; pg <= num_pops; pg++) { // Examine the allele over the populations
				if (n1[pg] != "NA" && Nc[pg] >= min_Nc) {	// Examine the population only when there are ML estimates with Nc equal to or greater than the specified value at the site
					if ( n1[pg] == alleles.at(ag) ) {
						id_pop_a.push_back(pg);
						freq_a.push_back( atof(s_best_p[pg].c_str()) );
						sum_freq_a = sum_freq_a + atof(s_best_p[pg].c_str());
					} else if ( n2[pg] == alleles.at(ag) ) {
						id_pop_a.push_back(pg);
						freq_a.push_back( atof(s_best_q[pg].c_str()) );
						sum_freq_a = sum_freq_a + atof(s_best_q[pg].c_str());
					}
				}
			}
			num_pops_a = id_pop_a.size();
			mean_freq_a = sum_freq_a/ne_pops;
			if (ag == 0) {
				maf_total = mean_freq_a;
			} else {
				if (mean_freq_a < maf_total) {
					maf_total = mean_freq_a;
				}
			}
			if (ne_pops >= 2 && num_pops_a == 1) {		// private allele
				private_allele.push_back(alleles.at(ag).c_str());
				id_pop_pa.push_back(id_pop_a.at(0));
				focal_paf.push_back(freq_a.at(0));	
				total_paf.push_back(mean_freq_a);
				Nc_focal = Nc[id_pop_a.at(0)];
				Nc_other = sum_Nc - Nc_focal;
				t_prob_pa = ( 1.0-pow(1.0-mean_freq_a,Nc_focal) )*pow(1.0-mean_freq_a,Nc_other);
				log_prob_pa.push_back( log10(t_prob_pa) );
			}
		}
		// print out the results
		num_pa = private_allele.size();
		if (num_pa >= 1) {
			for (ag = 0; ag < num_pa; ag++) {
				fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, ne_pops, num_alleles, private_allele.at(ag).c_str(), id_pop_pa.at(ag), focal_paf.at(ag), total_paf.at(ag), log_prob_pa.at(ag), maf_total);
				// printf("%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), tot_cov, ne_pops, num_alleles, private_allele.at(ag).c_str(), id_pop_pa.at(ag), focal_paf.at(ag), total_paf.at(ag), log_prob_pa.at(ag), maf_total);
			}
		}   						
	} 		

	return 0;
}
