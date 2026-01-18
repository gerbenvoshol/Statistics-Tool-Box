/* spearman.c - Spearman's Rank Correlation Calculator
 * 
 * A sample C program for calculating Spearman's rank correlation coefficient
 * between two data files using stb_stats.h.
 *
 * This program performs:
 * 1. Reads tab-separated data files (e.g., HTSeq count files)
 * 2. Extracts count values from the third column (gene_ID, gene_name, counts format)
 * 3. Optionally filters data based on minimum value threshold
 * 4. Calculates Spearman's rank correlation coefficient
 *
 * Input format: TAB-delimited files with gene_ID, gene_name, and counts
 * Expected format: gene_ID\tgene_name\tcount_value
 *
 * Usage: spearman [options] file1 file2
 */

#define STB_STATS_DEFINE
#include "stb_stats.h"
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_usage(const char *prog_name) {
	fprintf(stderr, "Usage: %s [options] file1 file2\n\n", prog_name);
	fprintf(stderr, "Spearman's Rank Correlation Calculator\n");
	fprintf(stderr, "Calculates Spearman's rank correlation coefficient between two data files.\n\n");
	fprintf(stderr, "Arguments:\n");
	fprintf(stderr, "  file1, file2            Input TAB-delimited files (e.g., HTSeq count files)\n");
	fprintf(stderr, "                          Expected format: gene_ID\\tgene_name\\tcount_value\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -h, --help              Show this help message\n");
	fprintf(stderr, "  -s, --skip-header       Skip the first line (header) in both files\n");
	fprintf(stderr, "  -m, --min-val VALUE     Minimum value threshold for filtering (keep pairs where\n");
	fprintf(stderr, "                          at least one value meets the threshold)\n\n");
	fprintf(stderr, "Output:\n");
	fprintf(stderr, "  file1\\tfile2\\tcorrelation_coefficient\n\n");
	fprintf(stderr, "Example:\n");
	fprintf(stderr, "  %s sample1.txt sample2.txt -s -m 10\n", prog_name);
}

int main(int argc, char *argv[])
{
	int skip_header = 0;
	int use_min_val = 0;
	double min_val = 0.0;
	
	static struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"skip-header", no_argument, 0, 's'},
		{"min-val", required_argument, 0, 'm'},
		{0, 0, 0, 0}
	};
	
	int opt;
	while ((opt = getopt_long(argc, argv, "hsm:", long_options, NULL)) != -1) {
		switch (opt) {
			case 'h':
				print_usage(argv[0]);
				return 0;
			case 's':
				skip_header = 1;
				break;
			case 'm':
				min_val = atof(optarg);
				use_min_val = 1;
				break;
			default:
				print_usage(argv[0]);
				return 1;
		}
	}
	
	if (optind + 2 != argc) {
		fprintf(stderr, "Error: Two input files required\n\n");
		print_usage(argv[0]);
		return 1;
	}
	
	const char *file1 = argv[optind];
	const char *file2 = argv[optind + 1];

	/* Reads an entire file into an array of strings, needs only a single call to free */
	size_t n1 = 0;
	int nrfields = 0;
	char **fdata1 = stb_fgetlns((char *)file1, &n1);
	if (!fdata1) {
		fprintf(stderr, "Failed to read file: %s\n", file1);
		return 1;
	}
	double *data1 = calloc(n1, sizeof(double));
	if (!data1) {
		fprintf(stderr, "Memory allocation failed\n");
		free(fdata1);
		return 1;
	}
	int valid1 = 0;
	size_t start_idx = skip_header ? 1 : 0;
	for (size_t i = start_idx; i < n1; i++) {
		if (!fdata1[i] || fdata1[i][0] == '\0') {
			continue;
		}
		nrfields = 0;
		char **line = stb_parse_csv(fdata1[i], '\t', &nrfields);
		if (!nrfields || nrfields == -1 || nrfields < 3) {
			if (line) stb_free_csv_line(line);
			continue;
		}
		data1[valid1] = atof(line[2]);
		valid1++;
		stb_free_csv_line(line);
	}
	free(fdata1);

	size_t n2 = 0;
	char **fdata2 = stb_fgetlns((char *)file2, &n2);
	if (!fdata2) {
		fprintf(stderr, "Failed to read file: %s\n", file2);
		free(data1);
		return 1;
	}
	double *data2 = calloc(n2, sizeof(double));
	if (!data2) {
		fprintf(stderr, "Memory allocation failed\n");
		free(fdata2);
		free(data1);
		return 1;
	}
	int valid2 = 0;
	for (size_t i = start_idx; i < n2; i++) {
		if (!fdata2[i] || fdata2[i][0] == '\0') {
			continue;
		}
		nrfields = 0;
		char **line = stb_parse_csv(fdata2[i], '\t', &nrfields);
		if (!nrfields || nrfields == -1 || nrfields < 3) {
			if (line) stb_free_csv_line(line);
			continue;
		}
		data2[valid2] = atof(line[2]);
		valid2++;
		stb_free_csv_line(line);
	}
	free(fdata2);
	
	// Use minimum length to avoid out-of-bounds access
	size_t n = (valid1 < valid2) ? valid1 : valid2;
	
	// Filter data if min_val is provided
	double *f1, *f2;
	size_t filtered_n;
	
	if (use_min_val) {
		// Allocate filtered arrays
		f1 = malloc(n * sizeof(double));
		f2 = malloc(n * sizeof(double));
		if (!f1 || !f2) {
			fprintf(stderr, "Memory allocation failed\n");
			if (f1) free(f1);
			if (f2) free(f2);
			free(data1);
			free(data2);
			return 1;
		}
		
		filtered_n = 0;
		for (size_t i = 0; i < n; i++) {
			// Keep row if at least one value meets the minimum threshold
			if (data1[i] >= min_val || data2[i] >= min_val) {
				f1[filtered_n] = data1[i];
				f2[filtered_n] = data2[i];
				filtered_n++;
			}
		}
		
		if (filtered_n == 0) {
			fprintf(stderr, "No data pairs meet the minimum value criterion\n");
			free(f1);
			free(f2);
			free(data1);
			free(data2);
			return 1;
		}
		
		double spearman = stb_spearman(f1, f2, filtered_n);
		printf("%s\t%s\t%lf\n", file1, file2, spearman);
		
		free(f1);
		free(f2);
	} else {
		double spearman = stb_spearman(data1, data2, n);
		printf("%s\t%s\t%lf\n", file1, file2, spearman);
	}

	free(data1);
	free(data2);

	return 0;
}
