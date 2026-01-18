/* deseq2_example.c - DESeq2-style Differential Expression Analysis
 * 
 * A sample C program demonstrating DESeq2-style differential expression analysis
 * using stb_stats.h for normalization and statistical testing.
 *
 * This program performs:
 * 1. RSE (Relative Log Expression) normalization using geometric means
 * 2. Dispersion estimation using stb_fit_f_dist
 * 3. Moderated t-test for differential expression
 * 4. Multiple testing correction using Benjamini-Hochberg FDR
 *
 * Usage: deseq2_example [counts_file] [options]
 */

#define STB_STATS_DEFINE
#include "stb_stats.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#define MAX_LINE_LENGTH 1048576
#define PSEUDOCOUNT 1.0

typedef struct {
    char *input_file;
    char *output_file;
    int group1_start;  /* Column index where group 1 starts */
    int group1_count;  /* Number of samples in group 1 */
    int group2_start;  /* Column index where group 2 starts */
    int group2_count;  /* Number of samples in group 2 */
    double fdr;        /* False discovery rate */
} Config;

/* Default configuration */
void init_config(Config *config) {
    config->input_file = NULL;
    config->output_file = NULL;
    config->group1_start = 1;  /* Skip first column (gene names) */
    config->group1_count = 3;
    config->group2_start = 4;
    config->group2_count = 3;
    config->fdr = 0.05;
}

void print_usage(const char *prog_name) {
    fprintf(stderr, "Usage: %s [input_file] [options]\n\n", prog_name);
    fprintf(stderr, "DESeq2-style Differential Expression Analysis\n");
    fprintf(stderr, "Performs normalization, dispersion estimation, and statistical testing.\n\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  input_file              Path to input count matrix (TAB-delimited)\n");
    fprintf(stderr, "                          Format: Gene\\tSample1\\tSample2\\t...\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o, --output FILE       Output filename (default: stdout)\n");
    fprintf(stderr, "  --g1-start N            Group 1 start column (default: 1)\n");
    fprintf(stderr, "  --g1-count N            Group 1 sample count (default: 3)\n");
    fprintf(stderr, "  --g2-start N            Group 2 start column (default: 4)\n");
    fprintf(stderr, "  --g2-count N            Group 2 sample count (default: 3)\n");
    fprintf(stderr, "  --fdr FLOAT             False discovery rate (default: 0.05)\n");
    fprintf(stderr, "  -h, --help              Show this help message\n\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "  %s counts.txt --g1-start 1 --g1-count 3 --g2-start 4 --g2-count 3 -o results.txt\n", prog_name);
}

/* Parse command line arguments */
int parse_args(int argc, char **argv, Config *config) {
    static struct option long_options[] = {
        {"output", required_argument, 0, 'o'},
        {"g1-start", required_argument, 0, '1'},
        {"g1-count", required_argument, 0, '2'},
        {"g2-start", required_argument, 0, '3'},
        {"g2-count", required_argument, 0, '4'},
        {"fdr", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;

    /* Check for help first */
    if (argc < 2) {
        return -1;
    }
    
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        return -1;
    }

    /* First non-option argument is input file */
    config->input_file = argv[1];

    /* Parse options */
    while ((opt = getopt_long(argc - 1, argv + 1, "o:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'o':
                config->output_file = optarg;
                break;
            case '1':
                config->group1_start = atoi(optarg);
                break;
            case '2':
                config->group1_count = atoi(optarg);
                break;
            case '3':
                config->group2_start = atoi(optarg);
                break;
            case '4':
                config->group2_count = atoi(optarg);
                break;
            case 'f':
                config->fdr = atof(optarg);
                break;
            case 'h':
                return -1;
            default:
                return -1;
        }
    }

    return 0;
}

/* Perform DESeq2-style differential expression analysis */
void perform_deseq2_analysis(STB_MAT *counts, Config *config, FILE *output) {
    int n_genes = counts->rows;
    int n_samples = counts->columns;
    int i, j;

    fprintf(stderr, "Processing %d genes across %d samples\n", n_genes, n_samples);
    fprintf(stderr, "Group 1: columns %d-%d (%d samples)\n", 
            config->group1_start, config->group1_start + config->group1_count - 1, config->group1_count);
    fprintf(stderr, "Group 2: columns %d-%d (%d samples)\n", 
            config->group2_start, config->group2_start + config->group2_count - 1, config->group2_count);

    /* Step 1: Calculate geometric scaling factors for normalization */
    fprintf(stderr, "Step 1: Calculating size factors (RSE normalization)...\n");
    double *scaling_factors;
    stb_calc_geometric_scaling_factors(counts, &scaling_factors);
    
    fprintf(stderr, "Size factors: ");
    for (i = 0; i < n_samples; i++) {
        fprintf(stderr, "%.4f ", scaling_factors[i]);
    }
    fprintf(stderr, "\n");

    /* Step 2: Normalize counts to common scale */
    fprintf(stderr, "Step 2: Normalizing counts...\n");
    double *common_means;
    double *common_vars;
    STB_MAT *normalized_counts = stb_dup_matrix(counts);
    stb_meanvar_counts_to_common_scale(normalized_counts, scaling_factors, &common_means, &common_vars);

    /* Step 3: Estimate dispersion parameters using variance data */
    fprintf(stderr, "Step 3: Estimating dispersion parameters...\n");
    double prior_variance;
    double prior_df;
    stb_fit_f_dist(common_vars, n_genes, n_samples - 2, &prior_variance, &prior_df);
    
    fprintf(stderr, "Prior variance: %.6f, Prior df: %.6f\n", prior_variance, prior_df);

    /* Step 4: Perform moderated t-test for each gene */
    fprintf(stderr, "Step 4: Performing statistical tests...\n");
    double *p_values = malloc(n_genes * sizeof(double));
    double *t_values = malloc(n_genes * sizeof(double));
    double *log2fc = malloc(n_genes * sizeof(double));
    double *baseMean = malloc(n_genes * sizeof(double));
    
    if (!p_values || !t_values || !log2fc || !baseMean) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free(p_values);
        free(t_values);
        free(log2fc);
        free(baseMean);
        free(scaling_factors);
        free(common_means);
        free(common_vars);
        stb_free_matrix(normalized_counts);
        return;
    }

    for (i = 0; i < n_genes; i++) {
        /* Extract group data from normalized counts */
        double *group1_data = malloc(config->group1_count * sizeof(double));
        double *group2_data = malloc(config->group2_count * sizeof(double));
        
        if (!group1_data || !group2_data) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            free(group1_data);
            free(group2_data);
            free(p_values);
            free(t_values);
            free(log2fc);
            free(baseMean);
            free(scaling_factors);
            free(common_means);
            free(common_vars);
            stb_free_matrix(normalized_counts);
            return;
        }
        
        for (j = 0; j < config->group1_count; j++) {
            group1_data[j] = normalized_counts->data[i][config->group1_start + j];
        }
        for (j = 0; j < config->group2_count; j++) {
            group2_data[j] = normalized_counts->data[i][config->group2_start + j];
        }

        /* Calculate means for log2 fold change */
        double mean1, var1, mean2, var2;
        stb_meanvar(group1_data, config->group1_count, &mean1, &var1);
        stb_meanvar(group2_data, config->group2_count, &mean2, &var2);
        
        baseMean[i] = (mean1 + mean2) / 2.0;
        log2fc[i] = stb_log2_fold_change(mean1, mean2, PSEUDOCOUNT);

        /* Perform moderated t-test */
        double t, p;
        stb_moderated_ttest(group1_data, config->group1_count, 
                           group2_data, config->group2_count,
                           prior_variance, prior_df, &t, &p);
        
        /* Convert to two-sided p-value */
        t_values[i] = t;
        p_values[i] = 2.0 * p;

        free(group1_data);
        free(group2_data);
    }

    /* Step 5: Apply Benjamini-Hochberg FDR correction */
    fprintf(stderr, "Step 5: Applying multiple testing correction (FDR)...\n");
    double *adjusted_p = malloc(n_genes * sizeof(double));
    if (!adjusted_p) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        free(p_values);
        free(t_values);
        free(log2fc);
        free(baseMean);
        free(scaling_factors);
        free(common_means);
        free(common_vars);
        stb_free_matrix(normalized_counts);
        return;
    }
    
    stb_adjust_pvalues_bh(p_values, n_genes, adjusted_p, config->fdr);

    /* Step 6: Output results */
    fprintf(stderr, "Step 6: Writing results...\n");
    fprintf(output, "Gene\tbaseMean\tlog2FoldChange\tstat\tpvalue\tpadj\n");
    
    int significant_count = 0;
    for (i = 0; i < n_genes; i++) {
        fprintf(output, "Gene_%d\t%.4f\t%.4f\t%.4f\t%.6e\t%.6e\n",
                i + 1, baseMean[i], log2fc[i], t_values[i], p_values[i], adjusted_p[i]);
        if (adjusted_p[i] < config->fdr) {
            significant_count++;
        }
    }

    fprintf(stderr, "Analysis complete!\n");
    fprintf(stderr, "Significant genes (FDR < %.2f): %d / %d (%.1f%%)\n", 
            config->fdr, significant_count, n_genes, 
            100.0 * significant_count / n_genes);

    /* Cleanup */
    free(scaling_factors);
    free(common_means);
    free(common_vars);
    free(p_values);
    free(t_values);
    free(log2fc);
    free(baseMean);
    free(adjusted_p);
    stb_free_matrix(normalized_counts);
}

int main(int argc, char **argv) {
    Config config;
    init_config(&config);

    /* Parse command line arguments */
    if (parse_args(argc, argv, &config) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    /* Load count matrix */
    fprintf(stderr, "Loading count matrix from: %s\n", config.input_file);
    STB_MAT *counts = stb_matrix_from_file(config.input_file);
    if (!counts) {
        fprintf(stderr, "Error: Could not load count matrix from %s\n", config.input_file);
        return 1;
    }

    /* Validate group configuration */
    int max_col = config.group1_start + config.group1_count;
    if (config.group2_start + config.group2_count > max_col) {
        max_col = config.group2_start + config.group2_count;
    }
    if (max_col > counts->columns) {
        fprintf(stderr, "Error: Group configuration exceeds available columns (%d)\n", counts->columns);
        stb_free_matrix(counts);
        return 1;
    }

    /* Open output file */
    FILE *output = stdout;
    if (config.output_file) {
        output = fopen(config.output_file, "w");
        if (!output) {
            fprintf(stderr, "Error: Could not open output file %s\n", config.output_file);
            stb_free_matrix(counts);
            return 1;
        }
    }

    /* Perform analysis */
    perform_deseq2_analysis(counts, &config, output);

    /* Cleanup */
    if (output != stdout) {
        fclose(output);
    }
    stb_free_matrix(counts);

    return 0;
}
