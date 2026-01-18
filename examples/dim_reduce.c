/* dim_reduce.c - Generic Dimension Reduction Tool
 * 
 * A flexible, robust C program for performing dimensionality reduction on tabular data.
 * Supports PCA, t-SNE, and UMAP with various input formats.
 *
 * Usage: dim_reduce [input_file] [options]
 */

#define STB_STATS_DEFINE
#include "stb_stats.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>

#define MAX_LINE_LENGTH 1048576
#define MAX_FIELDS 100000

typedef enum {
    ALGO_PCA,
    ALGO_TSNE,
    ALGO_UMAP
} Algorithm;

typedef struct {
    char *input_file;
    char *output_file;
    int col_major;
    int skip_cols;
    Algorithm algo;
    int n_components;
    int pre_pca;
    int no_scale;
    int no_eigen;
    double perplexity;
    int neighbors;
    double min_dist;
    int max_iter;
} Config;

/* Default configuration */
void init_config(Config *config) {
    config->input_file = NULL;
    config->output_file = NULL;
    config->col_major = 0;
    config->skip_cols = 1;
    config->algo = ALGO_PCA;
    config->n_components = 2;
    config->pre_pca = 0;
    config->no_scale = 0;
    config->no_eigen = 0;
    config->perplexity = 30.0;
    config->neighbors = 15;
    config->min_dist = 0.1;
    config->max_iter = 1000;
}

void print_usage(const char *prog_name) {
    fprintf(stderr, "Usage: %s [input_file] [options]\n\n", prog_name);
    fprintf(stderr, "Generic Dimension Reduction Tool\n");
    fprintf(stderr, "Supports PCA, t-SNE, and UMAP for dimensionality reduction.\n\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  input_file              Path to input TAB-delimited file (can be .gz)\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o, --output FILE       Output filename (default: stdout)\n");
    fprintf(stderr, "  -c, --col-major         Column-Major Mode (samples are columns)\n");
    fprintf(stderr, "  -s, --skip N            Number of description columns to skip (default: 1)\n");
    fprintf(stderr, "  -a, --algo ALGO         Algorithm: pca, tsne, umap (default: pca)\n");
    fprintf(stderr, "  -n, --n-components N    Number of output dimensions (default: 2)\n");
    fprintf(stderr, "  --pre-pca N             Run PCA first to reduce to N dimensions\n");
    fprintf(stderr, "  --no-scale              Disable Z-score normalization\n");
    fprintf(stderr, "  --no-eigen              (PCA) Output raw scores, not normalized\n");
    fprintf(stderr, "  --perplexity FLOAT      (t-SNE) Perplexity parameter (default: 30.0)\n");
    fprintf(stderr, "  --neighbors N           (UMAP) Number of neighbors (default: 15)\n");
    fprintf(stderr, "  --min-dist FLOAT        (UMAP) Minimum distance (default: 0.1)\n");
    fprintf(stderr, "  --max-iter N            Maximum iterations (default: 1000)\n");
    fprintf(stderr, "  -h, --help              Show this help message\n");
}

/* Read a line from regular or gzipped file */
char *read_line(FILE *fp, gzFile gzfp, char *buffer, size_t size) {
    if (gzfp) {
        return gzgets(gzfp, buffer, size);
    } else {
        return fgets(buffer, size, fp);
    }
}

/* Parse TAB-delimited line into fields */
int parse_line(char *line, char **fields, int max_fields) {
    int n = 0;
    char *token = strtok(line, "\t\n\r");
    while (token && n < max_fields) {
        fields[n++] = token;
        token = strtok(NULL, "\t\n\r");
    }
    return n;
}

/* Z-score normalization */
void zscore_normalize(double **data, int n, int p) {
    for (int j = 0; j < p; j++) {
        double mean = 0.0, var = 0.0;
        
        // Calculate mean
        for (int i = 0; i < n; i++) {
            mean += data[i][j];
        }
        mean /= n;
        
        // Calculate variance
        for (int i = 0; i < n; i++) {
            double diff = data[i][j] - mean;
            var += diff * diff;
        }
        var /= (n - 1);
        
        double std = sqrt(var);
        if (std < 1e-10) std = 1.0;  // Avoid division by zero
        
        // Normalize
        for (int i = 0; i < n; i++) {
            data[i][j] = (data[i][j] - mean) / std;
        }
    }
}

/* Load data from file */
int load_data(Config *config, double ***data_out, char ***sample_names_out, 
              int *n_samples, int *n_features) {
    FILE *fp = NULL;
    gzFile gzfp = NULL;
    
    // Check if gzipped
    size_t len = strlen(config->input_file);
    if (len > 3 && strcmp(config->input_file + len - 3, ".gz") == 0) {
        gzfp = gzopen(config->input_file, "r");
        if (!gzfp) {
            fprintf(stderr, "Error: Cannot open file %s\n", config->input_file);
            return -1;
        }
    } else {
        fp = fopen(config->input_file, "r");
        if (!fp) {
            fprintf(stderr, "Error: Cannot open file %s\n", config->input_file);
            return -1;
        }
    }
    
    char *buffer = (char *)malloc(MAX_LINE_LENGTH);
    char **fields = (char **)malloc(MAX_FIELDS * sizeof(char *));
    
    // Read header
    if (!read_line(fp, gzfp, buffer, MAX_LINE_LENGTH)) {
        fprintf(stderr, "Error: Empty file\n");
        free(buffer);
        free(fields);
        return -1;
    }
    
    int n_header_fields = parse_line(buffer, fields, MAX_FIELDS);
    
    // Read data lines to count rows
    int n_rows = 0;
    int n_data_fields = 0;
    long file_pos = fp ? ftell(fp) : gztell(gzfp);
    
    while (read_line(fp, gzfp, buffer, MAX_LINE_LENGTH)) {
        int nf = parse_line(buffer, fields, MAX_FIELDS);
        if (nf > 0) {
            n_rows++;
            if (n_data_fields == 0) n_data_fields = nf;
        }
    }
    
    // Reset file position
    if (fp) {
        fseek(fp, file_pos, SEEK_SET);
    } else {
        gzclose(gzfp);
        gzfp = gzopen(config->input_file, "r");
        read_line(fp, gzfp, buffer, MAX_LINE_LENGTH);  // Skip header again
    }
    
    // Determine format and dimensions
    int r_format = (n_header_fields == n_data_fields - 1);
    
    if (!config->col_major) {
        // Row-major: rows are samples
        *n_samples = n_rows;
        *n_features = n_data_fields - config->skip_cols;
        
        *sample_names_out = (char **)malloc(*n_samples * sizeof(char *));
        *data_out = (double **)stb_allocmat(*n_samples, *n_features, sizeof(double));
        
        int row_idx = 0;
        while (read_line(fp, gzfp, buffer, MAX_LINE_LENGTH)) {
            int nf = parse_line(buffer, fields, MAX_FIELDS);
            if (nf > 0 && row_idx < *n_samples) {
                (*sample_names_out)[row_idx] = strdup(fields[0]);
                
                for (int j = 0; j < *n_features; j++) {
                    (*data_out)[row_idx][j] = atof(fields[config->skip_cols + j]);
                }
                row_idx++;
            }
        }
    } else {
        // Column-major: columns are samples
        *n_samples = n_header_fields - config->skip_cols - (r_format ? 1 : 0);
        *n_features = n_rows;
        
        *sample_names_out = (char **)malloc(*n_samples * sizeof(char *));
        *data_out = (double **)stb_allocmat(*n_samples, *n_features, sizeof(double));
        
        // Store sample names from header
        int offset = config->skip_cols + (r_format ? 1 : 0);
        
        // Re-read header for sample names
        if (fp) {
            fseek(fp, 0, SEEK_SET);
        } else {
            gzclose(gzfp);
            gzfp = gzopen(config->input_file, "r");
        }
        read_line(fp, gzfp, buffer, MAX_LINE_LENGTH);
        parse_line(buffer, fields, MAX_FIELDS);
        
        for (int i = 0; i < *n_samples; i++) {
            (*sample_names_out)[i] = strdup(fields[offset + i]);
        }
        
        // Read data rows
        int row_idx = 0;
        while (read_line(fp, gzfp, buffer, MAX_LINE_LENGTH)) {
            int nf = parse_line(buffer, fields, MAX_FIELDS);
            if (nf > 0 && row_idx < *n_features) {
                for (int i = 0; i < *n_samples; i++) {
                    (*data_out)[i][row_idx] = atof(fields[offset + i]);
                }
                row_idx++;
            }
        }
    }
    
    free(buffer);
    free(fields);
    
    if (fp) fclose(fp);
    if (gzfp) gzclose(gzfp);
    
    return 0;
}

/* Write results to output */
void write_output(Config *config, double **result, char **sample_names, int n_samples, int n_components) {
    FILE *out = stdout;
    if (config->output_file) {
        out = fopen(config->output_file, "w");
        if (!out) {
            fprintf(stderr, "Error: Cannot open output file %s\n", config->output_file);
            out = stdout;
        }
    }
    
    // Write header
    fprintf(out, "SampleID");
    const char *prefix = (config->algo == ALGO_PCA) ? "PC" : 
                        (config->algo == ALGO_TSNE) ? "tSNE" : "UMAP";
    for (int i = 0; i < n_components; i++) {
        fprintf(out, "\t%s_%d", prefix, i + 1);
    }
    fprintf(out, "\n");
    
    // Write data
    for (int i = 0; i < n_samples; i++) {
        fprintf(out, "%s", sample_names[i]);
        for (int j = 0; j < n_components; j++) {
            fprintf(out, "\t%.6f", result[i][j]);
        }
        fprintf(out, "\n");
    }
    
    if (out != stdout) {
        fclose(out);
    }
}

int main(int argc, char **argv) {
    Config config;
    init_config(&config);
    
    // Parse command-line arguments
    static struct option long_options[] = {
        {"output", required_argument, 0, 'o'},
        {"col-major", no_argument, 0, 'c'},
        {"skip", required_argument, 0, 's'},
        {"algo", required_argument, 0, 'a'},
        {"n-components", required_argument, 0, 'n'},
        {"pre-pca", required_argument, 0, 1},
        {"no-scale", no_argument, 0, 2},
        {"no-eigen", no_argument, 0, 3},
        {"perplexity", required_argument, 0, 4},
        {"neighbors", required_argument, 0, 5},
        {"min-dist", required_argument, 0, 6},
        {"max-iter", required_argument, 0, 7},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "o:cs:a:n:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'o':
                config.output_file = optarg;
                break;
            case 'c':
                config.col_major = 1;
                break;
            case 's':
                config.skip_cols = atoi(optarg);
                break;
            case 'a':
                if (strcmp(optarg, "pca") == 0) {
                    config.algo = ALGO_PCA;
                } else if (strcmp(optarg, "tsne") == 0) {
                    config.algo = ALGO_TSNE;
                } else if (strcmp(optarg, "umap") == 0) {
                    config.algo = ALGO_UMAP;
                } else {
                    fprintf(stderr, "Error: Unknown algorithm '%s'\n", optarg);
                    return 1;
                }
                break;
            case 'n':
                config.n_components = atoi(optarg);
                break;
            case 1:  // --pre-pca
                config.pre_pca = atoi(optarg);
                break;
            case 2:  // --no-scale
                config.no_scale = 1;
                break;
            case 3:  // --no-eigen
                config.no_eigen = 1;
                break;
            case 4:  // --perplexity
                config.perplexity = atof(optarg);
                break;
            case 5:  // --neighbors
                config.neighbors = atoi(optarg);
                break;
            case 6:  // --min-dist
                config.min_dist = atof(optarg);
                break;
            case 7:  // --max-iter
                config.max_iter = atoi(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    
    // Get input file
    if (optind < argc) {
        config.input_file = argv[optind];
    } else {
        fprintf(stderr, "Error: Input file required\n\n");
        print_usage(argv[0]);
        return 1;
    }
    
    // Load data
    double **data;
    char **sample_names;
    int n_samples, n_features;
    
    fprintf(stderr, "Loading data from %s...\n", config.input_file);
    if (load_data(&config, &data, &sample_names, &n_samples, &n_features) != 0) {
        return 1;
    }
    
    fprintf(stderr, "Loaded %d samples with %d features\n", n_samples, n_features);
    
    // Apply Z-score normalization
    if (!config.no_scale) {
        fprintf(stderr, "Applying Z-score normalization...\n");
        zscore_normalize(data, n_samples, n_features);
    }
    
    // Apply pre-PCA if requested
    double **working_data = data;
    int working_features = n_features;
    double **pca_result = NULL;
    
    if (config.pre_pca > 0 && config.pre_pca < n_features && config.algo != ALGO_PCA) {
        fprintf(stderr, "Running pre-PCA to reduce to %d dimensions...\n", config.pre_pca);
        
        double *eigenvalues;
        double **eigenvectors;
        
        stb_pca(data, n_samples, n_features, NULL, NULL, config.pre_pca, 0, 
                &eigenvalues, &eigenvectors, &pca_result);
        
        working_data = pca_result;
        working_features = config.pre_pca;
        
        free(eigenvalues);
        free(eigenvectors);  // stb_allocmat allocated, single free
    }
    
    // Run dimensionality reduction
    double **result = NULL;
    
    switch (config.algo) {
        case ALGO_PCA: {
            fprintf(stderr, "Running PCA...\n");
            double *eigenvalues;
            double **eigenvectors;
            
            stb_pca(working_data, n_samples, working_features, NULL, NULL, 
                    config.n_components, 0, &eigenvalues, &eigenvectors, &result);
            
            // Print variance explained
            double total_var = 0.0;
            for (int i = 0; i < working_features; i++) {
                total_var += eigenvalues[i];
            }
            
            fprintf(stderr, "Variance explained by components:\n");
            for (int i = 0; i < config.n_components && i < working_features; i++) {
                fprintf(stderr, "  PC%d: %.2f%%\n", i + 1, 
                        (eigenvalues[i] / total_var) * 100.0);
            }
            
            free(eigenvalues);
            free(eigenvectors);  // stb_allocmat allocated, single free
            break;
        }
        
        case ALGO_TSNE: {
            fprintf(stderr, "Running t-SNE (perplexity=%.1f, max_iter=%d)...\n", 
                    config.perplexity, config.max_iter);
            stb_tsne(working_data, n_samples, working_features, config.n_components,
                    config.perplexity, config.max_iter, 200.0, &result);
            break;
        }
        
        case ALGO_UMAP: {
            fprintf(stderr, "Running UMAP (n_neighbors=%d, min_dist=%.2f)...\n",
                    config.neighbors, config.min_dist);
            stb_umap(working_data, n_samples, working_features, config.n_components,
                    config.neighbors, config.min_dist, 200, &result);
            break;
        }
    }
    
    // Write output
    fprintf(stderr, "Writing output...\n");
    write_output(&config, result, sample_names, n_samples, config.n_components);
    
    // Cleanup
    free(data);  // stb_allocmat allocated, single free
    for (int i = 0; i < n_samples; i++) {
        free(sample_names[i]);
    }
    free(sample_names);
    if (result) free(result);  // stb_allocmat allocated, single free
    
    if (pca_result) {
        free(pca_result);  // stb_allocmat allocated, single free
    }
    
    fprintf(stderr, "Done!\n");
    return 0;
}
