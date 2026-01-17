# dim_reduce - Generic Dimension Reduction Tool

A flexible, robust C program for performing dimensionality reduction on tabular data. It supports **PCA**, **t-SNE**, and **UMAP**, handling various input formats including gzipped files, row/column-major orientations, and automatic header detection.

## Features

*   **Algorithms**: PCA, t-SNE, and UMAP.
*   **Input Flexibility**: Reads TAB-delimited text (`.txt`, `.tsv`) and Gzip compressed files (`.gz`) seamlessly.
*   **Orientation Support**: 
    *   **Row-major** (Standard ML): Rows are samples, columns are features.
    *   **Column-major** (Genomics/Bioinformatics): Columns are samples, rows are features (e.g., gene expression matrices).
*   **Robust Parsing**: 
    *   Automatically detects headers.
    *   Handles "R-format" files (where the header is missing the ID column name).
    *   Skips non-data description columns.
*   **Preprocessing**:
    *   Automatic Z-score normalization (StandardScaler).
    *   **Pre-PCA**: Optional step to reduce dimensions via PCA *before* running UMAP/t-SNE to improve speed and denoising on high-dimensional data.

## Building

The program requires only standard C libraries and the STB Stats header file.

### Prerequisites

*   GCC or compatible C compiler
*   Standard C library
*   zlib (for gzipped file support)
*   Math library

### Compilation

```bash
make
```

Or manually:

```bash
gcc -O2 -Wall -Wextra -std=c99 -D_GNU_SOURCE -o dim_reduce dim_reduce.c -lm -lz
```

## Usage

```bash
./dim_reduce [input_file] [options]
```

### Arguments

| Argument | Description |
| :--- | :--- |
| **Input** | |
| `input_file` | Path to the input TAB-delimited file (can be `.gz`). |
| `-o`, `--output` | Output filename. Defaults to Standard Output (stdout). |
| **Parsing** | |
| `-c`, `--col-major` | **Column-Major Mode**. Switch to this if your samples are columns (e.g., single-cell matrices). Default is Row-Major (rows are samples). |
| `-s`, `--skip` | Number of description/ID columns to skip at the start of each row. Default: `1`. |
| **Algorithm** | |
| `-a`, `--algo` | Algorithm to use: `pca` (default), `umap`, or `tsne`. |
| `-n`, `--n-components`| Number of dimensions to output (e.g., 2 or 3). Default: `2`. |
| **Processing** | |
| `--pre-pca INT` | Run PCA first to reduce to `INT` dimensions before running UMAP/t-SNE. Useful for very high-dim inputs (e.g., `--pre-pca 50`). |
| `--no-scale` | Disable Z-score normalization (scaling). |
| `--no-eigen` | (PCA only) Output raw scores instead of normalized transformed coordinates. |
| **Hyperparameters** | |
| `--perplexity` | Perplexity for t-SNE. Default: `30.0`. |
| `--neighbors` | Number of neighbors for UMAP. Default: `15`. |
| `--min-dist` | Minimum distance for UMAP. Default: `0.1`. |
| `--max-iter` | Maximum iterations. Default: `1000`. |

## Input Formats

### 1. Row-Major (Default)
Each row is a data point (sample). The first `-s` columns are IDs/Descriptions.

```text
SampleID    Feature1    Feature2    Feature3
Sample_A    1.2         0.5         3.3
Sample_B    0.9         0.1         2.1
```

### 2. Column-Major (`-c`)
Each column is a data point (sample). Rows are features (genes, metrics). This is common in bioinformatics.

```text
GeneID      Sample_A    Sample_B    Sample_C
Gene_1      1.2         0.9         1.5
Gene_2      0.5         0.1         0.4
```

*Note: In Column-Major mode, the script expects sample names in the header. If the file is in "R-format" (header length = data length - 1), the script automatically adjusts alignment.*

## Examples

**1. Basic PCA on a standard CSV/TSV**
Rows are samples. Skip the first column (ID).
```bash
./dim_reduce data.txt -o result_pca.txt
```

**2. UMAP on Gene Expression Data (Column-Major)**
Samples are columns. The file is gzipped. We want to skip the first 2 columns (e.g., GeneID and GeneSymbol).
```bash
./dim_reduce expression.tab.gz -c -s 2 -a umap -o result_umap.txt
```

**3. High-Performance t-SNE**
For a very large dataset (e.g., 20k features), use `--pre-pca` to reduce to 50 dimensions before running t-SNE.
```bash
./dim_reduce big_data.txt -a tsne --pre-pca 50 --perplexity 50 -o result_tsne.txt
```

**4. 3D visualization with PCA**
```bash
./dim_reduce data.txt -a pca -n 3 -o result_pca_3d.txt
```

## Output Format

The output is a simple TAB-delimited file containing the Sample ID and the coordinates.

```text
SampleID    PC_1        PC_2
Sample_A    3.412301    -1.203910
Sample_B    1.902311    2.401293
...
```

For t-SNE, columns are named `tSNE_1`, `tSNE_2`, etc.
For UMAP, columns are named `UMAP_1`, `UMAP_2`, etc.

## Algorithms

### PCA (Principal Component Analysis)
- Linear dimensionality reduction technique
- Fast and deterministic
- Best for linearly separable data
- Provides variance explained information

### t-SNE (t-Distributed Stochastic Neighbor Embedding)
- Non-linear dimensionality reduction
- Excellent for visualization
- Preserves local structure
- Stochastic (results may vary between runs)
- Uses Barnes-Hut approximation for efficiency

### UMAP (Uniform Manifold Approximation and Projection)
- Non-linear dimensionality reduction
- Faster than t-SNE on large datasets
- Preserves both local and global structure
- Better for general-purpose dimensionality reduction

## Implementation Details

The implementation uses the STB Stats library (`stb_stats.h`) which provides:
- KD-tree for efficient nearest neighbor search
- Optimized matrix operations
- Barnes-Hut approximation for t-SNE
- Fast UMAP implementation with membership strength computation

## License

This software is dual-licensed to the public domain and under the following license: you are granted a perpetual, irrevocable license to copy, modify, publish, and distribute this file as you see fit.

## Citation

If you use this tool in a publication, please reference:

Voshol, G.P. (2024). STB: A simple Statistics Tool Box (Version 1.25) [Software]. 
Available from https://github.com/gerbenvoshol/Statistics-Tool-Box
