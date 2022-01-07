/* stb_stats.h - v1.21 - Statistics Tool Box -- public domain
					no warranty is offered or implied; use this code at your own risk

	 This is a single header file with a bunch of useful statistical functions

 ============================================================================
	 You MUST

			#define STB_STATS_DEFINE

	 in EXACTLY _one_ C or C++ file that includes this header, BEFORE the
	 include, like this:

			#define STB_STATS_DEFINE
			#include "stb_stats.h"

	 All other files should just #include "stb_stats.h" without the #define.
 ============================================================================

 Version History
        1.22  stb_shannon, stb_simpson, stb_jaccard, stb_bray_curtis, simple hash table
        1.21  stb_pdf_hypgeo hypergeometric distribution probability density function, speedup stb_log_factorial using lookup table
        1.20  stb_fisher2x2 simple fisher exact test for 2x2 contigency tables
        1.19  stb_pdf_binom and stb_pdf_pois, the binomial and poison probability density functions
        1.18  stb_polygamma, stb_trigamma_inverse gamme functions and stb_fit_f_dist for moment estimation of the scaled F-distribution
        1.17  stb_qnorm and stb_qnorm_with_reference (also matrix variants) quantile normalization between columns with and without a reference
        1.16  stb_neugas Neural gas clustering algorith
        1.15  stb_pca Principal Component Analysis
        1.14  stb_csm (confidence sequence method) for monte-carlo simulations
        1.13  stb_kmeans k-means++ classical data clustering
 		1.12  stb_qsort (Quicksort), could be used to replace current sorting method
 		1.11  stb_cdf_gumbel, stb_pdf_gumbel, stb_icdf_gumbel and stb_est_gumbel, the (inverse) cumulative/probability 
 		      density functions for the gumbel distribution and the ML estimator of the gumbel parameters
 		1.10  stb_kendall (Kendall's Rank correlation)
 		1.09  stb_jenks Initial port of O(k×n×log(n)) Jenks-Fisher algorithm originally created by Maarten Hilferink
 		1.08  stb_logistic_regression_L2 simple L2-regularized logistic regression
 		1.07  stb_spearman (Spearman's Rank correlation)
		1.06  stb_invert_matrix, stb_transpose_matrix, stb_matrix_multiply, etc., stb_multi_linear_regression
 		      stb_multi_logistic_regression 
		1.05  stb_ksample_anderson_darling, stb_2sample_anderson_darling, stb_expfit (Exponential fitting), 
		      stb_polyfit (Polynomial fitting), stb_powfit (Power curve fitting), stb_trap, stb_trapezoidal 
		      (returns the integral (area under the cruve) of a given function and interval) and 
		      stb_lagrange (polynomial interpolation), stb_sum (Neumaier summation algorithm)
 		1.04  stb_kruskal_wallis (Unfinished, needs a better way to handle Dunn's post-hoc test), stb_combinations, 
 		      stb_allocmat (simple allocation of 2d array, but might not work on all systems?!), stb_fgetln, stb_fgetlns
		1.03  stb_pcg32 (PCG-XSH-RR) and stb_xoshiro512 (xoshiro512**) Pseudo Random Number Generators
		1.02  stb_anova (One-Way Anova with Tukey HSD test and Scheffe T-statistics method (post-hoc) (Unfinished))
		1.01  stb_quartiles, stb_histogram (very simple histogram), stb_print_histogram, stb_free_histogram, 
		      stb_histogram_add, stb_factorial, stb_linfit (Liniear fitting)
		1.00  stb_meanvar, stb_ttest, stb_uttest, stb_ftest, stb_phi, stb_benjamini_hochberg
              stb_chisqr, stb_chisqr_matrix, stb_gtest, stb_gtest_matrix, stb_anderson_darling
		      stb_mann_whitney

 LICENSE

 This software is dual-licensed to the public domain and under the following
 license: you are granted a perpetual, irrevocable license to copy, modify,
 publish, and distribute this file as you see fit.

 CREDITS
 Jacob Wells (Gamma functions, Public domain), John D. Cook (Phi function, Public domain),
 David Blackman and Sebastiano Vigna (xoshiro512**, Public Domain), Sean T. Barrett (stb.h, Public Domain)
 Manas Sharma (Polynomial fitting, Public domain), Remi Dufour for the initial QuickSORT implementation (Public domain),
 Andy Allinger (k-means, Public domain)

 Written by Gerben Voshol.
*/

#ifndef STB__STATS__H
#define STB__STATS__H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <float.h>
#include <time.h>

#ifdef __cplusplus
#define STB_EXTERN   extern "C"
#else
#define STB_EXTERN   extern
#endif

#ifndef UINT64_MAX
#define UINT64_MAX (0xFFFFFFFFFFFFFFFFull)
#endif
#ifndef UINT32_MAX
#define UINT32_MAX (0xffffffff)
#endif

#define IS_ODD(n) (n & 1)
#define IS_EVEN(n) !IS_ODD(n)

#define CSV_QUOTE '\"'

/* Insertion sort threshold shift
 *
 * This macro defines the threshold shift (power of 2) at which the insertion
 * sort algorithm replaces the Quicksort.  A zero threshold shift disables the
 * insertion sort completely.
 *
 * The value is optimized for Linux and MacOS on the Intel x86 platform.
 */
#ifndef STB_INSERTION_SORT_THRESHOLD_SHIFT
# if defined (__APPLE__) && defined (__MACH__)
#  define STB_INSERTION_SORT_THRESHOLD_SHIFT (0)
# else
#  define STB_INSERTION_SORT_THRESHOLD_SHIFT (2)
# endif
#endif

/* Macro STB_DEFAULT_SWAP and STB_QSORT_SWAP
 *
 * Swaps the elements of two arrays.
 *
 * The length of the swap is determined by the value of "SIZE".  While both
 * arrays can't overlap, the case in which both pointers are the same works.
 */
#define STB_DEFAULT_SWAP(A,B,SIZE)               \
    do {                                             \
        register char       *a_byte = A;             \
        register char       *b_byte = B;             \
        register const char *a_end = a_byte + SIZE;  \
                                                     \
        while (a_byte < a_end)                       \
        {                                            \
            register const char swap_byte = *b_byte; \
            *b_byte++ = *a_byte;                     \
            *a_byte++ = swap_byte;                   \
        }                                            \
    } while (0)

#ifndef STB_QSORT_SWAP
#define STB_QSORT_SWAP(A,B,SIZE) STB_DEFAULT_SWAP(A, B, SIZE)
#endif

/* Performs a recursion to the left */
#define STB_RECURSE_LEFT             \
    if (left < store - size)             \
    {                                    \
        (++recursion)->left = left;      \
        recursion->right = store - size; \
    }

/* Performs a recursion to the right */
#define STB_RECURSE_RIGHT               \
    if (store + size < right)               \
    {                                       \
        (++recursion)->left = store + size; \
        recursion->right = right;           \
    }

/* Insertion sort inner-loop */
#define STB_INSERTION_SORT_LOOP(LEFT)                         \
    {                                                             \
        register char *trail = index - size;                      \
        while (trail >= LEFT && compare(trail, trail + size) > 0) \
        {                                                         \
            STB_QSORT_SWAP(trail, trail+size, size);                        \
            trail -= size;                                        \
        }                                                         \
    }

/* Performs insertion sort left of the pivot */
#define STB_INSERTION_STB_SORT_LEFT                \
    for (index = left + size; index < store; index +=size) \
        STB_INSERTION_SORT_LOOP(left)

/* Performs insertion sort right of the pivot */
#define STB_INSERTION_STB_SORT_RIGHT                        \
    for (index = store + (size << 1); index <= right; index +=size) \
        STB_INSERTION_SORT_LOOP(store + size)

/* Sorts to the left */
#if STB_INSERTION_SORT_THRESHOLD_SHIFT == 0
# define STB_SORT_LEFT STB_RECURSE_LEFT
#else
# define STB_SORT_LEFT                   \
    if (store - left <= threshold)           \
    {                                        \
        STB_INSERTION_STB_SORT_LEFT  \
    }                                        \
    else                                     \
    {                                        \
        STB_RECURSE_LEFT                 \
    }
#endif

/* Sorts to the right */
#if STB_INSERTION_SORT_THRESHOLD_SHIFT == 0
# define STB_SORT_RIGHT STB_RECURSE_RIGHT
#else
# define STB_SORT_RIGHT                  \
    if (right - store <= threshold)          \
    {                                        \
        STB_INSERTION_STB_SORT_RIGHT \
    }                                        \
    else                                     \
    {                                        \
        STB_RECURSE_RIGHT                \
    }
#endif

/* Function stb_qsort (Public Domain QuickSORT)
 *
 * This function performs a basic Quicksort.  This implementation is the
 * in-place version of the algorithm and is done in he following way:
 *
 * 1. In the middle of the array, we determine a pivot that we temporarily swap
 *    to the end.
 * 2. From the beginning to the end of the array, we swap any elements smaller
 *    than this pivot to the start, adjacent to other elements that were
 *    already moved.
 * 3. We swap the pivot next to these smaller elements.
 * 4. For both sub-arrays on sides of the pivot, we repeat this process
 *    recursively.
 * 5. For a sub-array smaller than a certain threshold, the insertion sort
 *    algorithm takes over.
 *
 * As an optimization, rather than performing a real recursion, we keep a
 * global stack to track boundaries for each recursion level.
 *
 * To ensure that at most O(log2 N) space is used, we recurse into the smaller
 * partition first.  The log2 of the highest unsigned value of an integer type
 * is the number of bits needed to store that integer.
 */
#define stb_qsort(A, L, S, C)                                             \
do {                                                                      \
    void *array = A;                                                      \
    size_t length = L;                                                    \
    size_t size = S;                                                      \
    int(*compare)(const void *, const void *) = C;                        \
                                                                          \
    /* Recursive stacks for array boundaries (both inclusive) */          \
    struct stackframe                                                     \
    {                                                                     \
        void *left;                                                       \
        void *right;                                                      \
    } stack[CHAR_BIT * sizeof(void *)];                                   \
                                                                          \
    /* Recursion level */                                                 \
    struct stackframe *recursion = stack;                                 \
                                                                          \
    /* Insertion sort threshold */                                        \
    const int threshold = size << STB_INSERTION_SORT_THRESHOLD_SHIFT;     \
                                                                          \
    /* Assign the first recursion level of the sorting */                 \
    recursion->left = array;                                              \
    recursion->right = (char *)array + size * (length - 1);               \
                                                                          \
    do                                                                    \
    {                                                                     \
        /* Partition the array */                                         \
        register char *index = recursion->left;                           \
        register char *right = recursion->right;                          \
        char          *left  = index;                                     \
                                                                          \
        /* Assigning store to the left */                                 \
        register char *store = index;                                     \
                                                                          \
        /* Pop the stack */                                               \
        --recursion;                                                      \
                                                                          \
        /* Determine a pivot (in the middle) and move it to the end */    \
        const size_t middle = (right - left) >> 1;                        \
        STB_QSORT_SWAP(left + middle - middle % size,right,size);           \
                                                                          \
        /* From left to right */                                          \
        while (index < right)                                             \
        {                                                                 \
            /* If item is smaller than pivot */                           \
            if (compare(right, index) > 0)                                \
            {                                                             \
                /* Swap item and store */                                 \
                STB_QSORT_SWAP(index,store,size);                           \
                                                                          \
                /* We increment store */                                  \
                store += size;                                            \
            }                                                             \
                                                                          \
            index += size;                                                \
        }                                                                 \
                                                                          \
        /* Move the pivot to its final place */                           \
        STB_QSORT_SWAP(right,store,size);                                   \
                                                                          \
        /* Recurse into the smaller partition first */                    \
        if (store - left < right - store)                                 \
        {                                                                 \
            /* Left side is smaller */                                    \
            STB_SORT_RIGHT                                            \
            STB_SORT_LEFT                                             \
                                                                          \
            continue;                                                     \
        }                                                                 \
                                                                          \
        /* Right side is smaller */                                       \
        STB_SORT_LEFT                                                 \
        STB_SORT_RIGHT                                                \
                                                                          \
    }                                                                     \
    while (recursion >= stack);                                           \
} while (0)

/* Function stb_creat_htable (Public Domain HASH Table)
 *
 * This function performs operations for a very basic hash table.
 *  
 * P is a prefix for all symbols of the table
 * K is the key type
 * V is the value type
 * H is an expression which, given `it` of type K returns the hash as a uint64_t
 * E is an expression which, given `l` and `r` of type K returns whether l and r are equal
 */
#define NB_START 8
#define SB_START 4
typedef unsigned long long _flag_t;
#define stb_create_htable(P, K, V, H, E)                                       \
  typedef struct {                                                             \
    K *keys;                                                                   \
    V *values;                                                                 \
    size_t *f_bucket;                                                          \
    size_t n_bucket;                                                           \
    size_t s_bucket;                                                           \
  } _##P##t;                                                                   \
  typedef _##P##t *P##t;                                                       \
  size_t P##bsize(P##t ptr) { return ptr->n_bucket * ptr->s_bucket; }          \
  size_t P##bstart(P##t ptr, size_t bucket) { return bucket * ptr->s_bucket; } \
  P##t _##P##new (size_t nb, size_t sb) {                                      \
    P##t ptr = malloc(sizeof(_##P##t));                                        \
    ptr->n_bucket = nb;                                                        \
    ptr->s_bucket = sb;                                                        \
    size_t total_size = P##bsize(ptr);                                         \
    ptr->keys = calloc(total_size, sizeof(K));                                 \
    ptr->values = calloc(total_size, sizeof(V));                               \
    ptr->f_bucket = calloc(ptr->n_bucket, sizeof(size_t));                     \
    return ptr;                                                                \
  }                                                                            \
                                                                               \
  P##t P##new (void) { return _##P##new (NB_START, SB_START); }                \
  void P##free(P##t ptr) {                                                     \
    free(ptr->keys);                                                           \
    free(ptr->values);                                                         \
    free(ptr->f_bucket);                                                       \
    free(ptr);                                                                 \
  }                                                                            \
  void P##put(P##t, K, V);                                                     \
  size_t P##nexti(P##t, size_t);                                               \
  size_t P##bufsize(P##t);                                                     \
  void _##P##resize(P##t table) {                                              \
    size_t ncount = table->n_bucket * 2;                                       \
    size_t nsize = table->s_bucket * 2;                                        \
    P##t n = _##P##new (ncount, nsize);                                        \
    for (size_t i = 0; i < P##bsize(table); i = P##nexti(table, i)) {          \
      P##put(n, table->keys[i], table->values[i]);                             \
    }                                                                          \
    free(table->keys);                                                         \
    free(table->values);                                                       \
    free(table->f_bucket);                                                     \
    *table = *n;                                                               \
    free(n);                                                                   \
  }                                                                            \
  static inline uint64_t P##hash(K it) { return H; }                                  \
  static inline _Bool P##eq(K l, K r) { return E; }                                   \
  void P##put(P##t table, K k, V v) {                                          \
    uint64_t hash = P##hash(k);                                                \
    size_t bucket = hash % table->n_bucket;                                    \
    size_t start = P##bstart(table, bucket);                                   \
    size_t insertion_point = start;                                            \
    for (;;) {                                                                 \
      K candidate = table->keys[insertion_point];                              \
      if (table->f_bucket[bucket] && !P##eq(candidate, k) &&                   \
          insertion_point - start < table->f_bucket[bucket]) {                 \
        if (insertion_point - start < table->s_bucket - 1) {                   \
          insertion_point++;                                                   \
        } else {                                                               \
          _##P##resize(table);                                                 \
          P##put(table, k, v);                                                 \
          break;                                                               \
        }                                                                      \
      } else {                                                                 \
        table->keys[insertion_point] = k;                                      \
        table->values[insertion_point] = v;                                    \
        table->f_bucket[bucket]++;                                             \
        break;                                                                 \
      }                                                                        \
    }                                                                          \
  }                                                                            \
  _Bool P##get(P##t table, K k, V *place) {                                    \
    uint64_t hash = P##hash(k);                                                \
    size_t bucket = hash % table->n_bucket;                                    \
    size_t start = P##bstart(table, bucket);                                   \
    for (size_t i = 0; i <= table->f_bucket[bucket]; ++i) {                    \
      if (P##eq(k, table->keys[start + i])) {                                  \
        *place = table->values[start + i];                                     \
        return true;                                                           \
      }                                                                        \
    }                                                                          \
    return false;                                                              \
  }                                                                            \
                                                                               \
  size_t P##nexti(P##t table, size_t idx) {                                    \
    if (idx >= P##bsize(table) - 1)                                            \
      return idx + 1;                                                          \
    idx++;                                                                     \
    if (table->keys[idx])                                                      \
      return idx;                                                              \
    return P##nexti(table, idx);                                               \
  }

STB_EXTERN int stb_strcmp(char * l, char *r);

/* Example usage
// Note the lack of ; !!!
stb_create_htable(stb_ii_, int, int, it, l == r)
stb_create_htable(stb_ss_, char const *, char const *, stb_murmer64(it), stb_strcmp(l, r) == 0)

int main() {}
    // create hashtable with int key and int value
    stb_ii_t t = stb_ii_new();
    stb_ii_put(t, 1, 3);
    stb_ii_put(t, 2, 3);
    stb_ii_put(t, 3, 3);
    stb_ii_put(t, 11, 3);
    stb_ii_put(t, 4, 3);
    stb_ii_put(t, 8, 3);
    stb_ii_put(t, 7, 3);
    stb_ii_put(t, 7, 3);
    
    for (int i = 10000; i < 10016; ++i) {
        stb_ii_put(t, i, i * 2);
    }

    for (size_t i = 0; i < stb_ii_bsize(t); i = stb_ii_nexti(t, i)) {
        if (t->keys[i] * 2 != t->values[i]) {
            fprintf(stderr, "Key ERROR\n");
        }
    }

    for (int i = 10000; i < 10016; ++i) {
        int u;
        stb_ii_get(t, i, &u);
        fprintf(stderr, "Key: %i Value: %i", i, u);
    }

    stb_ii_free(t);

    // create hash table with string key and string value
    stb_ss_t t = stb_ss_new();
    stb_ss_put(t, "1", "2");
    stb_ss_put(t, "2", "3");
    stb_ss_put(t, "6", "9");
    stb_ss_put(t, "3", "2");
    stb_ss_put(t, "1", "2");

    char const *ptr = NULL;
    stb_ss_get(t, "1", &ptr));

    stb_ss_free(t);
}
*/

struct stb_hist {
	int number_of_bins;
	double bin_width;
	double min;
	double max;
	double *count;
};

/* Given  an  array  of data (length n), this routine  returns  its  mean and sample variance calculated using the Welford’s method. */
STB_EXTERN void stb_meanvar(double *data, int n, double *mean, double *sample_variance);

/* Given  two arrays  of data (length n1 and n2), this routine  returns  its  t value and its significance (p) as probabiilty
 * NOTE: assumes equal variance 
 */
STB_EXTERN void stb_ttest(double *data1, int n1, double *data2, int n2, double *t, double *p);

/* Given  two arrays  of data (length n1 and n2), this routine  returns  its  t value and its significance (p) as probabiilty
 * NOTE: assumes unequal variance 
 */
STB_EXTERN void stb_uttest(double *data1, int n1, double *data2, int n2, double *t, double *p);

/* Given  two arrays  of data (length n1 and n2), this routine  returns  its  f value and its significance (p) as probabiilty */
STB_EXTERN void stb_ftest(double *data1, int n1, double *data2, int n2, double *f, double *p);

/* This function takes two arrays, one with the data (of total size n) and one specifying the group (0, 1, 2, ...) (with g groups) the data belongs to (total size n) and
 * performs an one-way anova and followed by the Tukey HSD test and Scheffe's method of mutliple comparisons 
 */
STB_EXTERN void stb_anova(double *data, int n, int *groups, int g, double *f, double *p);

/* This function takes two arrays, one with the data (of total size n) and one specifying the group (0, 1, 2, ...) (with g groups) the data belongs to (total size n) and
 * performs an one-way Kruskal-Wallis test.
 * Note:the number of samples per groups should be >> 5, otherwise the distribution is not similar enough to the ChiSqr
 */
STB_EXTERN void stb_kruskal_wallis(double *data, int n, int *groups, int g, double *H, double *p);

/* The function Φ(x) is the cumulative density function (CDF) of a standard normal (Gaussian) random variable */
STB_EXTERN double stb_phi(double x);

/* Given a rank (starting at 1 for the lowest p-value), the total number of comparisons performed and
 * the desired false discovery rate this function returns the Benjamini-Hochberg threshold. In other words,
 * if you put the individual P values in order, from smallest to largest. The smallest P value has a rank
 * 1, the next has 2, etc. Compare each individual P value to its Benjamini-Hochberg critical value (BHCV,
 * (rank/number_of_samples)*FDR. The largest P value that has P<BHCV is significant, and all of the P values
 * smaller than it are also significant, even the ones that aren't less than their Benjamini-Hochberg
 * critical value (see example).In general a FDR between 0.1 - 0.2 is reasonable.
 *
 * Dietary variable P value Rank BHCV
 * Total           <0.002   1    0.010
 * Carbohydrates    0.009   2    0.020
 * Fat              0.035   3    0.030
 * Sugars           0.044   4    0.040 <- Note that this value is higher than the BHCV but is
 * Proteins         0.046   5    0.050 <- with a FDR of 0.25 everyting BEFORE protein is siginificant also the sugars
 * Saturated        0.062   6    0.060
 */
STB_EXTERN double stb_benjamini_hochberg(int rank, int number_of_comparisons, double FDR);

/* Returns the adjusted Bonferroni P-value, the advantage is that you don't need to have a sorted list of P-values,
 * However, bonferroni adjustment is extremely conservative and leads to many false negatives!
 */
STB_EXTERN double stb_bonferroni(double p, int number_of_comparisons);

/* Perform the Mann-Whitney U test, the resulting corrected P-value should be the same as that of the R statistical program.
 * However, for tied data no exact P values is calculated. It is possible, but would require a more complicated program see
 * for example:
 * Alexander Marx, Christina Backes, Eckart Meese, Hans-Peter Lenhof, Andreas Keller,
 * EDISON-WMW: Exact Dynamic Programing Solution of the Wilcoxon–Mann–Whitney Test,
 * Genomics, Proteomics & Bioinformatics, Volume 14, Issue 1, 2016, Pages 55-61
 */
STB_EXTERN void stb_mann_whitney(double *data1, int n1, double *data2, int n2, double *U, double *p);

/* Calculates the two-tailed P-value for the Fisher Exact test.
 *                Men Women   Row total
 * Studying         1     9   10
 * Not-studying    11     3   14
 * Column total    12    12   24
 *
 * p = 2.759456e-03 two-sided
 */
STB_EXTERN void stb_fisher2x2(int a, int b, int c, int d, double *p);

/* Given a vector with observed and expected values with the amount of rows and columns, this function
 * performs the chisquare test and returns 
 */
STB_EXTERN void stb_chisqr(double *observed, double *expected, int rows, int columns, double *CV, double *p);

/* Given a matrix with observed and expected values with the amount of rows anc columns, this function
 * performs the chisquare test and returns 
 */
STB_EXTERN void stb_chisqr_matrix(double **observed, double **expected, int rows, int columns, double *CV, double *p);

/* Given a vector with observed and expected values with the amount of rows anc columns, this function
 * performs the G-test and returns the probability (based on the ChiSquared distribution)
 */
STB_EXTERN void stb_gtest(double *observed, double *expected, int rows, int columns, double *CV, double *p);

/* Given a matrix with observed and expected values with the amount of rows and columns, this function
 * performs the G-test and returns the probability (based on the ChiSquared distribution)
 */
STB_EXTERN void stb_gtest_matrix(double **observed, double **expected, int rows, int columns, double *CV, double *p);

/* Given an array of data, this function tests if the data is normally distributed
 * and returns the A value and the p-value. If the p-value is small (<0.05) the data is
 * not normally distributed! Care should be taken that the number of samples is not too
 * small, since this can lead to false positives/negatives (sample size > 7)
 */
STB_EXTERN void stb_anderson_darling(double *data, int n, double *A, double *p);

/* Anderson-Darling 2-sample test. */
STB_EXTERN void stb_2sample_anderson_darling(double *data1, int n1, double *data2, int n2, double *A, double *p);

/* Anderson-Darling k-sample test. */
STB_EXTERN void stb_ksample_anderson_darling(double *data, int n, int *groups, int g, double *A, double *p);

/* This function returns the min , q1, q2, q3 and max values according to https://en.wikipedia.org/wiki/Quartile */
STB_EXTERN void stb_quartiles(double *data, int n, double *min, double *q1, double *median, double *q3, double *max);

/* Free a histogram */
STB_EXTERN void stb_free_histogram(struct stb_hist *hist);

/* Keeps a score histogram, in which scores are counted into bins 0 <= x < 1, bin_width can be
 * a desired bin_width, or to calculate bin width using Freedman–Diaconis rule or -2 to calculate bin width using Rice rule
 * for this the amount of data points (n) should be > 3. To just initialize an empty histogram use n = 0
 * AND provide umin, umax and bin_width!
 */
STB_EXTERN struct stb_hist *stb_histogram(double *data, int n, double bin_width, double umin, double umax);

// print the histogram
STB_EXTERN void stb_print_histogram(struct stb_hist *hist);

// add some data to the histogram previously made using stb_histogram
STB_EXTERN void stb_histogram_add(double *data, int n, struct stb_hist *hist);

/* Calculate the factorial */
STB_EXTERN double stb_factorial(int n);

/* Calculate the factorial */
STB_EXTERN double stb_log_factorial(int n);

/* Compute (numerical stable) sum of the data using the Neumaier summation algorithm */
STB_EXTERN double stb_sum(double *data, int n);

/* Returns Spearman's Rank correlation of two vectors x and y, each of size n */
STB_EXTERN double stb_spearman(double *x, double *y, int n);

/* Returns Kendall's Rank correlation and the probability of two vectors x and y, each of size n */
STB_EXTERN double stb_kendall(double *x, double *y, int n, double *tau, double *z, double *prob);

// PDF
STB_EXTERN void stb_pdf_gumbel(double x, double mu, double sig, double *p);

// CDF
STB_EXTERN void stb_cdf_gumbel(double x, double mu, double sig, double *p);

// inverse CDF a.k.a quantile function
STB_EXTERN double stb_icdf_gumbel(double mu, double sig, double p);

/* Fit the Gumbel distribution using maximum likelihood estimation (MLE) */
STB_EXTERN void stb_est_gumbel(double *data, int n, double *mu, double *sig);

/* Calculate the linear regression
 * y = ax + b
 *
 * x,y  = arrays of data
 * n = number of data points
 * a = output slope
 * b = output intercept
 * r = output correlation coefficient (can be NULL if you don't want it)
 */
STB_EXTERN int stb_linfit(const double *x, const double *y, int n, double *a, double *b, double *r);

/* Calculate the exponential regression
 * Note: The data is tranformed before fitting it, so we might not get the best fit
 * y = a * e^(b * x)
 *
 * x,y  = arrays of data
 * n = number of data points
 * a = output base
 * b = output exp
 * r = output correlation coefficient (can be NULL if you don't want it)
 */
STB_EXTERN int stb_expfit(const double *x, const double *y, int n, double *a, double *b, double *r);

/* Calculate the polynomial regression
 * y = retval[2]* x^2 + retval[1]* x^1 + retval[0]* x^0
 *
 * x = x-axis values
 * y = y-axis values
 * N = number of data-points
 * n = degree of polynomial
 */
STB_EXTERN double *stb_polyfit(const double *x, const double *y, int N, int n);

/* Calculate the power regression
 * Note: The data is tranformed before fitting it, so we might not get the best fit
 * y = a * x^b
 *
 * x,y  = arrays of data
 * n = number of data points
 * a = output base
 * b = output exponent
 * r = output correlation coefficient (can be NULL if you don't want it)
 */
STB_EXTERN int stb_powfit(const double *x, const double *y, int n, double *a, double *b, double *r);

/* Lagrange interpolation 
 * x = x values
 * y = y values
 * n = number of data points
 * xp = value of x to interpolate y -> f(xp) = sum (the interpolated value)
 */
STB_EXTERN double stb_lagrange(double *x, double *y, int n, double xp);

typedef double (*stb_function)(double);

/* This function returns the area under a curve of a given function f
 * a = x[0]
 * b = x[n]
 * n = number of intervals
 * f = the curve function that is used to calculate y (the height of each slice)
 */
STB_EXTERN double stb_trap(double a, double b, int n, stb_function f);

/* This function perform integration by trapezoidal rule until converged to the given accuracy
 * or till LIMIT iterations has been reached and returns the area under a curve of a given function f
 * (Requires: stb_trap)
 * a = x[0]
 * b = x[n]
 * n = current number of intervals
 * accuracy = the desired accuracy
 * f = the curve function that is used to calculate y (the height of each slice)
 *
 * Example:
 * double func1(double x)
 * {
 *     return (1.0 / (1.0 + stb_sqr(x)));
 * }
 * 
 * int main()
 * {
 *     int intervals = 10;
 *     double accuracy = 0.000001;
 *     double x0 = 1;
 *     double xn = 5;
 *     area = stb_trapezoidal(x0, xn, &intervals, accuracy, func1);
 *     printf("Area = %lf (required: %i intervals)\n", area, intervals);
 *     return 0;
 * }
 */
STB_EXTERN double stb_trapezoidal(double a, double b, int *n, double accuracy, stb_function f);

/* stb_uintXX_to_double generates a real number in the interval [1..2), and then subtracts 1.0 to
 * obtain a real number in the interval [0..1) as follows:
 * First it sets the binary exponent of a floating point number to 1+bias and set the
 * mantissa to random bits. This will give a random number in the interval [1,2).
 * Then subtract 1.0 to get a random number in the interval [0,1).
 */
static inline double stb_uint64_to_double(uint64_t x)
{
	const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
	return u.d - 1.0;
}

static inline double stb_uint32_to_double(uint32_t x)
{
	const union { uint32_t i[2]; double d; } u = { .i[0] = x >> 20, .i[1] = (x >> 12) | UINT32_C(0x3FF00000)};
	return u.d - 1.0;
}

/* Helper functions for xoshiro512 and pcg32 */
static inline uint64_t stb_rotate64(const uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}

static inline uint32_t stb_rotate32(const uint32_t x, int k)
{
	return (x << k) | (x >> (32 - k));
}

static inline double stb_sqr(double x)
{
	return x * x;
}

/* Used for the initialization of xoshiro512**, but could also be used as a stand alone PRNG */
STB_EXTERN uint64_t stb_splitmix64(uint64_t seed);

/* Seed xoshiro512** */
STB_EXTERN uint64_t *stb_sxoshiro512(uint64_t seed);

/* xoshiro512** PRNG with 64-bit output and 512-bit state*/
STB_EXTERN uint64_t stb_xoshiro512(uint64_t *s);

/* stb_xoshiro512_bounded returns a uniformly distributed interger, r, in the range [0, n). */
STB_EXTERN uint64_t stb_xoshiro512_bounded(uint64_t *s, uint64_t n);

/* PCG-XSH-RR PRNG with 64-bit state and 32-bit output */
STB_EXTERN uint32_t stb_pcg32(uint64_t *s);

/* Seed the pcg32 PRNG */
STB_EXTERN uint64_t stb_spcg32(uint64_t seed);

/* stb_pcg32_bounded returns a uniformly distributed integer, r, in the range [0, n). */
STB_EXTERN uint32_t stb_pcg32_bounded(uint64_t *s, uint32_t n);

STB_EXTERN double stb_pcg32_uniform(uint64_t *seed);

// Gaussian (normal) random sample with mean 0 and standard deviation 1 from
// Knuth and Marsaglia and Bray, ``A Convenient Method for Generating Normal Variables''
STB_EXTERN double stb_pcg32_gauss(uint64_t *seed);

// Gaussian (normal) random sample with specified mean and standard deviation
STB_EXTERN double stb_pcg32_gauss_msd(uint64_t *seed, double mean, double stdev);

// Implementation based on "A Simple Method for Generating Gamma Variables"
// by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
// Vol 26, No 3, September 2000, pages 363-372.
// shape (alpha)  and scale (lambda)
STB_EXTERN double stb_pcg32_gamma(uint64_t *seed, double shape, double scale);

STB_EXTERN double stb_pcg32_exponential(uint64_t *seed);

// exponential random sample with specified mean
STB_EXTERN double stb_pcg32_exponential_m(uint64_t *seed, double mean);

// Knuth: mean (lambda)
STB_EXTERN double stb_pcg32_poisson(uint64_t *seed, const double mean);

STB_EXTERN double stb_pcg32_nbinom(uint64_t *seed, double size, double prob);

STB_EXTERN double stb_pcg32_chisquare(uint64_t *seed, double degrees_of_freedom);

STB_EXTERN double stb_pcg32_invgamma(uint64_t *seed, double shape, double scale);

STB_EXTERN double stb_pcg32_beta(uint64_t *seed, double a, double b);

STB_EXTERN double stb_pcg32_nbinom_mu(uint64_t *seed, double size, double mu);

/**
 * Confidence Sequence Method.
 *
 * See "A simple method for implementing Monte Carlo tests,"
 * Ding, Gandy, and Hahn, 2017 (https://arxiv.org/abs/1611.01675).
 *
 * Given n trials and s successes, can we conclude that the success rate
 * differs from alpha with exp(log_eps) false positive rate?
 *
 * Output the current log confidence level in OUT_log_level if non-NULL.
 */
STB_EXTERN int stb_csm(uint64_t n, double alpha, uint64_t s, double log_eps, double *OUT_log_level);

/* returns the number of combinations of t out of n and returns a matrix (2d array) of them
 * uses Knuth Algorithm T (section 2.7.1.3 from The Art of Computer Programming)
 */
STB_EXTERN int stb_combinations(int n, int t, int ***combinations);

/* Allocates a matric of rowsize*colsize with a single call to malloc, so
 * only a single call to free is necessary. This is a helper function for stb_combinations,
 * but might be useful as a general function. Please note that it might not work on all systems, 
 * since it assumes a void* is the same as int* etc and that there are no alignment offsets.
 * rowsize = the number of rows
 * colsize = the number of columns
 * item_size = the size of the elements (e.g sizeof(int) for a matrix with ints)
 */
STB_EXTERN void **stb_allocmat(int rowsize, int colsize, size_t item_size);

/* Reads an entire file into an array of strings, needs only a single call to free */
STB_EXTERN char **stb_fgetlns(char *filename, size_t *number_of_lines);

/* Dynamic allocation version of fgets(), capable of reading unlimited line lengths. */
STB_EXTERN char *stb_fgetln(char **buf, int *n, FILE *fp);

/*
 *  Given a string containing no linebreaks, or containing line breaks
 *  which are escaped by "double quotes", extract a NULL-terminated
 *  array of strings, one for every cell in the row.
 */
STB_EXTERN char **stb_parse_csv(const char *line, const char delim, int *nrfields);
STB_EXTERN int stb_count_fields(const char *line, const char delim);
STB_EXTERN void stb_free_csv_line(char **parsed);

/* Simple Matrix Function */

/* Matrix structure */
typedef struct {
	double **data;
	int rows;
	int columns;
} STB_MAT;

/* Print a matrix */
STB_EXTERN void stb_matrix_print(STB_MAT *matrix);

/* Free a matrix */
STB_EXTERN void stb_free_matrix(STB_MAT *matrix);

/* Return a new matrix of rowx column size initialized with 0 */
STB_EXTERN STB_MAT *stb_new_matrix(int rows, int columns);

/* Retrn an identity matrix of rowxcolumn size */
STB_EXTERN STB_MAT *stb_identity_matrix(int rows, int columns);

/* Return a copy of a matrix */
STB_EXTERN STB_MAT *stb_dup_matrix(STB_MAT *matrix);

/* Fill a matrix with a value */
STB_EXTERN void stb_fill_matrix(STB_MAT *matrix, double value);

/* Join two matrixes together and store the result in D. Example:
 * A |   B   =    D 
 *
 * 1 | 1 0 0   1 1 0 0
 * 1 | 0 1 0 = 1 0 1 0
 * 1 | 0 0 1   1 0 0 1
 */
STB_EXTERN void stb_join_matrix(STB_MAT *A, STB_MAT *B, STB_MAT **D);

/* Multiply two matrixes A and B, store the result in D */
STB_EXTERN void stb_matrix_multiply(STB_MAT *A, STB_MAT *B, STB_MAT **D);

/* Substract two equal sized matrixes A and B, store the result in D */
STB_EXTERN void stb_matrix_subs(STB_MAT *A, STB_MAT *B, STB_MAT **D);

/* Add two equal sized matrixes together A and B, store the result in D */
STB_EXTERN void stb_matrix_add(STB_MAT *A, STB_MAT *B, STB_MAT **D);

/* Transopse matrix A  and store it in Atransposed */
STB_EXTERN void stb_transpose_matrix(STB_MAT *A, STB_MAT **Atransposed);

/* Invert matrix A and store it in Ainverted*/
STB_EXTERN void stb_invert_matrix(STB_MAT *A, STB_MAT **Ainverted);

/* Multiply an entire matrix by a single value */
STB_EXTERN void stb_matrix_multiply_by_value(STB_MAT *A, double value);

/* Divide an entire matrix by a single value */
STB_EXTERN void stb_matrix_divide_by_value(STB_MAT *A, double value);

/* Returns the sum of a matrix */
double stb_matrix_sum(STB_MAT *A);

/* Returns a struct STB_MAT from a double matrix and the number of rows and columns */
STB_EXTERN STB_MAT *stb_matrix_from_double(double **data, int rows, int columns);

/* Reads a matrix from a tab seperated file.
 * The first line has two values the number of rows and the number of columns in the file.
 * The rest of the data are the actual values. Example:
 * 3	1 (3 rows and 1 column)
 * 62	
 * 66	
 * 61	
 */
STB_EXTERN STB_MAT *stb_matrix_from_file(char *filename);

/* Perform quantile normalization between columns without a reference
In statistics, quantile normalization is a technique for making two distributions identical in statistical properties.
see: https://en.wikipedia.org/wiki/Quantile_normalization
Arrays 1 to 3, genes A to D
A    5    4    3
B    2    1    4
C    3    4    6
D    4    2    8
Converted to STB_MAT file test.mat:
4 3
5 4 3
2 1 4
3 4 6
4 2 8

will become:
5.666667E+00 5.166667E+00 2.000000E+00 
2.000000E+00 2.000000E+00 3.000000E+00 
3.000000E+00 5.166667E+00 4.666667E+00 
4.666667E+00 3.000000E+00 5.666667E+00

statistics:
Min.   :2.000   Min.   :2.000   Min.   :2.000  
1st Qu.:2.750   1st Qu.:2.750   1st Qu.:2.750  
Median :3.833   Median :4.083   Median :3.833  
Mean   :3.833   Mean   :3.833   Mean   :3.833  
3rd Qu.:4.917   3rd Qu.:5.167   3rd Qu.:4.917  
Max.   :5.667   Max.   :5.167   Max.   :5.667  
 */
STB_EXTERN void stb_qnorm_matrix(STB_MAT *original);
STB_EXTERN void stb_qnorm(double **data, int rows, int columns);

/* Perform quantile normalization between columns with a reference, making the distribution of original equal to that of the reference
In statistics, quantile normalization is a technique for making two distributions identical in statistical properties.

NOTE:
Given a reference distribution, the original distribution is normalized by replacing each of its values by the value of the variable with the same rank 
in the reference distribution. If the reference distribution contains multiple samples, the original and reference distributions will only be identical if 
the reference distribution is first quantile normalized across all samples!

see: https://en.wikipedia.org/wiki/Quantile_normalization
Arrays 1 to 3, genes A to D
A    5    4    3
B    2    1    4
C    3    4    6
D    4    2    8
Converted to STB_MAT file test.mat:
4 3
5 4 3
2 1 4
3 4 6
4 2 8

Using the QUNATILE NORMALIZED reference:
4 3
5.666667E+00 5.166667E+00 2.000000E+00 
2.000000E+00 2.000000E+00 3.000000E+00 
3.000000E+00 5.166667E+00 4.666667E+00 
4.666667E+00 3.000000E+00 5.666667E+00

will become:
5.666667E+00 5.166667E+00 2.000000E+00 
2.000000E+00 2.000000E+00 3.000000E+00 
3.000000E+00 5.166667E+00 4.666667E+00 
4.666667E+00 3.000000E+00 5.666667E+00

statistics:
Min.   :2.000   Min.   :2.000   Min.   :2.000  
1st Qu.:2.750   1st Qu.:2.750   1st Qu.:2.750  
Median :3.833   Median :4.083   Median :3.833  
Mean   :3.833   Mean   :3.833   Mean   :3.833  
3rd Qu.:4.917   3rd Qu.:5.167   3rd Qu.:4.917  
Max.   :5.667   Max.   :5.167   Max.   :5.667  
 */
STB_EXTERN void stb_qnorm_matrix_with_reference(STB_MAT *original, STB_MAT *reference);
STB_EXTERN void stb_qnorm_with_reference(double **data, double **reference, int rows, int columns);

/* Perform a simple linear regression and return a vector containing the Beta values, the T-test values 
 * and the corresponding P-values. The formula determined using the least squared method is:
 * Y = Beta[0] + Beta[1] * X[0] + Beta[2] * X[1] + Beta[n] * X[n-1]
 * 
 * Note: This can also be calculated using a design matrix (1 on first column and X values for the rest)
 */
STB_EXTERN void stb_multi_linear_regression(STB_MAT *A, STB_MAT *Y, double **beta, double **tvalue, double **pvalue);

/* Perform a simple logistic regression and return a vector containing the Beta values, the Z-test values 
 * , the corresponding P-values and return the log-likelihood. The formula determined using the newton method is:
 * ln(1 / 1 - Y) = Beta[0] + Beta[1] * X[0] + Beta[2] * X[1] + Beta[n] * X[n-1]
 *
 * Note: This can also be calculated using a design matrix (1 on first column and X values for the rest)
 */
STB_EXTERN double stb_multi_logistic_regression(STB_MAT *A, STB_MAT *Y, double **beta, double **zvalue, double **pvalue);

/* L2-regularized logistic regression
 * NOTE: Unlike stb_multi_logistic_regression, to get the intercept (bias), add a column of 1.0s to matrix A
 * NOTE: Unlike stb_multi_logistic_regression, Y = 1 or -1 instead of 1 or 0!
 */
STB_EXTERN void stb_logistic_regression_L2(STB_MAT *A, STB_MAT *Y, double **beta, double **zvalue, double **pvalue);

/* simple logistic regression
 * NOTE: Unlike stb_multi_logistic_regression, to get the intercept (bias), add a column of 1.0s to matrix A
 * NOTE: Unlike stb_multi_logistic_regression, Y = 1 or -1 instead of 1 or 0!
 */
STB_EXTERN void stb_logistic_regression(STB_MAT *A, STB_MAT *Y, double **beta, double **zvalue, double **pvalue);

/* Filter items in place and return list of unique numbers, their count and the number of unique numbers.
 * NOTE: if a separate list is desired, duplicate it before calling this function 
 */
STB_EXTERN int stb_dunique(double *values, double **counts, int len);

/* Filter items in place and return list of unique numbers, their count and the tumber of unique numbers.
 * NOTE: if a separate list is desired, duplicate it before calling this function 
 */
STB_EXTERN int stb_iunique(int *values, int **counts, int len);

/**
 * Main entry point for creation of Jenks-Fisher natural breaks.
 * Port of Jenks/Fisher breaks originally created in C++ by Maarten Hilferink.
 * @param values array of the values, do not need to be sorted.
 * @param k number of breaks to create
 * @param len length of values array
 * @return Array with breaks
 */
STB_EXTERN double *stb_jenks(double *values, int len, int k);

/* calculate the euclidean distance between two points */
STB_EXTERN double stb_euclidean_distance(const double *a, const double *b, const int size);

/* Most of the time, this is sufficient and slightly faster calculate the sqr euclidean distance between two points */
STB_EXTERN double stb_euclidean_distance_sqr(const double *a, const double *b, const int size);

/* K-Means++ data clustering
 * x[n][d] = the data points
 * n       = Number of points
 * d       = Dimension of the data (e.g. color, weight, etc.) 
 * k       = # clusters
*  c[k][d] = Center points of clusters
*  z[n]    = What cluster a point is in
*  wss[k]  = The within-cluster sum of square of each cluster (optional)
*
* Note: This algorithm does not scale very well to a large number of points, consider using k-Means||
*/
STB_EXTERN void stb_kmeans(double **x, int n, int d, int k, double ***cret, int **zret, double **wssret);

/* Compute eigenvalues and eigenvectors of a symmetric matrix
 * a[n][n] = The matrix
 * n       = Order of a
 * w[n]    = Eigenvalues
 * z[n][n] = Eigenvectors
 */
STB_EXTERN int stb_eigenv (double **a, int n, double **wret, double ***zret);

/* This function returns the median */
STB_EXTERN double stb_median(double *data, int n);

/* stb_pca Principal Component Analysis
 * x[n][p]    = Data matrix
 * nx[n][p]   = 1 if data exists this point, 0 otherwise (optional (no missing data), can be NULL)
 * n          = Number of objects
 * p          = Variables each object
 * weights[p] = Weight of each variable (optional, can be NULL)
 * eret[p]    = Eigenvalues of covariance matrix
 * vret[p]    = Eigenvectors of covariance matrix
 * rret[n][m] = Projected data
 * m          = # of dimensions to project
 * level      = Level of robustness:
 *      -1 => flimsy statistics, Chebyshev codeviation
 *       0 => regular statistics, covariance matrix
 *       1 => semi-robust statistics, Manahattan codeviation
 *       2 => robust statistics, comedian matrix
 */
STB_EXTERN void stb_pca(double **x, int n, int p, int **nx, double *weights, int m, int level, double **eret, double ***vret, double ***rret);

/* Arrange the N elements of data in random order. */
STB_EXTERN void stb_shuffle(void *data, size_t n, size_t size, uint64_t *seed);

/* returns a with n non-repeating random integers in [0,high). This 
 * is useful when n is small, but the range of numbers is big. Otherwise
 * it might be best to use stb_shuffle
 */
STB_EXTERN int *stb_unique_random(int n, int high, uint64_t *seed);

/* stb_neugas, a data neural gas clustering algorithm
 * See: www.demogng.de/JavaPaper/node16.html
 *
 * x[n][p] = Data Matrix
 * n       = Number of objects
 * p       = Measurements per object
 * k       = Number of clusters
 * c[k][p] = Cluster centers
 * z[n]    = What cluster a point is in (optional)
 * wss[k]  = The within-cluster sum of square of each cluster (optional only possible in combination with z!)
 *
 * Note: Neural gas was developed with a focus on learning a representation of the data space, rather than partitioning a data set
 */
STB_EXTERN void stb_neugas(double **x, int n, int p, int k, double ***cret, int **zret, double **wssret);

static  double  horner[] = {1.6666666666666666e-01, 3.3333333333333333e-02,
                         2.3809523809523809e-02, 3.3333333333333333e-02, 7.5757575757575757e-02,
                         2.5311355311355311e-01, 1.1666666666666667e+00, 7.0921568627450980e+00,
                         5.4971177944862155e+01, 5.2912424242424242e+02
                      };

#define EHORNER (9) /* the elements number of `horner[]` */

const double stb_nan = 0.0/0.0;
const double stb_pos_inf = 1.0 /0.0;
const double stb_neg_inf = -1.0/0.0;

/* stb_polygamma return the polygamma function {\psi}^k(x)
 * of variable x If k=0,1,2,..., then returns digamma, trigamma,
 * tetragamma,... function values respectively.
 * Special cases:
 *  polygamma(k, x) is NaN with signal if k < 0 or x < 0;
 *  polygamma(k, x) is INF with signal if x = 0;
 *  polygamma(k, +-Inf) is NaN with signal;
 *  polygamma(k, NaN) is that NaN with no signal.
 * k
 * x
 */
STB_EXTERN double stb_polygamma(int k, double x);

// Solve trigamma(y) = x for y
// Newton's method
// 1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
// so iteration to solve 1/x = 1/trigamma is monotonically convergent
STB_EXTERN double stb_trigamma_inverse(double x);

/* Moment estimation of the parameters of a scaled F-distribution (the prior) 
 * The first degrees of freedom is given 
 */
STB_EXTERN void stb_fit_f_dist(double *var, int len, int df1, double *pvar, double *pdf2);

STB_EXTERN double stb_incgamma(double S, double Z);

STB_EXTERN long double stb_log_incgamma(long double S, long double Z);

STB_EXTERN long double stb_gamma(double N);

STB_EXTERN long double stb_beta(double a, double b);

STB_EXTERN long double stb_log_gamma(double N);

/* The normalized incomplete beta function. */
STB_EXTERN long double stb_incbeta(double a, double b, double x);

/* The function dbinom returns the value of the probability density function (pdf) 
 * of the binomial distribution given a certain random variable x, number of trials
 * (size) and probability of success on each trial (prob) 
 */
STB_EXTERN double stb_pdf_binom(double x, double size, double prob);

/* Returns the value of the Poisson probability density function. In other
 * words, the dpois function finds the probability that a certain number 
 * of successes (x) occur based on an average rate of success (lambda)
 */
STB_EXTERN double stb_pdf_pois(double x, double lambda);

/* Returns the probability of obtaining x from a hypergeometric distribution 
 * with parameters N, n, k 
 * [[https://stattrek.com/probability-distributions/hypergeometric.aspx]]
 * Suppose we randomly select 5 cards without replacement from an ordinary deck of 
 * playing cards. What is the probability of getting exactly 2 red cards (i.e., hearts or diamonds)?
 *
 * N = 52; since there are 52 cards in a deck.
 * k = 26; since there are 26 red cards in a deck.
 * n = 5; since we randomly select 5 cards from the deck.
 * x = 2; since 2 of the cards we select are red.
 */      
STB_EXTERN double stb_pdf_hypgeo(int x, int N, int n, int k);

// Hash function: Murmer One At A Time 32 bit
STB_EXTERN uint32_t inline stb_murmer32(const char *key);

// Hash function: Murmer One At A Time 64 bit
STB_EXTERN uint64_t inline stb_murmer64(const char *key);

/* Calculate the Shannon Index (a.k.a. Shannon's diversity index, the Shannon–Wiener index, 
 * the Shannon–Weaver index, the Shannon entropy). Pilou evenness compares the actual diversity
 * value (such as the Shannon Index, H′) to the maximum possible diversity value (when all species
 * are equally common, Hmax=ln S where S is the total number of species). 
 * The higher the Shannon index, the more diverse the species are in the habitat.
 * Evenness gives you a value between 0 and 1. The lower the evenness, the higher the diversity.
 */
STB_EXTERN void stb_shannon(double *data, size_t n, double *index, double *evenness);

/* Simpson's Diversity Index is a measure of diversity which takes into account the number of species
 * present, as well as the relative abundance of each species. As species richness and evenness 
 * increase, so diversity increases. The value ranges between 0 and 1. One represents infinite diversity
 * and 0, no diversity.
 */
STB_EXTERN void stb_simpson(double *data, size_t n, double *index);

/* The Jaccard similarity index (sometimes called the Jaccard similarity coefficient) compares members 
 * for two sets to see which members are shared and which are distinct. It’s a measure of similarity for 
 * the two sets of data, with a range from 0% to 100%. The higher the percentage, the more similar the 
 * two populations. Although it’s easy to interpret, it is extremely sensitive to small samples sizes and 
 * may give erroneous results, especially with very small samples or data sets with missing observations.
 */
STB_EXTERN double stb_jaccard(char **setA, size_t n_setA, char **setB, size_t n_setB);

/* The Bray–Curtis dissimilarity is a measure used to quantify the compositional dissimilarity between two 
 * different sites, based on counts (c_setA and c_setB) at each site. The Bray–Curtis dissimilarity is bounded between 0 and 1, 
 * where 0 means the two sites have the same composition and 1 means the two sites do not share any species.
 */
STB_EXTERN double stb_bray_curtis(char **setA, double *c_setA, size_t n_setA, char **setB, double *c_setB, size_t n_setB);

#ifdef STB_STATS_DEFINE

/* The incomplete gamma function (Thanks Jacob Wells) */
double stb_incgamma(double S, double Z)
{
	if(Z < 0.0) {
		return 0.0;
	}

	long double Sc = (1.0 / S);
	Sc *= powl(Z, S);
	Sc *= expl(-Z);

	long double Sum = 1.0;
	long double numerator = 1.0;
	long double denominator = 1.0;

	for(int I = 0; I < 200; I++) { // 200
		numerator *= Z;
		S++;
		denominator *= S;
		Sum += (numerator / denominator);
	}

	return Sum * Sc;
}

/* Returns the Natural Logarithm of the Incomplete Gamma Function */
long double stb_log_incgamma(long double S, long double Z)
{
	if(Z < 0.0) {
		return 0.0;
	}

	long double Sc;
	Sc = (logl(Z) * S) - Z - logl(S);

	long double Sum = 1.0;
	long double numerator = 1.0;
	long double denominator = 1.0;

	for (int i = 0; i < 1000; i++) { // Loops for 1000 iterations
		numerator *= Z;
		S++;
		denominator *= S;
		Sum += (numerator / denominator);
	}

	return logl(Sum) + Sc;
}

#define ACCURACY 15 // 15
/*
 * ACCURACY is the level of accuracy you wish to calculate.
 * Spouge's Approximation is slightly tricky, as you
 * can only reach the desired level of precision, if
 * you have EXTRA precision available so that it can
 * build up to the desired level.
 *
 * If you're using double (64 bit wide datatype), you
 * will need to set A to 11, as well as remember to
 * change the math functions to the regular
 * (i.e. pow() instead of powl())
 *
 * double A = 11
 * long double A = 15
 * 
 * !! IF YOU GO OVER OR UNDER THESE VALUES YOU WILL !!!
 *               !!! LOSE PRECISION !!!
 */

/* The gamma function using Spouge's Approximation */
long double stb_gamma(double N)
{
	const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;

	long double Z = (long double)N;
	long double Sc = powl((Z + ACCURACY), (Z + 0.5));
	Sc *= expl(-1.0 * (Z + ACCURACY));
	Sc /= Z;

	long double F = 1.0;
	long double Ck;
	long double Sum = SQRT2PI;


	for(int K = 1; K < ACCURACY; K++) {
		Z++;
		Ck = powl(ACCURACY - K, K - 0.5);
		Ck *= expl(ACCURACY - K);
		Ck /= F;

		Sum += (Ck / Z);

		F *= (-1.0 * K);
	}

	return (long double)(Sum * Sc);
}

long double stb_log_gamma(double N)
{
	/* The constant SQRT2PI is defined as sqrt(2.0 * PI);
	 * For speed the constant is already defined in decimal
	 * form.  However, if you wish to ensure that you achieve
	 *  maximum precision on your own machine, you can calculate
	 *  it yourself using (sqrt(atan(1.0) * 8.0)):
	 *
	 *  const long double SQRT2PI = sqrtl(atanl(1.0) * 8.0);
	 */

	const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;

	long double Z = (long double)N;
	long double Sc;

	Sc = (logl(Z + ACCURACY) * (Z + 0.5)) - (Z + ACCURACY) - logl(Z);

	long double F = 1.0;
	long double Ck;
	long double Sum = SQRT2PI;


	for(int K = 1; K < ACCURACY; K++) {
		Z++;
		Ck = powl(ACCURACY - K, K - 0.5);
		Ck *= expl(ACCURACY - K);
		Ck /= F;

		Sum += (Ck / Z);

		F *= (-1.0 * K);
	}

	return logl(Sum) + Sc;
}

long double stb_beta(double a, double b)
{
	return stb_gamma(a) * stb_gamma(b) / stb_gamma(a + b);
}

#define EPSILON 1.0e-8
#define MINIMAL 1.0e-30

/* Helper function to calulcate the regularized incomplete beta function. */
long double stb_beta_contfrac(double a, double b, double x)
{
	long double cf = 0.0;
	long double coeff = 0.0;
	long double delta = 0.0;

	long double m = 1.0;
	long double numerator = 1.0;
	long double denominator = 1.0 - (a + b) * x / (a + 1.0);

	if (fabsl(denominator) < MINIMAL) {
		denominator = MINIMAL;
	}

	denominator = 1.0 / denominator;
	cf = denominator;

	/* 1024 = 512 "real" iterations! */
	for (int i = 0; i <= 1024; i++) {

		if (!(i % 2)) {
			coeff = (m * (b - m) * x) / (((a - 1.0) + 2.0 * m) * (a + 2.0 * m));
		} else {
			coeff = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1.0));
			m++; /* m should only be increased once every two cycles! */
		}

		/* Do an iteration of Lentz's algorithm. */
		denominator = 1.0 + coeff * denominator;
		numerator = 1.0 + coeff / numerator;

		if (fabsl(denominator) < MINIMAL) {
			denominator = MINIMAL;
		}

		if (fabsl(numerator) < MINIMAL) {
			numerator = MINIMAL;
		}

		denominator = 1.0 / denominator;
		delta = denominator * numerator;
		cf *= delta;

		if (fabsl(delta - 1.0) < EPSILON) {
			/* Found it */
			break;
		}
	}

	if (fabsl(delta - 1.0) > EPSILON) {
		fprintf(stderr, "stb_incbeta: did not converge, more loops needed!\n");
		return 1.0 / 0.0;
	}

	return cf;
}

/* The normalized incomplete beta function. */
long double stb_incbeta(double a, double b, double x)
{
	long double preval;

	if (x < 0.0 || x > 1.0) {
		fprintf(stderr, "stb_incbeta: x should be between 0.0 and 1.0!\n");
		return 1.0 / 0.0; /* return NAN */
	}

	preval = -logl(stb_beta(a, b)) + a * logl(x) + b * logl(1.0 - x);

	if (x > (a + 1.0) / (a + b + 2.0)) {
		/* Apply continued fraction */
		return expl(preval) * stb_beta_contfrac(a, b, x) / a;
	} else {
		/* The incomplete beta is symmetrical, so:
		 * Ix(a,b) == I1-x(b,a)
		 */
		return 1.0 - stb_incbeta(b, a, 1.0 - x);
	}
}

/* stb_polygamma return the polygamma function {\psi}^k(x)
 * of variable x If k=0,1,2,..., then returns digamma, trigamma,
 * tetragamma,... function values respectively.
 * Special cases:
 *  polygamma(k, x) is NaN with signal if k < 0 or x < 0;
 *  polygamma(k, x) is INF with signal if x = 0;
 *  polygamma(k, +-Inf) is NaN with signal;
 *  polygamma(k, NaN) is that NaN with no signal.
 * k
 * x
 */
double stb_polygamma(int k, double x)
{
    double  s;  /* return value */
    double  y;  /* minimum value more than `large_value', adding `x' to integers */
    double  x_sqr;  /* x * x */
    double  fact_k; /* k! */
    double  pow_x_k;    /* pow_x_k = pow(x, k+1)    */
    double  large_value;    /* sufficient large value applied for asymptotic expansion */
    double  f;
    int n;  /* [large_value - x] */
    int i, j;
    int i2, isgn;

    if (k < 0 || x < 0.0 || !finite(x)) {   /* k < 0 or x < 0 or x is neither infinite nor a "not-a-number" (NaN) */
        /*
         * DOMAIN error: polygamma(x) return NaN
         */
        return stb_nan;
    } else if (k > 3) {
        /*
         * calculation of `large_value'
         */
        f = 1.0;
        for (i = k + 19; i > 20; i--) {
            f *= (double) i;
        }
        for (i = k + 1; i > 2; i--) {
            f /= (double) i;
        }
        f *= (174611.0 / 55.0); /* B_{20} / B_{2} */
        large_value = 6.812921 * pow(f, 1.0 / 18.0);
        if (large_value < 13.06) {
            large_value = 13.06;
        }
    } else {    /* 0 <= k <= 3 */
        large_value = 13.06;
    }

    /* fact_k = k! */
    fact_k = stb_factorial(k);

    if (x == 0.0) {
        /*
         * SING error: polygamma(x) return infinity
         */
        return stb_pos_inf;
    } else if (x >= large_value) {
        /* Adopted `x' to the asymptotic expansion. */
        s = 0.0;
        x_sqr = stb_sqr(x);
        isgn = k % 2 ? -1 : 1;
        if (k == 0) {
            /* digamma function */
            for (i = EHORNER; i >= 0; i--) {
                i2 = 2 * (i + 1);
                s += horner[i] / (double) i2 * isgn;
                s /= x_sqr;
                isgn *= -1;
            }
            s += log(x) - 0.5 / x;
        } else {
            /* k >= 1; trigamm, tetragamma, ... */
            for (i = EHORNER; i >= 0; i--) {
                f = 1.0;
                i2 = 2 * (i + 1);
                j = i2 + k - 1;
                while (j > i2) {
                    f *= (double) j--;
                }
                s += horner[i] * f * isgn;
                s /= x_sqr;
                isgn *= -1;
            }
            for (i = 0; i < k; i++) {
                s /= x;
            }
            pow_x_k = 1.0;
            for (i = 0; i < k; i++) {
                pow_x_k *= x;    /* pow_x_k = pow(x, k) */
            }

            s -= fact_k * 0.5 / pow_x_k / x * isgn;
            f = fact_k / (double) k;
            s -= f / pow_x_k * isgn;
        }
    } else {
        /*
         * x < large_value;
         * Adopted `y' instead of `x' to the asymptotic expansion,
         * we calculation the value.
         */
        n = (int)(large_value - x);
        y = (double) n + x + 1.0;
        s = stb_polygamma(k, y);
        isgn = k % 2 ? 1 : -1;
        for (i = 0; i <= n; i++) {
            y -= 1.0;
            if (fabs(y) < 1.e-3) {
                if (x > 0) {
                    y = x - (double)((int)(x + 0.5));
                } else {
                    y = x - (double)((int)(x - 0.5));
                }
            }
            pow_x_k = 1.0;
            for (j = 0; j < k; j++) {
                pow_x_k *= y;    /* pow_x_k = pow(y, k) */
            }
            s += isgn *  fact_k / pow_x_k / y;
        }
    }
    return (s);
}


// Solve trigamma(y) = x for y
// Newton's method
// 1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
// so iteration to solve 1/x = 1/trigamma is monotonically convergent
double stb_trigamma_inverse(double x)
{
    double y = 0.5 + 1.0 / x;
    int i = 0;
    double tri, dif;
    while(1) {
        i++;
        tri = stb_polygamma(1, y);
        dif = tri * ( 1 - tri / x) / stb_polygamma(2, y);
        y = y + dif;

        if ((-dif / y) < 1e-8) {
            break;
        }

        if (i > 50) {
            printf("stb_trigamma_inverse: Iteration limit exceeded\n");
            break;
        }
    }

    return y;
}

#define PI2 6.283185307179586476925286
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508 /* 1/1188 */

static double sfe[16] = {0, 0.081061466795327258219670264, 0.041340695955409294093822081, 0.0276779256849983391487892927,
                         0.020790672103765093111522771, 0.0166446911898211921631948653, 0.013876128823070747998745727,
                         0.0118967099458917700950557241, 0.010411265261972096497478567, 0.0092554621827127329177286366,
                         0.008330563433362871256469318, 0.0075736754879518407949720242, 0.006942840107209529865664152,
                         0.0064089941880042070684396310, 0.005951370112758847735624416, 0.00555473355196280137103868999
                        };

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n ) */
double stirlerr(double n)
{
    double nn;
    if (n < 16) {
        return (sfe[(int)n]);
    }
    nn = (double)n;
    nn = nn * nn;
    if (n > 500) {
        return ((S0 - S1 / nn) / n);
    }
    if (n > 80) {
        return ((S0 - (S1 - S2 / nn) / nn) / n);
    }
    if (n > 35) {
        return ((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
    }
    return ((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
}

/* Evaluate the deviance termbd0(x,np) = x log(x/np) + np - x*/
double bd0(double x,double np)
{
    double ej, s, s1, v;
    int j;
    if (fabs(x - np) < 0.1 * (x + np)) {
        s = (x - np) * (x - np) / (x + np);
        v = (x - np) / (x + np);
        ej = 2 * x * v;
        for (j = 1; ; j++) {
            ej *= v * v;
            s1 = s + ej / (2 * j + 1);
            if (s1 == s) {
                return (s1);
            }
            s = s1;
        }
    }
    return (x * log(x / np) + np - x);
}

/* The function dbinom returns the value of the probability density function (pdf) 
 * of the binomial distribution given a certain random variable x, number of trials
 * (size) and probability of success on each trial (prob) 
 */
double stb_pdf_binom(double x, double size, double prob)
{
    double lc;
    if (prob == 0.0) {
        return ((x == 0) ? 1.0 : 0.0);
    }
    if (prob == 1.0) {
        return ((x == size) ? 1.0 : 0.0);
    }
    if (x == 0) {
        return (exp(size * log(1 - prob)));
    }
    if (x == size) {
        return (exp(size * log(prob)));
    }
    lc = stirlerr(size) - stirlerr(x) - stirlerr(size - x)
         - bd0(x, size * prob) - bd0(size - x, size * (1.0 - prob));
    return (exp(lc) * sqrt(size / (PI2 * x * (size - x))));
}

/* Returns the value of the Poisson probability density function. In other
 * words, the dpois function finds the probability that a certain number 
 * of successes (x) occur based on an average rate of success (lambda)
 */
double stb_pdf_pois(double x, double lambda)
{
    if (lambda == 0) {
        return ((x == 0) ? 1.0 : 0.0);
    }
    if (x == 0) {
        return (exp(-lambda));
    }
    return (exp(-stirlerr(x) - bd0(x, lambda)) / sqrt(PI2 * x));
}

/* Returns the probability of obtaining x from a hypergeometric distribution 
 * with parameters N, n, k 
 * [[https://stattrek.com/probability-distributions/hypergeometric.aspx]]
 * Suppose we randomly select 5 cards without replacement from an ordinary deck of 
 * playing cards. What is the probability of getting exactly 2 red cards (i.e., hearts or diamonds)?
 *
 * n = 26; since there are 26 (52-26) NON red cards in a deck.
 * N = 26; since there are 26 red cards in a deck.
 * k = 5; since we randomly select 5 cards from the deck.
 * x = 2; since 2 of the cards we select are red.
 * Pval = 0.325130
 *
 * AND
 * Suppose we select 5 cards from an ordinary deck of playing cards. What is the probability of obtaining 2 or fewer hearts?
 * Solution: This is a hypergeometric experiment in which we know the following:
 * 
 * n = 52-13; since there are 52-13 NON heart cards in a deck.
 * N = 13; since there are 13 hearts in a deck.
 * k = 5; since we randomly select 5 cards from the deck.
 * x = 0 to 3; since our selection includes 0, 1, or 2 hearts.

int main(int argc, char const *argv[])
{
    // PDF
    double pval = sstb_pdf_hypgeo(2, 26, 52-26, 5);
    printf("%lf\n", pval);

    // CDF
    pval = 0;
    for (int i = 0; i <= 2; i++) {
        pval+= sstb_pdf_hypgeo(i, 13, 52-13, 5);
    }
    printf("%lf\n", pval);
    return 0;
}
 */  
double sstb_pdf_hypgeo(int x, int N, int n, int k)
{
    // The number of combinations of x out of N
    double c1 = stb_log_factorial(N) - (stb_log_factorial(x) + stb_log_factorial(N - x));
     // The number of combinations of k-x out of n
    double c2 = stb_log_factorial(n) - (stb_log_factorial(k-x) + stb_log_factorial(n - (k-x)));
    // The number of combinations of k out of N+n
    double c3 = stb_log_factorial(N+n) - (stb_log_factorial(k) + stb_log_factorial((N+n) - k));

    return exp(c1 + c2 - c3);
}

// PDF
void stb_pdf_gumbel(double x, double mu, double sig, double *p)
{
	*p = 1 / sig * exp(- (x - mu) / sig) * exp(-exp(- (x - mu) / sig));
}

// CDF
void stb_cdf_gumbel(double x, double mu, double sig, double *p)
{
	*p = exp(-exp(- (x - mu) / sig));
}

// inverse CDF a.k.a quantile function
double stb_icdf_gumbel(double mu, double sig, double p)
{
	return mu - sig * log( -log(p));
}

/* Fit the Gumbel distribution using maximum likelihood estimation (MLE)
	//Usage example:
	double test[53] = {312,590,248,670,365,770,465,545,315,115,232,260,655,675,
	                   455,1020,700,570,853,395,926,99,680,121,976,916,921,191,
	                   187,377,128,582,744,710,520,672,645,655,918,512,255,1126,
	                   1386,1394,600,950,731,700,1407,1284,165,1496,809};

	//estimate mean and var
	stb_est_gumbel(test, 53, &mu, &sig);
	// Use the estimate parameters to determine the probability of location 430
	stb_pdf_gumbel(430, mu, sig, &p);

	Depending on the epsilon used the values should be similar to the ones below 
	mu ~471
	s  ~298
*/ 
void stb_est_gumbel(double *data, int n, double *mu, double *sig)
{
	double m = 0.0, s = 0.0;
	double mean, var;
	stb_meanvar(data, n, &mean, &var);

	double threshold = EPSILON;

	double ll_old = 0.0;
	m = mean;
	s = var;
	for (;;) {
		// Estimate m and s
		double s_num = 0.0, s_denom = 0.0;
		for (int i = 0; i < n; i++) {
			s_num += data[i] * exp(- data[i] / s);
			s_denom += exp(- data[i] / s);
		}
		s = mean - s_num / s_denom;
		double m_sum = 0.0;
		for (int i = 0; i < n; i++) {
			m_sum += exp(- data[i] / s);
		}
		m = -s  * log( (1.0 / (double) n) * m_sum);

		// Log likelihood
		double ll_sum1 = 0.0, ll_sum2 = 0.0;
		for (int i = 0; i < n; i++) {
			ll_sum1 += ((data[i] - m) / s);
			ll_sum2 += exp(-((data[i] - m) / s));
		}
		double ll_new = -n * log(s) - ll_sum1 - ll_sum2;
		if (ll_old) {
			if (fabs(ll_old - ll_new) < threshold) {
				break;
			}
			//printf("delta: %lf\n", fabs(ll_old - ll_new));
			ll_old = ll_new;
		} else {
			ll_old = ll_new;
		}
	}

	*mu = m;
	*sig = s;
}

/* This function computes the cumulative distribution functions P(x)for the t-distribution with degrees_of_freedom degrees of freedom */
double stb_cdf_student_t(double t, double degrees_of_freedom)
{
	double x = (t + sqrt(t * t + degrees_of_freedom)) / (2.0 * sqrt(t * t + degrees_of_freedom));
	double probability = stb_incbeta(degrees_of_freedom / 2.0, degrees_of_freedom / 2.0, x);

	if(isnan(probability) || isinf(probability) || probability <= 1e-8) {
		return 1e-14;
	}

	return probability;
}

/*This function computes the cumulative distribution functions P(x)for the f-distribution with degrees_of_freedom1 and 2 degrees of freedom 
 * NOTE degrees_of_freedom1 = from the population with the largest variance !
 */
double stb_cdf_f_distribution(double f, double degrees_of_freedom1, double degrees_of_freedom2)
{
	/* The incomplete beta integral is used, according to the
	* formula
	*
	*	P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).
	*/
	double x = degrees_of_freedom2 / (degrees_of_freedom2 + degrees_of_freedom1 * f);
	double probability = stb_incbeta(degrees_of_freedom2 / 2.0, degrees_of_freedom1 / 2.0, x);

	if(isnan(probability) || isinf(probability) || probability <= 1e-8) {
		return 1e-14;
	}

	return probability;
}

/* This function computes the cumulative distribution functions P(x)for the chisqr-distribution with degrees_of_freedom degrees of freedom */
double stb_cdf_chisqr(double CV, double degrees_of_freedom)
{
	double probability = stb_incgamma(degrees_of_freedom / 2.0, CV / 2.0);

	if (degrees_of_freedom == 2) {
		return exp(-1.0 * CV / 2.0);
	}

	if(isnan(probability) || isinf(probability) || probability <= 1e-8) {
		return 1e-14;
	}

	probability /= stb_gamma(degrees_of_freedom / 2.0);

	return (1.0 - probability);

	// long double PValue, Gam;
	// long double ln_PV;
	// ln_PV = stb_log_incgamma(degrees_of_freedom/2.0, CV/2.0);

	// //Gam = approx_gamma(K);
	/* Using C99 lgammal */
	// Gam = lgammal(K);
	// //Gam = stb_log_gamma(K);

	// ln_PV -= Gam;
	// probability = 1.0 - expl(ln_PV);

	// return (double)probability;

}

/* This function evaluates the Gamma distribution CDF.
 * a	the gamma shape parameter \alpha
 * b	the gamma shape parameter \theta
 * x	the value at which to evaluate the distribution, must not be negative
 * upper return upper (if 1) or lower (if 0)
 */
double stb_cdf_gamma(double a, double b, double x, int upper)
{
	assert(x >= 0);

	double probability = stb_incgamma(b, a * x) / stb_gamma(b);
	
	return upper ? probability : 1 - probability;
}

/* Given  an  array  of data (length n), this routine  returns  its  mean and standard  deviation standard_dev (optional can be NULL) calculated using the Welford’s method. */
void stb_meanvar(double *data, int n, double *mean, double *sample_variance)
{
	int i;
	*mean = 0.0;
	double delta = 0.0;
	double delta2 = 0.0;
	double sq_sum = 0.0;

	if (n <= 1) {
		fprintf(stderr, "stb_stdev: At least 2 values are needed to calculate the mean and standard deviation\n");
		exit(1);
	}

	for (i = 0; i < n; i++) {
		delta = data[i] - *mean;
		*mean += delta / (i + 1);
		delta2 = data[i] - *mean;
		sq_sum += delta * delta2;
	}

	/* We could also return (population_variance = sq_sum / n) instead */
	if (sample_variance) {
        *sample_variance = sq_sum / (n - 1);
    }
}

/* Given  two arrays  of data (length n1 and n2), this routine  returns  its  t value and its significance (p) as probabiilty
 * NOTE: assumes equal variance 
 */
void stb_ttest(double *data1, int n1, double *data2, int n2, double *t, double *p)
{
	double mean1;
	double mean2;
	double var1;
	double var2;
	double degrees_of_freedom;
	double combined_variance;

	stb_meanvar(data1, n1, &mean1, &var1);
	stb_meanvar(data2, n2, &mean2, &var2);

	degrees_of_freedom = n1 + n2 - 2;

	combined_variance = ((n1 - 1) * var1 + (n2 - 1) * var2) / degrees_of_freedom;
	*t = (mean1 - mean2) / sqrt(combined_variance * (1.0 / n1 + 1.0 / n2));
	*t = fabs(*t);

	/* one sided ! multiply by two for two sided */
	*p = stb_cdf_student_t(-1 * (*t), degrees_of_freedom);
}

/* Given  two arrays  of data (length n1 and n2), this routine  returns  its  t value and its significance (p) as probabiilty
 * NOTE: assumes unequal variance 
 */
void stb_uttest(double *data1, int n1, double *data2, int n2, double *t, double *p)
{
	double mean1;
	double mean2;
	double var1;
	double var2;
	double degrees_of_freedom;

	stb_meanvar(data1, n1, &mean1, &var1);
	stb_meanvar(data2, n2, &mean2, &var2);

	degrees_of_freedom = ( ((var1 / n1) + (var2 / n2)) * ((var1 / n1) + (var2 / n2))) / ((var1 * var1) / (n1 * n1 * (n1 - 1)) + (var2 * var2) / (n2 * n2 * (n2 - 1)));
	*t = (mean1 - mean2) / sqrt((var1 / n1) + (var2 / n2));
	*t = fabs(*t);

	/* one sided! multiply by two for two sided */
	*p = stb_cdf_student_t(-1 * (*t), degrees_of_freedom);
}

/* Given  two arrays  of data (length n1 and n2), this routine  returns  its  f value and its significance (p) as probabiilty */
void stb_ftest(double *data1, int n1, double *data2, int n2, double *f, double *p)
{
	double mean1;
	double mean2;
	double var1;
	double var2;
	double degrees_of_freedom1;
	double degrees_of_freedom2;

	/* Step 1: If you are given standard deviations, go to Step 2. If you are given variances to compare, go to Step 3. */
	stb_meanvar(data1, n1, &mean1, &var1);
	stb_meanvar(data2, n2, &mean2, &var2);

	/* Step 3:Take the largest variance, and divide it by the smallest variance to get the f-value. For example, if your two variances were s1 = 2.5 and s2 = 9.4, divide 9.4 / 2.5 = 3.76.
	 Why? Placing the largest variance on top will force the F-test into a right tailed test, which is much easier to calculate than a left-tailed test.
	 */
	if (var1 > var2) {
		*f = var1 / var2;
		degrees_of_freedom1 = n1 - 1;
		degrees_of_freedom2 = n2 - 1;
	} else {
		*f = var2 / var1;
		/* This is swiched on purpose! See below */
		degrees_of_freedom1 = n2 - 1;
		degrees_of_freedom2 = n1 - 1;
	}

	/* two sided */
	/* NOTE degrees_of_freedom1 = from the population with the largest variance */
	*p = 2.0 * stb_cdf_f_distribution(*f, degrees_of_freedom1, degrees_of_freedom2);

	if (*p > 1.0)
	{ *p = 2.0 - *p; }

}

/* As an alternative to stb_phi, we could use erfc from C99! */
/*
double stb_phi(double x)
{
     return 0.5 * erfc(-x * M_SQRT1_2);
}
*/

/* This function takes two arrays, one with the data (of total size n) and one specifying the group (0, 1, 2, ...) (with g groups) the data belongs to (total size n) and
 * performs an one-way anova and followed by the Tukey HSD test and Scheffe's method of mutliple comparisons */
void stb_anova(double *data, int n, int *groups, int g, double *f, double *p)
{
	/* The overall mean */
	double Y = 0.0;

	/* The mean within each group */
	double *Yi = calloc(g, sizeof(double));
	/* The number of samples per group */
	int *a =  calloc(g, sizeof(int));

	/* Step 1: Calculate the mean within each group: */
	for (int i = 0; i < n; i++) {
		int index = groups[i];
		a[index]++;
		Yi[index] += data[i];
	}

	/* Step 2: Calculate the overall mean: */
	for (int i = 0; i < g; i++) {
		Yi[i] /= (double) a[i];
		Y += Yi[i];
	}

	Y /= (double) g;

	/* Step 3: Calculate the "between-group" sum of squared differences: */
	double sb = 0.0;

	for (int i = 0; i < g; i++) {
		sb += (double) a[i] * stb_sqr(Y - Yi[i]);
	}

	/* between-group degrees of freedom */
	double fb = (double) g - 1.0;

	/* between-group mean square value */
	double MSb = sb / fb;

	/* within-group sum of squares */
	double sw = 0.0;

	/* Step 4: Calculate the "within-group" sum of squares. Begin by centering the data in each group */
	for (int i = 0; i < n; i++) {
		int index = groups[i];
		sw += stb_sqr(data[i] - Yi[index]);
	}

	/* within-group degrees of freedom */
	double fw = (double) n - (double) g;
	/* within-group mean square */
	double MSw = sw / fw;

	*f = MSb / MSw;

	*p = stb_cdf_f_distribution(*f, fb, fw);

	double Q;
	/* Harmonic mean */
	double Hm;

	/* Tukey Honestly Significant Difference */
	for (int i = 0; i < g - 1; i++) {
		for (int j = i + 1; j < g; j++) {
			Hm = 2.0 / ((1.0 / (double) a[i]) + (1 / (double) a[j]));
			Q = fabs(Yi[i] - Yi[j]) / sqrt(MSw / Hm);
			printf("%i vs %i: %f", i, j, Q);
			/* Scheffé multiple comparison */
			double T;
			T = Q / sqrt(2);
			printf(" (Scheffe P: %f)\n", stb_cdf_f_distribution(stb_sqr(T) / ((double) g - 1.0), g - 1.0, fw));
		}
	}

	free(Yi);
	free(a);
}

/* The function Φ(x) is the cumulative density function (CDF) of a standard normal (Gaussian) random variable
 * Thanks John D. Cook 
 */
double stb_phi(double x)
{
	/* constants */
	double a1 =  0.254829592;
	double a2 = -0.284496736;
	double a3 =  1.421413741;
	double a4 = -1.453152027;
	double a5 =  1.061405429;
	double p  =  0.3275911;

	/* Save the sign of x */
	int sign = 1;

	if (x < 0)
	{ sign = -1; }

	x = fabs(x) / sqrt(2.0);

	/* A&S formula 7.1.26 */
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

	return 0.5 * (1.0 + sign * y);
}

/* Given an array of data, this function tests if the data is normally distributed
 * and returns the A value and the p-value. If the p-value is small (<0.05) the data is
 * not normally distributed! Care should be taken that the number of samples is not too
 * small, since this can lead to false positives/negatives (sample size > 7)
 * Example test data:
 * w should give the following results -> A = 0.1699 P-value = 0.88
 * double w[6]={294.2,308.5,313.1,317.7,322.7,338.7};
 * v should give the following results -> A = 0.3891 Adjusted_A = 0.4195 P-value = 0.3263
 * double v[12]={32.2,32.3,33.1,33.2,33.3,34.5,35.2,35.3,36.5,36.8,37.0,37.6};
 */
void stb_anderson_darling(double *data, int n, double *A, double *p)
{
	/* To perform the Anderson-Darling test, the array needs to be sorted */
	double *sorted_data = (double *) calloc(n, sizeof(double));
	double *standardized_data = (double *) calloc(n, sizeof(double));

	for(int i = 0; i < n; i++) {
		sorted_data[i] = data[i];
	}

	for(int i = 0; i < n; i++) {
		for(int j = i + 1; j < n; j++) {
			if(sorted_data[i] > sorted_data[j]) {
				double aux = sorted_data[i];
				sorted_data[i] = sorted_data[j];
				sorted_data[j] = aux;
			}
		}
	}

	/* Find the mean and variance of the data */
	double mean, var, stdev;
	stb_meanvar(sorted_data, n, &mean, &var);
	stdev = sqrt(var);

	/* The values are normalized/standardized */
	for(int i = 0; i < n; i++) {
		standardized_data[i] = (sorted_data[i] - mean) / stdev;
	}

	/* With a standard normal CDF (Phi), we calculate the Anderson_Darling Statistic */
	double temp_A = -n;

	for(int i = 0; i < n; i++) {
		temp_A +=  -1.0 / (double)n * (2 * (i + 1) - 1) * (log(stb_phi(standardized_data[i])) + log(1 - stb_phi(standardized_data[n - 1 - i])));
	}

	/* Adjust for small sample size (small <= 5) */
	*A = temp_A * (1 + .75 / (double)n + 2.25 / ((double)n * (double)n));

	/* Calculate the p-value according to https://www.spcforexcel.com/knowledge/basic-statistics/anderson-darling-test-for-normality */
	if (*A > 153.467) {
		*p = 0;
	} else if (*A > 0.6) {
		*p = exp(1.2937 - 5.709 * (*A) + 0.0186 * (*A) * (*A));
	} else if (*A > 0.34) {
		*p = exp(0.9177 - 4.279 * (*A) - 1.38 * (*A) * (*A));
	} else if (*A > 0.2) {
		*p = 1 - exp(-8.318 + 42.796 * (*A) - 59.938 * (*A) * (*A));
	} else {
		*p = 1 - exp(-13.436 + 101.14 * (*A) - 223.73 * (*A) * (*A));
	} 
	
	free(sorted_data);
	free(standardized_data);
}

void stb_2sample_anderson_darling(double *data1, int n1, double *data2, int n2, double *A, double *p)
{
	/* The number of groups */
	int g;
	/* The total number of values */
	int n;
	n = n1 + n2;
	g = 2;

	double *data;
	data = calloc(n , sizeof(double));
	int *groups;
	groups = calloc(n , sizeof(int));
	for (int i = 0; i < n1; i++) {
		data[i] = data1[i];
		groups[i] = 0;
	}
	for (int i = n1; i < n; i++) {
		data[i] = data2[i - n1];
		groups[i] = 1;
	}

	stb_ksample_anderson_darling(data, n, groups, g, A, p);

	free(data);
	free(groups);
}

/* Anderson-Darling k-sample test.
 *
 * Number of samples:  3
 * Sample sizes: 5 6 5
 * Total number of values: 16
 * Number of unique values: 11
 *
 * Mean of Anderson-Darling Criterion: 2
 * Standard deviation of Anderson-Darling Criterion: 0.92837
 *
 * T.AD = (Anderson-Darling Criterion - mean)/sigma
 *
 * Null Hypothesis: All samples come from a common population.
 *
 *                      T.AD P-value extrapolation
 * not adj. for ties 1.41756 0.08956             0
 * adj. for ties     1.62856 0.07084             0
 * double data[16] = {1, 3, 2, 5, 7, 2, 8, 1, 6, 9, 4, 12,  5,  7,  9, 11};
 * int groups[16] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
 * int n = 16;
 * int g = 3;
 *
 * double A, p;
 * stb_ksample_anderson_darling(data, n, groups, g, &A, &p);
 * printf("A: %lf p: %lf\n", A, p);
 *
 */
void stb_ksample_anderson_darling(double *data, int n, int *groups, int g, double *A, double *p)
{
	/* To perform the Anderson-Darling Test, the array needs to be sorted */
	double *Zstar = (double *) calloc(n + 1, sizeof(double));
	double *Z = (double *) calloc(n + 1, sizeof(double));
	double *observations = (double *) calloc(g, sizeof(double));

	*A = 0.0;
	*p = 0.0;

	/* Copy in the data */
	for (int i = 0; i < n; i++) {
		Zstar[i] = data[i];
		Z[i] = data[i];
		observations[groups[i]]++;
	}

	/* Sort the data */
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if(Zstar[i] > Zstar[j]) {
				/* Sort the data */
				double aux = Zstar[i];
				Zstar[i] = Zstar[j];
				Zstar[j] = aux;
				/* Sort the data */
				aux = Z[i];
				Z[i] = Z[j];
				Z[j] = aux;
			}
		}
	}

	/* Remove duplicates */
	/* L is the number of unique values */
	int L = 1;
	for (int i = 1; i < n; i++) {
		/* If Equal, continue */
        if (!(Zstar[L-1] < Zstar[i])) {
        	continue;
        }
        /* Store the unique value */
        Zstar[L++] = Zstar[i];
    }

	/* Determine the observations in the ith sample coinciding with Zstar */
	double **fij = (double **)stb_allocmat(g, L, sizeof(double));
	double *temp;
	for (int i = 0; i < g; i++) {
		temp = calloc(observations[i], sizeof(double));
		
		int index = 0;
		for (int j = 0; j < n; j++) {
			if (groups[j] == i) {
				temp[index++] = data[j];
			}
		}

		/* Sort the data */
		for (int k = 0; k < observations[i]; k++) {
			for (int j = k + 1; j < observations[i]; j++) {
				if(temp[k] > temp[j]) {
					/* Sort the data */
					double aux = temp[k];
					temp[k] = temp[j];
					temp[j] = aux;
				}
			}
		}

		for (int j = 0; j < observations[i]; j++) {
			for (int k = 0; k < L; k++) {
				fij[i][k] += (fabs(temp[j] - Zstar[k]) < EPSILON);
			}
		}
		free(temp);
	}

	double inner = 0, Mij = 0, Bj = 0, Maij = 0, Baj = 0, lj = 0, A2k = 0;
	for (int i = 0; i < g; i++) {
		inner = 0;
		Mij = 0;
		Bj = 0;
		for (int j = 0; j < L; j++) {
			lj = 0;
			for (int k = 0; k < g; k++) {
				lj += fij[k][j];
			}
			Mij += fij[i][j];
			Bj += lj;
			Maij = Mij - fij[i][j] / 2.0;
			Baj = Bj - lj / 2.0;
			inner += (lj / (double) n) * stb_sqr((double) n * Maij - observations[i] * Baj) / (Baj * ((double) n - Baj) - (double) n * lj / 4.0);			
		}
		A2k += inner / observations[i];
	}
	A2k *= ((double) n - 1.0) / (double) n;

	double H = 0.0;
	for (int i = 1; i <= g; i++) {
		H += 1.0 / observations[i - 1];
	}

	double h = 0.0;
	for (int i = 1; i <= (n - 1); i++) {
		h += 1.0 / (double) i;
	}

	double G = 0.0;
	for (int i = 1; i <= (n - 2); i++) {
		for (int j = i + 1; j <= (n - 1); j++) {
			G += 1.0 / (((double) n - (double) i) * (double) j);
		}
	}

	double a = (4.0 * G - 6.0) * ((double) g - 1.0) + (10.0 - 6.0 * G) * H;
	double b = (2.0 * G - 4.0) * stb_sqr((double) g) + 8.0 * h * (double) g + (2.0 * G - 14.0 * h - 4.0) * H - 8.0 * h + 4.0 * G - 6.0;
	double c = (6.0 * h + 2.0 * G - 2.0 ) * stb_sqr((double) g) + (4.0 * h - 4.0 * G + 6.0) * (double) g + (2.0 * h - 6.0) * H + 4.0 * h;
	double d = (2.0 * h + 6.0) * stb_sqr((double) g) - 4.0 * h * (double) g;
	double var = ( a * pow((double) n, 3.0) + b * pow((double) n, 2) + c * (double) n + d) / (((double) n - 1.0) * ((double) n - 2.0) * ((double) n - 3.0));

	double mean = (double) g - 1;
	*A = (A2k - mean) / sqrt(var);

    double sig[5] = {0.25, 0.1, 0.05, 0.025, 0.01};
    double b0[5] = {0.675, 1.281, 1.645, 1.96, 2.326};
    double b1[5] = {-0.245, 0.25, 0.678, 1.149, 1.822};
    double b2[5] = {-0.105, -0.305, -0.362, -0.391, -0.396};
    
    /* Populate the table with critical values */
    double tm[3][5];
    for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 5; j++) {
    		tm[i][j] = b0[j] + b1[j] / sqrt(mean) + b2[j] / mean;
    	}
    }

    double logsig[5];
    for (int i = 0; i < 5; i++) {
    	logsig[i] = log(sig[i] / (1.0 - sig[i]));
    }

    /* Apply a polynomial fit */
    double *fit;
    fit = stb_polyfit(tm[0], logsig, 5, 2);

    *p = fit[0] * pow(*A, 0) + fit[1] * pow(*A, 1) + fit[2] * pow(*A, 2);
    *p = exp(*p) / (1 + exp(*p));

    free(fit);
	free(observations);
	free(fij);
	free(Zstar);
}

/* Perform the Mann-Whitney U test, the resulting corrected P-value should be the same as that of the R statistical program.
 * However, for tied data no exact P values is calculated. It is possible, but would require a more complicated program see
 * for example:
 * Alexander Marx, Christina Backes, Eckart Meese, Hans-Peter Lenhof, Andreas Keller,
 * EDISON-WMW: Exact Dynamic Programing Solution of the Wilcoxon–Mann–Whitney Test,
 * Genomics, Proteomics & Bioinformatics, Volume 14, Issue 1, 2016, Pages 55-61
 */
void stb_mann_whitney(double *data1, int n1, double *data2, int n2, double *U, double *p)
{
	double n = (double) n1 + n2;

	/* To perform the Mann-Whitney Test, the array needs to be sorted */
	double *sorted_data = (double *) calloc(n1 + n2 + 1, sizeof(double));
	/* Keeps which group the sorted data belongs to */
	char *groups = (char *) calloc(n1 + n2, sizeof(char));
	/* Keeps the rank */
	double *rank = (double *) calloc(n1 + n2, sizeof(double));

	/* Copy in the data */
	for (int i = 0; i < n1; i++) {
		sorted_data[i] = data1[i];
		groups[i] = 1;
	}

	for (int i = 0; i < n2; i++) {
		sorted_data[i + n1] = data2[i];
		groups[i + n1] = 2;
	}

	/* Sort the data */
	for (int i = 0; i < n1 + n2; i++) {
		for (int j = i + 1; j < n1 + n2; j++) {
			if(sorted_data[i] > sorted_data[j]) {
				/* Sort the data */
				double aux = sorted_data[i];
				sorted_data[i] = sorted_data[j];
				sorted_data[j] = aux;
				/* Take along the groups */
				char aux_group = groups[i];
				groups[i] = groups[j];
				groups[j] = aux_group;
			}
		}
	}

	/* Adding an additional value bigger then the biggest value observed simplifies the ranking */
	sorted_data[n1 + n2] = sorted_data[n1 + n2 - 1] + 1;

	/* Initialize the Rank */
	for (int i = 0; i < (n1 + n2); i++) {
		rank[i] = i + 1;
	}

	/* Resolve ties */
	double sum = 0;
	double ties = 0;
	/* count frequencies of the ties to correct the standard deviation for them later */
	int total_number_of_ties_observed = 0;
	double *ties_freq = NULL;

	/* There will never be more than (n1+n2+1)/2 ties (= every number has 1 duplicate) */
	if ((n1 + n2) % 2) {
		ties_freq = (double *) calloc(n1 + n2 + 1, sizeof(double));
	} else {
		ties_freq = (double *) calloc(n1 + n2, sizeof(double));
	}

	for (int i = 0; i < (n1 + n2); i++) {
		if (!(sorted_data[i] < sorted_data[i + 1])) {
			/* Equal */
			sum = rank[i];
			ties = 1;

			for (int j = i; !(sorted_data[j] < sorted_data[j + 1]) && (j < n1 + n2); j++) {
				sum += rank[j + 1];
				ties++;
			}

			for (int j = i; j < i + ties; j++) {
				rank[j] = (double) sum / (double) ties;
			}

			/* Skip over the ties */
			i += ties - 1;
			ties_freq[total_number_of_ties_observed] = ties;
			total_number_of_ties_observed++;
		}
	}

	/* Sum the ranks of data1 */
	sum = 0;

	for(int i = 0; i < n1 + n2; i++) {
		if (groups[i] == 1) {
			sum += rank[i];
		}
	}

	double r1 = sum;
	double u1 = r1 - (double) n1 * ((double) n1 + 1) / 2.0;
	double u2 = (double) n1 * (double) n2 - u1;

	if (u1 < u2) {
		*U = u1;
	} else {
		*U = u2;
	}

	double sum_ties = 0.0;

	for(int i = 0; i < total_number_of_ties_observed; i++) {
		sum_ties += (pow(ties_freq[i], 3) - ties_freq[i]) / 12.0;
	}

	/* Calculate the normal approximation p-value (only valid for sample size > 20 */
	double E = (double) n1 * (double) n2 / 2.0;
	double stdev = 0.0;
	double z = 0.0;

	/* UNcorrected Stdev -> Not corrected for ties! */
	// stdev = sqrt( ((double) n1 * (double) n2 * (n + 1.0)) / 12.0);

	/* Corrected stdev according to https://secure.brightstat.com/index.php?p=c&d=1&c=2&i=5 */
	stdev = sqrt(((double) n1 * (double) n2 / (n * (n - 1.0))) * (((pow(n, 3) - n) / 12) - sum_ties));
	z =  fabs(*U - E) / stdev;

	/* Note one sided so if testing for not equal multiply by 2! And a continuity correction is applied (0.5)*/
	*p = 1 - stb_phi(z + 0.5);

	/* two sided without continuity correction */
	//*p = 2*(1 - stb_phi(z));

	free(ties_freq);
	free(sorted_data);
	free(groups);
}

/* This function takes two arrays, one with the data (of total size n) and one specifying the group (0, 1, 2, ...) (with g groups) the data belongs to (total size n) and
 * performs an one-way Kruskal-Wallis test.
 * Note:the number of samples per groups should be >> 5, otherwise the distribution is not similar enough to the ChiSqr
 */
void stb_kruskal_wallis(double *data, int n, int *groups, int g, double *H, double *p)
{
	/* To perform the Mann-Whitney Test, the array needs to be sorted */
	double *sorted_data = (double *) calloc(n + 1, sizeof(double));
	/* Keeps which group the sorted data belongs to */
	int *sorted_groups = (int *) calloc(n, sizeof(int));
	/* Keeps the rank */
	double *rank = (double *) calloc(n, sizeof(double));

	double *observations = (double *) calloc(g, sizeof(double));

	/* Copy in the data */
	for (int i = 0; i < n; i++) {
		sorted_data[i] = data[i];
		sorted_groups[i] = groups[i];
		observations[groups[i]]++;
	}

	/* Sort the data */
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if(sorted_data[i] > sorted_data[j]) {
				/* Sort the data */
				double aux = sorted_data[i];
				sorted_data[i] = sorted_data[j];
				sorted_data[j] = aux;
				/* Take along the groups */
				char aux_group = sorted_groups[i];
				sorted_groups[i] = sorted_groups[j];
				sorted_groups[j] = aux_group;
			}
		}
	}

	/* Adding an additional value bigger then the biggest value observed simplifies the ranking */
	sorted_data[n] = sorted_data[n - 1] + 1;

	/* Initialize the Rank */
	for (int i = 0; i < n; i++) {
		rank[i] = i + 1;
	}

	/* Resolve ties */
	double sum = 0;
	double ties = 0;
	/* count frequencies of the ties to correct the standard deviation for them later */
	int total_number_of_ties_observed = 0;
	double *ties_freq = NULL;

	/* There will never be more than (n+1)/2 ties (= every number has 1 duplicate) */
	if (n % 2) {
		ties_freq = (double *) calloc(n + 1, sizeof(double));
	} else {
		ties_freq = (double *) calloc(n, sizeof(double));
	}

	for (int i = 0; i < n; i++) {
		if (!(sorted_data[i] < sorted_data[i + 1])) {
			/* Equal */
			sum = rank[i];
			ties = 1;

			for (int j = i; !(sorted_data[j] < sorted_data[j + 1]) && (j < n); j++) {
				sum += rank[j + 1];
				ties++;
			}

			for (int j = i; j < i + ties; j++) {
				rank[j] = (double) sum / (double) ties;
			}

			/* Skip over the ties */
			i += ties - 1;
			ties_freq[total_number_of_ties_observed] = ties;
			total_number_of_ties_observed++;
		}
	}

	/* Sum the ranks of data1 */
	double *sum_rank = calloc(g, sizeof(double));

	for(int i = 0; i < n; i++) {
		sum_rank[sorted_groups[i]] += rank[i];
	}

	*H = 0;
	for (int i = 0; i < g; i++) {
		*H += observations[i] * stb_sqr(sum_rank[i] / observations[i]);
	}
	*H *= (12 / ((double) n * ((double)n + 1.0)));
	*H -= 3 * ((double)n + 1);

	double sum_ties = 0.0;

	for(int i = 0; i < total_number_of_ties_observed; i++) {
		sum_ties += (pow(ties_freq[i], 3) - ties_freq[i]);
	}

	/* Correct for ties */
	*H /= (1 - sum_ties / (pow(n, 3) - n));

	/* Note this is only valid for >> 5 samples per group */
	*p = stb_cdf_chisqr(fabs(*H), g - 1);

	/* Dunn's test */
	double sigma = 0.0;
	for (int i = 0; i < g - 1; i++) {
		for (int j = i + 1; j < g; j++) {
			double Ni = observations[i];
			double Nj = observations[j];
			sigma = sqrt( (((n * (n + 1)) - sum_ties / (n - 1)) / 12) * (1 / Ni + 1 / Nj));
			double z = fabs(((sum_rank[i] / Ni) - (sum_rank[j] / Nj)) / sigma);
			/* Two sided */
			printf("%i vs %i: %f (P = %lf) ", i, j, z, 2 * (1.0 - stb_phi(z)));
			printf("(Padj = %.4lf)\n", stb_bonferroni(2 * (1.0 - stb_phi(z)), ((double) g * ((double) g - 1.0)) / 2.0 ));
			/* One sided */
			// printf("%i vs %i: %f (P = %lf) ", i, j, z, (1.0 - stb_phi(z)));
			// printf("(Padj = %.4lf)\n", stb_bonferroni((1.0 - stb_phi(z)), ((double) g * ((double) g - 1.0)) / 2.0 ));
		}
	}

	free(ties_freq);
	free(sorted_data);
	free(sorted_groups);
	free(observations);
	free(sum_rank);
}

/* Given a rank (starting at 1 for the lowest p-value), the total number of comparisons performed and
 * the desired false discovery rate this function returns the Benjamini-Hochberg threshold. In other words,
 * if you put the individual P values in order, from smallest to largest. The smallest P value has a rank
 * 1, the next has 2, etc. Compare each individual P value to its Benjamini-Hochberg critical value (BHCV,
 * (rank/number_of_samples)*FDR. The largest P value that has P<BHCV is significant, and all of the P values
 * smaller than it are also significant, even the ones that aren't less than their Benjamini-Hochberg
 * critical value (see example).In general a FDR between 0.1 - 0.2 is reasonable.
 *
 * Dietary variable P value Rank BHCV
 * Total           <0.002   1    0.010
 * Carbohydrates    0.009   2    0.020
 * Fat              0.035   3    0.030
 * Sugars           0.044   4    0.040 <- Note that this value is higher than the BHCV but is
 * Proteins         0.046   5    0.050 <- with a FDR of 0.25 everyting BEFORE protein is siginificant also the sugars
 * Saturated        0.062   6    0.060
 */
double stb_benjamini_hochberg(int rank, int number_of_comparisons, double FDR)
{
	return ((double) rank / (double) number_of_comparisons) * FDR;
}

/* Returns the adjusted Bonferroni P-value */
double stb_bonferroni(double p, int number_of_comparisons)
{
	return (1.0 - pow(1.0 - p, (double) number_of_comparisons));
}

double stb_log_hypergeometric_prob(int a,int b,int c,int d) 
{
    return stb_log_factorial(a+b) + stb_log_factorial(c+d) + stb_log_factorial(a+c) + stb_log_factorial(b+d) - stb_log_factorial(a)- stb_log_factorial(b) - stb_log_factorial(c) - stb_log_factorial(d) - stb_log_factorial(a+b+c+d);
}

/* Calculates the two-tailed P-value for the Fisher Exact test.
 *                Men Women   Row total
 * Studying         1     9   10
 * Not-studying    11     3   14
 * Column total    12    12   24
 *
 * p = 2.759456e-03 two-sided
 */
void stb_fisher2x2(int a, int b, int c, int d, double *p)
{
    int n = a + b + c + d;
    double cutoff = stb_log_hypergeometric_prob(a,b,c,d);
    double Pvalue = 0;
    for(int x = 0; x <= n; x++) {
        if( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ) {
            double l = stb_log_hypergeometric_prob(x,a+b-x,a+c-x,d-a+x);
            if( l <= cutoff ) {
                Pvalue += exp(l);
            }
        }
    }
    
    *p = Pvalue;
}

/* Given a vector with observed and expected values with the amount of rows anc columns, this function
 * performs the chisquare test and returns */
void stb_chisqr(double *observed, double *expected, int rows, int columns, double *CV, double *p)
{
	*CV = 0.0;
	double degrees_of_freedom, sqr;

	for (int n = 0; n < rows; n++) {
		for (int i = 0; i < columns; i++) {
			sqr = observed[n * columns + i] - expected[n * columns + i];
			*CV += (sqr * sqr) / expected[n * columns + i];
		}
	}

	degrees_of_freedom = (rows - 1) * (columns - 1);

	if (degrees_of_freedom == 0) {
		degrees_of_freedom = 1;
	}

	*p = stb_cdf_chisqr(*CV, degrees_of_freedom);
}

/* Given a matrix with observed and expected values with the amount of rows anc columns, this function
 * performs the chisquare test and returns 
 */
void stb_chisqr_matrix(double **observed, double **expected, int rows, int columns, double *CV, double *p)
{
	*CV = 0.0;
	double degrees_of_freedom, sqr;

	for (int n = 0; n < rows; n++) {
		for (int i = 0; i < columns; i++) {
			sqr = observed[n][i] - expected[n][i];
			*CV += (sqr * sqr) / expected[n][i];
		}
	}

	degrees_of_freedom = (rows - 1) * (columns - 1);

	if (degrees_of_freedom == 0) {
		degrees_of_freedom = 1;
	}

	*p = stb_cdf_chisqr(*CV, degrees_of_freedom);
}

/* Given a vector with observed and expected values with the amount of rows anc columns, this function
 * performs the G-test and returns the probability (based on the ChiSquared distribution)
 */
void stb_gtest(double *observed, double *expected, int rows, int columns, double *CV, double *p)
{
	*CV = 0.0;
	double degrees_of_freedom;

	for (int n = 0; n < rows; n++) {
		for (int i = 0; i < columns; i++) {
			if (!observed[n * columns + i] || !expected[n * columns + i]) {
				continue;
			} else {
				*CV += observed[n * columns + i] * log(observed[n * columns + i] / expected[n * columns + i]);
			}
		}
	}

	*CV *= 2;

	degrees_of_freedom = (rows - 1) * (columns - 1);

	if (degrees_of_freedom == 0) {
		degrees_of_freedom = 1;
	}

	*p = stb_cdf_chisqr(*CV, degrees_of_freedom);
}

/* Given a matrix with observed and expected values with the amount of rows and columns, this function
 * performs the G-test and returns the probability (based on the ChiSquared distribution)
 */
void stb_gtest_matrix(double **observed, double **expected, int rows, int columns, double *CV, double *p)
{
	*CV = 0.0;
	double degrees_of_freedom;

	for (int n = 0; n < rows; n++) {
		for (int i = 0; i < columns; i++) {
			if (!observed[n][i] || !expected[n][i]) {
				continue;
			} else {
				*CV += observed[n][i] * log(observed[n][i] / expected[n][i]);
			}
		}
	}

	*CV *= 2;

	degrees_of_freedom = (rows - 1) * (columns - 1);

	if (degrees_of_freedom == 0) {
		degrees_of_freedom = 1;
	}

	*p = stb_cdf_chisqr(*CV, degrees_of_freedom);
}

/* This function returns the min , q1, q2, q3 and max values according to https://en.wikipedia.org/wiki/Quartile */
void stb_quartiles(double *data, int n, double *min, double *q1, double *median, double *q3, double *max)
{
	if (n < 3) {
		fprintf(stderr, "Need at least 3 elements\n");
		exit(1);
	}

	/* Copy the data */
	double *sorted_data = (double *) calloc(n, sizeof(double));

	for(int i = 0; i < n; i++) {
		sorted_data[i] = data[i];
	}

	/* Sort it */
	for(int i = 0; i < n; i++) {
		for(int j = i + 1; j < n; j++) {
			if(sorted_data[i] > sorted_data[j]) {
				double aux = sorted_data[i];
				sorted_data[i] = sorted_data[j];
				sorted_data[j] = aux;
			}
		}
	}

	if (IS_ODD(n)) {
		int split = (int) floor((double) n / 2.0);
		*median = sorted_data[split];

		if ((n - 1) % 4 == 0) {
			/* If there are (4n+1) data points, then the:
			 * lower quartile is 25% of the nth data value plus 75% of the (n+1)th data value;
			 * the upper quartile is 75% of the (3n+1)th data point plus 25% of the (3n+2)th data point.
			 */
			int quarter = (n - 1) / 4;
			*q1 = (sorted_data[quarter - 1] * .25) + (sorted_data[quarter] * .75);
			*q3 = (sorted_data[3 * quarter] * .75) + (sorted_data[3 * quarter + 1] * .25);
		} else if ((n - 3) % 4 == 0) {
			/* If there are (4n+3) data points, then the:
			 * lower quartile is 75% of the (n+1)th data value plus 25% of the (n+2)th data value;
			 * the upper quartile is 25% of the (3n+2)th data point plus 75% of the (3n+3)th data point.
			 */
			int quarter = (n - 3) / 4;
			*q1 = (sorted_data[quarter] * .75) + (sorted_data[quarter + 1] * .25);
			*q3 = (sorted_data[3 * quarter + 1] * .25) + (sorted_data[3 * quarter + 2] * .75);
		}
	} else {
		/* If there are an even number of data points in the original ordered data set, split this data set exactly in half.
		 * The lower quartile value is the median of the lower half of the data. The upper quartile value is the median of the upper half of the data.
		 */
		int split = n / 2;
		*median = (sorted_data[split - 1] + sorted_data[split]) / 2.0;

		int middle = floor(split / 2);

		if (IS_ODD(middle)) {
			*q1 = sorted_data[middle];
			*q3 = sorted_data[split + middle];
		} else {
			*q1 = (sorted_data[middle - 1] + sorted_data[middle]) / 2;
			*q3 = (sorted_data[split + middle - 1] + sorted_data[split + middle]) / 2;
		}
	}

	*min = sorted_data[0];
	*max = sorted_data[n - 1];

	free(sorted_data);
}

/* Free a histogram */
void stb_free_histogram(struct stb_hist *hist)
{
	free(hist->count);
	free(hist);
}

/* Keeps a score histogram, in which scores are counted into bins 0 <= x < 1 */
struct stb_hist *stb_histogram(double *data, int n, double bin_width, double umin, double umax)
{
	struct stb_hist *hist = NULL;
	hist = calloc(1, sizeof(struct stb_hist));

	int bin;

	/* If there is no data yet, we need a umin, umax and bin_width provided by the user */
	if (n == 0) {
		if (umin < 0) {
			fprintf(stderr, "stb_histogram: cannot initialize an empty histogram without a minimal value\n");
			exit(1);
		}

		if (umax < 0) {
			fprintf(stderr, "stb_histogram: cannot initialize an empty histogram without a maximum value\n");
			exit(1);
		}

		if (bin_width < 0) {
			fprintf(stderr, "stb_histogram: cannot initialize an empty histogram without a bin width\n");
			exit(1);
		}

		hist->min = umin;
		hist->max = umax;
		hist->bin_width = bin_width;
		hist->number_of_bins = (int) ceil((hist->max - hist->min) / (hist->bin_width) + 2.0);
		hist->count = (double *)calloc(hist->number_of_bins, sizeof(double));
	} else if (n < 4) {
		if (bin_width < 0) {
			fprintf(stderr, "stb_histogram: need atleast 3 data points to initialize the histogram\n");
			exit(1);
		}
	} else {
		double min, max, q1, q2, q3;

		stb_quartiles(data, n, &min, &q1, &q2, &q3, &max);

		if (umin < 0) {
			hist->min = min;
		} else {
			hist->min = umin;
		}

		if (umax < 0) {
			hist->max = max;
		} else {
			hist->max = umax;
		}

		if (bin_width > 0) {
			hist->bin_width = bin_width;
			hist->number_of_bins = (int) ceil((hist->max - hist->min) / (hist->bin_width) + 2.0);
		} else if ((bin_width < 0) && (bin_width > -2)) {
			/* Calculate bin width using Freedman–Diaconis rule */
			hist->bin_width = 2 * (q3 - q1) / pow(n, 1 / 3);
			hist->number_of_bins = (int) ceil((hist->max - hist->min) / (hist->bin_width) + 2.0);
		} else if ((bin_width < -1) && (bin_width > -3)) {
			/* Calculate bin width using Rice rule */
			hist->bin_width = 2 * pow(n, 1 / 3);
			hist->number_of_bins = (int) ceil((hist->max - hist->min) / (hist->bin_width) + 2.0);
		}

		hist->count = (double *)calloc(hist->number_of_bins, sizeof(double));

		if (hist == NULL) {
			return hist;
		}

		for (int i = 0; i < n; ++i) {
			bin = (int) ceil( ( (data[i] - hist->min) / hist->bin_width )) + 1;

			if (bin < 1) {
				/* less */
				hist->count[0]++;
			} else if (bin > hist->number_of_bins - 1) {
				/* more */
				hist->count[hist->number_of_bins]++;
			} else {
				hist->count[bin]++;
			}
		}
	}

	return hist;
}

/* print the histogram */
void stb_print_histogram(struct stb_hist *hist)
{
	int i;
	printf("Less %.2lf\n", hist->count[0]);

	for (i = 1; i < hist->number_of_bins; i++) {
		printf("%.2lf %.2lf\n", hist->min + (i - 1) * hist->bin_width, hist->count[i]);
	}

	printf("More %.2lf\n", hist->count[i]);
}

/* add some data to the histogram */
void stb_histogram_add(double *data, int n, struct stb_hist *hist)
{
	int bin;

	for (int i = 0; i < n; ++i) {
		bin = (int) ceil( ( (data[i] - hist->min) / hist->bin_width )) + 1;

		if (bin < 1) {
			/* less */
			hist->count[0]++;
		} else if (bin > hist->number_of_bins - 1) {
			/* more */
			hist->count[hist->number_of_bins]++;
		} else {
			hist->count[bin]++;
		}
	}
}

/* Compute (numerical stable) sum of the data using the Neumaier summation algorithm */
double stb_sum(double *data, int n)
{
	double c = 0.0; /* A running compensation for lost low-order bits. */

	double sum = data[0];
	for (int i = 1; i < n; i++) {
		double t = sum + data[i];
		if (fabs(sum) >= fabs(data[i])) {
			/* If sum is bigger, low-order digits of input[i] are lost. */
			c += (sum - t) + data[i];
		} else {
			/* Else low-order digits of sum are lost */
			c += (data[i] - t) + sum;
		}
		sum = t;
	}

	return sum + c;
}

// Maximum factorial value present in the precalculated lookup table
#define STB_LOG_FACT_MAX (20000)
// Maximum factorial value in the precalculated lookup table that does not overflow
#define STB_FACT_MAX (170)
/* Precalculated lookup table for stb_factorial and stb_log_factorial generated using:

#ifdef DBL_DECIMAL_DIG
#define OP_DBL_Digs (DBL_DECIMAL_DIG)
#else
#ifdef DECIMAL_DIG
#define OP_DBL_Digs (DECIMAL_DIG)
#else
#define OP_DBL_Digs (DBL_DIG + 3)
#endif
#endif

double stb_log_factorial(int n)
{
    double factorial = 0;

    while (n > 1) {
        factorial = factorial + log(n);
        n = n - 1;
    }

    return factorial;
}

int main(int argc, char const *argv[])
{
    double logfac;
    printf("double log_fact_table[%i] = { 0.0000000000000000e+00", STB_LOG_FACT_MAX + 1);
    for (int i = 1; i <= STB_LOG_FACT_MAX; i++) {
        logfac = stb_log_factorial(i);
        printf(", %.*e", OP_DBL_Digs - 1, logfac);
    }
    printf("};\n");
    return 0;
}

*/
double log_fact_table[20001] = {0.0000000000000000e+00, 0.0000000000000000e+00, 6.9314718055994529e-01, 1.7917594692280550e+00, 3.1780538303479458e+00, 4.7874917427820458e+00, 6.5792512120101012e+00, 8.5251613610654147e+00, 1.0604602902745249e+01, 1.2801827480081469e+01, 1.5104412573075514e+01, 1.7502307845873883e+01, 1.9987214495661885e+01, 2.2552163853123425e+01, 2.5191221182738680e+01, 2.7899271383840887e+01, 3.0671860106080665e+01, 3.3505073450136891e+01, 3.6395445208033046e+01, 3.9339884187199495e+01, 4.2335616460753485e+01, 4.5380138898476908e+01, 4.8471181351835220e+01, 5.1606675567764370e+01, 5.4784729398112319e+01, 5.8003605222980518e+01, 6.1261701761002001e+01, 6.4557538627006323e+01, 6.7889743137181540e+01, 7.1257038967168015e+01, 7.4658236348830172e+01, 7.8092223553315321e+01, 8.1557959456115029e+01, 8.5054467017581530e+01, 8.8580827542197696e+01, 9.2136175603687093e+01, 9.5719694542143202e+01, 9.9330612454787428e+01, 1.0296819861451380e+02, 1.0663176026064347e+02, 1.1032063971475741e+02, 1.1403421178146172e+02, 1.1777188139974508e+02, 1.2153308151543864e+02, 1.2531727114935690e+02, 1.2912393363912724e+02, 1.3295257503561632e+02, 1.3680272263732635e+02, 1.4067392364823428e+02, 1.4456574394634487e+02, 1.4847776695177302e+02, 1.5240959258449737e+02, 1.5636083630307880e+02, 1.6033112821663090e+02, 1.6432011226319520e+02, 1.6832744544842765e+02, 1.7235279713916282e+02, 1.7639584840699735e+02, 1.8045629141754375e+02, 1.8453382886144948e+02, 1.8862817342367160e+02, 1.9273904728784490e+02, 1.9686618167288998e+02, 2.0100931639928152e+02, 2.0516819948264117e+02, 2.0934258675253682e+02, 2.1353224149456321e+02, 2.1773693411395419e+02, 2.2195644181913030e+02, 2.2619054832372754e+02, 2.3043904356577690e+02, 2.3470172344281821e+02, 2.3897838956183426e+02, 2.4326884900298265e+02, 2.4757291409618685e+02, 2.5189040220972316e+02, 2.5622113555000948e+02, 2.6056494097186322e+02, 2.6492164979855272e+02, 2.6929109765101981e+02, 2.7367312428569369e+02, 2.7806757344036612e+02, 2.8247429268763040e+02, 2.8689313329542699e+02, 2.9132395009427034e+02, 2.9576660135076065e+02, 3.0022094864701415e+02, 3.0468685676566872e+02, 3.0916419358014690e+02, 3.1365282994987905e+02, 3.1815263962020930e+02, 3.2266349912672615e+02, 3.2718528770377515e+02, 3.3171788719692842e+02, 3.3626118197919845e+02, 3.4081505887079902e+02, 3.4537940706226686e+02, 3.4995411804077025e+02, 3.5453908551944085e+02, 3.5913420536957534e+02, 3.6373937555556353e+02, 3.6835449607240480e+02, 3.7297946888568902e+02, 3.7761419787391867e+02, 3.8225858877306001e+02, 3.8691254912321750e+02, 3.9157598821732955e+02, 3.9624881705179149e+02, 4.0093094827891571e+02, 4.0562229616114485e+02, 4.1032277652693728e+02, 4.1503230672824964e+02, 4.1975080559954472e+02, 4.2447819341825704e+02, 4.2921439186665151e+02, 4.3395932399501476e+02, 4.3871291418612111e+02, 4.4347508812091888e+02, 4.4824577274538456e+02, 4.5302489623849607e+02, 4.5781238798127811e+02, 4.6260817852687489e+02, 4.6741219957160814e+02, 4.7222438392698058e+02, 4.7704466549258564e+02, 4.8187297922988790e+02, 4.8670926113683947e+02, 4.9155344822329801e+02, 4.9640547848721769e+02, 5.0126529089157935e+02, 5.0613282534203483e+02, 5.1100802266523601e+02, 5.1589082458782241e+02, 5.2078117371604412e+02, 5.2567901351599505e+02, 5.3058428829443346e+02, 5.3549694318016941e+02, 5.4041692410599762e+02, 5.4534417779115472e+02, 5.5027865172428551e+02, 5.5522029414689484e+02, 5.6016905403727299e+02, 5.6512488109487424e+02, 5.7008772572513408e+02, 5.7505753902471020e+02, 5.8003427276713080e+02, 5.8501787938883922e+02, 5.9000831197561786e+02, 5.9500552424938189e+02, 6.0000947055532731e+02, 6.0502010584942366e+02, 6.1003738568623862e+02, 6.1506126620708494e+02, 6.2009170412847743e+02, 6.2512865673089095e+02, 6.3017208184781020e+02, 6.3522193785505988e+02, 6.4027818366040822e+02, 6.4534077869343514e+02, 6.5040968289565535e+02, 6.5548485671088906e+02, 6.6056626107587351e+02, 6.6565385741110583e+02, 6.7074760761191260e+02, 6.7584747403973699e+02, 6.8095341951363753e+02, 6.8606540730199413e+02, 6.9118340111441091e+02, 6.9630736509381416e+02, 7.0143726380873727e+02, 7.0657306224578758e+02, 7.1171472580229022e+02, 7.1686222027910367e+02, 7.2201551187360144e+02, 7.2717456717281584e+02, 7.3233935314673931e+02, 7.3750983714177744e+02, 7.4268598687435122e+02, 7.4786777042464337e+02, 7.5305515623048393e+02, 7.5824811308137430e+02, 7.6344661011264009e+02, 7.6865061679971700e+02, 7.7386010295255835e+02, 7.7907503871016741e+02, 7.8429539453524569e+02, 7.8952114120895897e+02, 7.9475224982581346e+02, 7.9998869178864345e+02, 8.0523043880370312e+02, 8.1047746287586358e+02, 8.1572973630391027e+02, 8.2098723167593801e+02, 8.2624992186484303e+02, 8.3151778002390620e+02, 8.3679077958246990e+02, 8.4206889424170049e+02, 8.4735209797043854e+02, 8.5264036500113298e+02, 8.5793366982585746e+02, 8.6323198719240543e+02, 8.6853529210046452e+02, 8.7384355979786574e+02, 8.7915676577690749e+02, 8.8447488577075171e+02, 8.8979789574989002e+02, 8.9512577191867956e+02, 9.0045849071194493e+02, 9.0579602879164622e+02, 9.1113836304361109e+02, 9.1648547057432859e+02, 9.2183732870780455e+02, 9.2719391498247660e+02, 9.3255520714818601e+02, 9.3792118316320796e+02, 9.4329182119133554e+02, 9.4866709959901971e+02, 9.5404699695256022e+02, 9.5943149201534925e+02, 9.6482056374516571e+02, 9.7021419129151809e+02, 9.7561235399303598e+02, 9.8101503137490818e+02, 9.8642220314636836e+02, 9.9183384919822333e+02, 9.9724994960042784e+02, 1.0026704845997001e+03, 1.0080954346171816e+03, 1.0135247802461359e+03, 1.0189585022496901e+03, 1.0243965815586134e+03, 1.0298389992691350e+03, 1.0352857366408014e+03, 1.0407367750943674e+03, 1.0461920962097252e+03, 1.0516516817238694e+03, 1.0571155135288950e+03, 1.0625835736700303e+03, 1.0680558443437019e+03, 1.0735323078956333e+03, 1.0790129468189753e+03, 1.0844977437524660e+03, 1.0899866814786226e+03, 1.0954797429219630e+03, 1.1009769111472563e+03, 1.1064781693578009e+03, 1.1119835008937332e+03, 1.1174928892303612e+03, 1.1230063179765261e+03, 1.1285237708729910e+03, 1.1340452317908534e+03, 1.1395706847299853e+03, 1.1451001138174968e+03, 1.1506335033062242e+03, 1.1561708375732428e+03, 1.1617121011184013e+03, 1.1672572785628809e+03, 1.1728063546477761e+03, 1.1783593142326977e+03, 1.1839161422943971e+03, 1.1894768239254124e+03, 1.1950413443327352e+03, 1.2006096888364962e+03, 1.2061818428686738e+03, 1.2117577919718201e+03, 1.2173375217978064e+03, 1.2229210181065882e+03, 1.2285082667649885e+03, 1.2340992537454990e+03, 1.2396939651251012e+03, 1.2452923870840996e+03, 1.2508945059049790e+03, 1.2565003079712751e+03, 1.2621097797664600e+03, 1.2677229078728481e+03, 1.2733396789705146e+03, 1.2789600798362319e+03, 1.2845840973424192e+03, 1.2902117184561096e+03, 1.2958429302379309e+03, 1.3014777198411000e+03, 1.3071160745104335e+03, 1.3127579815813717e+03, 1.3184034284790152e+03, 1.3240524027171762e+03, 1.3297048918974449e+03, 1.3353608837082643e+03, 1.3410203659240242e+03, 1.3466833264041600e+03, 1.3523497530922721e+03, 1.3580196340152527e+03, 1.3636929572824242e+03, 1.3693697110846927e+03, 1.3750498836937097e+03, 1.3807334634610484e+03, 1.3864204388173882e+03, 1.3921107982717122e+03, 1.3978045304105149e+03, 1.4035016238970202e+03, 1.4092020674704113e+03, 1.4149058499450673e+03, 1.4206129602098163e+03, 1.4263233872271910e+03, 1.4320371200327004e+03, 1.4377541477341067e+03, 1.4434744595107143e+03, 1.4491980446126668e+03, 1.4549248923602538e+03, 1.4606549921432277e+03, 1.4663883334201253e+03, 1.4721249057176049e+03, 1.4778646986297838e+03, 1.4836077018175934e+03, 1.4893539050081336e+03, 1.4951032979940417e+03, 1.5008558706328674e+03, 1.5066116128464541e+03, 1.5123705146203315e+03, 1.5181325660031114e+03, 1.5238977571058963e+03, 1.5296660781016901e+03, 1.5354375192248201e+03, 1.5412120707703646e+03, 1.5469897230935871e+03, 1.5527704666093796e+03, 1.5585542917917094e+03, 1.5643411891730759e+03, 1.5701311493439732e+03, 1.5759241629523574e+03, 1.5817202207031228e+03, 1.5875193133575833e+03, 1.5933214317329600e+03, 1.5991265667018765e+03, 1.6049347091918571e+03, 1.6107458501848337e+03, 1.6165599807166586e+03, 1.6223770918766218e+03, 1.6281971748069741e+03, 1.6340202207024570e+03, 1.6398462208098374e+03, 1.6456751664274477e+03, 1.6515070489047307e+03, 1.6573418596417935e+03, 1.6631795900889595e+03, 1.6690202317463329e+03, 1.6748637761633643e+03, 1.6807102149384223e+03, 1.6865595397183690e+03, 1.6924117421981432e+03, 1.6982668141203458e+03, 1.7041247472748294e+03, 1.7099855334982951e+03, 1.7158491646738933e+03, 1.7217156327308264e+03, 1.7275849296439601e+03, 1.7334570474334355e+03, 1.7393319781642874e+03, 1.7452097139460673e+03, 1.7510902469324681e+03, 1.7569735693209561e+03, 1.7628596733524066e+03, 1.7687485513107395e+03, 1.7746401955225649e+03, 1.7805345983568297e+03, 1.7864317522244664e+03, 1.7923316495780489e+03, 1.7982342829114502e+03, 1.8041396447595050e+03, 1.8100477276976740e+03, 1.8159585243417143e+03, 1.8218720273473525e+03, 1.8277882294099600e+03, 1.8337071232642334e+03, 1.8396287016838771e+03, 1.8455529574812917e+03, 1.8514798835072620e+03, 1.8574094726506519e+03, 1.8633417178380998e+03, 1.8692766120337194e+03, 1.8752141482388017e+03, 1.8811543194915223e+03, 1.8870971188666490e+03, 1.8930425394752554e+03, 1.8989905744644359e+03, 1.9049412170170235e+03, 1.9108944603513114e+03, 1.9168502977207761e+03, 1.9228087224138060e+03, 1.9287697277534294e+03, 1.9347333070970478e+03, 1.9406994538361712e+03, 1.9466681613961566e+03, 1.9526394232359473e+03, 1.9586132328478166e+03, 1.9645895837571145e+03, 1.9705684695220154e+03, 1.9765498837332698e+03, 1.9825338200139568e+03, 1.9885202720192412e+03, 1.9945092334361311e+03, 2.0005006979832392e+03, 2.0064946594105456e+03, 2.0124911114991646e+03, 2.0184900480611113e+03, 2.0244914629390726e+03, 2.0304953500061793e+03, 2.0365017031657810e+03, 2.0425105163512237e+03, 2.0485217835256281e+03, 2.0545354986816710e+03, 2.0605516558413692e+03, 2.0665702490558651e+03, 2.0725912724052150e+03, 2.0786147199981760e+03, 2.0846405859720012e+03, 2.0906688644922319e+03, 2.0966995497524931e+03, 2.1027326359742924e+03, 2.1087681174068166e+03, 2.1148059883267388e+03, 2.1208462430380173e+03, 2.1268888758716998e+03, 2.1329338811857356e+03, 2.1389812533647814e+03, 2.1450309868200134e+03, 2.1510830759889373e+03, 2.1571375153352069e+03, 2.1631942993484354e+03, 2.1692534225440177e+03, 2.1753148794629456e+03, 2.1813786646716330e+03, 2.1874447727617367e+03, 2.1935131983499809e+03, 2.1995839360779842e+03, 2.2056569806120838e+03, 2.2117323266431727e+03, 2.2178099688865218e+03, 2.2238899020816175e+03, 2.2299721209919940e+03, 2.2360566204050692e+03, 2.2421433951319814e+03, 2.2482324400074281e+03, 2.2543237498895064e+03, 2.2604173196595516e+03, 2.2665131442219845e+03, 2.2726112185041502e+03, 2.2787115374561699e+03, 2.2848140960507840e+03, 2.2909188892831990e+03, 2.2970259121709414e+03, 2.3031351597537064e+03, 2.3092466270932086e+03, 2.3153603092730409e+03, 2.3214762013985242e+03, 2.3275942985965662e+03, 2.3337145960155167e+03, 2.3398370888250311e+03, 2.3459617722159255e+03, 2.3520886414000393e+03, 2.3582176916101007e+03, 2.3643489180995834e+03, 2.3704823161425797e+03, 2.3766178810336614e+03, 2.3827556080877480e+03, 2.3888954926399738e+03, 2.3950375300455612e+03, 2.4011817156796869e+03, 2.4073280449373560e+03, 2.4134765132332741e+03, 2.4196271160017200e+03, 2.4257798486964243e+03, 2.4319347067904405e+03, 2.4380916857760253e+03, 2.4442507811645182e+03, 2.4504119884862130e+03, 2.4565753032902480e+03, 2.4627407211444793e+03, 2.4689082376353676e+03, 2.4750778483678587e+03, 2.4812495489652702e+03, 2.4874233350691720e+03, 2.4935992023392778e+03, 2.4997771464533284e+03, 2.5059571631069807e+03, 2.5121392480136979e+03, 2.5183233969046355e+03, 2.5245096055285362e+03, 2.5306978696516185e+03, 2.5368881850574717e+03, 2.5430805475469469e+03, 2.5492749529380508e+03, 2.5554713970658454e+03, 2.5616698757823374e+03, 2.5678703849563799e+03, 2.5740729204735680e+03, 2.5802774782361366e+03, 2.5864840541628619e+03, 2.5926926441889582e+03, 2.5989032442659832e+03, 2.6051158503617348e+03, 2.6113304584601574e+03, 2.6175470645612418e+03, 2.6237656646809337e+03, 2.6299862548510332e+03, 2.6362088311191046e+03, 2.6424333895483805e+03, 2.6486599262176678e+03, 2.6548884372212592e+03, 2.6611189186688375e+03, 2.6673513666853883e+03, 2.6735857774111073e+03, 2.6798221470013109e+03, 2.6860604716263501e+03, 2.6923007474715209e+03, 2.6985429707369763e+03, 2.7047871376376397e+03, 2.7110332444031214e+03, 2.7172812872776294e+03, 2.7235312625198890e+03, 2.7297831664030541e+03, 2.7360369952146302e+03, 2.7422927452563836e+03, 2.7485504128442667e+03, 2.7548099943083312e+03, 2.7610714859926525e+03, 2.7673348842552437e+03, 2.7736001854679812e+03, 2.7798673860165231e+03, 2.7861364823002295e+03, 2.7924074707320879e+03, 2.7986803477386338e+03, 2.8049551097598760e+03, 2.8112317532492175e+03, 2.8175102746733833e+03, 2.8237906705123437e+03, 2.8300729372592396e+03, 2.8363570714203106e+03, 2.8426430695148197e+03, 2.8489309280749812e+03, 2.8552206436458900e+03, 2.8615122127854488e+03, 2.8678056320642950e+03, 2.8741008980657348e+03, 2.8803980073856692e+03, 2.8866969566325247e+03, 2.8929977424271879e+03, 2.8993003614029326e+03, 2.9056048102053546e+03, 2.9119110854923028e+03, 2.9182191839338120e+03, 2.9245291022120382e+03, 2.9308408370211914e+03, 2.9371543850674684e+03, 2.9434697430689912e+03, 2.9497869077557393e+03, 2.9561058758694853e+03, 2.9624266441637360e+03, 2.9687492094036634e+03, 2.9750735683660446e+03, 2.9813997178391996e+03, 2.9877276546229291e+03, 2.9940573755284518e+03, 3.0003888773783456e+03, 3.0067221570064853e+03, 3.0130572112579835e+03, 3.0193940369891297e+03, 3.0257326310673334e+03, 3.0320729903710612e+03, 3.0384151117897823e+03, 3.0447589922239085e+03, 3.0511046285847369e+03, 3.0574520177943928e+03, 3.0638011567857720e+03, 3.0701520425024873e+03, 3.0765046718988060e+03, 3.0828590419396032e+03, 3.0892151496002989e+03, 3.0955729918668067e+03, 3.1019325657354793e+03, 3.1082938682130525e+03, 3.1146568963165928e+03, 3.1210216470734445e+03, 3.1273881175211759e+03, 3.1337563047075264e+03, 3.1401262056903543e+03, 3.1464978175375859e+03, 3.1528711373271631e+03, 3.1592461621469911e+03, 3.1656228890948896e+03, 3.1720013152785414e+03, 3.1783814378154411e+03, 3.1847632538328476e+03, 3.1911467604677314e+03, 3.1975319548667294e+03, 3.2039188341860918e+03, 3.2103073955916379e+03, 3.2166976362587038e+03, 3.2230895533720964e+03, 3.2294831441260471e+03, 3.2358784057241619e+03, 3.2422753353793782e+03, 3.2486739303139134e+03, 3.2550741877592222e+03, 3.2614761049559497e+03, 3.2678796791538848e+03, 3.2742849076119155e+03, 3.2806917875979848e+03, 3.2871003163890446e+03, 3.2935104912710103e+03, 3.2999223095387201e+03, 3.3063357684958869e+03, 3.3127508654550588e+03, 3.3191675977375717e+03, 3.3255859626735078e+03, 3.3320059576016547e+03, 3.3384275798694616e+03, 3.3448508268329952e+03, 3.3512756958569007e+03, 3.3577021843143584e+03, 3.3641302895870426e+03, 3.3705600090650823e+03, 3.3769913401470153e+03, 3.3834242802397548e+03, 3.3898588267585424e+03, 3.3962949771269118e+03, 3.4027327287766484e+03, 3.4091720791477483e+03, 3.4156130256883812e+03, 3.4220555658548496e+03, 3.4284996971115502e+03, 3.4349454169309356e+03, 3.4413927227934769e+03, 3.4478416121876235e+03, 3.4542920826097679e+03, 3.4607441315642050e+03, 3.4671977565630978e+03, 3.4736529551264375e+03, 3.4801097247820098e+03, 3.4865680630653546e+03, 3.4930279675197321e+03, 3.4994894356960858e+03, 3.5059524651530064e+03, 3.5124170534566961e+03, 3.5188831981809335e+03, 3.5253508969070381e+03, 3.5318201472238334e+03, 3.5382909467276163e+03, 3.5447632930221171e+03, 3.5512371837184692e+03, 3.5577126164351730e+03, 3.5641895887980636e+03, 3.5706680984402715e+03, 3.5771481430021981e+03, 3.5836297201314742e+03, 3.5901128274829321e+03, 3.5965974627185669e+03, 3.6030836235075112e+03, 3.6095713075259955e+03, 3.6160605124573208e+03, 3.6225512359918230e+03, 3.6290434758268439e+03, 3.6355372296666947e+03, 3.6420324952226320e+03, 3.6485292702128177e+03, 3.6550275523622940e+03, 3.6615273394029500e+03, 3.6680286290734903e+03, 3.6745314191194061e+03, 3.6810357072929423e+03, 3.6875414913530708e+03, 3.6940487690654559e+03, 3.7005575382024281e+03, 3.7070677965429509e+03, 3.7135795418725952e+03, 3.7200927719835081e+03, 3.7266074846743800e+03, 3.7331236777504232e+03, 3.7396413490233358e+03, 3.7461604963112759e+03, 3.7526811174388349e+03, 3.7592032102370049e+03, 3.7657267725431543e+03, 3.7722518022009981e+03, 3.7787782970605690e+03, 3.7853062549781916e+03, 3.7918356738164539e+03, 3.7983665514441795e+03, 3.8048988857364020e+03, 3.8114326745743351e+03, 3.8179679158453491e+03, 3.8245046074429401e+03, 3.8310427472667079e+03, 3.8375823332223258e+03, 3.8441233632215158e+03, 3.8506658351820229e+03, 3.8572097470275876e+03, 3.8637550966879217e+03, 3.8703018820986822e+03, 3.8768501012014444e+03, 3.8833997519436784e+03, 3.8899508322787215e+03, 3.8965033401657556e+03, 3.9030572735697815e+03, 3.9096126304615923e+03, 3.9161694088177505e+03, 3.9227276066205627e+03, 3.9292872218580565e+03, 3.9358482525239529e+03, 3.9424106966176469e+03, 3.9489745521441787e+03, 3.9555398171142142e+03, 3.9621064895440172e+03, 3.9686745674554295e+03, 3.9752440488758439e+03, 3.9818149318381834e+03, 3.9883872143808776e+03, 3.9949608945478381e+03, 4.0015359703884378e+03, 4.0081124399574860e+03, 4.0146903013152069e+03, 4.0212695525272165e+03, 4.0278501916645014e+03, 4.0344322168033941e+03, 4.0410156260255526e+03, 4.0476004174179384e+03, 4.0541865890727931e+03, 4.0607741390876181e+03, 4.0673630655651518e+03, 4.0739533666133484e+03, 4.0805450403453569e+03, 4.0871380848794993e+03, 4.0937324983392491e+03, 4.1003282788532106e+03, 4.1069254245550974e+03, 4.1135239335837114e+03, 4.1201238040829248e+03, 4.1267250342016532e+03, 4.1333276220938424e+03, 4.1399315659184422e+03, 4.1465368638393911e+03, 4.1531435140255899e+03, 4.1597515146508858e+03, 4.1663608638940523e+03, 4.1729715599387700e+03, 4.1795836009736040e+03, 4.1861969851919830e+03, 4.1928117107921862e+03, 4.1994277759773195e+03, 4.2060451789552935e+03, 4.2126639179388112e+03, 4.2192839911453430e+03, 4.2259053967971067e+03, 4.2325281331210563e+03, 4.2391521983488556e+03, 4.2457775907168625e+03, 4.2524043084661116e+03, 4.2590323498422913e+03, 4.2656617130957293e+03, 4.2722923964813717e+03, 4.2789243982587677e+03, 4.2855577166920484e+03, 4.2921923500499097e+03, 4.2988282966055958e+03, 4.3054655546368804e+03, 4.3121041224260471e+03, 4.3187439982598744e+03, 4.3253851804296146e+03, 4.3320276672309819e+03, 4.3386714569641299e+03, 4.3453165479336367e+03, 4.3519629384484842e+03, 4.3586106268220483e+03, 4.3652596113720729e+03, 4.3719098904206603e+03, 4.3785614622942503e+03, 4.3852143253236027e+03, 4.3918684778437855e+03, 4.3985239181941524e+03, 4.4051806447183308e+03, 4.4118386557642016e+03, 4.4184979496838860e+03, 4.4251585248337251e+03, 4.4318203795742702e+03, 4.4384835122702616e+03, 4.4451479212906115e+03, 4.4518136050083940e+03, 4.4584805618008231e+03, 4.4651487900492411e+03, 4.4718182881390985e+03, 4.4784890544599448e+03, 4.4851610874054068e+03, 4.4918343853731740e+03, 4.4985089467649887e+03, 4.5051847699866230e+03, 4.5118618534478701e+03, 4.5185401955625248e+03, 4.5252197947483683e+03, 4.5319006494271580e+03, 4.5385827580246078e+03, 4.5452661189703740e+03, 4.5519507306980422e+03, 4.5586365916451105e+03, 4.5653237002529777e+03, 4.5720120549669255e+03, 4.5787016542361034e+03, 4.5853924965135229e+03, 4.5920845802560298e+03, 4.5987779039243005e+03, 4.6054724659828216e+03, 4.6121682648998813e+03, 4.6188652991475474e+03, 4.6255635672016624e+03, 4.6322630675418231e+03, 4.6389637986513708e+03, 4.6456657590173727e+03, 4.6523689471306134e+03, 4.6590733614855772e+03, 4.6657790005804372e+03, 4.6724858629170412e+03, 4.6791939470008938e+03, 4.6859032513411521e+03, 4.6926137744506050e+03, 4.6993255148456610e+03, 4.7060384710463386e+03, 4.7127526415762477e+03, 4.7194680249625826e+03, 4.7261846197361037e+03, 4.7329024244311267e+03, 4.7396214375855125e+03, 4.7463416577406470e+03, 4.7530630834414387e+03, 4.7597857132362942e+03, 4.7665095456771160e+03, 4.7732345793192826e+03, 4.7799608127216416e+03, 4.7866882444464918e+03, 4.7934168730595766e+03, 4.8001466971300670e+03, 4.8068777152305483e+03, 4.8136099259370158e+03, 4.8203433278288530e+03, 4.8270779194888264e+03, 4.8338136995030691e+03, 4.8405506664610712e+03, 4.8472888189556670e+03, 4.8540281555830234e+03, 4.8607686749426302e+03, 4.8675103756372819e+03, 4.8742532562730739e+03, 4.8809973154593845e+03, 4.8877425518088685e+03, 4.8944889639374423e+03, 4.9012365504642721e+03, 4.9079853100117643e+03, 4.9147352412055525e+03, 4.9214863426744896e+03, 4.9282386130506320e+03, 4.9349920509692292e+03, 4.9417466550687177e+03, 4.9485024239907016e+03, 4.9552593563799492e+03, 4.9620174508843766e+03, 4.9687767061550394e+03, 4.9755371208461229e+03, 4.9822986936149264e+03, 4.9890614231218578e+03, 4.9958253080304203e+03, 5.0025903470072008e+03, 5.0093565387218605e+03, 5.0161238818471256e+03, 5.0228923750587746e+03, 5.0296620170356264e+03, 5.0364328064595356e+03, 5.0432047420153749e+03, 5.0499778223910307e+03, 5.0567520462773882e+03, 5.0635274123683239e+03, 5.0703039193606955e+03, 5.0770815659543314e+03, 5.0838603508520164e+03, 5.0906402727594896e+03, 5.0974213303854249e+03, 5.1042035224414321e+03, 5.1109868476420361e+03, 5.1177713047046736e+03, 5.1245568923496812e+03, 5.1313436093002865e+03, 5.1381314542825958e+03, 5.1449204260255883e+03, 5.1517105232611020e+03, 5.1585017447238288e+03, 5.1652940891512999e+03, 5.1720875552838797e+03, 5.1788821418647558e+03, 5.1856778476399295e+03, 5.1924746713582053e+03, 5.1992726117711809e+03, 5.2060716676332395e+03, 5.2128718377015412e+03, 5.2196731207360126e+03, 5.2264755154993372e+03, 5.2332790207569442e+03, 5.2400836352770066e+03, 5.2468893578304233e+03, 5.2536961871908152e+03, 5.2605041221345155e+03, 5.2673131614405584e+03, 5.2741233038906730e+03, 5.2809345482692752e+03, 5.2877468933634536e+03, 5.2945603379629647e+03, 5.3013748808602240e+03, 5.3081905208502985e+03, 5.3150072567308925e+03, 5.3218250873023471e+03, 5.3286440113676217e+03, 5.3354640277322960e+03, 5.3422851352045527e+03, 5.3491073325951738e+03, 5.3559306187175289e+03, 5.3627549923875713e+03, 5.3695804524238274e+03, 5.3764069976473829e+03, 5.3832346268818865e+03, 5.3900633389535278e+03, 5.3968931326910406e+03, 5.4037240069256877e+03, 5.4105559604912532e+03, 5.4173889922240396e+03, 5.4242231009628531e+03, 5.4310582855490002e+03, 5.4378945448262766e+03, 5.4447318776409620e+03, 5.4515702828418098e+03, 5.4584097592800381e+03, 5.4652503058093271e+03, 5.4720919212858043e+03, 5.4789346045680431e+03, 5.4857783545170487e+03, 5.4926231699962564e+03, 5.4994690498715199e+03, 5.5063159930111051e+03, 5.5131639982856814e+03, 5.5200130645683139e+03, 5.5268631907344607e+03, 5.5337143756619544e+03, 5.5405666182310069e+03, 5.5474199173241932e+03, 5.5542742718264481e+03, 5.5611296806250584e+03, 5.5679861426096531e+03, 5.5748436566721975e+03, 5.5817022217069889e+03, 5.5885618366106419e+03, 5.5954225002820904e+03, 5.6022842116225711e+03, 5.6091469695356227e+03, 5.6160107729270767e+03, 5.6228756207050474e+03, 5.6297415117799319e+03, 5.6366084450643939e+03, 5.6434764194733634e+03, 5.6503454339240297e+03, 5.6572154873358277e+03, 5.6640865786304375e+03, 5.6709587067317771e+03, 5.6778318705659894e+03, 5.6847060690614426e+03, 5.6915813011487198e+03, 5.6984575657606110e+03, 5.7053348618321079e+03, 5.7122131883003995e+03, 5.7190925441048603e+03, 5.7259729281870459e+03, 5.7328543394906892e+03, 5.7397367769616867e+03, 5.7466202395480987e+03, 5.7535047262001417e+03, 5.7603902358701762e+03, 5.7672767675127079e+03, 5.7741643200843728e+03, 5.7810528925439376e+03, 5.7879424838522918e+03, 5.7948330929724389e+03, 5.8017247188694910e+03, 5.8086173605106624e+03, 5.8155110168652654e+03, 5.8224056869046990e+03, 5.8293013696024473e+03, 5.8361980639340691e+03, 5.8430957688771978e+03, 5.8499944834115277e+03, 5.8568942065188130e+03, 5.8637949371828590e+03, 5.8706966743895164e+03, 5.8775994171266739e+03, 5.8845031643842594e+03, 5.8914079151542219e+03, 5.8983136684305318e+03, 5.9052204232091817e+03, 5.9121281784881639e+03, 5.9190369332674800e+03, 5.9259466865491249e+03, 5.9328574373370857e+03, 5.9397691846373382e+03, 5.9466819274578320e+03, 5.9535956648084921e+03, 5.9605103957012107e+03, 5.9674261191498426e+03, 5.9743428341701956e+03, 5.9812605397800307e+03, 5.9881792349990510e+03, 5.9950989188488993e+03, 6.0020195903531485e+03, 6.0089412485372977e+03, 6.0158638924287743e+03, 6.0227875210569118e+03, 6.0297121334529620e+03, 6.0366377286500710e+03, 6.0435643056832942e+03, 6.0504918635895729e+03, 6.0574204014077368e+03, 6.0643499181785000e+03, 6.0712804129444521e+03, 6.0782118847500515e+03, 6.0851443326416238e+03, 6.0920777556673547e+03, 6.0990121528772825e+03, 6.1059475233232970e+03, 6.1128838660591318e+03, 6.1198211801403559e+03, 6.1267594646243733e+03, 6.1336987185704138e+03, 6.1406389410395341e+03, 6.1475801310946026e+03, 6.1545222878003024e+03, 6.1614654102231225e+03, 6.1684094974313521e+03, 6.1753545484950773e+03, 6.1823005624861762e+03, 6.1892475384783129e+03, 6.1961954755469269e+03, 6.2031443727692404e+03, 6.2100942292242416e+03, 6.2170450439926826e+03, 6.2239968161570823e+03, 6.2309495448017069e+03, 6.2379032290125779e+03, 6.2448578678774584e+03, 6.2518134604858551e+03, 6.2587700059290073e+03, 6.2657275032998841e+03, 6.2726859516931818e+03, 6.2796453502053155e+03, 6.2866056979344166e+03, 6.2935669939803274e+03, 6.3005292374445935e+03, 6.3074924274304640e+03, 6.3144565630428815e+03, 6.3214216433884831e+03, 6.3283876675755891e+03, 6.3353546347142028e+03, 6.3423225439160051e+03, 6.3492913942943469e+03, 6.3562611849642471e+03, 6.3632319150423918e+03, 6.3702035836471177e+03, 6.3771761898984187e+03, 6.3841497329179392e+03, 6.3911242118289638e+03, 6.3980996257564202e+03, 6.4050759738268671e+03, 6.4120532551684964e+03, 6.4190314689111274e+03, 6.4260106141861961e+03, 6.4329906901267577e+03, 6.4399716958674790e+03, 6.4469536305446354e+03, 6.4539364932961043e+03, 6.4609202832613628e+03, 6.4679049995814812e+03, 6.4748906413991208e+03, 6.4818772078585271e+03, 6.4888646981055272e+03, 6.4958531112875271e+03, 6.5028424465535008e+03, 6.5098327030539958e+03, 6.5168238799411165e+03, 6.5238159763685326e+03, 6.5308089914914653e+03, 6.5378029244666886e+03, 6.5447977744525206e+03, 6.5517935406088254e+03, 6.5587902220970027e+03, 6.5657878180799844e+03, 6.5727863277222350e+03, 6.5797857501897433e+03, 6.5867860846500180e+03, 6.5937873302720873e+03, 6.6007894862264911e+03, 6.6077925516852756e+03, 6.6147965258219983e+03, 6.6218014078117112e+03, 6.6288071968309641e+03, 6.6358138920578012e+03, 6.6428214926717528e+03, 6.6498299978538353e+03, 6.6568394067865438e+03, 6.6638497186538516e+03, 6.6708609326412025e+03, 6.6778730479355090e+03, 6.6848860637251482e+03, 6.6918999791999586e+03, 6.6989147935512337e+03, 6.7059305059717217e+03, 6.7129471156556156e+03, 6.7199646217985564e+03, 6.7269830235976260e+03, 6.7340023202513403e+03, 6.7410225109596522e+03, 6.7480435949239400e+03, 6.7550655713470123e+03, 6.7620884394330960e+03, 6.7691121983878347e+03, 6.7761368474182882e+03, 6.7831623857329269e+03, 6.7901888125416262e+03, 6.7972161270556662e+03, 6.8042443284877236e+03, 6.8112734160518739e+03, 6.8183033889635790e+03, 6.8253342464396956e+03, 6.8323659876984584e+03, 6.8393986119594865e+03, 6.8464321184437740e+03, 6.8534665063736911e+03, 6.8605017749729723e+03, 6.8675379234667225e+03, 6.8745749510814085e+03, 6.8816128570448554e+03, 6.8886516405862440e+03, 6.8956913009361051e+03, 6.9027318373263215e+03, 6.9097732489901173e+03, 6.9168155351620571e+03, 6.9238586950780455e+03, 6.9309027279753209e+03, 6.9379476330924508e+03, 6.9449934096693305e+03, 6.9520400569471785e+03, 6.9590875741685359e+03, 6.9661359605772568e+03, 6.9731852154185126e+03, 6.9802353379387823e+03, 6.9872863273858511e+03, 6.9943381830088065e+03, 7.0013909040580384e+03, 7.0084444897852318e+03, 7.0154989394433642e+03, 7.0225542522867036e+03, 7.0296104275708049e+03, 7.0366674645525018e+03, 7.0437253624899140e+03, 7.0507841206424328e+03, 7.0578437382707243e+03, 7.0649042146367237e+03, 7.0719655490036348e+03, 7.0790277406359210e+03, 7.0860907887993098e+03, 7.0931546927607815e+03, 7.1002194517885728e+03, 7.1072850651521694e+03, 7.1143515321223076e+03, 7.1214188519709614e+03, 7.1284870239713482e+03, 7.1355560473979267e+03, 7.1426259215263863e+03, 7.1496966456336459e+03, 7.1567682189978577e+03, 7.1638406408983956e+03, 7.1709139106158536e+03, 7.1779880274320521e+03, 7.1850629906300192e+03, 7.1921387994939960e+03, 7.1992154533094399e+03, 7.2062929513630106e+03, 7.2133712929425665e+03, 7.2204504773371764e+03, 7.2275305038370989e+03, 7.2346113717337894e+03, 7.2416930803198957e+03, 7.2487756288892506e+03, 7.2558590167368775e+03, 7.2629432431589739e+03, 7.2700283074529270e+03, 7.2771142089172927e+03, 7.2842009468518036e+03, 7.2912885205573612e+03, 7.2983769293360365e+03, 7.3054661724910638e+03, 7.3125562493268390e+03, 7.3196471591489189e+03, 7.3267389012640142e+03, 7.3338314749799893e+03, 7.3409248796058582e+03, 7.3480191144517821e+03, 7.3551141788290697e+03, 7.3622100720501676e+03, 7.3693067934286628e+03, 7.3764043422792784e+03, 7.3835027179178687e+03, 7.3906019196614216e+03, 7.3977019468280523e+03, 7.4048027987369960e+03, 7.4119044747086155e+03, 7.4190069740643912e+03, 7.4261102961269171e+03, 7.4332144402199056e+03, 7.4403194056681750e+03, 7.4474251917976562e+03, 7.4545317979353840e+03, 7.4616392234094947e+03, 7.4687474675492258e+03, 7.4758565296849138e+03, 7.4829664091479854e+03, 7.4900771052709633e+03, 7.4971886173874591e+03, 7.5043009448321700e+03, 7.5114140869408757e+03, 7.5185280430504426e+03, 7.5256428124988097e+03, 7.5327583946249943e+03, 7.5398747887690870e+03, 7.5469919942722518e+03, 7.5541100104767174e+03, 7.5612288367257797e+03, 7.5683484723637966e+03, 7.5754689167361885e+03, 7.5825901691894323e+03, 7.5897122290710622e+03, 7.5968350957296607e+03, 7.6039587685148645e+03, 7.6110832467773589e+03, 7.6182085298688698e+03, 7.6253346171421690e+03, 7.6324615079510686e+03, 7.6395892016504167e+03, 7.6467176975960965e+03, 7.6538469951450252e+03, 7.6609770936551513e+03, 7.6681079924854475e+03, 7.6752396909959143e+03, 7.6823721885475734e+03, 7.6895054845024697e+03, 7.6966395782236632e+03, 7.7037744690752288e+03, 7.7109101564222565e+03, 7.7180466396308466e+03, 7.7251839180681072e+03, 7.7323219911021524e+03, 7.7394608581020984e+03, 7.7466005184380638e+03, 7.7537409714811647e+03, 7.7608822166035152e+03, 7.7680242531782214e+03, 7.7751670805793819e+03, 7.7823106981820865e+03, 7.7894551053624073e+03, 7.7966003014974049e+03, 7.8037462859651196e+03, 7.8108930581445720e+03, 7.8180406174157606e+03, 7.8251889631596605e+03, 7.8323380947582191e+03, 7.8394880115943506e+03, 7.8466387130519433e+03, 7.8537901985158487e+03, 7.8609424673718813e+03, 7.8680955190068189e+03, 7.8752493528083969e+03, 7.8824039681653112e+03, 7.8895593644672072e+03, 7.8967155411046879e+03, 7.9038724974693032e+03, 7.9110302329535543e+03, 7.9181887469508829e+03, 7.9253480388556800e+03, 7.9325081080632754e+03, 7.9396689539699391e+03, 7.9468305759728801e+03, 7.9539929734702355e+03, 7.9611561458610822e+03, 7.9683200925454248e+03, 7.9754848129241964e+03, 7.9826503063992577e+03, 7.9898165723733910e+03, 7.9969836102503032e+03, 8.0041514194346191e+03, 8.0113199993318831e+03, 8.0184893493485542e+03, 8.0256594688920031e+03, 8.0328303573705152e+03, 8.0400020141932837e+03, 8.0471744387704075e+03, 8.0543476305128943e+03, 8.0615215888326511e+03, 8.0686963131424873e+03, 8.0758718028561116e+03, 8.0830480573881287e+03, 8.0902250761540390e+03, 8.0974028585702345e+03, 8.1045814040539981e+03, 8.1117607120235025e+03, 8.1189407818978052e+03, 8.1261216130968505e+03, 8.1333032050414622e+03, 8.1404855571533471e+03, 8.1476686688550899e+03, 8.1548525395701517e+03, 8.1620371687228690e+03, 8.1692225557384490e+03, 8.1764087000429718e+03, 8.1835956010633836e+03, 8.1907832582274978e+03, 8.1979716709639943e+03, 8.2051608387024153e+03, 8.2123507608731579e+03, 8.2195414369074933e+03, 8.2267328662375312e+03, 8.2339250482962452e+03, 8.2411179825174604e+03, 8.2483116683358530e+03, 8.2555061051869543e+03, 8.2627012925071340e+03, 8.2698972297336095e+03, 8.2770939163044459e+03, 8.2842913516585431e+03, 8.2914895352356452e+03, 8.2986884664763347e+03, 8.3058881448220236e+03, 8.3130885697149679e+03, 8.3202897405982494e+03, 8.3274916569157813e+03, 8.3346943181123061e+03, 8.3418977236333903e+03, 8.3491018729254265e+03, 8.3563067654356328e+03, 8.3635124006120423e+03, 8.3707187779035139e+03, 8.3779258967597234e+03, 8.3851337566311558e+03, 8.3923423569691167e+03, 8.3995516972257192e+03, 8.4067617768538894e+03, 8.4139725953073612e+03, 8.4211841520406742e+03, 8.4283964465091740e+03, 8.4356094781690081e+03, 8.4428232464771263e+03, 8.4500377508912788e+03, 8.4572529908700089e+03, 8.4644689658726602e+03, 8.4716856753593693e+03, 8.4789031187910659e+03, 8.4861212956294694e+03, 8.4933402053370883e+03, 8.5005598473772206e+03, 8.5077802212139432e+03, 8.5150013263121273e+03, 8.5222231621374158e+03, 8.5294457281562391e+03, 8.5366690238358024e+03, 8.5438930486440877e+03, 8.5511178020498573e+03, 8.5583432835226395e+03, 8.5655694925327389e+03, 8.5727964285512317e+03, 8.5800240910499597e+03, 8.5872524795015343e+03, 8.5944815933793270e+03, 8.6017114321574791e+03, 8.6089419953108900e+03, 8.6161732823152179e+03, 8.6234052926468812e+03, 8.6306380257830606e+03, 8.6378714812016842e+03, 8.6451056583814334e+03, 8.6523405568017461e+03, 8.6595761759428133e+03, 8.6668125152855682e+03, 8.6740495743116935e+03, 8.6812873525036157e+03, 8.6885258493445090e+03, 8.6957650643182897e+03, 8.7030049969096108e+03, 8.7102456466038657e+03, 8.7174870128871862e+03, 8.7247290952464446e+03, 8.7319718931692405e+03, 8.7392154061439051e+03, 8.7464596336595077e+03, 8.7537045752058439e+03, 8.7609502302734381e+03, 8.7681965983535410e+03, 8.7754436789381252e+03, 8.7826914715198945e+03, 8.7899399755922659e+03, 8.7971891906493802e+03, 8.8044391161860967e+03, 8.8116897516979952e+03, 8.8189410966813684e+03, 8.8261931506332203e+03, 8.8334459130512751e+03, 8.8406993834339592e+03, 8.8479535612804157e+03, 8.8552084460904916e+03, 8.8624640373647453e+03, 8.8697203346044371e+03, 8.8769773373115295e+03, 8.8842350449886908e+03, 8.8914934571392860e+03, 8.8987525732673821e+03, 8.9060123928777448e+03, 8.9132729154758363e+03, 8.9205341405678082e+03, 8.9277960676605126e+03, 8.9350586962614852e+03, 8.9423220258789624e+03, 8.9495860560218607e+03, 8.9568507861997896e+03, 8.9641162159230425e+03, 8.9713823447025989e+03, 8.9786491720501199e+03, 8.9859166974779473e+03, 8.9931849204991067e+03, 9.0004538406273005e+03, 9.0077234573769092e+03, 9.0149937702629886e+03, 9.0222647788012691e+03, 9.0295364825081560e+03, 9.0368088809007240e+03, 9.0440819734967245e+03, 9.0513557598145690e+03, 9.0586302393733440e+03, 9.0659054116927964e+03, 9.0731812762933441e+03, 9.0804578326960636e+03, 9.0877350804226935e+03, 9.0950130189956380e+03, 9.1022916479379583e+03, 9.1095709667733736e+03, 9.1168509750262583e+03, 9.1241316722216434e+03, 9.1314130578852164e+03, 9.1386951315433089e+03, 9.1459778927229145e+03, 9.1532613409516725e+03, 9.1605454757578682e+03, 9.1678302966704377e+03, 9.1751158032189596e+03, 9.1824019949336598e+03, 9.1896888713454100e+03, 9.1969764319857204e+03, 9.2042646763867397e+03, 9.2115536040812603e+03, 9.2188432146027117e+03, 9.2261335074851577e+03, 9.2334244822633009e+03, 9.2407161384724750e+03, 9.2480084756486467e+03, 9.2553014933284194e+03, 9.2625951910490239e+03, 9.2698895683483115e+03, 9.2771846247647736e+03, 9.2844803598375238e+03, 9.2917767731062959e+03, 9.2990738641114549e+03, 9.3063716323939861e+03, 9.3136700774954952e+03, 9.3209691989582061e+03, 9.3282689963249650e+03, 9.3355694691392328e+03, 9.3428706169450888e+03, 9.3501724392872275e+03, 9.3574749357109540e+03, 9.3647781057621905e+03, 9.3720819489874684e+03, 9.3793864649339266e+03, 9.3866916531493207e+03, 9.3939975131820065e+03, 9.4013040445809456e+03, 9.4086112468957090e+03, 9.4159191196764714e+03, 9.4232276624740116e+03, 9.4305368748397032e+03, 9.4378467563255290e+03, 9.4451573064840650e+03, 9.4524685248684855e+03, 9.4597804110325633e+03, 9.4670929645306660e+03, 9.4744061849177560e+03, 9.4817200717493870e+03, 9.4890346245817109e+03, 9.4963498429714637e+03, 9.5036657264759742e+03, 9.5109822746531554e+03, 9.5182994870615166e+03, 9.5256173632601440e+03, 9.5329359028087110e+03, 9.5402551052674808e+03, 9.5475749701972909e+03, 9.5548954971595649e+03, 9.5622166857163047e+03, 9.5695385354300943e+03, 9.5768610458640924e+03, 9.5841842165820381e+03, 9.5915080471482397e+03, 9.5988325371275914e+03, 9.6061576860855457e+03, 9.6134834935881427e+03, 9.6208099592019844e+03, 9.6281370824942423e+03, 9.6354648630326628e+03, 9.6427933003855578e+03, 9.6501223941218050e+03, 9.6574521438108459e+03, 9.6647825490226933e+03, 9.6721136093279129e+03, 9.6794453242976379e+03, 9.6867776935035672e+03, 9.6941107165179528e+03, 9.7014443929136105e+03, 9.7087787222639108e+03, 9.7161137041427810e+03, 9.7234493381247103e+03, 9.7307856237847300e+03, 9.7381225606984390e+03, 9.7454601484419763e+03, 9.7527983865920432e+03, 9.7601372747258829e+03, 9.7674768124212915e+03, 9.7748169992566109e+03, 9.7821578348107323e+03, 9.7894993186630963e+03, 9.7968414503936801e+03, 9.8041842295830138e+03, 9.8115276558121623e+03, 9.8188717286627361e+03, 9.8262164477168863e+03, 9.8335618125573037e+03, 9.8409078227672162e+03, 9.8482544779303917e+03, 9.8556017776311328e+03, 9.8629497214542807e+03, 9.8702983089852096e+03, 9.8776475398098228e+03, 9.8849974135145621e+03, 9.8923479296863952e+03, 9.8996990879128280e+03, 9.9070508877818847e+03, 9.9144033288821283e+03, 9.9217564108026436e+03, 9.9291101331330428e+03, 9.9364644954634641e+03, 9.9438194973845693e+03, 9.9511751384875442e+03, 9.9585314183640930e+03, 9.9658883366064474e+03, 9.9732458928073593e+03, 9.9806040865600899e+03, 9.9879629174584315e+03, 9.9953223850966897e+03, 1.0002682489069677e+04, 1.0010043228972736e+04, 1.0017404604401714e+04, 1.0024766614952974e+04, 1.0032129260223392e+04, 1.0039492539810353e+04, 1.0046856453311759e+04, 1.0054221000326015e+04, 1.0061586180452037e+04, 1.0068951993289247e+04, 1.0076318438437575e+04, 1.0083685515497456e+04, 1.0091053224069829e+04, 1.0098421563756139e+04, 1.0105790534158336e+04, 1.0113160134878861e+04, 1.0120530365520668e+04, 1.0127901225687207e+04, 1.0135272714982420e+04, 1.0142644833010758e+04, 1.0150017579377163e+04, 1.0157390953687072e+04, 1.0164764955546421e+04, 1.0172139584561643e+04, 1.0179514840339651e+04, 1.0186890722487866e+04, 1.0194267230614192e+04, 1.0201644364327027e+04, 1.0209022123235254e+04, 1.0216400506948250e+04, 1.0223779515075879e+04, 1.0231159147228489e+04, 1.0238539403016914e+04, 1.0245920282052479e+04, 1.0253301783946987e+04, 1.0260683908312723e+04, 1.0268066654762462e+04, 1.0275450022909454e+04, 1.0282834012367432e+04, 1.0290218622750608e+04, 1.0297603853673672e+04, 1.0304989704751799e+04, 1.0312376175600628e+04, 1.0319763265836284e+04, 1.0327150975075365e+04, 1.0334539302934943e+04, 1.0341928249032562e+04, 1.0349317812986241e+04, 1.0356707994414466e+04, 1.0364098792936202e+04, 1.0371490208170877e+04, 1.0378882239738390e+04, 1.0386274887259113e+04, 1.0393668150353877e+04, 1.0401062028643984e+04, 1.0408456521751203e+04, 1.0415851629297764e+04, 1.0423247350906366e+04, 1.0430643686200166e+04, 1.0438040634802788e+04, 1.0445438196338313e+04, 1.0452836370431281e+04, 1.0460235156706702e+04, 1.0467634554790033e+04, 1.0475034564307196e+04, 1.0482435184884567e+04, 1.0489836416148979e+04, 1.0497238257727722e+04, 1.0504640709248542e+04, 1.0512043770339631e+04, 1.0519447440629645e+04, 1.0526851719747681e+04, 1.0534256607323297e+04, 1.0541662102986496e+04, 1.0549068206367732e+04, 1.0556474917097908e+04, 1.0563882234808378e+04, 1.0571290159130938e+04, 1.0578698689697834e+04, 1.0586107826141755e+04, 1.0593517568095836e+04, 1.0600927915193657e+04, 1.0608338867069240e+04, 1.0615750423357053e+04, 1.0623162583691997e+04, 1.0630575347709424e+04, 1.0637988715045120e+04, 1.0645402685335308e+04, 1.0652817258216659e+04, 1.0660232433326271e+04, 1.0667648210301688e+04, 1.0675064588780882e+04, 1.0682481568402263e+04, 1.0689899148804678e+04, 1.0697317329627404e+04, 1.0704736110510155e+04, 1.0712155491093074e+04, 1.0719575471016737e+04, 1.0726996049922147e+04, 1.0734417227450742e+04, 1.0741839003244386e+04, 1.0749261376945373e+04, 1.0756684348196422e+04, 1.0764107916640682e+04, 1.0771532081921723e+04, 1.0778956843683547e+04, 1.0786382201570574e+04, 1.0793808155227651e+04, 1.0801234704300048e+04, 1.0808661848433456e+04, 1.0816089587273987e+04, 1.0823517920468179e+04, 1.0830946847662983e+04, 1.0838376368505771e+04, 1.0845806482644333e+04, 1.0853237189726880e+04, 1.0860668489402035e+04, 1.0868100381318842e+04, 1.0875532865126759e+04, 1.0882965940475658e+04, 1.0890399607015823e+04, 1.0897833864397957e+04, 1.0905268712273170e+04, 1.0912704150292986e+04, 1.0920140178109335e+04, 1.0927576795374569e+04, 1.0935014001741440e+04, 1.0942451796863112e+04, 1.0949890180393155e+04, 1.0957329151985552e+04, 1.0964768711294684e+04, 1.0972208857975349e+04, 1.0979649591682737e+04, 1.0987090912072454e+04, 1.0994532818800504e+04, 1.1001975311523298e+04, 1.1009418389897648e+04, 1.1016862053580762e+04, 1.1024306302230259e+04, 1.1031751135504152e+04, 1.1039196553060852e+04, 1.1046642554559176e+04, 1.1054089139658334e+04, 1.1061536308017936e+04, 1.1068984059297984e+04, 1.1076432393158881e+04, 1.1083881309261424e+04, 1.1091330807266808e+04, 1.1098780886836616e+04, 1.1106231547632828e+04, 1.1113682789317816e+04, 1.1121134611554344e+04, 1.1128587014005567e+04, 1.1136039996335032e+04, 1.1143493558206674e+04, 1.1150947699284821e+04, 1.1158402419234186e+04, 1.1165857717719868e+04, 1.1173313594407360e+04, 1.1180770048962537e+04, 1.1188227081051660e+04, 1.1195684690341375e+04, 1.1203142876498716e+04, 1.1210601639191096e+04, 1.1218060978086318e+04, 1.1225520892852559e+04, 1.1232981383158385e+04, 1.1240442448672738e+04, 1.1247904089064948e+04, 1.1255366304004716e+04, 1.1262829093162127e+04, 1.1270292456207648e+04, 1.1277756392812116e+04, 1.1285220902646754e+04, 1.1292685985383154e+04, 1.1300151640693288e+04, 1.1307617868249503e+04, 1.1315084667724521e+04, 1.1322552038791438e+04, 1.1330019981123725e+04, 1.1337488494395222e+04, 1.1344957578280140e+04, 1.1352427232453074e+04, 1.1359897456588973e+04, 1.1367368250363168e+04, 1.1374839613451357e+04, 1.1382311545529603e+04, 1.1389784046274341e+04, 1.1397257115362374e+04, 1.1404730752470870e+04, 1.1412204957277367e+04, 1.1419679729459764e+04, 1.1427155068696331e+04, 1.1434630974665701e+04, 1.1442107447046863e+04, 1.1449584485519183e+04, 1.1457062089762379e+04, 1.1464540259456540e+04, 1.1472018994282107e+04, 1.1479498293919891e+04, 1.1486978158051055e+04, 1.1494458586357128e+04, 1.1501939578519998e+04, 1.1509421134221908e+04, 1.1516903253145461e+04, 1.1524385934973614e+04, 1.1531869179389687e+04, 1.1539352986077352e+04, 1.1546837354720637e+04, 1.1554322285003927e+04, 1.1561807776611955e+04, 1.1569293829229820e+04, 1.1576780442542962e+04, 1.1584267616237174e+04, 1.1591755349998612e+04, 1.1599243643513772e+04, 1.1606732496469504e+04, 1.1614221908553012e+04, 1.1621711879451846e+04, 1.1629202408853907e+04, 1.1636693496447440e+04, 1.1644185141921047e+04, 1.1651677344963664e+04, 1.1659170105264588e+04, 1.1666663422513449e+04, 1.1674157296400233e+04, 1.1681651726615266e+04, 1.1689146712849217e+04, 1.1696642254793102e+04, 1.1704138352138278e+04, 1.1711635004576445e+04, 1.1719132211799648e+04, 1.1726629973500272e+04, 1.1734128289371038e+04, 1.1741627159105014e+04, 1.1749126582395606e+04, 1.1756626558936559e+04, 1.1764127088421954e+04, 1.1771628170546215e+04, 1.1779129805004099e+04, 1.1786631991490702e+04, 1.1794134729701453e+04, 1.1801638019332126e+04, 1.1809141860078826e+04, 1.1816646251637989e+04, 1.1824151193706382e+04, 1.1831656685981123e+04, 1.1839162728159639e+04, 1.1846669319939710e+04, 1.1854176461019435e+04, 1.1861684151097257e+04, 1.1869192389871938e+04, 1.1876701177042569e+04, 1.1884210512308588e+04, 1.1891720395369741e+04, 1.1899230825926117e+04, 1.1906741803678133e+04, 1.1914253328326524e+04, 1.1921765399572358e+04, 1.1929278017117033e+04, 1.1936791180662263e+04, 1.1944304889910103e+04, 1.1951819144562918e+04, 1.1959333944323409e+04, 1.1966849288894589e+04, 1.1974365177979804e+04, 1.1981881611282717e+04, 1.1989398588507322e+04, 1.1996916109357922e+04, 1.2004434173539155e+04, 1.2011952780755970e+04, 1.2019471930713642e+04, 1.2026991623117761e+04, 1.2034511857674233e+04, 1.2042032634089295e+04, 1.2049553952069497e+04, 1.2057075811321696e+04, 1.2064598211553084e+04, 1.2072121152471156e+04, 1.2079644633783731e+04, 1.2087168655198935e+04, 1.2094693216425219e+04, 1.2102218317171346e+04, 1.2109743957146389e+04, 1.2117270136059737e+04, 1.2124796853621088e+04, 1.2132324109540461e+04, 1.2139851903528184e+04, 1.2147380235294890e+04, 1.2154909104551532e+04, 1.2162438511009370e+04, 1.2169968454379969e+04, 1.2177498934375215e+04, 1.2185029950707294e+04, 1.2192561503088704e+04, 1.2200093591232244e+04, 1.2207626214851034e+04, 1.2215159373658489e+04, 1.2222693067368336e+04, 1.2230227295694609e+04, 1.2237762058351649e+04, 1.2245297355054094e+04, 1.2252833185516893e+04, 1.2260369549455298e+04, 1.2267906446584866e+04, 1.2275443876621450e+04, 1.2282981839281221e+04, 1.2290520334280633e+04, 1.2298059361336456e+04, 1.2305598920165756e+04, 1.2313139010485900e+04, 1.2320679632014557e+04, 1.2328220784469695e+04, 1.2335762467569575e+04, 1.2343304681032771e+04, 1.2350847424578142e+04, 1.2358390697924848e+04, 1.2365934500792349e+04, 1.2373478832900402e+04, 1.2381023693969060e+04, 1.2388569083718672e+04, 1.2396115001869885e+04, 1.2403661448143628e+04, 1.2411208422261147e+04, 1.2418755923943963e+04, 1.2426303952913897e+04, 1.2433852508893067e+04, 1.2441401591603879e+04, 1.2448951200769034e+04, 1.2456501336111522e+04, 1.2464051997354627e+04, 1.2471603184221925e+04, 1.2479154896437276e+04, 1.2486707133724836e+04, 1.2494259895809048e+04, 1.2501813182414649e+04, 1.2509366993266658e+04, 1.2516921328090384e+04, 1.2524476186611422e+04, 1.2532031568555662e+04, 1.2539587473649275e+04, 1.2547143901618716e+04, 1.2554700852190728e+04, 1.2562258325092344e+04, 1.2569816320050873e+04, 1.2577374836793921e+04, 1.2584933875049361e+04, 1.2592493434545371e+04, 1.2600053515010391e+04, 1.2607614116173163e+04, 1.2615175237762689e+04, 1.2622736879508278e+04, 1.2630299041139506e+04, 1.2637861722386226e+04, 1.2645424922978584e+04, 1.2652988642646998e+04, 1.2660552881122168e+04, 1.2668117638135072e+04, 1.2675682913416973e+04, 1.2683248706699402e+04, 1.2690815017714171e+04, 1.2698381846193379e+04, 1.2705949191869391e+04, 1.2713517054474856e+04, 1.2721085433742695e+04, 1.2728654329406103e+04, 1.2736223741198553e+04, 1.2743793668853794e+04, 1.2751364112105852e+04, 1.2758935070689020e+04, 1.2766506544337872e+04, 1.2774078532787245e+04, 1.2781651035772267e+04, 1.2789224053028318e+04, 1.2796797584291066e+04, 1.2804371629296436e+04, 1.2811946187780641e+04, 1.2819521259480147e+04, 1.2827096844131707e+04, 1.2834672941472329e+04, 1.2842249551239305e+04, 1.2849826673170179e+04, 1.2857404307002782e+04, 1.2864982452475202e+04, 1.2872561109325796e+04, 1.2880140277293192e+04, 1.2887719956116283e+04, 1.2895300145534229e+04, 1.2902880845286456e+04, 1.2910462055112648e+04, 1.2918043774752774e+04, 1.2925626003947049e+04, 1.2933208742435965e+04, 1.2940791989960266e+04, 1.2948375746260974e+04, 1.2955960011079364e+04, 1.2963544784156975e+04, 1.2971130065235613e+04, 1.2978715854057345e+04, 1.2986302150364498e+04, 1.2993888953899661e+04, 1.3001476264405681e+04, 1.3009064081625676e+04, 1.3016652405303010e+04, 1.3024241235181316e+04, 1.3031830571004486e+04, 1.3039420412516671e+04, 1.3047010759462273e+04, 1.3054601611585962e+04, 1.3062192968632660e+04, 1.3069784830347549e+04, 1.3077377196476069e+04, 1.3084970066763917e+04, 1.3092563440957036e+04, 1.3100157318801643e+04, 1.3107751700044193e+04, 1.3115346584431414e+04, 1.3122941971710266e+04, 1.3130537861627987e+04, 1.3138134253932052e+04, 1.3145731148370194e+04, 1.3153328544690410e+04, 1.3160926442640930e+04, 1.3168524841970253e+04, 1.3176123742427126e+04, 1.3183723143760542e+04, 1.3191323045719751e+04, 1.3198923448054253e+04, 1.3206524350513793e+04, 1.3214125752848378e+04, 1.3221727654808254e+04, 1.3229330056143919e+04, 1.3236932956606122e+04, 1.3244536355945866e+04, 1.3252140253914384e+04, 1.3259744650263181e+04, 1.3267349544743991e+04, 1.3274954937108809e+04, 1.3282560827109863e+04, 1.3290167214499632e+04, 1.3297774099030854e+04, 1.3305381480456495e+04, 1.3312989358529774e+04, 1.3320597733004155e+04, 1.3328206603633347e+04, 1.3335815970171299e+04, 1.3343425832372212e+04, 1.3351036189990524e+04, 1.3358647042780920e+04, 1.3366258390498324e+04, 1.3373870232897903e+04, 1.3381482569735073e+04, 1.3389095400765482e+04, 1.3396708725745020e+04, 1.3404322544429828e+04, 1.3411936856576280e+04, 1.3419551661940992e+04, 1.3427166960280820e+04, 1.3434782751352855e+04, 1.3442399034914435e+04, 1.3450015810723135e+04, 1.3457633078536763e+04, 1.3465250838113368e+04, 1.3472869089211248e+04, 1.3480487831588916e+04, 1.3488107065005142e+04, 1.3495726789218925e+04, 1.3503347003989498e+04, 1.3510967709076338e+04, 1.3518588904239148e+04, 1.3526210589237869e+04, 1.3533832763832686e+04, 1.3541455427784012e+04, 1.3549078580852487e+04, 1.3556702222798996e+04, 1.3564326353384660e+04, 1.3571950972370820e+04, 1.3579576079519060e+04, 1.3587201674591191e+04, 1.3594827757349261e+04, 1.3602454327555553e+04, 1.3610081384972573e+04, 1.3617708929363062e+04, 1.3625336960489991e+04, 1.3632965478116565e+04, 1.3640594482006220e+04, 1.3648223971922615e+04, 1.3655853947629641e+04, 1.3663484408891427e+04, 1.3671115355472319e+04, 1.3678746787136895e+04, 1.3686378703649965e+04, 1.3694011104776566e+04, 1.3701643990281964e+04, 1.3709277359931641e+04, 1.3716911213491321e+04, 1.3724545550726953e+04, 1.3732180371404700e+04, 1.3739815675290954e+04, 1.3747451462152350e+04, 1.3755087731755730e+04, 1.3762724483868167e+04, 1.3770361718256956e+04, 1.3777999434689626e+04, 1.3785637632933909e+04, 1.3793276312757787e+04, 1.3800915473929444e+04, 1.3808555116217305e+04, 1.3816195239390001e+04, 1.3823835843216393e+04, 1.3831476927465568e+04, 1.3839118491906827e+04, 1.3846760536309701e+04, 1.3854403060443934e+04, 1.3862046064079497e+04, 1.3869689546986572e+04, 1.3877333508935577e+04, 1.3884977949697133e+04, 1.3892622869042092e+04, 1.3900268266741519e+04, 1.3907914142566706e+04, 1.3915560496289150e+04, 1.3923207327680580e+04, 1.3930854636512939e+04, 1.3938502422558378e+04, 1.3946150685589282e+04, 1.3953799425378236e+04, 1.3961448641698060e+04, 1.3969098334321772e+04, 1.3976748503022618e+04, 1.3984399147574053e+04, 1.3992050267749755e+04, 1.3999701863323617e+04, 1.4007353934069730e+04, 1.4015006479762427e+04, 1.4022659500176231e+04, 1.4030312995085891e+04, 1.4037966964266370e+04, 1.4045621407492841e+04, 1.4053276324540691e+04, 1.4060931715185516e+04, 1.4068587579203131e+04, 1.4076243916369562e+04, 1.4083900726461043e+04, 1.4091558009254020e+04, 1.4099215764525152e+04, 1.4106873992051316e+04, 1.4114532691609582e+04, 1.4122191862977250e+04, 1.4129851505931814e+04, 1.4137511620250991e+04, 1.4145172205712692e+04, 1.4152833262095055e+04, 1.4160494789176415e+04, 1.4168156786735317e+04, 1.4175819254550517e+04, 1.4183482192400979e+04, 1.4191145600065875e+04, 1.4198809477324578e+04, 1.4206473823956678e+04, 1.4214138639741961e+04, 1.4221803924460435e+04, 1.4229469677892293e+04, 1.4237135899817955e+04, 1.4244802590018035e+04, 1.4252469748273355e+04, 1.4260137374364940e+04, 1.4267805468074022e+04, 1.4275474029182040e+04, 1.4283143057470626e+04, 1.4290812552721634e+04, 1.4298482514717109e+04, 1.4306152943239298e+04, 1.4313823838070663e+04, 1.4321495198993855e+04, 1.4329167025791734e+04, 1.4336839318247361e+04, 1.4344512076144003e+04, 1.4352185299265127e+04, 1.4359858987394395e+04, 1.4367533140315676e+04, 1.4375207757813036e+04, 1.4382882839670756e+04, 1.4390558385673290e+04, 1.4398234395605321e+04, 1.4405910869251709e+04, 1.4413587806397531e+04, 1.4421265206828046e+04, 1.4428943070328723e+04, 1.4436621396685228e+04, 1.4444300185683429e+04, 1.4451979437109383e+04, 1.4459659150749349e+04, 1.4467339326389783e+04, 1.4475019963817345e+04, 1.4482701062818885e+04, 1.4490382623181444e+04, 1.4498064644692271e+04, 1.4505747127138804e+04, 1.4513430070308683e+04, 1.4521113473989733e+04, 1.4528797337969992e+04, 1.4536481662037671e+04, 1.4544166445981196e+04, 1.4551851689589172e+04, 1.4559537392650405e+04, 1.4567223554953898e+04, 1.4574910176288840e+04, 1.4582597256444626e+04, 1.4590284795210822e+04, 1.4597972792377217e+04, 1.4605661247733768e+04, 1.4613350161070632e+04, 1.4621039532178163e+04, 1.4628729360846899e+04, 1.4636419646867575e+04, 1.4644110390031119e+04, 1.4651801590128640e+04, 1.4659493246951452e+04, 1.4667185360291043e+04, 1.4674877929939112e+04, 1.4682570955687528e+04, 1.4690264437328364e+04, 1.4697958374653872e+04, 1.4705652767456502e+04, 1.4713347615528886e+04, 1.4721042918663852e+04, 1.4728738676654402e+04, 1.4736434889293751e+04, 1.4744131556375278e+04, 1.4751828677692562e+04, 1.4759526253039365e+04, 1.4767224282209636e+04, 1.4774922764997516e+04, 1.4782621701197329e+04, 1.4790321090603587e+04, 1.4798020933010985e+04, 1.4805721228214405e+04, 1.4813421976008916e+04, 1.4821123176189774e+04, 1.4828824828552415e+04, 1.4836526932892468e+04, 1.4844229489005736e+04, 1.4851932496688214e+04, 1.4859635955736079e+04, 1.4867339865945700e+04, 1.4875044227113611e+04, 1.4882749039036544e+04, 1.4890454301511409e+04, 1.4898160014335304e+04, 1.4905866177305505e+04, 1.4913572790219469e+04, 1.4921279852874837e+04, 1.4928987365069437e+04, 1.4936695326601273e+04, 1.4944403737268529e+04, 1.4952112596869576e+04, 1.4959821905202962e+04, 1.4967531662067415e+04, 1.4975241867261850e+04, 1.4982952520585352e+04, 1.4990663621837191e+04, 1.4998375170816818e+04, 1.5006087167323867e+04, 1.5013799611158143e+04, 1.5021512502119635e+04, 1.5029225840008507e+04, 1.5036939624625104e+04, 1.5044653855769955e+04, 1.5052368533243754e+04, 1.5060083656847388e+04, 1.5067799226381909e+04, 1.5075515241648553e+04, 1.5083231702448729e+04, 1.5090948608584027e+04, 1.5098665959856213e+04, 1.5106383756067227e+04, 1.5114101997019186e+04, 1.5121820682514384e+04, 1.5129539812355293e+04, 1.5137259386344549e+04, 1.5144979404284981e+04, 1.5152699865979581e+04, 1.5160420771231516e+04, 1.5168142119844137e+04, 1.5175863911620954e+04, 1.5183586146365664e+04, 1.5191308823882131e+04, 1.5199031943974400e+04, 1.5206755506446678e+04, 1.5214479511103355e+04, 1.5222203957748987e+04, 1.5229928846188310e+04, 1.5237654176226226e+04, 1.5245379947667816e+04, 1.5253106160318321e+04, 1.5260832813983170e+04, 1.5268559908467947e+04, 1.5276287443578423e+04, 1.5284015419120531e+04, 1.5291743834900373e+04, 1.5299472690724228e+04, 1.5307201986398535e+04, 1.5314931721729919e+04, 1.5322661896525167e+04, 1.5330392510591229e+04, 1.5338123563735235e+04, 1.5345855055764479e+04, 1.5353586986486431e+04, 1.5361319355708712e+04, 1.5369052163239137e+04, 1.5376785408885664e+04, 1.5384519092456439e+04, 1.5392253213759768e+04, 1.5399987772604121e+04, 1.5407722768798145e+04, 1.5415458202150645e+04, 1.5423194072470598e+04, 1.5430930379567149e+04, 1.5438667123249603e+04, 1.5446404303327437e+04, 1.5454141919610296e+04, 1.5461879971907985e+04, 1.5469618460030479e+04, 1.5477357383787918e+04, 1.5485096742990607e+04, 1.5492836537449019e+04, 1.5500576766973782e+04, 1.5508317431375699e+04, 1.5516058530465734e+04, 1.5523800064055014e+04, 1.5531542031954834e+04, 1.5539284433976651e+04, 1.5547027269932079e+04, 1.5554770539632909e+04, 1.5562514242891082e+04, 1.5570258379518715e+04, 1.5578002949328069e+04, 1.5585747952131582e+04, 1.5593493387741857e+04, 1.5601239255971648e+04, 1.5608985556633879e+04, 1.5616732289541635e+04, 1.5624479454508155e+04, 1.5632227051346848e+04, 1.5639975079871283e+04, 1.5647723539895183e+04, 1.5655472431232438e+04, 1.5663221753697098e+04, 1.5670971507103370e+04, 1.5678721691265628e+04, 1.5686472305998401e+04, 1.5694223351116374e+04, 1.5701974826434393e+04, 1.5709726731767472e+04, 1.5717479066930775e+04, 1.5725231831739622e+04, 1.5732985026009506e+04, 1.5740738649556066e+04, 1.5748492702195103e+04, 1.5756247183742571e+04, 1.5764002094014593e+04, 1.5771757432827440e+04, 1.5779513199997544e+04, 1.5787269395341491e+04, 1.5795026018676026e+04, 1.5802783069818061e+04, 1.5810540548584644e+04, 1.5818298454792994e+04, 1.5826056788260486e+04, 1.5833815548804645e+04, 1.5841574736243152e+04, 1.5849334350393847e+04, 1.5857094391074728e+04, 1.5864854858103945e+04, 1.5872615751299794e+04, 1.5880377070480741e+04, 1.5888138815465400e+04, 1.5895900986072540e+04, 1.5903663582121082e+04, 1.5911426603430098e+04, 1.5919190049818828e+04, 1.5926953921106646e+04, 1.5934718217113097e+04, 1.5942482937657869e+04, 1.5950248082560805e+04, 1.5958013651641902e+04, 1.5965779644721310e+04, 1.5973546061619329e+04, 1.5981312902156415e+04, 1.5989080166153175e+04, 1.5996847853430359e+04, 1.6004615963808885e+04, 1.6012384497109813e+04, 1.6020153453154349e+04, 1.6027922831763864e+04, 1.6035692632759869e+04, 1.6043462855964028e+04, 1.6051233501198158e+04, 1.6059004568284223e+04, 1.6066776057044341e+04, 1.6074547967300778e+04, 1.6082320298875942e+04, 1.6090093051592414e+04, 1.6097866225272895e+04, 1.6105639819740254e+04, 1.6113413834817506e+04, 1.6121188270327808e+04, 1.6128963126094472e+04, 1.6136738401940960e+04, 1.6144514097690875e+04, 1.6152290213167975e+04, 1.6160066748196159e+04, 1.6167843702599484e+04, 1.6175621076202142e+04, 1.6183398868828479e+04, 1.6191177080302992e+04, 1.6198955710450316e+04, 1.6206734759095247e+04, 1.6214514226062705e+04, 1.6222294111177775e+04, 1.6230074414265682e+04, 1.6237855135151802e+04, 1.6245636273661648e+04, 1.6253417829620879e+04, 1.6261199802855315e+04, 1.6268982193190901e+04, 1.6276765000453741e+04, 1.6284548224470078e+04, 1.6292331865066297e+04, 1.6300115922068939e+04, 1.6307900395304676e+04, 1.6315685284600329e+04, 1.6323470589782870e+04, 1.6331256310679406e+04, 1.6339042447117188e+04, 1.6346828998923616e+04, 1.6354615965926228e+04, 1.6362403347952715e+04, 1.6370191144830895e+04, 1.6377979356388742e+04, 1.6385767982454367e+04, 1.6393557022856025e+04, 1.6401346477422110e+04, 1.6409136345981169e+04, 1.6416926628361867e+04, 1.6424717324393052e+04, 1.6432508433903655e+04, 1.6440299956722807e+04, 1.6448091892679742e+04, 1.6455884241603850e+04, 1.6463677003324668e+04, 1.6471470177671854e+04, 1.6479263764475229e+04, 1.6487057763564728e+04, 1.6494852174770454e+04, 1.6502646997922631e+04, 1.6510442232851634e+04, 1.6518237879387969e+04, 1.6526033937362285e+04, 1.6533830406605372e+04, 1.6541627286948151e+04, 1.6549424578221700e+04, 1.6557222280257218e+04, 1.6565020392886050e+04, 1.6572818915939675e+04, 1.6580617849249717e+04, 1.6588417192647936e+04, 1.6596216945966222e+04, 1.6604017109036617e+04, 1.6611817681691293e+04, 1.6619618663762543e+04, 1.6627420055082835e+04, 1.6635221855484742e+04, 1.6643024064800989e+04, 1.6650826682864434e+04, 1.6658629709508066e+04, 1.6666433144565017e+04, 1.6674236987868549e+04, 1.6682041239252081e+04, 1.6689845898549138e+04, 1.6697650965593395e+04, 1.6705456440218670e+04, 1.6713262322258899e+04, 1.6721068611548162e+04, 1.6728875307920680e+04, 1.6736682411210808e+04, 1.6744489921253022e+04, 1.6752297837881950e+04, 1.6760106160932341e+04, 1.6767914890239088e+04, 1.6775724025637206e+04, 1.6783533566961858e+04, 1.6791343514048331e+04, 1.6799153866732053e+04, 1.6806964624848581e+04, 1.6814775788233608e+04, 1.6822587356722954e+04, 1.6830399330152573e+04, 1.6838211708358558e+04, 1.6846024491177133e+04, 1.6853837678444655e+04, 1.6861651269997608e+04, 1.6869465265672610e+04, 1.6877279665306422e+04, 1.6885094468735904e+04, 1.6892909675798095e+04, 1.6900725286330133e+04, 1.6908541300169290e+04, 1.6916357717152987e+04, 1.6924174537118750e+04, 1.6931991759904256e+04, 1.6939809385347311e+04, 1.6947627413285845e+04, 1.6955445843557914e+04, 1.6963264676001720e+04, 1.6971083910455574e+04, 1.6978903546757945e+04, 1.6986723584747400e+04, 1.6994544024262661e+04, 1.7002364865142572e+04, 1.7010186107226094e+04, 1.7018007750352335e+04, 1.7025829794360525e+04, 1.7033652239090014e+04, 1.7041475084380290e+04, 1.7049298330070978e+04, 1.7057121976001810e+04, 1.7064946022012664e+04, 1.7072770467943545e+04, 1.7080595313634571e+04, 1.7088420558926005e+04, 1.7096246203658226e+04, 1.7104072247671746e+04, 1.7111898690807197e+04, 1.7119725532905359e+04, 1.7127552773807114e+04, 1.7135380413353480e+04, 1.7143208451385606e+04, 1.7151036887744769e+04, 1.7158865722272352e+04, 1.7166694954809896e+04, 1.7174524585199048e+04, 1.7182354613281579e+04, 1.7190185038899399e+04, 1.7198015861894535e+04, 1.7205847082109140e+04, 1.7213678699385491e+04, 1.7221510713566000e+04, 1.7229343124493189e+04, 1.7237175932009712e+04, 1.7245009135958353e+04, 1.7252842736182018e+04, 1.7260676732523723e+04, 1.7268511124826633e+04, 1.7276345912934026e+04, 1.7284181096689292e+04, 1.7292016675935960e+04, 1.7299852650517681e+04, 1.7307689020278231e+04, 1.7315525785061494e+04, 1.7323362944711494e+04, 1.7331200499072376e+04, 1.7339038447988401e+04, 1.7346876791303963e+04, 1.7354715528863559e+04, 1.7362554660511829e+04, 1.7370394186093537e+04, 1.7378234105453554e+04, 1.7386074418436874e+04, 1.7393915124888626e+04, 1.7401756224654047e+04, 1.7409597717578508e+04, 1.7417439603507493e+04, 1.7425281882286603e+04, 1.7433124553761587e+04, 1.7440967617778279e+04, 1.7448811074182653e+04, 1.7456654922820806e+04, 1.7464499163538949e+04, 1.7472343796183413e+04, 1.7480188820600655e+04, 1.7488034236637250e+04, 1.7495880044139885e+04, 1.7503726242955385e+04, 1.7511572832930677e+04, 1.7519419813912813e+04, 1.7527267185748973e+04, 1.7535114948286442e+04, 1.7542963101372643e+04, 1.7550811644855097e+04, 1.7558660578581461e+04, 1.7566509902399506e+04, 1.7574359616157111e+04, 1.7582209719702281e+04, 1.7590060212883156e+04, 1.7597911095547966e+04, 1.7605762367545078e+04, 1.7613614028722965e+04, 1.7621466078930232e+04, 1.7629318518015592e+04, 1.7637171345827872e+04, 1.7645024562216029e+04, 1.7652878167029128e+04, 1.7660732160116349e+04, 1.7668586541327008e+04, 1.7676441310510505e+04, 1.7684296467516386e+04, 1.7692152012194299e+04, 1.7700007944394019e+04, 1.7707864263965428e+04, 1.7715720970758524e+04, 1.7723578064623423e+04, 1.7731435545410364e+04, 1.7739293412969699e+04, 1.7747151667151884e+04, 1.7755010307807504e+04, 1.7762869334787250e+04, 1.7770728747941943e+04, 1.7778588547122508e+04, 1.7786448732179975e+04, 1.7794309302965514e+04, 1.7802170259330393e+04, 1.7810031601125993e+04, 1.7817893328203816e+04, 1.7825755440415476e+04, 1.7833617937612711e+04, 1.7841480819647357e+04, 1.7849344086371362e+04, 1.7857207737636811e+04, 1.7865071773295880e+04, 1.7872936193200876e+04, 1.7880800997204209e+04, 1.7888666185158396e+04, 1.7896531756916076e+04, 1.7904397712330010e+04, 1.7912264051253056e+04, 1.7920130773538196e+04, 1.7927997879038514e+04, 1.7935865367607214e+04, 1.7943733239097612e+04, 1.7951601493363130e+04, 1.7959470130257316e+04, 1.7967339149633815e+04, 1.7975208551346390e+04, 1.7983078335248923e+04, 1.7990948501195391e+04, 1.7998819049039899e+04, 1.8006689978636648e+04, 1.8014561289839974e+04, 1.8022432982504295e+04, 1.8030305056484165e+04, 1.8038177511634225e+04, 1.8046050347809254e+04, 1.8053923564864115e+04, 1.8061797162653802e+04, 1.8069671141033403e+04, 1.8077545499858134e+04, 1.8085420238983304e+04, 1.8093295358264346e+04, 1.8101170857556794e+04, 1.8109046736716293e+04, 1.8116922995598597e+04, 1.8124799634059571e+04, 1.8132676651955197e+04, 1.8140554049141549e+04, 1.8148431825474829e+04, 1.8156309980811329e+04, 1.8164188515007470e+04, 1.8172067427919770e+04, 1.8179946719404848e+04, 1.8187826389319453e+04, 1.8195706437520421e+04, 1.8203586863864715e+04, 1.8211467668209392e+04, 1.8219348850411614e+04, 1.8227230410328670e+04, 1.8235112347817943e+04, 1.8242994662736925e+04, 1.8250877354943212e+04, 1.8258760424294513e+04, 1.8266643870648652e+04, 1.8274527693863543e+04, 1.8282411893797220e+04, 1.8290296470307814e+04, 1.8298181423253576e+04, 1.8306066752492843e+04, 1.8313952457884086e+04, 1.8321838539285865e+04, 1.8329724996556841e+04, 1.8337611829555801e+04, 1.8345499038141610e+04, 1.8353386622173271e+04, 1.8361274581509875e+04, 1.8369162916010613e+04, 1.8377051625534794e+04, 1.8384940709941831e+04, 1.8392830169091238e+04, 1.8400720002842630e+04, 1.8408610211055744e+04, 1.8416500793590392e+04, 1.8424391750306535e+04, 1.8432283081064194e+04, 1.8440174785723524e+04, 1.8448066864144770e+04, 1.8455959316188291e+04, 1.8463852141714542e+04, 1.8471745340584082e+04, 1.8479638912657589e+04, 1.8487532857795828e+04, 1.8495427175859670e+04, 1.8503321866710095e+04, 1.8511216930208189e+04, 1.8519112366215137e+04, 1.8527008174592218e+04, 1.8534904355200830e+04, 1.8542800907902474e+04, 1.8550697832558741e+04, 1.8558595129031339e+04, 1.8566492797182062e+04, 1.8574390836872826e+04, 1.8582289247965637e+04, 1.8590188030322606e+04, 1.8598087183805950e+04, 1.8605986708277982e+04, 1.8613886603601124e+04, 1.8621786869637886e+04, 1.8629687506250906e+04, 1.8637588513302901e+04, 1.8645489890656692e+04, 1.8653391638175213e+04, 1.8661293755721490e+04, 1.8669196243158651e+04, 1.8677099100349933e+04, 1.8685002327158665e+04, 1.8692905923448277e+04, 1.8700809889082309e+04, 1.8708714223924395e+04, 1.8716618927838274e+04, 1.8724524000687772e+04, 1.8732429442336834e+04, 1.8740335252649489e+04, 1.8748241431489885e+04, 1.8756147978722252e+04, 1.8764054894210931e+04, 1.8771962177820358e+04, 1.8779869829415074e+04, 1.8787777848859707e+04, 1.8795686236018995e+04, 1.8803594990757778e+04, 1.8811504112940987e+04, 1.8819413602433662e+04, 1.8827323459100931e+04, 1.8835233682808033e+04, 1.8843144273420286e+04, 1.8851055230803133e+04, 1.8858966554822095e+04, 1.8866878245342796e+04, 1.8874790302230980e+04, 1.8882702725352454e+04, 1.8890615514573143e+04, 1.8898528669759071e+04, 1.8906442190776354e+04, 1.8914356077491211e+04, 1.8922270329769952e+04, 1.8930184947478992e+04, 1.8938099930484841e+04, 1.8946015278654100e+04, 1.8953930991853482e+04, 1.8961847069949785e+04, 1.8969763512809906e+04, 1.8977680320300846e+04, 1.8985597492289689e+04, 1.8993515028643636e+04, 1.9001432929229966e+04, 1.9009351193916067e+04, 1.9017269822569404e+04, 1.9025188815057569e+04, 1.9033108171248230e+04, 1.9041027891009155e+04, 1.9048947974208208e+04, 1.9056868420713352e+04, 1.9064789230392638e+04, 1.9072710403114226e+04, 1.9080631938746363e+04, 1.9088553837157382e+04, 1.9096476098215731e+04, 1.9104398721789948e+04, 1.9112321707748662e+04, 1.9120245055960590e+04, 1.9128168766294559e+04, 1.9136092838619483e+04, 1.9144017272804369e+04, 1.9151942068718326e+04, 1.9159867226230548e+04, 1.9167792745210336e+04, 1.9175718625527079e+04, 1.9183644867050250e+04, 1.9191571469649432e+04, 1.9199498433194298e+04, 1.9207425757554611e+04, 1.9215353442600226e+04, 1.9223281488201101e+04, 1.9231209894227282e+04, 1.9239138660548910e+04, 1.9247067787036216e+04, 1.9254997273559529e+04, 1.9262927119989272e+04, 1.9270857326195957e+04, 1.9278787892050190e+04, 1.9286718817422672e+04, 1.9294650102184198e+04, 1.9302581746205655e+04, 1.9310513749358015e+04, 1.9318446111512359e+04, 1.9326378832539838e+04, 1.9334311912311721e+04, 1.9342245350699348e+04, 1.9350179147574170e+04, 1.9358113302807702e+04, 1.9366047816271584e+04, 1.9373982687837528e+04, 1.9381917917377348e+04, 1.9389853504762934e+04, 1.9397789449866283e+04, 1.9405725752559487e+04, 1.9413662412714712e+04, 1.9421599430204227e+04, 1.9429536804900396e+04, 1.9437474536675654e+04, 1.9445412625402554e+04, 1.9453351070953715e+04, 1.9461289873201869e+04, 1.9469229032019823e+04, 1.9477168547280486e+04, 1.9485108418856849e+04, 1.9493048646621995e+04, 1.9500989230449097e+04, 1.9508930170211428e+04, 1.9516871465782337e+04, 1.9524813117035264e+04, 1.9532755123843755e+04, 1.9540697486081430e+04, 1.9548640203622006e+04, 1.9556583276339283e+04, 1.9564526704107153e+04, 1.9572470486799612e+04, 1.9580414624290726e+04, 1.9588359116454660e+04, 1.9596303963165661e+04, 1.9604249164298071e+04, 1.9612194719726325e+04, 1.9620140629324935e+04, 1.9628086892968520e+04, 1.9636033510531764e+04, 1.9643980481889455e+04, 1.9651927806916472e+04, 1.9659875485487777e+04, 1.9667823517478413e+04, 1.9675771902763525e+04, 1.9683720641218337e+04, 1.9691669732718165e+04, 1.9699619177138415e+04, 1.9707568974354577e+04, 1.9715519124242233e+04, 1.9723469626677041e+04, 1.9731420481534758e+04, 1.9739371688691230e+04, 1.9747323248022392e+04, 1.9755275159404242e+04, 1.9763227422712898e+04, 1.9771180037824550e+04, 1.9779133004615473e+04, 1.9787086322962034e+04, 1.9795039992740680e+04, 1.9802994013827956e+04, 1.9810948386100488e+04, 1.9818903109434988e+04, 1.9826858183708246e+04, 1.9834813608797162e+04, 1.9842769384578693e+04, 1.9850725510929908e+04, 1.9858681987727945e+04, 1.9866638814850034e+04, 1.9874595992173494e+04, 1.9882553519575726e+04, 1.9890511396934216e+04, 1.9898469624126537e+04, 1.9906428201030350e+04, 1.9914387127523401e+04, 1.9922346403483516e+04, 1.9930306028788618e+04, 1.9938266003316698e+04, 1.9946226326945845e+04, 1.9954186999554229e+04, 1.9962148021020115e+04, 1.9970109391221835e+04, 1.9978071110037818e+04, 1.9986033177346577e+04, 1.9993995593026691e+04, 2.0001958356956860e+04, 2.0009921469015841e+04, 2.0017884929082476e+04, 2.0025848737035707e+04, 2.0033812892754544e+04, 2.0041777396118101e+04, 2.0049742247005546e+04, 2.0057707445296161e+04, 2.0065672990869291e+04, 2.0073638883604377e+04, 2.0081605123380938e+04, 2.0089571710078577e+04, 2.0097538643576983e+04, 2.0105505923755929e+04, 2.0113473550495262e+04, 2.0121441523674930e+04, 2.0129409843174941e+04, 2.0137378508875408e+04, 2.0145347520656516e+04, 2.0153316878398535e+04, 2.0161286581981811e+04, 2.0169256631286789e+04, 2.0177227026193981e+04, 2.0185197766583991e+04, 2.0193168852337498e+04, 2.0201140283335266e+04, 2.0209112059458144e+04, 2.0217084180587066e+04, 2.0225056646603040e+04, 2.0233029457387162e+04, 2.0241002612820605e+04, 2.0248976112784625e+04, 2.0256949957160574e+04, 2.0264924145829857e+04, 2.0272898678673988e+04, 2.0280873555574544e+04, 2.0288848776413201e+04, 2.0296824341071697e+04, 2.0304800249431857e+04, 2.0312776501375607e+04, 2.0320753096784923e+04, 2.0328730035541885e+04, 2.0336707317528642e+04, 2.0344684942627424e+04, 2.0352662910720555e+04, 2.0360641221690421e+04, 2.0368619875419507e+04, 2.0376598871790356e+04, 2.0384578210685620e+04, 2.0392557891988010e+04, 2.0400537915580317e+04, 2.0408518281345430e+04, 2.0416498989166303e+04, 2.0424480038925969e+04, 2.0432461430507548e+04, 2.0440443163794240e+04, 2.0448425238669322e+04, 2.0456407655016148e+04, 2.0464390412718159e+04, 2.0472373511658872e+04, 2.0480356951721878e+04, 2.0488340732790854e+04, 2.0496324854749557e+04, 2.0504309317481820e+04, 2.0512294120871553e+04, 2.0520279264802753e+04, 2.0528264749159483e+04, 2.0536250573825906e+04, 2.0544236738686239e+04, 2.0552223243624790e+04, 2.0560210088525953e+04, 2.0568197273274185e+04, 2.0576184797754035e+04, 2.0584172661850120e+04, 2.0592160865447142e+04, 2.0600149408429883e+04, 2.0608138290683193e+04, 2.0616127512092004e+04, 2.0624117072541339e+04, 2.0632106971916281e+04, 2.0640097210102002e+04, 2.0648087786983742e+04, 2.0656078702446834e+04, 2.0664069956376672e+04, 2.0672061548658741e+04, 2.0680053479178594e+04, 2.0688045747821863e+04, 2.0696038354474262e+04, 2.0704031299021579e+04, 2.0712024581349677e+04, 2.0720018201344512e+04, 2.0728012158892085e+04, 2.0736006453878501e+04, 2.0744001086189932e+04, 2.0751996055712629e+04, 2.0759991362332919e+04, 2.0767987005937208e+04, 2.0775982986411971e+04, 2.0783979303643770e+04, 2.0791975957519229e+04, 2.0799972947925071e+04, 2.0807970274748070e+04, 2.0815967937875088e+04, 2.0823965937193068e+04, 2.0831964272589026e+04, 2.0839962943950042e+04, 2.0847961951163285e+04, 2.0855961294115998e+04, 2.0863960972695502e+04, 2.0871960986789181e+04, 2.0879961336284508e+04, 2.0887962021069023e+04, 2.0895963041030347e+04, 2.0903964396056170e+04, 2.0911966086034266e+04, 2.0919968110852489e+04, 2.0927970470398741e+04, 2.0935973164561030e+04, 2.0943976193227416e+04, 2.0951979556286049e+04, 2.0959983253625145e+04, 2.0967987285132996e+04, 2.0975991650697975e+04, 2.0983996350208527e+04, 2.0992001383553161e+04, 2.1000006750620480e+04, 2.1008012451299142e+04, 2.1016018485477893e+04, 2.1024024853045539e+04, 2.1032031553890982e+04, 2.1040038587903175e+04, 2.1048045954971160e+04, 2.1056053654984044e+04, 2.1064061687831014e+04, 2.1072070053401327e+04, 2.1080078751584319e+04, 2.1088087782269387e+04, 2.1096097145346019e+04, 2.1104106840703764e+04, 2.1112116868232246e+04, 2.1120127227821169e+04, 2.1128137919360299e+04, 2.1136148942739481e+04, 2.1144160297848648e+04, 2.1152171984577770e+04, 2.1160184002816932e+04, 2.1168196352456260e+04, 2.1176209033385967e+04, 2.1184222045496335e+04, 2.1192235388677724e+04, 2.1200249062820560e+04, 2.1208263067815336e+04, 2.1216277403552638e+04, 2.1224292069923100e+04, 2.1232307066817451e+04, 2.1240322394126470e+04, 2.1248338051741030e+04, 2.1256354039552054e+04, 2.1264370357450560e+04, 2.1272387005327619e+04, 2.1280403983074382e+04, 2.1288421290582068e+04, 2.1296438927741980e+04, 2.1304456894445473e+04, 2.1312475190583988e+04, 2.1320493816049035e+04, 2.1328512770732188e+04, 2.1336532054525102e+04, 2.1344551667319502e+04, 2.1352571609007176e+04, 2.1360591879479995e+04, 2.1368612478629890e+04, 2.1376633406348872e+04, 2.1384654662529018e+04, 2.1392676247062474e+04, 2.1400698159841464e+04, 2.1408720400758270e+04, 2.1416742969705258e+04, 2.1424765866574860e+04, 2.1432789091259576e+04, 2.1440812643651981e+04, 2.1448836523644717e+04, 2.1456860731130495e+04, 2.1464885266002100e+04, 2.1472910128152387e+04, 2.1480935317474279e+04, 2.1488960833860769e+04, 2.1496986677204921e+04, 2.1505012847399867e+04, 2.1513039344338813e+04, 2.1521066167915036e+04, 2.1529093318021871e+04, 2.1537120794552728e+04, 2.1545148597401094e+04, 2.1553176726460530e+04, 2.1561205181624646e+04, 2.1569233962787133e+04, 2.1577263069841756e+04, 2.1585292502682336e+04, 2.1593322261202778e+04, 2.1601352345297048e+04, 2.1609382754859173e+04, 2.1617413489783270e+04, 2.1625444549963515e+04, 2.1633475935294140e+04, 2.1641507645669462e+04, 2.1649539680983860e+04, 2.1657572041131789e+04, 2.1665604726007758e+04, 2.1673637735506352e+04, 2.1681671069522228e+04, 2.1689704727950117e+04, 2.1697738710684800e+04, 2.1705773017621137e+04, 2.1713807648654059e+04, 2.1721842603678560e+04, 2.1729877882589706e+04, 2.1737913485282628e+04, 2.1745949411652520e+04, 2.1753985661594652e+04, 2.1762022235004359e+04, 2.1770059131777045e+04, 2.1778096351808177e+04, 2.1786133894993298e+04, 2.1794171761228004e+04, 2.1802209950407974e+04, 2.1810248462428954e+04, 2.1818287297186744e+04, 2.1826326454577211e+04, 2.1834365934496313e+04, 2.1842405736840050e+04, 2.1850445861504497e+04, 2.1858486308385800e+04, 2.1866527077380168e+04, 2.1874568168383877e+04, 2.1882609581293269e+04, 2.1890651316004758e+04, 2.1898693372414818e+04, 2.1906735750419990e+04, 2.1914778449916885e+04, 2.1922821470802184e+04, 2.1930864812972628e+04, 2.1938908476325021e+04, 2.1946952460756242e+04, 2.1954996766163233e+04, 2.1963041392442999e+04, 2.1971086339492616e+04, 2.1979131607209220e+04, 2.1987177195490021e+04, 2.1995223104232293e+04, 2.2003269333333370e+04, 2.2011315882690655e+04, 2.2019362752201607e+04, 2.2027409941763784e+04, 2.2035457451274764e+04, 2.2043505280632227e+04, 2.2051553429733893e+04, 2.2059601898477562e+04, 2.2067650686761095e+04, 2.2075699794482425e+04, 2.2083749221539536e+04, 2.2091798967830484e+04, 2.2099849033253402e+04, 2.2107899417706467e+04, 2.2115950121087939e+04, 2.2124001143296129e+04, 2.2132052484229422e+04, 2.2140104143786266e+04, 2.2148156121865166e+04, 2.2156208418364702e+04, 2.2164261033183517e+04, 2.2172313966220314e+04, 2.2180367217373860e+04, 2.2188420786542996e+04, 2.2196474673626613e+04, 2.2204528878523673e+04, 2.2212583401133212e+04, 2.2220638241354311e+04, 2.2228693399086129e+04, 2.2236748874227887e+04, 2.2244804666678869e+04, 2.2252860776338410e+04, 2.2260917203105932e+04, 2.2268973946880909e+04, 2.2277031007562873e+04, 2.2285088385051433e+04, 2.2293146079246242e+04, 2.2301204090047049e+04, 2.2309262417353628e+04, 2.2317321061065850e+04, 2.2325380021083613e+04, 2.2333439297306919e+04, 2.2341498889635808e+04, 2.2349558797970381e+04, 2.2357619022210820e+04, 2.2365679562257355e+04, 2.2373740418010293e+04, 2.2381801589369985e+04, 2.2389863076236856e+04, 2.2397924878511392e+04, 2.2405986996094147e+04, 2.2414049428885730e+04, 2.2422112176786817e+04, 2.2430175239698143e+04, 2.2438238617520514e+04, 2.2446302310154777e+04, 2.2454366317501877e+04, 2.2462430639462789e+04, 2.2470495275938563e+04, 2.2478560226830312e+04, 2.2486625492039209e+04, 2.2494691071466492e+04, 2.2502756965013457e+04, 2.2510823172581462e+04, 2.2518889694071931e+04, 2.2526956529386349e+04, 2.2535023678426260e+04, 2.2543091141093271e+04, 2.2551158917289049e+04, 2.2559227006915331e+04, 2.2567295409873899e+04, 2.2575364126066612e+04, 2.2583433155395389e+04, 2.2591502497762198e+04, 2.2599572153069086e+04, 2.2607642121218145e+04, 2.2615712402111545e+04, 2.2623782995651494e+04, 2.2631853901740280e+04, 2.2639925120280252e+04, 2.2647996651173806e+04, 2.2656068494323412e+04, 2.2664140649631601e+04, 2.2672213117000952e+04, 2.2680285896334124e+04, 2.2688358987533818e+04, 2.2696432390502803e+04, 2.2704506105143912e+04, 2.2712580131360039e+04, 2.2720654469054130e+04, 2.2728729118129198e+04, 2.2736804078488312e+04, 2.2744879350034611e+04, 2.2752954932671280e+04, 2.2761030826301580e+04, 2.2769107030828818e+04, 2.2777183546156375e+04, 2.2785260372187669e+04, 2.2793337508826207e+04, 2.2801414955975539e+04, 2.2809492713539275e+04, 2.2817570781421095e+04, 2.2825649159524721e+04, 2.2833727847753951e+04, 2.2841806846012638e+04, 2.2849886154204691e+04, 2.2857965772234082e+04, 2.2866045700004837e+04, 2.2874125937421057e+04, 2.2882206484386879e+04, 2.2890287340806524e+04, 2.2898368506584247e+04, 2.2906449981624384e+04, 2.2914531765831322e+04, 2.2922613859109504e+04, 2.2930696261363431e+04, 2.2938778972497668e+04, 2.2946861992416838e+04, 2.2954945321025625e+04, 2.2963028958228770e+04, 2.2971112903931058e+04, 2.2979197158037361e+04, 2.2987281720452600e+04, 2.2995366591081740e+04, 2.3003451769829811e+04, 2.3011537256601918e+04, 2.3019623051303199e+04, 2.3027709153838867e+04, 2.3035795564114193e+04, 2.3043882282034498e+04, 2.3051969307505162e+04, 2.3060056640431638e+04, 2.3068144280719418e+04, 2.3076232228274057e+04, 2.3084320483001178e+04, 2.3092409044806453e+04, 2.3100497913595616e+04, 2.3108587089274457e+04, 2.3116676571748816e+04, 2.3124766360924605e+04, 2.3132856456707785e+04, 2.3140946859004380e+04, 2.3149037567720465e+04, 2.3157128582762176e+04, 2.3165219904035705e+04, 2.3173311531447307e+04, 2.3181403464903291e+04, 2.3189495704310015e+04, 2.3197588249573902e+04, 2.3205681100601440e+04, 2.3213774257299163e+04, 2.3221867719573664e+04, 2.3229961487331599e+04, 2.3238055560479665e+04, 2.3246149938924638e+04, 2.3254244622573337e+04, 2.3262339611332642e+04, 2.3270434905109480e+04, 2.3278530503810856e+04, 2.3286626407343821e+04, 2.3294722615615468e+04, 2.3302819128532967e+04, 2.3310915946003544e+04, 2.3319013067934462e+04, 2.3327110494233064e+04, 2.3335208224806724e+04, 2.3343306259562902e+04, 2.3351404598409092e+04, 2.3359503241252845e+04, 2.3367602188001791e+04, 2.3375701438563585e+04, 2.3383800992845958e+04, 2.3391900850756698e+04, 2.3400001012203629e+04, 2.3408101477094660e+04, 2.3416202245337736e+04, 2.3424303316840858e+04, 2.3432404691512085e+04, 2.3440506369259539e+04, 2.3448608349991398e+04, 2.3456710633615879e+04, 2.3464813220041266e+04, 2.3472916109175905e+04, 2.3481019300928194e+04, 2.3489122795206571e+04, 2.3497226591919552e+04, 2.3505330690975698e+04, 2.3513435092283620e+04, 2.3521539795751993e+04, 2.3529644801289542e+04, 2.3537750108805045e+04, 2.3545855718207345e+04, 2.3553961629405334e+04, 2.3562067842307955e+04, 2.3570174356824213e+04, 2.3578281172863157e+04, 2.3586388290333907e+04, 2.3594495709145631e+04, 2.3602603429207538e+04, 2.3610711450428913e+04, 2.3618819772719085e+04, 2.3626928395987441e+04, 2.3635037320143416e+04, 2.3643146545096512e+04, 2.3651256070756266e+04, 2.3659365897032279e+04, 2.3667476023834217e+04, 2.3675586451071795e+04, 2.3683697178654773e+04, 2.3691808206492962e+04, 2.3699919534496246e+04, 2.3708031162574556e+04, 2.3716143090637866e+04, 2.3724255318596217e+04, 2.3732367846359699e+04, 2.3740480673838450e+04, 2.3748593800942668e+04, 2.3756707227582610e+04, 2.3764820953668583e+04, 2.3772934979110938e+04, 2.3781049303820095e+04, 2.3789163927706515e+04, 2.3797278850680719e+04, 2.3805394072653278e+04, 2.3813509593534825e+04, 2.3821625413236037e+04, 2.3829741531667645e+04, 2.3837857948740442e+04, 2.3845974664365262e+04, 2.3854091678452998e+04, 2.3862208990914602e+04, 2.3870326601661069e+04, 2.3878444510603451e+04, 2.3886562717652858e+04, 2.3894681222720446e+04, 2.3902800025717428e+04, 2.3910919126555065e+04, 2.3919038525144679e+04, 2.3927158221397633e+04, 2.3935278215225357e+04, 2.3943398506539328e+04, 2.3951519095251067e+04, 2.3959639981272161e+04, 2.3967761164514239e+04, 2.3975882644888989e+04, 2.3984004422308153e+04, 2.3992126496683515e+04, 2.4000248867926923e+04, 2.4008371535950268e+04, 2.4016494500665503e+04, 2.4024617761984628e+04, 2.4032741319819685e+04, 2.4040865174082792e+04, 2.4048989324686099e+04, 2.4057113771541815e+04, 2.4065238514562199e+04, 2.4073363553659568e+04, 2.4081488888746280e+04, 2.4089614519734758e+04, 2.4097740446537464e+04, 2.4105866669066923e+04, 2.4113993187235705e+04, 2.4122120000956435e+04, 2.4130247110141780e+04, 2.4138374514704472e+04, 2.4146502214557291e+04, 2.4154630209613068e+04, 2.4162758499784672e+04, 2.4170887084985046e+04, 2.4179015965127172e+04, 2.4187145140124081e+04, 2.4195274609888871e+04, 2.4203404374334663e+04, 2.4211534433374658e+04, 2.4219664786922087e+04, 2.4227795434890249e+04, 2.4235926377192482e+04, 2.4244057613742178e+04, 2.4252189144452783e+04, 2.4260320969237786e+04, 2.4268453088010745e+04, 2.4276585500685247e+04, 2.4284718207174941e+04, 2.4292851207393520e+04, 2.4300984501254741e+04, 2.4309118088672403e+04, 2.4317251969560351e+04, 2.4325386143832489e+04, 2.4333520611402768e+04, 2.4341655372185181e+04, 2.4349790426093794e+04, 2.4357925773042702e+04, 2.4366061412946055e+04, 2.4374197345718061e+04, 2.4382333571272968e+04, 2.4390470089525083e+04, 2.4398606900388757e+04, 2.4406744003778400e+04, 2.4414881399608457e+04, 2.4423019087793433e+04, 2.4431157068247881e+04, 2.4439295340886412e+04, 2.4447433905623675e+04, 2.4455572762374370e+04, 2.4463711911053255e+04, 2.4471851351575129e+04, 2.4479991083854846e+04, 2.4488131107807305e+04, 2.4496271423347469e+04, 2.4504412030390326e+04, 2.4512552928850935e+04, 2.4520694118644395e+04, 2.4528835599685852e+04, 2.4536977371890505e+04, 2.4545119435173609e+04, 2.4553261789450458e+04, 2.4561404434636403e+04, 2.4569547370646833e+04, 2.4577690597397199e+04, 2.4585834114803001e+04, 2.4593977922779770e+04, 2.4602122021243111e+04, 2.4610266410108659e+04, 2.4618411089292105e+04, 2.4626556058709193e+04, 2.4634701318275711e+04, 2.4642846867907494e+04, 2.4650992707520436e+04, 2.4659138837030459e+04, 2.4667285256353556e+04, 2.4675431965405758e+04, 2.4683578964103148e+04, 2.4691726252361859e+04, 2.4699873830098062e+04, 2.4708021697227985e+04, 2.4716169853667907e+04, 2.4724318299334151e+04, 2.4732467034143086e+04, 2.4740616058011143e+04, 2.4748765370854777e+04, 2.4756914972590519e+04, 2.4765064863134921e+04, 2.4773215042404601e+04, 2.4781365510316224e+04, 2.4789516266786501e+04, 2.4797667311732184e+04, 2.4805818645070089e+04, 2.4813970266717060e+04, 2.4822122176590001e+04, 2.4830274374605862e+04, 2.4838426860681644e+04, 2.4846579634734389e+04, 2.4854732696681189e+04, 2.4862886046439187e+04, 2.4871039683925574e+04, 2.4879193609057584e+04, 2.4887347821752497e+04, 2.4895502321927648e+04, 2.4903657109500418e+04, 2.4911812184388229e+04, 2.4919967546508557e+04, 2.4928123195778928e+04, 2.4936279132116902e+04, 2.4944435355440091e+04, 2.4952591865666171e+04, 2.4960748662712849e+04, 2.4968905746497876e+04, 2.4977063116939062e+04, 2.4985220773954261e+04, 2.4993378717461368e+04, 2.5001536947378328e+04, 2.5009695463623135e+04, 2.5017854266113827e+04, 2.5026013354768493e+04, 2.5034172729505273e+04, 2.5042332390242336e+04, 2.5050492336897911e+04, 2.5058652569390280e+04, 2.5066813087637755e+04, 2.5074973891558711e+04, 2.5083134981071555e+04, 2.5091296356094754e+04, 2.5099458016546811e+04, 2.5107619962346282e+04, 2.5115782193411764e+04, 2.5123944709661901e+04, 2.5132107511015394e+04, 2.5140270597390980e+04, 2.5148433968707443e+04, 2.5156597624883605e+04, 2.5164761565838362e+04, 2.5172925791490627e+04, 2.5181090301759374e+04, 2.5189255096563622e+04, 2.5197420175822426e+04, 2.5205585539454900e+04, 2.5213751187380200e+04, 2.5221917119517519e+04, 2.5230083335786116e+04, 2.5238249836105271e+04, 2.5246416620394331e+04, 2.5254583688572671e+04, 2.5262751040559724e+04, 2.5270918676274970e+04, 2.5279086595637928e+04, 2.5287254798568167e+04, 2.5295423284985296e+04, 2.5303592054808967e+04, 2.5311761107958897e+04, 2.5319930444354821e+04, 2.5328100063916547e+04, 2.5336269966563901e+04, 2.5344440152216783e+04, 2.5352610620795112e+04, 2.5360781372218869e+04, 2.5368952406408076e+04, 2.5377123723282799e+04, 2.5385295322763141e+04, 2.5393467204769266e+04, 2.5401639369221379e+04, 2.5409811816039724e+04, 2.5417984545144591e+04, 2.5426157556456317e+04, 2.5434330849895283e+04, 2.5442504425381914e+04, 2.5450678282836689e+04, 2.5458852422180120e+04, 2.5467026843332762e+04, 2.5475201546215230e+04, 2.5483376530748174e+04, 2.5491551796852284e+04, 2.5499727344448303e+04, 2.5507903173457020e+04, 2.5516079283799256e+04, 2.5524255675395889e+04, 2.5532432348167840e+04, 2.5540609302036064e+04, 2.5548786536921580e+04, 2.5556964052745425e+04, 2.5565141849428706e+04, 2.5573319926892556e+04, 2.5581498285058162e+04, 2.5589676923846753e+04, 2.5597855843179601e+04, 2.5606035042978023e+04, 2.5614214523163380e+04, 2.5622394283657082e+04, 2.5630574324380574e+04, 2.5638754645255347e+04, 2.5646935246202942e+04, 2.5655116127144938e+04, 2.5663297288002963e+04, 2.5671478728698683e+04, 2.5679660449153809e+04, 2.5687842449290103e+04, 2.5696024729029366e+04, 2.5704207288293430e+04, 2.5712390127004201e+04, 2.5720573245083593e+04, 2.5728756642453591e+04, 2.5736940319036214e+04, 2.5745124274753518e+04, 2.5753308509527607e+04, 2.5761493023280640e+04, 2.5769677815934807e+04, 2.5777862887412339e+04, 2.5786048237635518e+04, 2.5794233866526665e+04, 2.5802419774008147e+04, 2.5810605960002369e+04, 2.5818792424431795e+04, 2.5826979167218909e+04, 2.5835166188286250e+04, 2.5843353487556404e+04, 2.5851541064951994e+04, 2.5859728920395690e+04, 2.5867917053810204e+04, 2.5876105465118282e+04, 2.5884294154242729e+04, 2.5892483121106376e+04, 2.5900672365632116e+04, 2.5908861887742863e+04, 2.5917051687361589e+04, 2.5925241764411312e+04, 2.5933432118815075e+04, 2.5941622750495979e+04, 2.5949813659377167e+04, 2.5958004845381809e+04, 2.5966196308433136e+04, 2.5974388048454413e+04, 2.5982580065368951e+04, 2.5990772359100101e+04, 2.5998964929571252e+04, 2.6007157776705844e+04, 2.6015350900427358e+04, 2.6023544300659312e+04, 2.6031737977325261e+04, 2.6039931930348830e+04, 2.6048126159653650e+04, 2.6056320665163414e+04, 2.6064515446801855e+04, 2.6072710504492752e+04, 2.6080905838159913e+04, 2.6089101447727204e+04, 2.6097297333118517e+04, 2.6105493494257797e+04, 2.6113689931069031e+04, 2.6121886643476242e+04, 2.6130083631403504e+04, 2.6138280894774918e+04, 2.6146478433514640e+04, 2.6154676247546860e+04, 2.6162874336795820e+04, 2.6171072701185785e+04, 2.6179271340641084e+04, 2.6187470255086067e+04, 2.6195669444445146e+04, 2.6203868908642758e+04, 2.6212068647603388e+04, 2.6220268661251564e+04, 2.6228468949511851e+04, 2.6236669512308861e+04, 2.6244870349567238e+04, 2.6253071461211683e+04, 2.6261272847166918e+04, 2.6269474507357729e+04, 2.6277676441708925e+04, 2.6285878650145358e+04, 2.6294081132591935e+04, 2.6302283888973594e+04, 2.6310486919215309e+04, 2.6318690223242102e+04, 2.6326893800979040e+04, 2.6335097652351225e+04, 2.6343301777283799e+04, 2.6351506175701950e+04, 2.6359710847530903e+04, 2.6367915792695923e+04, 2.6376121011122315e+04, 2.6384326502735435e+04, 2.6392532267460669e+04, 2.6400738305223447e+04, 2.6408944615949244e+04, 2.6417151199563563e+04, 2.6425358055991965e+04, 2.6433565185160034e+04, 2.6441772586993411e+04, 2.6449980261417764e+04, 2.6458188208358813e+04, 2.6466396427742311e+04, 2.6474604919494050e+04, 2.6482813683539873e+04, 2.6491022719805649e+04, 2.6499232028217291e+04, 2.6507441608700770e+04, 2.6515651461182071e+04, 2.6523861585587238e+04, 2.6532071981842342e+04, 2.6540282649873508e+04, 2.6548493589606882e+04, 2.6556704800968673e+04, 2.6564916283885119e+04, 2.6573128038282495e+04, 2.6581340064087122e+04, 2.6589552361225349e+04, 2.6597764929623587e+04, 2.6605977769208264e+04, 2.6614190879905858e+04, 2.6622404261642896e+04, 2.6630617914345927e+04, 2.6638831837941550e+04, 2.6647046032356404e+04, 2.6655260497517160e+04, 2.6663475233350542e+04, 2.6671690239783304e+04, 2.6679905516742241e+04, 2.6688121064154191e+04, 2.6696336881946023e+04, 2.6704552970044653e+04, 2.6712769328377042e+04, 2.6720985956870176e+04, 2.6729202855451087e+04, 2.6737420024046853e+04, 2.6745637462584582e+04, 2.6753855170991428e+04, 2.6762073149194581e+04, 2.6770291397121266e+04, 2.6778509914698760e+04, 2.6786728701854361e+04, 2.6794947758515420e+04, 2.6803167084609329e+04, 2.6811386680063508e+04, 2.6819606544805418e+04, 2.6827826678762569e+04, 2.6836047081862504e+04, 2.6844267754032797e+04, 2.6852488695201078e+04, 2.6860709905295003e+04, 2.6868931384242271e+04, 2.6877153131970616e+04, 2.6885375148407817e+04, 2.6893597433481689e+04, 2.6901819987120085e+04, 2.6910042809250899e+04, 2.6918265899802060e+04, 2.6926489258701535e+04, 2.6934712885877339e+04, 2.6942936781257518e+04, 2.6951160944770156e+04, 2.6959385376343380e+04, 2.6967610075905344e+04, 2.6975835043384257e+04, 2.6984060278708359e+04, 2.6992285781805931e+04, 2.7000511552605280e+04, 2.7008737591034766e+04, 2.7016963897022779e+04, 2.7025190470497753e+04, 2.7033417311388162e+04, 2.7041644419622513e+04, 2.7049871795129347e+04, 2.7058099437837252e+04, 2.7066327347674855e+04, 2.7074555524570806e+04, 2.7082783968453808e+04, 2.7091012679252606e+04, 2.7099241656895967e+04, 2.7107470901312703e+04, 2.7115700412431670e+04, 2.7123930190181749e+04, 2.7132160234491876e+04, 2.7140390545291008e+04, 2.7148621122508157e+04, 2.7156851966072350e+04, 2.7165083075912680e+04, 2.7173314451958257e+04, 2.7181546094138230e+04, 2.7189778002381790e+04, 2.7198010176618176e+04, 2.7206242616776643e+04, 2.7214475322786508e+04, 2.7222708294577100e+04, 2.7230941532077806e+04, 2.7239175035218042e+04, 2.7247408803927257e+04, 2.7255642838134947e+04, 2.7263877137770640e+04, 2.7272111702763905e+04, 2.7280346533044351e+04, 2.7288581628541608e+04, 2.7296816989185358e+04, 2.7305052614905322e+04, 2.7313288505631252e+04, 2.7321524661292933e+04, 2.7329761081820197e+04, 2.7337997767142911e+04, 2.7346234717190971e+04, 2.7354471931894321e+04, 2.7362709411182932e+04, 2.7370947154986825e+04, 2.7379185163236041e+04, 2.7387423435860681e+04, 2.7395661972790851e+04, 2.7403900773956721e+04, 2.7412139839288491e+04, 2.7420379168716390e+04, 2.7428618762170696e+04, 2.7436858619581715e+04, 2.7445098740879792e+04, 2.7453339125995306e+04, 2.7461579774858681e+04, 2.7469820687400374e+04, 2.7478061863550869e+04, 2.7486303303240697e+04, 2.7494545006400425e+04, 2.7502786972960661e+04, 2.7511029202852031e+04, 2.7519271696005224e+04, 2.7527514452350937e+04, 2.7535757471819925e+04, 2.7544000754342971e+04, 2.7552244299850896e+04, 2.7560488108274560e+04, 2.7568732179544859e+04, 2.7576976513592712e+04, 2.7585221110349092e+04, 2.7593465969745004e+04, 2.7601711091711484e+04, 2.7609956476179603e+04, 2.7618202123080475e+04, 2.7626448032345248e+04, 2.7634694203905106e+04, 2.7642940637691267e+04, 2.7651187333634989e+04, 2.7659434291667560e+04, 2.7667681511720300e+04, 2.7675928993724589e+04, 2.7684176737611811e+04, 2.7692424743313415e+04, 2.7700673010760860e+04, 2.7708921539885661e+04, 2.7717170330619352e+04, 2.7725419382893524e+04, 2.7733668696639783e+04, 2.7741918271789786e+04, 2.7750168108275215e+04, 2.7758418206027785e+04, 2.7766668564979263e+04, 2.7774919185061441e+04, 2.7783170066206138e+04, 2.7791421208345226e+04, 2.7799672611410606e+04, 2.7807924275334211e+04, 2.7816176200048012e+04, 2.7824428385484018e+04, 2.7832680831574264e+04, 2.7840933538250829e+04, 2.7849186505445836e+04, 2.7857439733091414e+04, 2.7865693221119764e+04, 2.7873946969463090e+04, 2.7882200978053654e+04, 2.7890455246823742e+04, 2.7898709775705684e+04, 2.7906964564631835e+04, 2.7915219613534588e+04, 2.7923474922346373e+04, 2.7931730490999656e+04, 2.7939986319426935e+04, 2.7948242407560749e+04, 2.7956498755333669e+04, 2.7964755362678294e+04, 2.7973012229527274e+04, 2.7981269355813271e+04, 2.7989526741468995e+04, 2.7997784386427204e+04, 2.8006042290620670e+04, 2.8014300453982207e+04, 2.8022558876444666e+04, 2.8030817557940933e+04, 2.8039076498403923e+04, 2.8047335697766590e+04, 2.8055595155961924e+04, 2.8063854872922944e+04, 2.8072114848582714e+04, 2.8080375082874318e+04, 2.8088635575730888e+04, 2.8096896327085589e+04, 2.8105157336871616e+04, 2.8113418605022191e+04, 2.8121680131470584e+04, 2.8129941916150103e+04, 2.8138203958994069e+04, 2.8146466259935856e+04, 2.8154728818908869e+04, 2.8162991635846542e+04, 2.8171254710682340e+04, 2.8179518043349784e+04, 2.8187781633782404e+04, 2.8196045481913774e+04, 2.8204309587677504e+04, 2.8212573951007234e+04, 2.8220838571836648e+04, 2.8229103450099443e+04, 2.8237368585729382e+04, 2.8245633978660233e+04, 2.8253899628825809e+04, 2.8262165536159966e+04, 2.8270431700596579e+04, 2.8278698122069560e+04, 2.8286964800512866e+04, 2.8295231735860478e+04, 2.8303498928046411e+04, 2.8311766377004711e+04, 2.8320034082669474e+04, 2.8328302044974811e+04, 2.8336570263854876e+04, 2.8344838739243856e+04, 2.8353107471075975e+04, 2.8361376459285482e+04, 2.8369645703806666e+04, 2.8377915204573848e+04, 2.8386184961521383e+04, 2.8394454974583656e+04, 2.8402725243695091e+04, 2.8410995768790148e+04, 2.8419266549803313e+04, 2.8427537586669103e+04, 2.8435808879322081e+04, 2.8444080427696837e+04, 2.8452352231727993e+04, 2.8460624291350206e+04, 2.8468896606498161e+04, 2.8477169177106585e+04, 2.8485442003110235e+04, 2.8493715084443898e+04, 2.8501988421042406e+04, 2.8510262012840605e+04, 2.8518535859773387e+04, 2.8526809961775682e+04, 2.8535084318782439e+04, 2.8543358930728649e+04, 2.8551633797549333e+04, 2.8559908919179550e+04, 2.8568184295554383e+04, 2.8576459926608964e+04, 2.8584735812278439e+04, 2.8593011952497996e+04, 2.8601288347202859e+04, 2.8609564996328278e+04, 2.8617841899809544e+04, 2.8626119057581971e+04, 2.8634396469580923e+04, 2.8642674135741774e+04, 2.8650952055999944e+04, 2.8659230230290887e+04, 2.8667508658550090e+04, 2.8675787340713061e+04, 2.8684066276715355e+04, 2.8692345466492545e+04, 2.8700624909980259e+04, 2.8708904607114135e+04, 2.8717184557829860e+04, 2.8725464762063137e+04, 2.8733745219749715e+04, 2.8742025930825381e+04, 2.8750306895225935e+04, 2.8758588112887221e+04, 2.8766869583745116e+04, 2.8775151307735527e+04, 2.8783433284794395e+04, 2.8791715514857689e+04, 2.8799997997861421e+04, 2.8808280733741620e+04, 2.8816563722434363e+04, 2.8824846963875749e+04, 2.8833130458001910e+04, 2.8841414204749017e+04, 2.8849698204053264e+04, 2.8857982455850884e+04, 2.8866266960078145e+04, 2.8874551716671329e+04, 2.8882836725566784e+04, 2.8891121986700851e+04, 2.8899407500009933e+04, 2.8907693265430447e+04, 2.8915979282898850e+04, 2.8924265552351633e+04, 2.8932552073725317e+04, 2.8940838846956445e+04, 2.8949125871981611e+04, 2.8957413148737425e+04, 2.8965700677160534e+04, 2.8973988457187628e+04, 2.8982276488755400e+04, 2.8990564771800608e+04, 2.8998853306260025e+04, 2.9007142092070451e+04, 2.9015431129168734e+04, 2.9023720417491731e+04, 2.9032009956976355e+04, 2.9040299747559537e+04, 2.9048589789178241e+04, 2.9056880081769465e+04, 2.9065170625270239e+04, 2.9073461419617619e+04, 2.9081752464748701e+04, 2.9090043760600609e+04, 2.9098335307110490e+04, 2.9106627104215539e+04, 2.9114919151852973e+04, 2.9123211449960036e+04, 2.9131503998474011e+04, 2.9139796797332212e+04, 2.9148089846471976e+04, 2.9156383145830692e+04, 2.9164676695345752e+04, 2.9172970494954599e+04, 2.9181264544594706e+04, 2.9189558844203559e+04, 2.9197853393718702e+04, 2.9206148193077690e+04, 2.9214443242218127e+04, 2.9222738541077626e+04, 2.9231034089593850e+04, 2.9239329887704487e+04, 2.9247625935347252e+04, 2.9255922232459892e+04, 2.9264218778980194e+04, 2.9272515574845962e+04, 2.9280812619995046e+04, 2.9289109914365312e+04, 2.9297407457894667e+04, 2.9305705250521049e+04, 2.9314003292182420e+04, 2.9322301582816777e+04, 2.9330600122362153e+04, 2.9338898910756601e+04, 2.9347197947938217e+04, 2.9355497233845112e+04, 2.9363796768415446e+04, 2.9372096551587394e+04, 2.9380396583299178e+04, 2.9388696863489029e+04, 2.9396997392095229e+04, 2.9405298169056077e+04, 2.9413599194309914e+04, 2.9421900467795105e+04, 2.9430201989450048e+04, 2.9438503759213159e+04, 2.9446805777022917e+04, 2.9455108042817792e+04, 2.9463410556536306e+04, 2.9471713318117007e+04, 2.9480016327498481e+04, 2.9488319584619330e+04, 2.9496623089418204e+04, 2.9504926841833767e+04, 2.9513230841804725e+04, 2.9521535089269801e+04, 2.9529839584167770e+04, 2.9538144326437410e+04, 2.9546449316017548e+04, 2.9554754552847044e+04, 2.9563060036864772e+04, 2.9571365768009648e+04, 2.9579671746220614e+04, 2.9587977971436645e+04, 2.9596284443596745e+04, 2.9604591162639950e+04, 2.9612898128505323e+04, 2.9621205341131950e+04, 2.9629512800458961e+04, 2.9637820506425513e+04, 2.9646128458970783e+04, 2.9654436658033988e+04, 2.9662745103554371e+04, 2.9671053795471209e+04, 2.9679362733723803e+04, 2.9687671918251490e+04, 2.9695981348993628e+04, 2.9704291025889615e+04, 2.9712600948878877e+04, 2.9720911117900858e+04, 2.9729221532895048e+04, 2.9737532193800955e+04, 2.9745843100558122e+04, 2.9754154253106124e+04, 2.9762465651384558e+04, 2.9770777295333060e+04, 2.9779089184891294e+04, 2.9787401319998942e+04, 2.9795713700595727e+04, 2.9804026326621402e+04, 2.9812339198015747e+04, 2.9820652314718569e+04, 2.9828965676669701e+04, 2.9837279283809021e+04, 2.9845593136076415e+04, 2.9853907233411825e+04, 2.9862221575755193e+04, 2.9870536163046516e+04, 2.9878850995225799e+04, 2.9887166072233089e+04, 2.9895481394008468e+04, 2.9903796960492029e+04, 2.9912112771623917e+04, 2.9920428827344280e+04, 2.9928745127593316e+04, 2.9937061672311251e+04, 2.9945378461438322e+04, 2.9953695494914809e+04, 2.9962012772681031e+04, 2.9970330294677318e+04, 2.9978648060844036e+04, 2.9986966071121587e+04, 2.9995284325450386e+04, 3.0003602823770892e+04, 3.0011921566023582e+04, 3.0020240552148971e+04, 3.0028559782087603e+04, 3.0036879255780044e+04, 3.0045198973166895e+04, 3.0053518934188778e+04, 3.0061839138786363e+04, 3.0070159586900314e+04, 3.0078480278471365e+04, 3.0086801213440245e+04, 3.0095122391747736e+04, 3.0103443813334634e+04, 3.0111765478141766e+04, 3.0120087386109997e+04, 3.0128409537180207e+04, 3.0136731931293321e+04, 3.0145054568390275e+04, 3.0153377448412048e+04, 3.0161700571299632e+04, 3.0170023936994072e+04, 3.0178347545436416e+04, 3.0186671396567755e+04, 3.0194995490329202e+04, 3.0203319826661907e+04, 3.0211644405507046e+04, 3.0219969226805813e+04, 3.0228294290499445e+04, 3.0236619596529199e+04, 3.0244945144836358e+04, 3.0253270935362245e+04, 3.0261596968048201e+04, 3.0269923242835597e+04, 3.0278249759665836e+04, 3.0286576518480349e+04, 3.0294903519220592e+04, 3.0303230761828050e+04, 3.0311558246244233e+04, 3.0319885972410695e+04, 3.0328213940268997e+04, 3.0336542149760746e+04, 3.0344870600827566e+04, 3.0353199293411108e+04, 3.0361528227453065e+04, 3.0369857402895141e+04, 3.0378186819679082e+04, 3.0386516477746653e+04, 3.0394846377039648e+04, 3.0403176517499891e+04, 3.0411506899069242e+04, 3.0419837521689577e+04, 3.0428168385302801e+04, 3.0436499489850852e+04, 3.0444830835275694e+04, 3.0453162421519326e+04, 3.0461494248523763e+04, 3.0469826316231054e+04, 3.0478158624583273e+04, 3.0486491173522525e+04, 3.0494823962990940e+04, 3.0503156992930682e+04, 3.0511490263283940e+04, 3.0519823773992925e+04, 3.0528157524999879e+04, 3.0536491516247072e+04, 3.0544825747676812e+04, 3.0553160219231409e+04, 3.0561494930853234e+04, 3.0569829882484657e+04, 3.0578165074068093e+04, 3.0586500505545973e+04, 3.0594836176860768e+04, 3.0603172087954968e+04, 3.0611508238771083e+04, 3.0619844629251671e+04, 3.0628181259339312e+04, 3.0636518128976601e+04, 3.0644855238106164e+04, 3.0653192586670659e+04, 3.0661530174612777e+04, 3.0669868001875224e+04, 3.0678206068400741e+04, 3.0686544374132100e+04, 3.0694882919012081e+04, 3.0703221702983523e+04, 3.0711560725989268e+04, 3.0719899987972192e+04, 3.0728239488875199e+04, 3.0736579228641218e+04, 3.0744919207213210e+04, 3.0753259424534157e+04, 3.0761599880547070e+04, 3.0769940575194996e+04, 3.0778281508420994e+04, 3.0786622680168166e+04, 3.0794964090379628e+04, 3.0803305738998530e+04, 3.0811647625968049e+04, 3.0819989751231384e+04, 3.0828332114731766e+04, 3.0836674716412454e+04, 3.0845017556216724e+04, 3.0853360634087894e+04, 3.0861703949969298e+04, 3.0870047503804304e+04, 3.0878391295536301e+04, 3.0886735325108708e+04, 3.0895079592464965e+04, 3.0903424097548555e+04, 3.0911768840302972e+04, 3.0920113820671744e+04, 3.0928459038598419e+04, 3.0936804494026583e+04, 3.0945150186899831e+04, 3.0953496117161812e+04, 3.0961842284756178e+04, 3.0970188689626611e+04, 3.0978535331716834e+04, 3.0986882210970580e+04, 3.0995229327331617e+04, 3.1003576680743739e+04, 3.1011924271150769e+04, 3.1020272098496549e+04, 3.1028620162724954e+04, 3.1036968463779893e+04, 3.1045317001605279e+04, 3.1053665776145070e+04, 3.1062014787343243e+04, 3.1070364035143808e+04, 3.1078713519490797e+04, 3.1087063240328269e+04, 3.1095413197600312e+04, 3.1103763391251032e+04, 3.1112113821224564e+04, 3.1120464487465084e+04, 3.1128815389916781e+04, 3.1137166528523863e+04, 3.1145517903230586e+04, 3.1153869513981212e+04, 3.1162221360720039e+04, 3.1170573443391393e+04, 3.1178925761939619e+04, 3.1187278316309093e+04, 3.1195631106444220e+04, 3.1203984132289421e+04, 3.1212337393789156e+04, 3.1220690890887898e+04, 3.1229044623530161e+04, 3.1237398591660476e+04, 3.1245752795223398e+04, 3.1254107234163512e+04, 3.1262461908425430e+04, 3.1270816817953786e+04, 3.1279171962693250e+04, 3.1287527342588503e+04, 3.1295882957584265e+04, 3.1304238807625272e+04, 3.1312594892656296e+04, 3.1320951212622123e+04, 3.1329307767467577e+04, 3.1337664557137497e+04, 3.1346021581576762e+04, 3.1354378840730260e+04, 3.1362736334542918e+04, 3.1371094062959681e+04, 3.1379452025925526e+04, 3.1387810223385455e+04, 3.1396168655284488e+04, 3.1404527321567683e+04, 3.1412886222180103e+04, 3.1421245357066859e+04, 3.1429604726173082e+04, 3.1437964329443923e+04, 3.1446324166824565e+04, 3.1454684238260208e+04, 3.1463044543696087e+04, 3.1471405083077461e+04, 3.1479765856349601e+04, 3.1488126863457834e+04, 3.1496488104347478e+04, 3.1504849578963895e+04, 3.1513211287252467e+04, 3.1521573229158614e+04, 3.1529935404627759e+04, 3.1538297813605375e+04, 3.1546660456036941e+04, 3.1555023331867967e+04, 3.1563386441044007e+04, 3.1571749783510604e+04, 3.1580113359213352e+04, 3.1588477168097870e+04, 3.1596841210109793e+04, 3.1605205485194783e+04, 3.1613569993298534e+04, 3.1621934734366754e+04, 3.1630299708345188e+04, 3.1638664915179612e+04, 3.1647030354815801e+04, 3.1655396027199575e+04, 3.1663761932276775e+04, 3.1672128069993272e+04, 3.1680494440294955e+04, 3.1688861043127738e+04, 3.1697227878437570e+04, 3.1705594946170408e+04, 3.1713962246272247e+04, 3.1722329778689109e+04, 3.1730697543367034e+04, 3.1739065540252090e+04, 3.1747433769290365e+04, 3.1755802230427980e+04, 3.1764170923611080e+04, 3.1772539848785826e+04, 3.1780909005898415e+04, 3.1789278394895060e+04, 3.1797648015722010e+04, 3.1806017868325529e+04, 3.1814387952651909e+04, 3.1822758268647460e+04, 3.1831128816258533e+04, 3.1839499595431498e+04, 3.1847870606112730e+04, 3.1856241848248665e+04, 3.1864613321785731e+04, 3.1872985026670394e+04, 3.1881356962849153e+04, 3.1889729130268519e+04, 3.1898101528875031e+04, 3.1906474158615256e+04, 3.1914847019435780e+04, 3.1923220111283223e+04, 3.1931593434104219e+04, 3.1939966987845437e+04, 3.1948340772453554e+04, 3.1956714787875298e+04, 3.1965089034057393e+04, 3.1973463510946611e+04, 3.1981838218489727e+04, 3.1990213156633563e+04, 3.1998588325324952e+04, 3.2006963724510748e+04, 3.2015339354137843e+04, 3.2023715214153144e+04, 3.2032091304503578e+04, 3.2040467625136112e+04, 3.2048844175997725e+04, 3.2057220957035424e+04, 3.2065597968196242e+04, 3.2073975209427230e+04, 3.2082352680675471e+04, 3.2090730381888068e+04, 3.2099108313012151e+04, 3.2107486473994875e+04, 3.2115864864783409e+04, 3.2124243485324965e+04, 3.2132622335566757e+04, 3.2141001415456045e+04, 3.2149380724940096e+04, 3.2157760263966215e+04, 3.2166140032481722e+04, 3.2174520030433960e+04, 3.2182900257770303e+04, 3.2191280714438144e+04, 3.2199661400384906e+04, 3.2208042315558032e+04, 3.2216423459904985e+04, 3.2224804833373262e+04, 3.2233186435910371e+04, 3.2241568267463859e+04, 3.2249950327981285e+04, 3.2258332617410233e+04, 3.2266715135698323e+04, 3.2275097882793183e+04, 3.2283480858642477e+04, 3.2291864063193891e+04, 3.2300247496395128e+04, 3.2308631158193919e+04, 3.2317015048538022e+04, 3.2325399167375213e+04, 3.2333783514653296e+04, 3.2342168090320101e+04, 3.2350552894323471e+04, 3.2358937926611281e+04, 3.2367323187131442e+04, 3.2375708675831862e+04, 3.2384094392660489e+04, 3.2392480337565295e+04, 3.2400866510494270e+04, 3.2409252911395437e+04, 3.2417639540216835e+04, 3.2426026396906524e+04, 3.2434413481412590e+04, 3.2442800793683156e+04, 3.2451188333666349e+04, 3.2459576101310322e+04, 3.2467964096563268e+04, 3.2476352319373385e+04, 3.2484740769688909e+04, 3.2493129447458086e+04, 3.2501518352629202e+04, 3.2509907485150554e+04, 3.2518296844970460e+04, 3.2526686432037273e+04, 3.2535076246299359e+04, 3.2543466287705116e+04, 3.2551856556202954e+04, 3.2560247051741328e+04, 3.2568637774268689e+04, 3.2577028723733529e+04, 3.2585419900084362e+04, 3.2593811303269718e+04, 3.2602202933238161e+04, 3.2610594789938266e+04, 3.2618986873318638e+04, 3.2627379183327914e+04, 3.2635771719914730e+04, 3.2644164483027766e+04, 3.2652557472615721e+04, 3.2660950688627319e+04, 3.2669344131011305e+04, 3.2677737799716437e+04, 3.2686131694691503e+04, 3.2694525815885328e+04, 3.2702920163246748e+04, 3.2711314736724617e+04, 3.2719709536267819e+04, 3.2728104561825265e+04, 3.2736499813345872e+04, 3.2744895290778604e+04, 3.2753290994072431e+04, 3.2761686923176356e+04, 3.2770083078039388e+04, 3.2778479458610585e+04, 3.2786876064839009e+04, 3.2795272896673756e+04, 3.2803669954063924e+04, 3.2812067236958675e+04, 3.2820464745307145e+04, 3.2828862479058524e+04, 3.2837260438162019e+04, 3.2845658622566851e+04, 3.2854057032222285e+04, 3.2862455667077571e+04, 3.2870854527082025e+04, 3.2879253612184963e+04, 3.2887652922335728e+04, 3.2896052457483675e+04, 3.2904452217578197e+04, 3.2912852202568705e+04, 3.2921252412404639e+04, 3.2929652847035439e+04, 3.2938053506410601e+04, 3.2946454390479616e+04, 3.2954855499192010e+04, 3.2963256832497325e+04, 3.2971658390345139e+04, 3.2980060172685044e+04, 3.2988462179466645e+04, 3.2996864410639595e+04, 3.3005266866153543e+04, 3.3013669545958168e+04, 3.3022072450003179e+04, 3.3030475578238307e+04, 3.3038878930613304e+04, 3.3047282507077936e+04, 3.3055686307581993e+04, 3.3064090332075306e+04, 3.3072494580507700e+04, 3.3080899052829060e+04, 3.3089303748989238e+04, 3.3097708668938176e+04, 3.3106113812625772e+04, 3.3114519180002011e+04, 3.3122924771016842e+04, 3.3131330585620286e+04, 3.3139736623762343e+04, 3.3148142885393048e+04, 3.3156549370462482e+04, 3.3164956078920724e+04, 3.3173363010717883e+04, 3.3181770165804097e+04, 3.3190177544129510e+04, 3.3198585145644291e+04, 3.3206992970298648e+04, 3.3215401018042809e+04, 3.3223809288827004e+04, 3.3232217782601503e+04, 3.3240626499316590e+04, 3.3249035438922561e+04, 3.3257444601369760e+04, 3.3265853986608541e+04, 3.3274263594589291e+04, 3.3282673425262372e+04, 3.3291083478578228e+04, 3.3299493754487303e+04, 3.3307904252940047e+04, 3.3316314973886954e+04, 3.3324725917278520e+04, 3.3333137083065296e+04, 3.3341548471197806e+04, 3.3349960081626654e+04, 3.3358371914302414e+04, 3.3366783969175711e+04, 3.3375196246197171e+04, 3.3383608745317470e+04, 3.3392021466487284e+04, 3.3400434409657326e+04, 3.3408847574778316e+04, 3.3417260961801010e+04, 3.3425674570676172e+04, 3.3434088401354587e+04, 3.3442502453787092e+04, 3.3450916727924494e+04, 3.3459331223717672e+04, 3.3467745941117500e+04, 3.3476160880074865e+04, 3.3484576040540720e+04, 3.3492991422465988e+04, 3.3501407025801644e+04, 3.3509822850498669e+04, 3.3518238896508075e+04, 3.3526655163780903e+04, 3.3535071652268198e+04, 3.3543488361921030e+04, 3.3551905292690506e+04, 3.3560322444527745e+04, 3.3568739817383881e+04, 3.3577157411210072e+04, 3.3585575225957509e+04, 3.3593993261577387e+04, 3.3602411518020941e+04, 3.3610829995239423e+04, 3.3619248693184090e+04, 3.3627667611806231e+04, 3.3636086751057177e+04, 3.3644506110888244e+04, 3.3652925691250792e+04, 3.3661345492096203e+04, 3.3669765513375867e+04, 3.3678185755041210e+04, 3.3686606217043663e+04, 3.3695026899334705e+04, 3.3703447801865797e+04, 3.3711868924588460e+04, 3.3720290267454220e+04, 3.3728711830414621e+04, 3.3737133613421232e+04, 3.3745555616425641e+04, 3.3753977839379462e+04, 3.3762400282234332e+04, 3.3770822944941894e+04, 3.3779245827453844e+04, 3.3787668929721855e+04, 3.3796092251697657e+04, 3.3804515793333005e+04, 3.3812939554579629e+04, 3.3821363535389326e+04, 3.3829787735713893e+04, 3.3838212155505156e+04, 3.3846636794714963e+04, 3.3855061653295175e+04, 3.3863486731197678e+04, 3.3871912028374398e+04, 3.3880337544777241e+04, 3.3888763280358173e+04, 3.3897189235069149e+04, 3.3905615408862184e+04, 3.3914041801689273e+04, 3.3922468413502458e+04, 3.3930895244253799e+04, 3.3939322293895370e+04, 3.3947749562379257e+04, 3.3956177049657592e+04, 3.3964604755682507e+04, 3.3973032680406162e+04, 3.3981460823780748e+04, 3.3989889185758460e+04, 3.3998317766291519e+04, 3.4006746565332171e+04, 3.4015175582832679e+04, 3.4023604818745334e+04, 3.4032034273022444e+04, 3.4040463945616335e+04, 3.4048893836479343e+04, 3.4057323945563854e+04, 3.4065754272822254e+04, 3.4074184818206944e+04, 3.4082615581670361e+04, 3.4091046563164964e+04, 3.4099477762643204e+04, 3.4107909180057599e+04, 3.4116340815360651e+04, 3.4124772668504898e+04, 3.4133204739442888e+04, 3.4141637028127217e+04, 3.4150069534510469e+04, 3.4158502258545253e+04, 3.4166935200184220e+04, 3.4175368359380031e+04, 3.4183801736085348e+04, 3.4192235330252886e+04, 3.4200669141835366e+04, 3.4209103170785522e+04, 3.4217537417056104e+04, 3.4225971880599929e+04, 3.4234406561369775e+04, 3.4242841459318457e+04, 3.4251276574398842e+04, 3.4259711906563774e+04, 3.4268147455766150e+04, 3.4276583221958877e+04, 3.4285019205094868e+04, 3.4293455405127075e+04, 3.4301891822008460e+04, 3.4310328455692012e+04, 3.4318765306130750e+04, 3.4327202373277687e+04, 3.4335639657085878e+04, 3.4344077157508385e+04, 3.4352514874498302e+04, 3.4360952808008726e+04, 3.4369390957992808e+04, 3.4377829324403669e+04, 3.4386267907194517e+04, 3.4394706706318502e+04, 3.4403145721728863e+04, 3.4411584953378813e+04, 3.4420024401221614e+04, 3.4428464065210515e+04, 3.4436903945298829e+04, 3.4445344041439857e+04, 3.4453784353586940e+04, 3.4462224881693415e+04, 3.4470665625712674e+04, 3.4479106585598085e+04, 3.4487547761303074e+04, 3.4495989152781083e+04, 3.4504430759985538e+04, 3.4512872582869939e+04, 3.4521314621387755e+04, 3.4529756875492501e+04, 3.4538199345137728e+04, 3.4546642030276969e+04, 3.4555084930863806e+04, 3.4563528046851832e+04, 3.4571971378194648e+04, 3.4580414924845900e+04, 3.4588858686759231e+04, 3.4597302663888317e+04, 3.4605746856186845e+04, 3.4614191263608533e+04, 3.4622635886107113e+04, 3.4631080723636340e+04, 3.4639525776149982e+04, 3.4647971043601829e+04, 3.4656416525945686e+04, 3.4664862223135395e+04, 3.4673308135124818e+04, 3.4681754261867798e+04, 3.4690200603318248e+04, 3.4698647159430075e+04, 3.4707093930157193e+04, 3.4715540915453566e+04, 3.4723988115273161e+04, 3.4732435529569964e+04, 3.4740883158297998e+04, 3.4749331001411279e+04, 3.4757779058863860e+04, 3.4766227330609814e+04, 3.4774675816603223e+04, 3.4783124516798191e+04, 3.4791573431148856e+04, 3.4800022559609351e+04, 3.4808471902133861e+04, 3.4816921458676567e+04, 3.4825371229191667e+04, 3.4833821213633390e+04, 3.4842271411955982e+04, 3.4850721824113709e+04, 3.4859172450060847e+04, 3.4867623289751711e+04, 3.4876074343140623e+04, 3.4884525610181918e+04, 3.4892977090829969e+04, 3.4901428785039148e+04, 3.4909880692763873e+04, 3.4918332813958550e+04, 3.4926785148577612e+04, 3.4935237696575539e+04, 3.4943690457906800e+04, 3.4952143432525889e+04, 3.4960596620387325e+04, 3.4969050021445662e+04, 3.4977503635655434e+04, 3.4985957462971230e+04, 3.4994411503347641e+04, 3.5002865756739287e+04, 3.5011320223100796e+04, 3.5019774902386816e+04, 3.5028229794552040e+04, 3.5036684899551139e+04, 3.5045140217338841e+04, 3.5053595747869869e+04, 3.5062051491098966e+04, 3.5070507446980919e+04, 3.5078963615470500e+04, 3.5087419996522520e+04, 3.5095876590091808e+04, 3.5104333396133217e+04, 3.5112790414601594e+04, 3.5121247645451847e+04, 3.5129705088638853e+04, 3.5138162744117559e+04, 3.5146620611842889e+04, 3.5155078691769806e+04, 3.5163536983853308e+04, 3.5171995488048371e+04, 3.5180454204310030e+04, 3.5188913132593319e+04, 3.5197372272853281e+04, 3.5205831625045008e+04, 3.5214291189123585e+04, 3.5222750965044135e+04, 3.5231210952761779e+04, 3.5239671152231684e+04, 3.5248131563408999e+04, 3.5256592186248927e+04, 3.5265053020706677e+04, 3.5273514066737473e+04, 3.5281975324296553e+04, 3.5290436793339199e+04, 3.5298898473820678e+04, 3.5307360365696317e+04, 3.5315822468921418e+04, 3.5324284783451330e+04, 3.5332747309241400e+04, 3.5341210046247026e+04, 3.5349672994423592e+04, 3.5358136153726518e+04, 3.5366599524111232e+04, 3.5375063105533198e+04, 3.5383526897947886e+04, 3.5391990901310783e+04, 3.5400455115577410e+04, 3.5408919540703289e+04, 3.5417384176643965e+04, 3.5425849023355004e+04, 3.5434314080791999e+04, 3.5442779348910553e+04, 3.5451244827666276e+04, 3.5459710517014828e+04, 3.5468176416911861e+04, 3.5476642527313059e+04, 3.5485108848174103e+04, 3.5493575379450711e+04, 3.5502042121098631e+04, 3.5510509073073619e+04, 3.5518976235331422e+04, 3.5527443607827860e+04, 3.5535911190518724e+04, 3.5544378983359849e+04, 3.5552846986307071e+04, 3.5561315199316268e+04, 3.5569783622343311e+04, 3.5578252255344109e+04, 3.5586721098274589e+04, 3.5595190151090675e+04, 3.5603659413748333e+04, 3.5612128886203536e+04, 3.5620598568412286e+04, 3.5629068460330578e+04, 3.5637538561914465e+04, 3.5646008873119979e+04, 3.5654479393903195e+04, 3.5662950124220202e+04, 3.5671421064027105e+04, 3.5679892213280014e+04, 3.5688363571935079e+04, 3.5696835139948467e+04, 3.5705306917276357e+04, 3.5713778903874932e+04, 3.5722251099700421e+04, 3.5730723504709043e+04, 3.5739196118857064e+04, 3.5747668942100740e+04, 3.5756141974396371e+04, 3.5764615215700258e+04, 3.5773088665968724e+04, 3.5781562325158120e+04, 3.5790036193224791e+04, 3.5798510270125131e+04, 3.5806984555815528e+04, 3.5815459050252415e+04, 3.5823933753392215e+04, 3.5832408665191382e+04, 3.5840883785606376e+04, 3.5849359114593703e+04, 3.5857834652109843e+04, 3.5866310398111345e+04, 3.5874786352554744e+04, 3.5883262515396600e+04, 3.5891738886593506e+04, 3.5900215466102032e+04, 3.5908692253878813e+04, 3.5917169249880477e+04, 3.5925646454063681e+04, 3.5934123866385075e+04, 3.5942601486801374e+04, 3.5951079315269271e+04, 3.5959557351745490e+04, 3.5968035596186775e+04, 3.5976514048549870e+04, 3.5984992708791571e+04, 3.5993471576868666e+04, 3.6001950652737971e+04, 3.6010429936356311e+04, 3.6018909427680555e+04, 3.6027389126667542e+04, 3.6035869033274175e+04, 3.6044349147457353e+04, 3.6052829469173987e+04, 3.6061309998381032e+04, 3.6069790735035436e+04, 3.6078271679094185e+04, 3.6086752830514248e+04, 3.6095234189252653e+04, 3.6103715755266421e+04, 3.6112197528512610e+04, 3.6120679508948269e+04, 3.6129161696530493e+04, 3.6137644091216367e+04, 3.6146126692963015e+04, 3.6154609501727566e+04, 3.6163092517467179e+04, 3.6171575740139022e+04, 3.6180059169700289e+04, 3.6188542806108177e+04, 3.6197026649319916e+04, 3.6205510699292739e+04, 3.6213994955983915e+04, 3.6222479419350711e+04, 3.6230964089350426e+04, 3.6239448965940363e+04, 3.6247934049077863e+04, 3.6256419338720261e+04, 3.6264904834824942e+04, 3.6273390537349267e+04, 3.6281876446250644e+04, 3.6290362561486487e+04, 3.6298848883014238e+04, 3.6307335410791347e+04, 3.6315822144775273e+04, 3.6324309084923523e+04, 3.6332796231193584e+04, 3.6341283583542987e+04, 3.6349771141929268e+04, 3.6358258906310002e+04, 3.6366746876642741e+04, 3.6375235052885087e+04, 3.6383723434994652e+04, 3.6392212022929059e+04, 3.6400700816645949e+04, 3.6409189816103004e+04, 3.6417679021257878e+04, 3.6426168432068283e+04, 3.6434658048491925e+04, 3.6443147870486551e+04, 3.6451637898009889e+04, 3.6460128131019723e+04, 3.6468618569473831e+04, 3.6477109213330019e+04, 3.6485600062546102e+04, 3.6494091117079901e+04, 3.6502582376889295e+04, 3.6511073841932135e+04, 3.6519565512166322e+04, 3.6528057387549758e+04, 3.6536549468040357e+04, 3.6545041753596066e+04, 3.6553534244174843e+04, 3.6562026939734656e+04, 3.6570519840233501e+04, 3.6579012945629387e+04, 3.6587506255880347e+04, 3.6595999770944407e+04, 3.6604493490779634e+04, 3.6612987415344105e+04, 3.6621481544595925e+04, 3.6629975878493198e+04, 3.6638470416994060e+04, 3.6646965160056636e+04, 3.6655460107639105e+04, 3.6663955259699651e+04, 3.6672450616196445e+04, 3.6680946177087746e+04, 3.6689441942331745e+04, 3.6697937911886707e+04, 3.6706434085710898e+04, 3.6714930463762605e+04, 3.6723427046000113e+04, 3.6731923832381755e+04, 3.6740420822865854e+04, 3.6748918017410768e+04, 3.6757415415974858e+04, 3.6765913018516512e+04, 3.6774410824994120e+04, 3.6782908835366128e+04, 3.6791407049590940e+04, 3.6799905467627024e+04, 3.6808404089432850e+04, 3.6816902914966915e+04, 3.6825401944187695e+04, 3.6833901177053740e+04, 3.6842400613523569e+04, 3.6850900253555737e+04, 3.6859400097108817e+04, 3.6867900144141400e+04, 3.6876400394612094e+04, 3.6884900848479505e+04, 3.6893401505702284e+04, 3.6901902366239075e+04, 3.6910403430048565e+04, 3.6918904697089434e+04, 3.6927406167320391e+04, 3.6935907840700143e+04, 3.6944409717187445e+04, 3.6952911796741042e+04, 3.6961414079319722e+04, 3.6969916564882260e+04, 3.6978419253387474e+04, 3.6986922144794182e+04, 3.6995425239061216e+04, 3.7003928536147447e+04, 3.7012432036011727e+04, 3.7020935738612956e+04, 3.7029439643910046e+04, 3.7037943751861909e+04, 3.7046448062427495e+04, 3.7054952575565760e+04, 3.7063457291235667e+04, 3.7071962209396210e+04, 3.7080467330006388e+04, 3.7088972653025230e+04, 3.7097478178411780e+04, 3.7105983906125082e+04, 3.7114489836124216e+04, 3.7122995968368275e+04, 3.7131502302816356e+04, 3.7140008839427588e+04, 3.7148515578161096e+04, 3.7157022518976060e+04, 3.7165529661831621e+04, 3.7174037006686973e+04, 3.7182544553501335e+04, 3.7191052302233926e+04, 3.7199560252843978e+04, 3.7208068405290738e+04, 3.7216576759533484e+04, 3.7225085315531505e+04, 3.7233594073244101e+04, 3.7242103032630592e+04, 3.7250612193650311e+04, 3.7259121556262609e+04, 3.7267631120426871e+04, 3.7276140886102461e+04, 3.7284650853248786e+04, 3.7293161021825261e+04, 3.7301671391791329e+04, 3.7310181963106428e+04, 3.7318692735730045e+04, 3.7327203709621637e+04, 3.7335714884740730e+04, 3.7344226261046824e+04, 3.7352737838499459e+04, 3.7361249617058173e+04, 3.7369761596682532e+04, 3.7378273777332121e+04, 3.7386786158966534e+04, 3.7395298741545397e+04, 3.7403811525028330e+04, 3.7412324509374972e+04, 3.7420837694544985e+04, 3.7429351080498054e+04, 3.7437864667193877e+04, 3.7446378454592159e+04, 3.7454892442652628e+04, 3.7463406631335019e+04, 3.7471921020599097e+04, 3.7480435610404646e+04, 3.7488950400711445e+04, 3.7497465391479294e+04, 3.7505980582668039e+04, 3.7514495974237507e+04, 3.7523011566147550e+04, 3.7531527358358057e+04, 3.7540043350828899e+04, 3.7548559543519979e+04, 3.7557075936391229e+04, 3.7565592529402573e+04, 3.7574109322513970e+04, 3.7582626315685382e+04, 3.7591143508876805e+04, 3.7599660902048228e+04, 3.7608178495159664e+04, 3.7616696288171144e+04, 3.7625214281042732e+04, 3.7633732473734475e+04, 3.7642250866206465e+04, 3.7650769458418792e+04, 3.7659288250331570e+04, 3.7667807241904928e+04, 3.7676326433099006e+04, 3.7684845823873962e+04, 3.7693365414189975e+04, 3.7701885204007245e+04, 3.7710405193285958e+04, 3.7718925381986359e+04, 3.7727445770068669e+04, 3.7735966357493147e+04, 3.7744487144220075e+04, 3.7753008130209724e+04, 3.7761529315422420e+04, 3.7770050699818457e+04, 3.7778572283358168e+04, 3.7787094066001911e+04, 3.7795616047710064e+04, 3.7804138228442993e+04, 3.7812660608161088e+04, 3.7821183186824783e+04, 3.7829705964394489e+04, 3.7838228940830668e+04, 3.7846752116093754e+04, 3.7855275490144246e+04, 3.7863799062942628e+04, 3.7872322834449413e+04, 3.7880846804625107e+04, 3.7889370973430261e+04, 3.7897895340825416e+04, 3.7906419906771160e+04, 3.7914944671228084e+04, 3.7923469634156761e+04, 3.7931994795517821e+04, 3.7940520155271908e+04, 3.7949045713379652e+04, 3.7957571469801725e+04, 3.7966097424498803e+04, 3.7974623577431594e+04, 3.7983149928560793e+04, 3.7991676477847141e+04, 3.8000203225251367e+04, 3.8008730170734219e+04, 3.8017257314256487e+04, 3.8025784655778960e+04, 3.8034312195262428e+04, 3.8042839932667717e+04, 3.8051367867955669e+04, 3.8059896001087131e+04, 3.8068424332022951e+04, 3.8076952860724028e+04, 3.8085481587151262e+04, 3.8094010511265558e+04, 3.8102539633027845e+04, 3.8111068952399059e+04, 3.8119598469340162e+04, 3.8128128183812136e+04, 3.8136658095775958e+04, 3.8145188205192644e+04, 3.8153718512023210e+04, 3.8162249016228685e+04, 3.8170779717770121e+04, 3.8179310616608600e+04, 3.8187841712705187e+04, 3.8196373006020978e+04, 3.8204904496517091e+04, 3.8213436184154663e+04, 3.8221968068894821e+04, 3.8230500150698732e+04, 3.8239032429527564e+04, 3.8247564905342508e+04, 3.8256097578104775e+04, 3.8264630447775577e+04, 3.8273163514316140e+04, 3.8281696777687743e+04, 3.8290230237851625e+04, 3.8298763894769072e+04, 3.8307297748401383e+04, 3.8315831798709871e+04, 3.8324366045655850e+04, 3.8332900489200671e+04, 3.8341435129305690e+04, 3.8349969965932280e+04, 3.8358504999041827e+04, 3.8367040228595724e+04, 3.8375575654555403e+04, 3.8384111276882286e+04, 3.8392647095537832e+04, 3.8401183110483493e+04, 3.8409719321680735e+04, 3.8418255729091077e+04, 3.8426792332676014e+04, 3.8435329132397062e+04, 3.8443866128215777e+04, 3.8452403320093705e+04, 3.8460940707992399e+04, 3.8469478291873464e+04, 3.8478016071698490e+04, 3.8486554047429097e+04, 3.8495092219026898e+04, 3.8503630586453546e+04, 3.8512169149670692e+04, 3.8520707908640026e+04, 3.8529246863323227e+04, 3.8537786013681995e+04, 3.8546325359678056e+04, 3.8554864901273140e+04, 3.8563404638428983e+04, 3.8571944571107371e+04, 3.8580484699270070e+04, 3.8589025022878872e+04, 3.8597565541895587e+04, 3.8606106256282044e+04, 3.8614647166000075e+04, 3.8623188271011539e+04, 3.8631729571278302e+04, 3.8640271066762245e+04, 3.8648812757425265e+04, 3.8657354643229271e+04, 3.8665896724136204e+04, 3.8674439000107981e+04, 3.8682981471106585e+04, 3.8691524137093977e+04, 3.8700066998032140e+04, 3.8708610053883080e+04, 3.8717153304608823e+04, 3.8725696750171388e+04, 3.8734240390532825e+04, 3.8742784225655188e+04, 3.8751328255500550e+04, 3.8759872480031023e+04, 3.8768416899208685e+04, 3.8776961512995680e+04, 3.8785506321354136e+04, 3.8794051324246189e+04, 3.8802596521634012e+04, 3.8811141913479791e+04, 3.8819687499745705e+04, 3.8828233280393979e+04, 3.8836779255386820e+04, 3.8845325424686474e+04, 3.8853871788255194e+04, 3.8862418346055238e+04, 3.8870965098048895e+04, 3.8879512044198469e+04, 3.8888059184466249e+04, 3.8896606518814573e+04, 3.8905154047205782e+04, 3.8913701769602230e+04, 3.8922249685966293e+04, 3.8930797796260340e+04, 3.8939346100446783e+04, 3.8947894598488027e+04, 3.8956443290346506e+04, 3.8964992175984655e+04, 3.8973541255364929e+04, 3.8982090528449808e+04, 3.8990639995201767e+04, 3.8999189655583323e+04, 3.9007739509556974e+04, 3.9016289557085263e+04, 3.9024839798130735e+04, 3.9033390232655933e+04, 3.9041940860623436e+04, 3.9050491681995831e+04, 3.9059042696735727e+04, 3.9067593904805733e+04, 3.9076145306168481e+04, 3.9084696900786614e+04, 3.9093248688622800e+04, 3.9101800669639691e+04, 3.9110352843800007e+04, 3.9118905211066434e+04, 3.9127457771401685e+04, 3.9136010524768506e+04, 3.9144563471129630e+04, 3.9153116610447818e+04, 3.9161669942685854e+04, 3.9170223467806522e+04, 3.9178777185772618e+04, 3.9187331096546972e+04, 3.9195885200092402e+04, 3.9204439496371779e+04, 3.9212993985347937e+04, 3.9221548666983763e+04, 3.9230103541242148e+04, 3.9238658608085992e+04, 3.9247213867478218e+04, 3.9255769319381754e+04, 3.9264324963759544e+04, 3.9272880800574552e+04, 3.9281436829789760e+04, 3.9289993051368147e+04, 3.9298549465272707e+04, 3.9307106071466485e+04, 3.9315662869912499e+04, 3.9324219860573787e+04, 3.9332777043413415e+04, 3.9341334418394472e+04, 3.9349891985480026e+04, 3.9358449744633188e+04, 3.9367007695817076e+04, 3.9375565838994829e+04, 3.9384124174129574e+04, 3.9392682701184487e+04, 3.9401241420122737e+04, 3.9409800330907499e+04, 3.9418359433501988e+04, 3.9426918727869423e+04, 3.9435478213973031e+04, 3.9444037891776054e+04, 3.9452597761241756e+04, 3.9461157822333393e+04, 3.9469718075014272e+04, 3.9478278519247680e+04, 3.9486839154996938e+04, 3.9495399982225375e+04, 3.9503961000896335e+04, 3.9512522210973162e+04, 3.9521083612419243e+04, 3.9529645205197958e+04, 3.9538206989272701e+04, 3.9546768964606890e+04, 3.9555331131163948e+04, 3.9563893488907313e+04, 3.9572456037800446e+04, 3.9581018777806814e+04, 3.9589581708889906e+04, 3.9598144831013211e+04, 3.9606708144140233e+04, 3.9615271648234513e+04, 3.9623835343259583e+04, 3.9632399229178998e+04, 3.9640963305956313e+04, 3.9649527573555119e+04, 3.9658092031938999e+04, 3.9666656681071567e+04, 3.9675221520916450e+04, 3.9683786551437275e+04, 3.9692351772597707e+04, 3.9700917184361395e+04, 3.9709482786692017e+04, 3.9718048579553266e+04, 3.9726614562908857e+04, 3.9735180736722490e+04, 3.9743747100957924e+04, 3.9752313655578873e+04, 3.9760880400549118e+04, 3.9769447335832432e+04, 3.9778014461392602e+04, 3.9786581777193423e+04, 3.9795149283198705e+04, 3.9803716979372301e+04, 3.9812284865678033e+04, 3.9820852942079764e+04, 3.9829421208541360e+04, 3.9837989665026718e+04, 3.9846558311499721e+04, 3.9855127147924286e+04, 3.9863696174264340e+04, 3.9872265390483830e+04, 3.9880834796546696e+04, 3.9889404392416902e+04, 3.9897974178058437e+04, 3.9906544153435294e+04, 3.9915114318511478e+04, 3.9923684673251017e+04, 3.9932255217617923e+04, 3.9940825951576269e+04, 3.9949396875090104e+04, 3.9957967988123506e+04, 3.9966539290640569e+04, 3.9975110782605392e+04, 3.9983682463982092e+04, 3.9992254334734796e+04, 4.0000826394827658e+04, 4.0009398644224821e+04, 4.0017971082890464e+04, 4.0026543710788770e+04, 4.0035116527883933e+04, 4.0043689534140176e+04, 4.0052262729521703e+04, 4.0060836113992773e+04, 4.0069409687517626e+04, 4.0077983450060536e+04, 4.0086557401585764e+04, 4.0095131542057628e+04, 4.0103705871440412e+04, 4.0112280389698448e+04, 4.0120855096796062e+04, 4.0129429992697609e+04, 4.0138005077367437e+04, 4.0146580350769938e+04, 4.0155155812869474e+04, 4.0163731463630465e+04, 4.0172307303017311e+04, 4.0180883330994453e+04, 4.0189459547526320e+04, 4.0198035952577367e+04, 4.0206612546112061e+04, 4.0215189328094886e+04, 4.0223766298490344e+04, 4.0232343457262927e+04, 4.0240920804377158e+04, 4.0249498339797581e+04, 4.0258076063488741e+04, 4.0266653975415196e+04, 4.0275232075541513e+04, 4.0283810363832294e+04, 4.0292388840252119e+04, 4.0300967504765627e+04, 4.0309546357337422e+04, 4.0318125397932163e+04, 4.0326704626514504e+04, 4.0335284043049098e+04, 4.0343863647500628e+04, 4.0352443439833805e+04, 4.0361023420013313e+04, 4.0369603588003883e+04, 4.0378183943770258e+04, 4.0386764487277171e+04, 4.0395345218489390e+04, 4.0403926137371680e+04, 4.0412507243888846e+04, 4.0421088538005664e+04, 4.0429670019686964e+04, 4.0438251688897566e+04, 4.0446833545602305e+04, 4.0455415589766038e+04, 4.0463997821353638e+04, 4.0472580240329975e+04, 4.0481162846659936e+04, 4.0489745640308436e+04, 4.0498328621240384e+04, 4.0506911789420723e+04, 4.0515495144814398e+04, 4.0524078687386354e+04, 4.0532662417101572e+04, 4.0541246333925024e+04, 4.0549830437821722e+04, 4.0558414728756681e+04, 4.0566999206694898e+04, 4.0575583871601426e+04, 4.0584168723441318e+04, 4.0592753762179629e+04, 4.0601338987781441e+04, 4.0609924400211836e+04, 4.0618509999435912e+04, 4.0627095785418795e+04, 4.0635681758125596e+04, 4.0644267917521480e+04, 4.0652854263571593e+04, 4.0661440796241084e+04, 4.0670027515495152e+04, 4.0678614421298975e+04, 4.0687201513617780e+04, 4.0695788792416759e+04, 4.0704376257661163e+04, 4.0712963909316226e+04, 4.0721551747347214e+04, 4.0730139771719383e+04, 4.0738727982398035e+04, 4.0747316379348456e+04, 4.0755904962535955e+04, 4.0764493731925868e+04, 4.0773082687483504e+04, 4.0781671829174236e+04, 4.0790261156963403e+04, 4.0798850670816406e+04, 4.0807440370698612e+04, 4.0816030256575417e+04, 4.0824620328412246e+04, 4.0833210586174522e+04, 4.0841801029827671e+04, 4.0850391659337176e+04, 4.0858982474668461e+04, 4.0867573475787023e+04, 4.0876164662658346e+04, 4.0884756035247941e+04, 4.0893347593521314e+04, 4.0901939337443990e+04, 4.0910531266981525e+04, 4.0919123382099453e+04, 4.0927715682763352e+04, 4.0936308168938805e+04, 4.0944900840591399e+04, 4.0953493697686739e+04, 4.0962086740190432e+04, 4.0970679968068129e+04, 4.0979273381285457e+04, 4.0987866979808074e+04, 4.0996460763601652e+04, 4.1005054732631870e+04, 4.1013648886864430e+04, 4.1022243226265025e+04, 4.1030837750799379e+04, 4.1039432460433221e+04, 4.1048027355132304e+04, 4.1056622434862373e+04, 4.1065217699589208e+04, 4.1073813149278583e+04, 4.1082408783896317e+04, 4.1091004603408182e+04, 4.1099600607780019e+04, 4.1108196796977660e+04, 4.1116793170966950e+04, 4.1125389729713745e+04, 4.1133986473183919e+04, 4.1142583401343356e+04, 4.1151180514157946e+04, 4.1159777811593609e+04, 4.1168375293616256e+04, 4.1176972960191815e+04, 4.1185570811286256e+04, 4.1194168846865519e+04, 4.1202767066895583e+04, 4.1211365471342418e+04, 4.1219964060172038e+04, 4.1228562833350450e+04, 4.1237161790843667e+04, 4.1245760932617733e+04, 4.1254360258638691e+04, 4.1262959768872599e+04, 4.1271559463285521e+04, 4.1280159341843559e+04, 4.1288759404512806e+04, 4.1297359651259350e+04, 4.1305960082049336e+04, 4.1314560696848901e+04, 4.1323161495624168e+04, 4.1331762478341319e+04, 4.1340363644966506e+04, 4.1348964995465925e+04, 4.1357566529805772e+04, 4.1366168247952257e+04, 4.1374770149871598e+04, 4.1383372235530034e+04, 4.1391974504893813e+04, 4.1400576957929181e+04, 4.1409179594602414e+04, 4.1417782414879795e+04, 4.1426385418727630e+04, 4.1434988606112209e+04, 4.1443591976999865e+04, 4.1452195531356920e+04, 4.1460799269149742e+04, 4.1469403190344667e+04, 4.1478007294908071e+04, 4.1486611582806341e+04, 4.1495216054005861e+04, 4.1503820708473046e+04, 4.1512425546174316e+04, 4.1521030567076101e+04, 4.1529635771144836e+04, 4.1538241158346988e+04, 4.1546846728649019e+04, 4.1555452482017419e+04, 4.1564058418418659e+04, 4.1572664537819270e+04, 4.1581270840185753e+04, 4.1589877325484646e+04, 4.1598483993682486e+04, 4.1607090844745835e+04, 4.1615697878641244e+04, 4.1624305095335309e+04, 4.1632912494794611e+04, 4.1641520076985762e+04, 4.1650127841875357e+04, 4.1658735789430051e+04, 4.1667343919616447e+04, 4.1675952232401229e+04, 4.1684560727751050e+04, 4.1693169405632594e+04, 4.1701778266012530e+04, 4.1710387308857578e+04, 4.1718996534134443e+04, 4.1727605941809838e+04, 4.1736215531850517e+04, 4.1744825304223225e+04, 4.1753435258894722e+04, 4.1762045395831781e+04, 4.1770655715001194e+04, 4.1779266216369746e+04, 4.1787876899904259e+04, 4.1796487765571532e+04, 4.1805098813338416e+04, 4.1813710043171763e+04, 4.1822321455038415e+04, 4.1830933048905252e+04, 4.1839544824739154e+04, 4.1848156782507016e+04, 4.1856768922175746e+04, 4.1865381243712254e+04, 4.1873993747083478e+04, 4.1882606432256347e+04, 4.1891219299197837e+04, 4.1899832347874893e+04, 4.1908445578254505e+04, 4.1917058990303660e+04, 4.1925672583989362e+04, 4.1934286359278623e+04, 4.1942900316138461e+04, 4.1951514454535929e+04, 4.1960128774438082e+04, 4.1968743275811961e+04, 4.1977357958624656e+04, 4.1985972822843243e+04, 4.1994587868434828e+04, 4.2003203095366516e+04, 4.2011818503605427e+04, 4.2020434093118703e+04, 4.2029049863873479e+04, 4.2037665815836925e+04, 4.2046281948976197e+04, 4.2054898263258481e+04, 4.2063514758650977e+04, 4.2072131435120871e+04, 4.2080748292635406e+04, 4.2089365331161782e+04, 4.2097982550667271e+04, 4.2106599951119104e+04, 4.2115217532484552e+04, 4.2123835294730896e+04, 4.2132453237825408e+04, 4.2141071361735405e+04, 4.2149689666428189e+04, 4.2158308151871090e+04, 4.2166926818031439e+04, 4.2175545664876583e+04, 4.2184164692373881e+04, 4.2192783900490700e+04, 4.2201403289194430e+04, 4.2210022858452467e+04, 4.2218642608232207e+04, 4.2227262538501076e+04, 4.2235882649226493e+04, 4.2244502940375911e+04, 4.2253123411916778e+04, 4.2261744063816557e+04, 4.2270364896042731e+04, 4.2278985908562783e+04, 4.2287607101344220e+04, 4.2296228474354546e+04, 4.2304850027561290e+04, 4.2313471760931978e+04, 4.2322093674434160e+04, 4.2330715768035399e+04, 4.2339338041703260e+04, 4.2347960495405328e+04, 4.2356583129109211e+04, 4.2365205942782486e+04, 4.2373828936392783e+04, 4.2382452109907739e+04, 4.2391075463294990e+04, 4.2399698996522173e+04, 4.2408322709556967e+04, 4.2416946602367039e+04, 4.2425570674920091e+04, 4.2434194927183795e+04, 4.2442819359125882e+04, 4.2451443970714070e+04, 4.2460068761916082e+04, 4.2468693732699670e+04, 4.2477318883032589e+04, 4.2485944212882605e+04, 4.2494569722217508e+04, 4.2503195411005079e+04, 4.2511821279213123e+04, 4.2520447326809444e+04, 4.2529073553761889e+04, 4.2537699960038277e+04, 4.2546326545606462e+04, 4.2554953310434306e+04, 4.2563580254489687e+04, 4.2572207377740473e+04, 4.2580834680154570e+04, 4.2589462161699885e+04, 4.2598089822344322e+04, 4.2606717662055831e+04, 4.2615345680802333e+04, 4.2623973878551791e+04, 4.2632602255272170e+04, 4.2641230810931440e+04, 4.2649859545497588e+04, 4.2658488458938613e+04, 4.2667117551222531e+04, 4.2675746822317349e+04, 4.2684376272191119e+04, 4.2693005900811862e+04, 4.2701635708147653e+04, 4.2710265694166541e+04, 4.2718895858836615e+04, 4.2727526202125955e+04, 4.2736156724002678e+04, 4.2744787424434893e+04, 4.2753418303390717e+04, 4.2762049360838282e+04, 4.2770680596745740e+04, 4.2779312011081252e+04, 4.2787943603812972e+04, 4.2796575374909095e+04, 4.2805207324337804e+04, 4.2813839452067310e+04, 4.2822471758065833e+04, 4.2831104242301590e+04, 4.2839736904742815e+04, 4.2848369745357762e+04, 4.2857002764114681e+04, 4.2865635960981854e+04, 4.2874269335927558e+04, 4.2882902888920100e+04, 4.2891536619927763e+04, 4.2900170528918861e+04, 4.2908804615861751e+04, 4.2917438880724752e+04, 4.2926073323476216e+04, 4.2934707944084512e+04, 4.2943342742517998e+04, 4.2951977718745067e+04, 4.2960612872734120e+04, 4.2969248204453550e+04, 4.2977883713871794e+04, 4.2986519400957259e+04, 4.2995155265678397e+04, 4.3003791308003652e+04, 4.3012427527901484e+04, 4.3021063925340379e+04, 4.3029700500288818e+04, 4.3038337252715297e+04, 4.3046974182588310e+04, 4.3055611289876389e+04, 4.3064248574548066e+04, 4.3072886036571872e+04, 4.3081523675916367e+04, 4.3090161492550113e+04, 4.3098799486441676e+04, 4.3107437657559647e+04, 4.3116076005872615e+04, 4.3124714531349207e+04, 4.3133353233958012e+04, 4.3141992113667693e+04, 4.3150631170446861e+04, 4.3159270404264185e+04, 4.3167909815088322e+04, 4.3176549402887955e+04, 4.3185189167631761e+04, 4.3193829109288439e+04, 4.3202469227826688e+04, 4.3211109523215244e+04, 4.3219749995422819e+04, 4.3228390644418163e+04, 4.3237031470170034e+04, 4.3245672472647173e+04, 4.3254313651818375e+04, 4.3262955007652417e+04, 4.3271596540118080e+04, 4.3280238249184200e+04, 4.3288880134819570e+04, 4.3297522196993028e+04, 4.3306164435673425e+04, 4.3314806850829591e+04, 4.3323449442430399e+04, 4.3332092210444731e+04, 4.3340735154841452e+04, 4.3349378275589457e+04, 4.3358021572657679e+04, 4.3366665046014998e+04, 4.3375308695630374e+04, 4.3383952521472733e+04, 4.3392596523511012e+04, 4.3401240701714181e+04, 4.3409885056051215e+04, 4.3418529586491088e+04, 4.3427174293002805e+04, 4.3435819175555363e+04, 4.3444464234117775e+04, 4.3453109468659080e+04, 4.3461754879148291e+04, 4.3470400465554478e+04, 4.3479046227846688e+04, 4.3487692165993991e+04, 4.3496338279965472e+04, 4.3504984569730223e+04, 4.3513631035257342e+04, 4.3522277676515943e+04, 4.3530924493475155e+04, 4.3539571486104112e+04, 4.3548218654371951e+04, 4.3556865998247828e+04, 4.3565513517700922e+04, 4.3574161212700405e+04, 4.3582809083215463e+04, 4.3591457129215290e+04, 4.3600105350669110e+04, 4.3608753747546143e+04, 4.3617402319815614e+04, 4.3626051067446773e+04, 4.3634699990408866e+04, 4.3643349088671166e+04, 4.3651998362202939e+04, 4.3660647810973474e+04, 4.3669297434952074e+04, 4.3677947234108040e+04, 4.3686597208410683e+04, 4.3695247357829357e+04, 4.3703897682333380e+04, 4.3712548181892103e+04, 4.3721198856474890e+04, 4.3729849706051122e+04, 4.3738500730590167e+04, 4.3747151930061431e+04, 4.3755803304434317e+04, 4.3764454853678239e+04, 4.3773106577762606e+04, 4.3781758476656876e+04, 4.3790410550330489e+04, 4.3799062798752901e+04, 4.3807715221893573e+04, 4.3816367819722000e+04, 4.3825020592207657e+04, 4.3833673539320051e+04, 4.3842326661028688e+04, 4.3850979957303098e+04, 4.3859633428112807e+04, 4.3868287073427360e+04, 4.3876940893216306e+04, 4.3885594887449217e+04, 4.3894249056095658e+04, 4.3902903399125229e+04, 4.3911557916507496e+04, 4.3920212608212096e+04, 4.3928867474208637e+04, 4.3937522514466749e+04, 4.3946177728956063e+04, 4.3954833117646231e+04, 4.3963488680506911e+04, 4.3972144417507778e+04, 4.3980800328618512e+04, 4.3989456413808788e+04, 4.3998112673048330e+04, 4.4006769106306834e+04, 4.4015425713554039e+04, 4.4024082494759663e+04, 4.4032739449893452e+04, 4.4041396578925167e+04, 4.4050053881824570e+04, 4.4058711358561435e+04, 4.4067369009105547e+04, 4.4076026833426709e+04, 4.4084684831494720e+04, 4.4093343003279399e+04, 4.4102001348750571e+04, 4.4110659867878079e+04, 4.4119318560631771e+04, 4.4127977426981503e+04, 4.4136636466897144e+04, 4.4145295680348579e+04, 4.4153955067305702e+04, 4.4162614627738403e+04, 4.4171274361616604e+04, 4.4179934268910220e+04, 4.4188594349589177e+04, 4.4197254603623434e+04, 4.4205915030982935e+04, 4.4214575631637650e+04, 4.4223236405557538e+04, 4.4231897352712600e+04, 4.4240558473072822e+04, 4.4249219766608214e+04, 4.4257881233288790e+04, 4.4266542873084567e+04, 4.4275204685965589e+04, 4.4283866671901909e+04, 4.4292528830863572e+04, 4.4301191162820658e+04, 4.4309853667743235e+04, 4.4318516345601398e+04, 4.4327179196365236e+04, 4.4335842220004859e+04, 4.4344505416490392e+04, 4.4353168785791961e+04, 4.4361832327879711e+04, 4.4370496042723797e+04, 4.4379159930294358e+04, 4.4387823990561592e+04, 4.4396488223495660e+04, 4.4405152629066753e+04, 4.4413817207245083e+04, 4.4422481958000855e+04, 4.4431146881304288e+04, 4.4439811977125624e+04, 4.4448477245435104e+04, 4.4457142686202984e+04, 4.4465808299399519e+04, 4.4474474084994989e+04, 4.4483140042959669e+04, 4.4491806173263860e+04, 4.4500472475877868e+04, 4.4509138950772001e+04, 4.4517805597916595e+04, 4.4526472417281962e+04, 4.4535139408838470e+04, 4.4543806572556467e+04, 4.4552473908406311e+04, 4.4561141416358383e+04, 4.4569809096383069e+04, 4.4578476948450778e+04, 4.4587144972531889e+04, 4.4595813168596840e+04, 4.4604481536616055e+04, 4.4613150076559963e+04, 4.4621818788399010e+04, 4.4630487672103671e+04, 4.4639156727644397e+04, 4.4647825954991669e+04, 4.4656495354115978e+04, 4.4665164924987817e+04, 4.4673834667577692e+04, 4.4682504581856127e+04, 4.4691174667793646e+04, 4.4699844925360791e+04, 4.4708515354528099e+04, 4.4717185955266148e+04, 4.4725856727545499e+04, 4.4734527671336720e+04, 4.4743198786610410e+04, 4.4751870073337166e+04, 4.4760541531487586e+04, 4.4769213161032312e+04, 4.4777884961941956e+04, 4.4786556934187152e+04, 4.4795229077738564e+04, 4.4803901392566855e+04, 4.4812573878642674e+04, 4.4821246535936712e+04, 4.4829919364419657e+04, 4.4838592364062213e+04, 4.4847265534835082e+04, 4.4855938876708991e+04, 4.4864612389654656e+04, 4.4873286073642834e+04, 4.4881959928644268e+04, 4.4890633954629709e+04, 4.4899308151569938e+04, 4.4907982519435725e+04, 4.4916657058197859e+04, 4.4925331767827149e+04, 4.4934006648294395e+04, 4.4942681699570421e+04, 4.4951356921626058e+04, 4.4960032314432152e+04, 4.4968707877959539e+04, 4.4977383612179088e+04, 4.4986059517061658e+04, 4.4994735592578130e+04, 4.5003411838699401e+04, 4.5012088255396360e+04, 4.5020764842639932e+04, 4.5029441600401013e+04, 4.5038118528650550e+04, 4.5046795627359483e+04, 4.5055472896498744e+04, 4.5064150336039304e+04, 4.5072827945952129e+04, 4.5081505726208197e+04, 4.5090183676778492e+04, 4.5098861797634010e+04, 4.5107540088745773e+04, 4.5116218550084792e+04, 4.5124897181622080e+04, 4.5133575983328687e+04, 4.5142254955175667e+04, 4.5150934097134072e+04, 4.5159613409174955e+04, 4.5168292891269419e+04, 4.5176972543388540e+04, 4.5185652365503403e+04, 4.5194332357585125e+04, 4.5203012519604825e+04, 4.5211692851533611e+04, 4.5220373353342642e+04, 4.5229054025003054e+04, 4.5237734866485996e+04, 4.5246415877762651e+04, 4.5255097058804175e+04, 4.5263778409581755e+04, 4.5272459930066594e+04, 4.5281141620229886e+04, 4.5289823480042854e+04, 4.5298505509476730e+04, 4.5307187708502737e+04, 4.5315870077092113e+04, 4.5324552615216111e+04, 4.5333235322846005e+04, 4.5341918199953063e+04, 4.5350601246508566e+04, 4.5359284462483811e+04, 4.5367967847850094e+04, 4.5376651402578726e+04, 4.5385335126641032e+04, 4.5394019020008338e+04, 4.5402703082651991e+04, 4.5411387314543339e+04, 4.5420071715653743e+04, 4.5428756285954565e+04, 4.5437441025417196e+04, 4.5446125934013013e+04, 4.5454811011713427e+04, 4.5463496258489839e+04, 4.5472181674313673e+04, 4.5480867259156352e+04, 4.5489553012989309e+04, 4.5498238935784000e+04, 4.5506925027511876e+04, 4.5515611288144406e+04, 4.5524297717653069e+04, 4.5532984316009351e+04, 4.5541671083184745e+04, 4.5550358019150743e+04, 4.5559045123878874e+04, 4.5567732397340667e+04, 4.5576419839507646e+04, 4.5585107450351345e+04, 4.5593795229843345e+04, 4.5602483177955182e+04, 4.5611171294658445e+04, 4.5619859579924698e+04, 4.5628548033725550e+04, 4.5637236656032597e+04, 4.5645925446817440e+04, 4.5654614406051711e+04, 4.5663303533707040e+04, 4.5671992829755050e+04, 4.5680682294167404e+04, 4.5689371926915759e+04, 4.5698061727971784e+04, 4.5706751697307147e+04, 4.5715441834893543e+04, 4.5724132140702670e+04, 4.5732822614706231e+04, 4.5741513256875936e+04, 4.5750204067183513e+04, 4.5758895045600701e+04, 4.5767586192099239e+04, 4.5776277506650884e+04, 4.5784968989227404e+04, 4.5793660639800561e+04, 4.5802352458342131e+04, 4.5811044444823921e+04, 4.5819736599217722e+04, 4.5828428921495353e+04, 4.5837121411628628e+04, 4.5845814069589382e+04, 4.5854506895349441e+04, 4.5863199888880663e+04, 4.5871893050154904e+04, 4.5880586379144028e+04, 4.5889279875819913e+04, 4.5897973540154439e+04, 4.5906667372119511e+04, 4.5915361371687039e+04, 4.5924055538828921e+04, 4.5932749873517088e+04, 4.5941444375723477e+04, 4.5950139045420023e+04, 4.5958833882578678e+04, 4.5967528887171407e+04, 4.5976224059170177e+04, 4.5984919398546976e+04, 4.5993614905273796e+04, 4.6002310579322620e+04, 4.6011006420665457e+04, 4.6019702429274344e+04, 4.6028398605121292e+04, 4.6037094948178339e+04, 4.6045791458417523e+04, 4.6054488135810912e+04, 4.6063184980330574e+04, 4.6071881991948561e+04, 4.6080579170636971e+04, 4.6089276516367892e+04, 4.6097974029113429e+04, 4.6106671708845694e+04, 4.6115369555536810e+04, 4.6124067569158891e+04, 4.6132765749684084e+04, 4.6141464097084543e+04, 4.6150162611332424e+04, 4.6158861292399881e+04, 4.6167560140259098e+04, 4.6176259154882275e+04, 4.6184958336241580e+04, 4.6193657684309241e+04, 4.6202357199057456e+04, 4.6211056880458447e+04, 4.6219756728484448e+04, 4.6228456743107694e+04, 4.6237156924300449e+04, 4.6245857272034955e+04, 4.6254557786283498e+04, 4.6263258467018350e+04, 4.6271959314211796e+04, 4.6280660327836136e+04, 4.6289361507863658e+04, 4.6298062854266704e+04, 4.6306764367017582e+04, 4.6315466046088623e+04, 4.6324167891452169e+04, 4.6332869903080580e+04, 4.6341572080946200e+04, 4.6350274425021424e+04, 4.6358976935278610e+04, 4.6367679611690153e+04, 4.6376382454228456e+04, 4.6385085462865914e+04, 4.6393788637574966e+04, 4.6402491978328006e+04, 4.6411195485097494e+04, 4.6419899157855849e+04, 4.6428602996575544e+04, 4.6437307001229019e+04, 4.6446011171788770e+04, 4.6454715508227258e+04, 4.6463420010516988e+04, 4.6472124678630433e+04, 4.6480829512540127e+04, 4.6489534512218568e+04, 4.6498239677638281e+04, 4.6506945008771821e+04, 4.6515650505591708e+04, 4.6524356168070502e+04, 4.6533061996180768e+04, 4.6541767989895074e+04, 4.6550474149186004e+04, 4.6559180474026136e+04, 4.6567886964388083e+04, 4.6576593620244450e+04, 4.6585300441567837e+04, 4.6594007428330886e+04, 4.6602714580506225e+04, 4.6611421898066503e+04, 4.6620129380984356e+04, 4.6628837029232469e+04, 4.6637544842783493e+04, 4.6646252821610127e+04, 4.6654960965685030e+04, 4.6663669274980915e+04, 4.6672377749470499e+04, 4.6681086389126489e+04, 4.6689795193921607e+04, 4.6698504163828591e+04, 4.6707213298820177e+04, 4.6715922598869125e+04, 4.6724632063948193e+04, 4.6733341694030139e+04, 4.6742051489087753e+04, 4.6750761449093829e+04, 4.6759471574021154e+04, 4.6768181863842525e+04, 4.6776892318530779e+04, 4.6785602938058728e+04, 4.6794313722399194e+04, 4.6803024671525025e+04, 4.6811735785409081e+04, 4.6820447064024214e+04, 4.6829158507343287e+04, 4.6837870115339188e+04, 4.6846581887984787e+04, 4.6855293825253000e+04, 4.6864005927116712e+04, 4.6872718193548848e+04, 4.6881430624522327e+04, 4.6890143220010075e+04, 4.6898855979985041e+04, 4.6907568904420157e+04, 4.6916281993288394e+04, 4.6924995246562714e+04, 4.6933708664216087e+04, 4.6942422246221511e+04, 4.6951135992551965e+04, 4.6959849903180453e+04, 4.6968563978079990e+04, 4.6977278217223597e+04, 4.6985992620584300e+04, 4.6994707188135137e+04, 4.7003421919849148e+04, 4.7012136815699399e+04, 4.7020851875658947e+04, 4.7029567099700856e+04, 4.7038282487798228e+04, 4.7046998039924139e+04, 4.7055713756051686e+04, 4.7064429636153975e+04, 4.7073145680204136e+04, 4.7081861888175299e+04, 4.7090578260040566e+04, 4.7099294795773116e+04, 4.7108011495346087e+04, 4.7116728358732638e+04, 4.7125445385905929e+04, 4.7134162576839160e+04, 4.7142879931505493e+04, 4.7151597449878143e+04, 4.7160315131930314e+04, 4.7169032977635208e+04, 4.7177750986966057e+04, 4.7186469159896085e+04, 4.7195187496398539e+04, 4.7203905996446658e+04, 4.7212624660013709e+04, 4.7221343487072954e+04, 4.7230062477597661e+04, 4.7238781631561127e+04, 4.7247500948936635e+04, 4.7256220429697481e+04, 4.7264940073816986e+04, 4.7273659881268461e+04, 4.7282379852025239e+04, 4.7291099986060661e+04, 4.7299820283348054e+04, 4.7308540743860773e+04, 4.7317261367572195e+04, 4.7325982154455691e+04, 4.7334703104484623e+04, 4.7343424217632382e+04, 4.7352145493872369e+04, 4.7360866933177997e+04, 4.7369588535522671e+04, 4.7378310300879813e+04, 4.7387032229222867e+04, 4.7395754320525259e+04, 4.7404476574760440e+04, 4.7413198991901867e+04, 4.7421921571923005e+04, 4.7430644314797340e+04, 4.7439367220498338e+04, 4.7448090288999505e+04, 4.7456813520274336e+04, 4.7465536914296339e+04, 4.7474260471039022e+04, 4.7482984190475938e+04, 4.7491708072580594e+04, 4.7500432117326549e+04, 4.7509156324687348e+04, 4.7517880694636558e+04, 4.7526605227147738e+04, 4.7535329922194476e+04, 4.7544054779750360e+04, 4.7552779799788979e+04, 4.7561504982283936e+04, 4.7570230327208854e+04, 4.7578955834537337e+04, 4.7587681504243024e+04, 4.7596407336299548e+04, 4.7605133330680568e+04, 4.7613859487359725e+04, 4.7622585806310693e+04, 4.7631312287507128e+04, 4.7640038930922725e+04, 4.7648765736531175e+04, 4.7657492704306162e+04, 4.7666219834221403e+04, 4.7674947126250612e+04, 4.7683674580367508e+04, 4.7692402196545831e+04, 4.7701129974759308e+04, 4.7709857914981701e+04, 4.7718586017186768e+04, 4.7727314281348270e+04, 4.7736042707439978e+04, 4.7744771295435676e+04, 4.7753500045309149e+04, 4.7762228957034211e+04, 4.7770958030584661e+04, 4.7779687265934313e+04, 4.7788416663057003e+04, 4.7797146221926567e+04, 4.7805875942516832e+04, 4.7814605824801663e+04, 4.7823335868754904e+04, 4.7832066074350434e+04, 4.7840796441562125e+04, 4.7849526970363862e+04, 4.7858257660729549e+04, 4.7866988512633063e+04, 4.7875719526048335e+04, 4.7884450700949274e+04, 4.7893182037309809e+04, 4.7901913535103871e+04, 4.7910645194305413e+04, 4.7919377014888370e+04, 4.7928108996826719e+04, 4.7936841140094424e+04, 4.7945573444665461e+04, 4.7954305910513809e+04, 4.7963038537613465e+04, 4.7971771325938440e+04, 4.7980504275462736e+04, 4.7989237386160370e+04, 4.7997970658005383e+04, 4.8006704090971791e+04, 4.8015437685033656e+04, 4.8024171440165017e+04, 4.8032905356339943e+04, 4.8041639433532502e+04, 4.8050373671716770e+04, 4.8059108070866838e+04, 4.8067842630956788e+04, 4.8076577351960732e+04, 4.8085312233852776e+04, 4.8094047276607045e+04, 4.8102782480197668e+04, 4.8111517844598769e+04, 4.8120253369784507e+04, 4.8128989055729027e+04, 4.8137724902406488e+04, 4.8146460909791058e+04, 4.8155197077856908e+04, 4.8163933406578246e+04, 4.8172669895929248e+04, 4.8181406545884114e+04, 4.8190143356417073e+04, 4.8198880327502324e+04, 4.8207617459114110e+04, 4.8216354751226652e+04, 4.8225092203814209e+04, 4.8233829816851015e+04, 4.8242567590311344e+04, 4.8251305524169453e+04, 4.8260043618399635e+04, 4.8268781872976157e+04, 4.8277520287873325e+04, 4.8286258863065443e+04, 4.8294997598526803e+04, 4.8303736494231736e+04, 4.8312475550154573e+04, 4.8321214766269630e+04, 4.8329954142551265e+04, 4.8338693678973825e+04, 4.8347433375511668e+04, 4.8356173232139161e+04, 4.8364913248830679e+04, 4.8373653425560617e+04, 4.8382393762303349e+04, 4.8391134259033279e+04, 4.8399874915724817e+04, 4.8408615732352380e+04, 4.8417356708890395e+04, 4.8426097845313299e+04, 4.8434839141595519e+04, 4.8443580597711516e+04, 4.8452322213635751e+04, 4.8461063989342671e+04, 4.8469805924806766e+04, 4.8478548020002505e+04, 4.8487290274904393e+04, 4.8496032689486914e+04, 4.8504775263724587e+04, 4.8513517997591916e+04, 4.8522260891063430e+04, 4.8531003944113654e+04, 4.8539747156717123e+04, 4.8548490528848401e+04, 4.8557234060482027e+04, 4.8565977751592574e+04, 4.8574721602154605e+04, 4.8583465612142703e+04, 4.8592209781531456e+04, 4.8600954110295461e+04, 4.8609698598409312e+04, 4.8618443245847622e+04, 4.8627188052585028e+04, 4.8635933018596137e+04, 4.8644678143855606e+04, 4.8653423428338065e+04, 4.8662168872018163e+04, 4.8670914474870566e+04, 4.8679660236869939e+04, 4.8688406157990961e+04, 4.8697152238208320e+04, 4.8705898477496696e+04, 4.8714644875830811e+04, 4.8723391433185352e+04, 4.8732138149535058e+04, 4.8740885024854630e+04, 4.8749632059118805e+04, 4.8758379252302329e+04, 4.8767126604379955e+04, 4.8775874115326435e+04, 4.8784621785116535e+04, 4.8793369613725023e+04, 4.8802117601126673e+04, 4.8810865747296288e+04, 4.8819614052208672e+04, 4.8828362515838613e+04, 4.8837111138160937e+04, 4.8845859919150455e+04, 4.8854608858781998e+04, 4.8863357957030392e+04, 4.8872107213870506e+04, 4.8880856629277172e+04, 4.8889606203225252e+04, 4.8898355935689622e+04, 4.8907105826645158e+04, 4.8915855876066737e+04, 4.8924606083929270e+04, 4.8933356450207641e+04, 4.8942106974876755e+04, 4.8950857657911540e+04, 4.8959608499286915e+04, 4.8968359498977807e+04, 4.8977110656959172e+04, 4.8985861973205945e+04, 4.8994613447693089e+04, 4.9003365080395561e+04, 4.9012116871288337e+04, 4.9020868820346397e+04, 4.9029620927544725e+04, 4.9038373192858322e+04, 4.9047125616262194e+04, 4.9055878197731334e+04, 4.9064630937240785e+04, 4.9073383834765555e+04, 4.9082136890280693e+04, 4.9090890103761238e+04, 4.9099643475182231e+04, 4.9108397004518753e+04, 4.9117150691745839e+04, 4.9125904536838600e+04, 4.9134658539772099e+04, 4.9143412700521425e+04, 4.9152167019061671e+04, 4.9160921495367955e+04, 4.9169676129415391e+04, 4.9178430921179097e+04, 4.9187185870634188e+04, 4.9195940977755825e+04, 4.9204696242519138e+04, 4.9213451664899279e+04, 4.9222207244871417e+04, 4.9230962982410725e+04, 4.9239718877492385e+04, 4.9248474930091557e+04, 4.9257231140183445e+04, 4.9265987507743244e+04, 4.9274744032746174e+04, 4.9283500715167436e+04, 4.9292257554982265e+04, 4.9301014552165885e+04, 4.9309771706693544e+04, 4.9318529018540488e+04, 4.9327286487681959e+04, 4.9336044114093223e+04, 4.9344801897749559e+04, 4.9353559838626243e+04, 4.9362317936698557e+04, 4.9371076191941786e+04, 4.9379834604331241e+04, 4.9388593173842230e+04, 4.9397351900450063e+04, 4.9406110784130084e+04, 4.9414869824857604e+04, 4.9423629022607973e+04, 4.9432388377356539e+04, 4.9441147889078653e+04, 4.9449907557749684e+04, 4.9458667383345004e+04, 4.9467427365839983e+04, 4.9476187505210000e+04, 4.9484947801430484e+04, 4.9493708254476798e+04, 4.9502468864324372e+04, 4.9511229630948612e+04, 4.9519990554324948e+04, 4.9528751634428816e+04, 4.9537512871235653e+04, 4.9546274264720909e+04, 4.9555035814860050e+04, 4.9563797521628512e+04, 4.9572559385001783e+04, 4.9581321404955350e+04, 4.9590083581464685e+04, 4.9598845914505284e+04, 4.9607608404052655e+04, 4.9616371050082300e+04, 4.9625133852569743e+04, 4.9633896811490515e+04, 4.9642659926820132e+04, 4.9651423198534139e+04, 4.9660186626608091e+04, 4.9668950211017545e+04, 4.9677713951738056e+04, 4.9686477848745191e+04, 4.9695241902014546e+04, 4.9704006111521689e+04, 4.9712770477242222e+04, 4.9721534999151743e+04, 4.9730299677225856e+04, 4.9739064511440185e+04, 4.9747829501770357e+04, 4.9756594648191996e+04, 4.9765359950680744e+04, 4.9774125409212254e+04, 4.9782891023762168e+04, 4.9791656794306153e+04, 4.9800422720819894e+04, 4.9809188803279038e+04, 4.9817955041659283e+04, 4.9826721435936328e+04, 4.9835487986085878e+04, 4.9844254692083632e+04, 4.9853021553905310e+04, 4.9861788571526617e+04, 4.9870555744923295e+04, 4.9879323074071093e+04, 4.9888090558945740e+04, 4.9896858199522998e+04, 4.9905625995778624e+04, 4.9914393947688382e+04, 4.9923162055228058e+04, 4.9931930318373430e+04, 4.9940698737100283e+04, 4.9949467311384426e+04, 4.9958236041201657e+04, 4.9967004926527792e+04, 4.9975773967338653e+04, 4.9984543163610069e+04, 4.9993312515317863e+04, 5.0002082022437891e+04, 5.0010851684946007e+04, 5.0019621502818067e+04, 5.0028391476029923e+04, 5.0037161604557456e+04, 5.0045931888376552e+04, 5.0054702327463099e+04, 5.0063472921792993e+04, 5.0072243671342134e+04, 5.0081014576086432e+04, 5.0089785636001805e+04, 5.0098556851064175e+04, 5.0107328221249489e+04, 5.0116099746533670e+04, 5.0124871426892685e+04, 5.0133643262302474e+04, 5.0142415252739003e+04, 5.0151187398178248e+04, 5.0159959698596183e+04, 5.0168732153968798e+04, 5.0177504764272067e+04, 5.0186277529482017e+04, 5.0195050449574643e+04, 5.0203823524525971e+04, 5.0212596754311999e+04, 5.0221370138908780e+04, 5.0230143678292334e+04, 5.0238917372438722e+04, 5.0247691221323985e+04, 5.0256465224924192e+04, 5.0265239383215398e+04, 5.0274013696173679e+04, 5.0282788163775127e+04, 5.0291562785995826e+04, 5.0300337562811867e+04, 5.0309112494199362e+04, 5.0317887580134418e+04, 5.0326662820593148e+04, 5.0335438215551694e+04, 5.0344213764986176e+04, 5.0352989468872751e+04, 5.0361765327187546e+04, 5.0370541339906726e+04, 5.0379317507006454e+04, 5.0388093828462908e+04, 5.0396870304252254e+04, 5.0405646934350676e+04, 5.0414423718734375e+04, 5.0423200657379552e+04, 5.0431977750262406e+04, 5.0440754997359159e+04, 5.0449532398646021e+04, 5.0458309954099233e+04, 5.0467087663695027e+04, 5.0475865527409645e+04, 5.0484643545219340e+04, 5.0493421717100377e+04, 5.0502200043029014e+04, 5.0510978522981517e+04, 5.0519757156934189e+04, 5.0528535944863288e+04, 5.0537314886745124e+04, 5.0546093982556005e+04, 5.0554873232272235e+04, 5.0563652635870138e+04, 5.0572432193326022e+04, 5.0581211904616233e+04, 5.0589991769717097e+04, 5.0598771788604965e+04, 5.0607551961256184e+04, 5.0616332287647136e+04, 5.0625112767754166e+04, 5.0633893401553658e+04, 5.0642674189021993e+04, 5.0651455130135575e+04, 5.0660236224870772e+04, 5.0669017473204010e+04, 5.0677798875111694e+04, 5.0686580430570248e+04, 5.0695362139556077e+04, 5.0704144002045636e+04, 5.0712926018015358e+04, 5.0721708187441684e+04, 5.0730490510301082e+04, 5.0739272986570009e+04, 5.0748055616224927e+04, 5.0756838399242326e+04, 5.0765621335598669e+04, 5.0774404425270463e+04, 5.0783187668234204e+04, 5.0791971064466401e+04, 5.0800754613943558e+04, 5.0809538316642196e+04, 5.0818322172538836e+04, 5.0827106181610019e+04, 5.0835890343832289e+04, 5.0844674659182194e+04, 5.0853459127636284e+04, 5.0862243749171124e+04, 5.0871028523763285e+04, 5.0879813451389346e+04, 5.0888598532025877e+04, 5.0897383765649494e+04, 5.0906169152236784e+04, 5.0914954691764346e+04, 5.0923740384208795e+04, 5.0932526229546755e+04, 5.0941312227754854e+04, 5.0950098378809729e+04, 5.0958884682688004e+04, 5.0967671139366350e+04, 5.0976457748821413e+04, 5.0985244511029858e+04, 5.0994031425968351e+04, 5.1002818493613566e+04, 5.1011605713942183e+04, 5.1020393086930912e+04, 5.1029180612556447e+04, 5.1037968290795485e+04, 5.1046756121624741e+04, 5.1055544105020941e+04, 5.1064332240960801e+04, 5.1073120529421059e+04, 5.1081908970378470e+04, 5.1090697563809757e+04, 5.1099486309691696e+04, 5.1108275208001047e+04, 5.1117064258714563e+04, 5.1125853461809042e+04, 5.1134642817261250e+04, 5.1143432325047994e+04, 5.1152221985146054e+04, 5.1161011797532236e+04, 5.1169801762183371e+04, 5.1178591879076266e+04, 5.1187382148187738e+04, 5.1196172569494644e+04, 5.1204963142973793e+04, 5.1213753868602056e+04, 5.1222544746356267e+04, 5.1231335776213316e+04, 5.1240126958150053e+04, 5.1248918292143353e+04, 5.1257709778170109e+04, 5.1266501416207189e+04, 5.1275293206231509e+04, 5.1284085148219965e+04, 5.1292877242149465e+04, 5.1301669487996936e+04, 5.1310461885739292e+04, 5.1319254435353469e+04, 5.1328047136816414e+04, 5.1336839990105051e+04, 5.1345632995196349e+04, 5.1354426152067266e+04, 5.1363219460694767e+04, 5.1372012921055815e+04, 5.1380806533127412e+04, 5.1389600296886529e+04, 5.1398394212310159e+04, 5.1407188279375310e+04, 5.1415982498058998e+04, 5.1424776868338209e+04, 5.1433571390190002e+04, 5.1442366063591384e+04, 5.1451160888519393e+04, 5.1459955864951087e+04, 5.1468750992863505e+04, 5.1477546272233689e+04, 5.1486341703038735e+04, 5.1495137285255689e+04, 5.1503933018861637e+04, 5.1512728903833668e+04, 5.1521524940148869e+04, 5.1530321127784344e+04, 5.1539117466717187e+04, 5.1547913956924524e+04, 5.1556710598383463e+04, 5.1565507391071136e+04, 5.1574304334964683e+04, 5.1583101430041228e+04, 5.1591898676277939e+04, 5.1600696073651947e+04, 5.1609493622140428e+04, 5.1618291321720542e+04, 5.1627089172369480e+04, 5.1635887174064395e+04, 5.1644685326782506e+04, 5.1653483630500989e+04, 5.1662282085197054e+04, 5.1671080690847906e+04, 5.1679879447430758e+04, 5.1688678354922849e+04, 5.1697477413301400e+04, 5.1706276622543643e+04, 5.1715075982626819e+04, 5.1723875493528190e+04, 5.1732675155225006e+04, 5.1741474967694521e+04, 5.1750274930914027e+04, 5.1759075044860801e+04, 5.1767875309512114e+04, 5.1776675724845263e+04, 5.1785476290837541e+04, 5.1794277007466255e+04, 5.1803077874708724e+04, 5.1811878892542270e+04, 5.1820680060944200e+04, 5.1829481379891862e+04, 5.1838282849362593e+04, 5.1847084469333742e+04, 5.1855886239782667e+04, 5.1864688160686717e+04, 5.1873490232023250e+04, 5.1882292453769653e+04, 5.1891094825903310e+04, 5.1899897348401588e+04, 5.1908700021241908e+04, 5.1917502844401650e+04, 5.1926305817858221e+04, 5.1935108941589053e+04, 5.1943912215571552e+04, 5.1952715639783150e+04, 5.1961519214201297e+04, 5.1970322938803409e+04, 5.1979126813566938e+04, 5.1987930838469358e+04, 5.1996735013488113e+04, 5.2005539338600669e+04, 5.2014343813784522e+04, 5.2023148439017139e+04, 5.2031953214276007e+04, 5.2040758139538630e+04, 5.2049563214782502e+04, 5.2058368439985141e+04, 5.2067173815124050e+04, 5.2075979340176760e+04, 5.2084785015120789e+04, 5.2093590839933691e+04, 5.2102396814593005e+04, 5.2111202939076276e+04, 5.2120009213361052e+04, 5.2128815637424908e+04, 5.2137622211245412e+04, 5.2146428934800133e+04, 5.2155235808066667e+04, 5.2164042831022598e+04, 5.2172850003645515e+04, 5.2181657325913031e+04, 5.2190464797802742e+04, 5.2199272419292283e+04, 5.2208080190359266e+04, 5.2216888110981323e+04, 5.2225696181136089e+04, 5.2234504400801212e+04, 5.2243312769954333e+04, 5.2252121288573107e+04, 5.2260929956635213e+04, 5.2269738774118319e+04, 5.2278547741000089e+04, 5.2287356857258215e+04, 5.2296166122870382e+04, 5.2304975537814280e+04, 5.2313785102067639e+04, 5.2322594815608150e+04, 5.2331404678413528e+04, 5.2340214690461493e+04, 5.2349024851729780e+04, 5.2357835162196141e+04, 5.2366645621838303e+04, 5.2375456230634016e+04, 5.2384266988561045e+04, 5.2393077895597147e+04, 5.2401888951720088e+04, 5.2410700156907653e+04, 5.2419511511137614e+04, 5.2428323014387781e+04, 5.2437134666635924e+04, 5.2445946467859867e+04, 5.2454758418037396e+04, 5.2463570517146363e+04, 5.2472382765164555e+04, 5.2481195162069824e+04, 5.2490007707839992e+04, 5.2498820402452911e+04, 5.2507633245886434e+04, 5.2516446238118398e+04, 5.2525259379126692e+04, 5.2534072668889159e+04, 5.2542886107383689e+04, 5.2551699694588162e+04, 5.2560513430480460e+04, 5.2569327315038492e+04, 5.2578141348240140e+04, 5.2586955530063329e+04, 5.2595769860485969e+04, 5.2604584339485977e+04, 5.2613398967041292e+04, 5.2622213743129832e+04, 5.2631028667729552e+04, 5.2639843740818396e+04, 5.2648658962374313e+04, 5.2657474332375270e+04, 5.2666289850799236e+04, 5.2675105517624179e+04, 5.2683921332828082e+04, 5.2692737296388936e+04, 5.2701553408284737e+04, 5.2710369668493469e+04, 5.2719186076993159e+04, 5.2728002633761796e+04, 5.2736819338777423e+04, 5.2745636192018042e+04, 5.2754453193461712e+04, 5.2763270343086457e+04, 5.2772087640870319e+04, 5.2780905086791361e+04, 5.2789722680827639e+04, 5.2798540422957209e+04, 5.2807358313158162e+04, 5.2816176351408554e+04, 5.2824994537686478e+04, 5.2833812871970040e+04, 5.2842631354237310e+04, 5.2851449984466410e+04, 5.2860268762635445e+04, 5.2869087688722539e+04, 5.2877906762705810e+04, 5.2886725984563382e+04, 5.2895545354273396e+04, 5.2904364871814003e+04, 5.2913184537163339e+04, 5.2922004350299561e+04, 5.2930824311200849e+04, 5.2939644419845354e+04, 5.2948464676211268e+04, 5.2957285080276757e+04, 5.2966105632020008e+04, 5.2974926331419221e+04, 5.2983747178452599e+04, 5.2992568173098356e+04, 5.3001389315334687e+04, 5.3010210605139822e+04, 5.3019032042491985e+04, 5.3027853627369419e+04, 5.3036675359750356e+04, 5.3045497239613032e+04, 5.3054319266935716e+04, 5.3063141441696665e+04, 5.3071963763874141e+04, 5.3080786233446408e+04, 5.3089608850391749e+04, 5.3098431614688459e+04, 5.3107254526314813e+04, 5.3116077585249121e+04, 5.3124900791469670e+04, 5.3133724144954787e+04, 5.3142547645682782e+04, 5.3151371293631972e+04, 5.3160195088780689e+04, 5.3169019031107273e+04, 5.3177843120590071e+04, 5.3186667357207414e+04, 5.3195491740937672e+04, 5.3204316271759199e+04, 5.3213140949650362e+04, 5.3221965774589538e+04, 5.3230790746555103e+04, 5.3239615865525455e+04, 5.3248441131478969e+04, 5.3257266544394050e+04, 5.3266092104249117e+04, 5.3274917811022562e+04, 5.3283743664692818e+04, 5.3292569665238298e+04, 5.3301395812637442e+04, 5.3310222106868678e+04, 5.3319048547910461e+04, 5.3327875135741233e+04, 5.3336701870339457e+04, 5.3345528751683589e+04, 5.3354355779752106e+04, 5.3363182954523465e+04, 5.3372010275976165e+04, 5.3380837744088683e+04, 5.3389665358839527e+04, 5.3398493120207189e+04, 5.3407321028170161e+04, 5.3416149082706972e+04, 5.3424977283796143e+04, 5.3433805631416195e+04, 5.3442634125545657e+04, 5.3451462766163066e+04, 5.3460291553246985e+04, 5.3469120486775944e+04, 5.3477949566728508e+04, 5.3486778793083242e+04, 5.3495608165818710e+04, 5.3504437684913486e+04, 5.3513267350346170e+04, 5.3522097162095328e+04, 5.3530927120139560e+04, 5.3539757224457477e+04, 5.3548587475027671e+04, 5.3557417871828766e+04, 5.3566248414839385e+04, 5.3575079104038145e+04, 5.3583909939403682e+04, 5.3592740920914628e+04, 5.3601572048549642e+04, 5.3610403322287370e+04, 5.3619234742106455e+04, 5.3628066307985573e+04, 5.3636898019903398e+04, 5.3645729877838596e+04, 5.3654561881769849e+04, 5.3663394031675853e+04, 5.3672226327535303e+04, 5.3681058769326883e+04, 5.3689891357029323e+04, 5.3698724090621319e+04, 5.3707556970081598e+04, 5.3716389995388883e+04, 5.3725223166521908e+04, 5.3734056483459404e+04, 5.3742889946180127e+04, 5.3751723554662814e+04, 5.3760557308886229e+04, 5.3769391208829133e+04, 5.3778225254470301e+04, 5.3787059445788502e+04, 5.3795893782762520e+04, 5.3804728265371145e+04, 5.3813562893593160e+04, 5.3822397667407378e+04, 5.3831232586792597e+04, 5.3840067651727630e+04, 5.3848902862191295e+04, 5.3857738218162413e+04, 5.3866573719619824e+04, 5.3875409366542357e+04, 5.3884245158908860e+04, 5.3893081096698181e+04, 5.3901917179889169e+04, 5.3910753408460696e+04, 5.3919589782391624e+04, 5.3928426301660831e+04, 5.3937262966247188e+04, 5.3946099776129588e+04, 5.3954936731286922e+04, 5.3963773831698083e+04, 5.3972611077341986e+04, 5.3981448468197530e+04, 5.3990286004243644e+04, 5.3999123685459243e+04, 5.4007961511823247e+04, 5.4016799483314615e+04, 5.4025637599912261e+04, 5.4034475861595150e+04, 5.4043314268342219e+04, 5.4052152820132440e+04, 5.4060991516944785e+04, 5.4069830358758220e+04, 5.4078669345551723e+04, 5.4087508477304262e+04, 5.4096347753994851e+04, 5.4105187175602470e+04, 5.4114026742106129e+04, 5.4122866453484836e+04, 5.4131706309717600e+04, 5.4140546310783437e+04, 5.4149386456661385e+04, 5.4158226747330467e+04, 5.4167067182769737e+04, 5.4175907762958224e+04, 5.4184748487874982e+04, 5.4193589357499077e+04, 5.4202430371809554e+04, 5.4211271530785503e+04, 5.4220112834405976e+04, 5.4228954282650084e+04, 5.4237795875496879e+04, 5.4246637612925479e+04, 5.4255479494914973e+04, 5.4264321521444479e+04, 5.4273163692493086e+04, 5.4282006008039934e+04, 5.4290848468064134e+04, 5.4299691072544811e+04, 5.4308533821461118e+04, 5.4317376714792175e+04, 5.4326219752517136e+04, 5.4335062934615155e+04, 5.4343906261065400e+04, 5.4352749731847034e+04, 5.4361593346939211e+04, 5.4370437106321122e+04, 5.4379281009971957e+04, 5.4388125057870900e+04, 5.4396969249997150e+04, 5.4405813586329903e+04, 5.4414658066848358e+04, 5.4423502691531743e+04, 5.4432347460359277e+04, 5.4441192373310172e+04, 5.4450037430363678e+04, 5.4458882631499015e+04, 5.4467727976695438e+04, 5.4476573465932190e+04, 5.4485419099188533e+04, 5.4494264876443718e+04, 5.4503110797677022e+04, 5.4511956862867715e+04, 5.4520803071995084e+04, 5.4529649425038398e+04, 5.4538495921976952e+04, 5.4547342562790051e+04, 5.4556189347456988e+04, 5.4565036275957085e+04, 5.4573883348269650e+04, 5.4582730564374011e+04, 5.4591577924249483e+04, 5.4600425427875402e+04, 5.4609273075231111e+04, 5.4618120866295962e+04, 5.4626968801049290e+04, 5.4635816879470454e+04, 5.4644665101538827e+04, 5.4653513467233766e+04, 5.4662361976534652e+04, 5.4671210629420864e+04, 5.4680059425871783e+04, 5.4688908365866817e+04, 5.4697757449385354e+04, 5.4706606676406795e+04, 5.4715456046910549e+04, 5.4724305560876041e+04, 5.4733155218282685e+04, 5.4742005019109907e+04, 5.4750854963337144e+04, 5.4759705050943834e+04, 5.4768555281909423e+04, 5.4777405656213363e+04, 5.4786256173835114e+04, 5.4795106834754130e+04, 5.4803957638949891e+04, 5.4812808586401858e+04, 5.4821659677089512e+04, 5.4830510910992365e+04, 5.4839362288089884e+04, 5.4848213808361572e+04, 5.4857065471786940e+04, 5.4865917278345492e+04, 5.4874769228016747e+04, 5.4883621320780228e+04, 5.4892473556615456e+04, 5.4901325935501969e+04, 5.4910178457419308e+04, 5.4919031122347013e+04, 5.4927883930264630e+04, 5.4936736881151730e+04, 5.4945589974987874e+04, 5.4954443211752623e+04, 5.4963296591425547e+04, 5.4972150113986238e+04, 5.4981003779414270e+04, 5.4989857587689250e+04, 5.4998711538790762e+04, 5.5007565632698424e+04, 5.5016419869391830e+04, 5.5025274248850605e+04, 5.5034128771054362e+04, 5.5042983435982729e+04, 5.5051838243615348e+04, 5.5060693193931846e+04, 5.5069548286911871e+04, 5.5078403522535074e+04, 5.5087258900781118e+04, 5.5096114421629652e+04, 5.5104970085060355e+04, 5.5113825891052889e+04, 5.5122681839586941e+04, 5.5131537930642196e+04, 5.5140394164198340e+04, 5.5149250540235065e+04, 5.5158107058732079e+04, 5.5166963719669104e+04, 5.5175820523025832e+04, 5.5184677468781985e+04, 5.5193534556917301e+04, 5.5202391787411507e+04, 5.5211249160244326e+04, 5.5220106675395524e+04, 5.5228964332844822e+04, 5.5237822132572008e+04, 5.5246680074556818e+04, 5.5255538158779011e+04, 5.5264396385218373e+04, 5.5273254753854686e+04, 5.5282113264667722e+04, 5.5290971917637275e+04, 5.5299830712743125e+04, 5.5308689649965090e+04, 5.5317548729282971e+04, 5.5326407950676585e+04, 5.5335267314125733e+04, 5.5344126819610254e+04, 5.5352986467109971e+04, 5.5361846256604716e+04, 5.5370706188074335e+04, 5.5379566261498665e+04, 5.5388426476857567e+04, 5.5397286834130893e+04, 5.5406147333298504e+04, 5.5415007974340275e+04, 5.5423868757236080e+04, 5.5432729681965793e+04, 5.5441590748509319e+04, 5.5450451956846518e+04, 5.5459313306957323e+04, 5.5468174798821608e+04, 5.5477036432419300e+04, 5.5485898207730301e+04, 5.5494760124734537e+04, 5.5503622183411928e+04, 5.5512484383742420e+04, 5.5521346725705931e+04, 5.5530209209282417e+04, 5.5539071834451832e+04, 5.5547934601194116e+04, 5.5556797509489232e+04, 5.5565660559317141e+04, 5.5574523750657834e+04, 5.5583387083491274e+04, 5.5592250557797437e+04, 5.5601114173556329e+04, 5.5609977930747933e+04, 5.5618841829352248e+04, 5.5627705869349287e+04, 5.5636570050719049e+04, 5.5645434373441552e+04, 5.5654298837496834e+04, 5.5663163442864905e+04, 5.5672028189525809e+04, 5.5680893077459586e+04, 5.5689758106646274e+04, 5.5698623277065926e+04, 5.5707488588698601e+04, 5.5716354041524348e+04, 5.5725219635523259e+04, 5.5734085370675377e+04, 5.5742951246960802e+04, 5.5751817264359619e+04, 5.5760683422851900e+04, 5.5769549722417752e+04, 5.5778416163037284e+04, 5.5787282744690587e+04, 5.5796149467357784e+04, 5.5805016331018989e+04, 5.5813883335654318e+04, 5.5822750481243922e+04, 5.5831617767767908e+04, 5.5840485195206435e+04, 5.5849352763539646e+04, 5.5858220472747686e+04, 5.5867088322810720e+04, 5.5875956313708906e+04, 5.5884824445422411e+04, 5.5893692717931415e+04, 5.5902561131216084e+04, 5.5911429685256619e+04, 5.5920298380033193e+04, 5.5929167215526017e+04, 5.5938036191715291e+04, 5.5946905308581219e+04, 5.5955774566104024e+04, 5.5964643964263909e+04, 5.5973513503041097e+04, 5.5982383182415833e+04, 5.5991253002368343e+04, 5.6000122962878864e+04, 5.6008993063927650e+04, 5.6017863305494946e+04, 5.6026733687561013e+04, 5.6035604210106110e+04, 5.6044474873110521e+04, 5.6053345676554505e+04, 5.6062216620418338e+04, 5.6071087704682315e+04, 5.6079958929326727e+04, 5.6088830294331870e+04, 5.6097701799678034e+04, 5.6106573445345544e+04, 5.6115445231314690e+04, 5.6124317157565813e+04, 5.6133189224079222e+04, 5.6142061430835245e+04, 5.6150933777814229e+04, 5.6159806264996507e+04, 5.6168678892362419e+04, 5.6177551659892335e+04, 5.6186424567566595e+04, 5.6195297615365569e+04, 5.6204170803269626e+04, 5.6213044131259128e+04, 5.6221917599314460e+04, 5.6230791207416012e+04, 5.6239664955544169e+04, 5.6248538843679329e+04, 5.6257412871801884e+04, 5.6266287039892246e+04, 5.6275161347930829e+04, 5.6284035795898053e+04, 5.6292910383774331e+04, 5.6301785111540099e+04, 5.6310659979175784e+04, 5.6319534986661834e+04, 5.6328410133978687e+04, 5.6337285421106797e+04, 5.6346160848026615e+04, 5.6355036414718605e+04, 5.6363912121163237e+04, 5.6372787967340977e+04, 5.6381663953232310e+04, 5.6390540078817699e+04, 5.6399416344077654e+04, 5.6408292748992659e+04, 5.6417169293543222e+04, 5.6426045977709844e+04, 5.6434922801473025e+04, 5.6443799764813280e+04, 5.6452676867711147e+04, 5.6461554110147139e+04, 5.6470431492101794e+04, 5.6479309013555649e+04, 5.6488186674489247e+04, 5.6497064474883133e+04, 5.6505942414717851e+04, 5.6514820493973981e+04, 5.6523698712632075e+04, 5.6532577070672705e+04, 5.6541455568076446e+04, 5.6550334204823870e+04, 5.6559212980895580e+04, 5.6568091896272155e+04, 5.6576970950934199e+04, 5.6585850144862306e+04, 5.6594729478037101e+04, 5.6603608950439171e+04, 5.6612488562049155e+04, 5.6621368312847670e+04, 5.6630248202815339e+04, 5.6639128231932809e+04, 5.6648008400180712e+04, 5.6656888707539692e+04, 5.6665769153990404e+04, 5.6674649739513508e+04, 5.6683530464089657e+04, 5.6692411327699519e+04, 5.6701292330323777e+04, 5.6710173471943104e+04, 5.6719054752538177e+04, 5.6727936172089685e+04, 5.6736817730578332e+04, 5.6745699427984800e+04, 5.6754581264289802e+04, 5.6763463239474055e+04, 5.6772345353518263e+04, 5.6781227606403147e+04, 5.6790109998109445e+04, 5.6798992528617884e+04, 5.6807875197909190e+04, 5.6816758005964111e+04, 5.6825640952763395e+04, 5.6834524038287804e+04, 5.6843407262518078e+04, 5.6852290625434995e+04, 5.6861174127019316e+04, 5.6870057767251819e+04, 5.6878941546113274e+04, 5.6887825463584486e+04, 5.6896709519646225e+04, 5.6905593714279297e+04, 5.6914478047464501e+04, 5.6923362519182636e+04, 5.6932247129414522e+04, 5.6941131878140972e+04, 5.6950016765342800e+04, 5.6958901791000848e+04, 5.6967786955095944e+04, 5.6976672257608923e+04, 5.6985557698520628e+04, 5.6994443277811908e+04, 5.7003328995463620e+04, 5.7012214851456629e+04, 5.7021100845771783e+04, 5.7029986978389963e+04, 5.7038873249292032e+04, 5.7047759658458875e+04, 5.7056646205871395e+04, 5.7065532891510462e+04, 5.7074419715356977e+04, 5.7083306677391847e+04, 5.7092193777595967e+04, 5.7101081015950258e+04, 5.7109968392435643e+04, 5.7118855907033030e+04, 5.7127743559723356e+04, 5.7136631350487543e+04, 5.7145519279306551e+04, 5.7154407346161308e+04, 5.7163295551032767e+04, 5.7172183893901871e+04, 5.7181072374749587e+04, 5.7189960993556888e+04, 5.7198849750304733e+04, 5.7207738644974095e+04, 5.7216627677545970e+04, 5.7225516848001331e+04, 5.7234406156321180e+04, 5.7243295602486498e+04, 5.7252185186478295e+04, 5.7261074908277580e+04, 5.7269964767865349e+04, 5.7278854765222626e+04, 5.7287744900330443e+04, 5.7296635173169831e+04, 5.7305525583721799e+04, 5.7314416131967409e+04, 5.7323306817887686e+04, 5.7332197641463696e+04, 5.7341088602676471e+04, 5.7349979701507087e+04, 5.7358870937936605e+04, 5.7367762311946091e+04, 5.7376653823516615e+04, 5.7385545472629266e+04, 5.7394437259265120e+04, 5.7403329183405280e+04, 5.7412221245030829e+04, 5.7421113444122864e+04, 5.7430005780662505e+04, 5.7438898254630847e+04, 5.7447790866009025e+04, 5.7456683614778136e+04, 5.7465576500919327e+04, 5.7474469524413726e+04, 5.7483362685242457e+04, 5.7492255983386676e+04, 5.7501149418827525e+04, 5.7510042991546157e+04, 5.7518936701523722e+04, 5.7527830548741396e+04, 5.7536724533180328e+04, 5.7545618654821716e+04, 5.7554512913646715e+04, 5.7563407309636525e+04, 5.7572301842772322e+04, 5.7581196513035306e+04, 5.7590091320406675e+04, 5.7598986264867628e+04, 5.7607881346399379e+04, 5.7616776564983149e+04, 5.7625671920600151e+04, 5.7634567413231598e+04, 5.7643463042858741e+04, 5.7652358809462792e+04, 5.7661254713025010e+04, 5.7670150753526628e+04, 5.7679046930948905e+04, 5.7687943245273098e+04, 5.7696839696480449e+04, 5.7705736284552237e+04, 5.7714633009469741e+04, 5.7723529871214218e+04, 5.7732426869766961e+04, 5.7741324005109258e+04, 5.7750221277222394e+04, 5.7759118686087662e+04, 5.7768016231686364e+04, 5.7776913913999822e+04, 5.7785811733009345e+04, 5.7794709688696224e+04, 5.7803607781041806e+04, 5.7812506010027413e+04, 5.7821404375634367e+04, 5.7830302877844013e+04, 5.7839201516637688e+04, 5.7848100291996743e+04, 5.7856999203902538e+04, 5.7865898252336425e+04, 5.7874797437279762e+04, 5.7883696758713915e+04, 5.7892596216620259e+04, 5.7901495810980188e+04, 5.7910395541775077e+04, 5.7919295408986298e+04, 5.7928195412595254e+04, 5.7937095552583356e+04, 5.7945995828931984e+04, 5.7954896241622570e+04, 5.7963796790636501e+04, 5.7972697475955218e+04, 5.7981598297560136e+04, 5.7990499255432682e+04, 5.7999400349554293e+04, 5.8008301579906401e+04, 5.8017202946470454e+04, 5.8026104449227903e+04, 5.8035006088160197e+04, 5.8043907863248809e+04, 5.8052809774475194e+04, 5.8061711821820820e+04, 5.8070614005267154e+04, 5.8079516324795688e+04, 5.8088418780387896e+04, 5.8097321372025275e+04, 5.8106224099689309e+04, 5.8115126963361508e+04, 5.8124029963023379e+04, 5.8132933098656409e+04, 5.8141836370242134e+04, 5.8150739777762072e+04, 5.8159643321197735e+04, 5.8168547000530656e+04, 5.8177450815742384e+04, 5.8186354766814438e+04, 5.8195258853728381e+04, 5.8204163076465746e+04, 5.8213067435008095e+04, 5.8221971929336985e+04, 5.8230876559433993e+04, 5.8239781325280681e+04, 5.8248686226858612e+04, 5.8257591264149378e+04, 5.8266496437134556e+04, 5.8275401745795745e+04, 5.8284307190114545e+04, 5.8293212770072540e+04, 5.8302118485651328e+04, 5.8311024336832539e+04, 5.8319930323597779e+04, 5.8328836445928660e+04, 5.8337742703806813e+04, 5.8346649097213871e+04, 5.8355555626131456e+04, 5.8364462290541225e+04, 5.8373369090424821e+04, 5.8382276025763887e+04, 5.8391183096540073e+04, 5.8400090302735043e+04, 5.8408997644330448e+04, 5.8417905121307980e+04, 5.8426812733649291e+04, 5.8435720481336080e+04, 5.8444628364350021e+04, 5.8453536382672806e+04, 5.8462444536286122e+04, 5.8471352825171685e+04, 5.8480261249311181e+04, 5.8489169808686325e+04, 5.8498078503278826e+04, 5.8506987333070421e+04, 5.8515896298042819e+04, 5.8524805398177748e+04, 5.8533714633456933e+04, 5.8542624003862133e+04, 5.8551533509375076e+04, 5.8560443149977516e+04, 5.8569352925651212e+04, 5.8578262836377915e+04, 5.8587172882139392e+04, 5.8596083062917402e+04, 5.8604993378693725e+04, 5.8613903829450144e+04, 5.8622814415168432e+04, 5.8631725135830384e+04, 5.8640635991417788e+04, 5.8649546981912456e+04, 5.8658458107296159e+04, 5.8667369367550731e+04, 5.8676280762657974e+04, 5.8685192292599706e+04, 5.8694103957357758e+04, 5.8703015756913941e+04, 5.8711927691250108e+04, 5.8720839760348077e+04, 5.8729751964189702e+04, 5.8738664302756821e+04, 5.8747576776031288e+04, 5.8756489383994958e+04, 5.8765402126629691e+04, 5.8774315003917363e+04, 5.8783228015839835e+04, 5.8792141162378990e+04, 5.8801054443516703e+04, 5.8809967859234857e+04, 5.8818881409515358e+04, 5.8827795094340087e+04, 5.8836708913690942e+04, 5.8845622867549835e+04, 5.8854536955898671e+04, 5.8863451178719377e+04, 5.8872365535993857e+04, 5.8881280027704051e+04, 5.8890194653831873e+04, 5.8899109414359271e+04, 5.8908024309268185e+04, 5.8916939338540542e+04, 5.8925854502158312e+04, 5.8934769800103430e+04, 5.8943685232357857e+04, 5.8952600798903564e+04, 5.8961516499722522e+04, 5.8970432334796693e+04, 5.8979348304108069e+04, 5.8988264407638620e+04, 5.8997180645370339e+04, 5.9006097017285218e+04, 5.9015013523365262e+04, 5.9023930163592457e+04, 5.9032846937948831e+04, 5.9041763846416376e+04, 5.9050680888977113e+04, 5.9059598065613063e+04, 5.9068515376306263e+04, 5.9077432821038739e+04, 5.9086350399792515e+04, 5.9095268112549646e+04, 5.9104185959292176e+04, 5.9113103940002147e+04, 5.9122022054661626e+04, 5.9130940303252653e+04, 5.9139858685757317e+04, 5.9148777202157675e+04, 5.9157695852435812e+04, 5.9166614636573788e+04, 5.9175533554553702e+04, 5.9184452606357627e+04, 5.9193371791967678e+04, 5.9202291111365928e+04, 5.9211210564534500e+04, 5.9220130151455494e+04, 5.9229049872111027e+04, 5.9237969726483221e+04, 5.9246889714554192e+04, 5.9255809836306056e+04, 5.9264730091720972e+04, 5.9273650480781056e+04, 5.9282571003468445e+04, 5.9291491659765306e+04, 5.9300412449653777e+04, 5.9309333373116002e+04, 5.9318254430134162e+04, 5.9327175620690410e+04, 5.9336096944766912e+04, 5.9345018402345864e+04, 5.9353939993409425e+04, 5.9362861717939792e+04, 5.9371783575919144e+04, 5.9380705567329678e+04, 5.9389627692153590e+04, 5.9398549950373090e+04, 5.9407472341970388e+04, 5.9416394866927687e+04, 5.9425317525227212e+04, 5.9434240316851181e+04, 5.9443163241781818e+04, 5.9452086300001363e+04, 5.9461009491492056e+04, 5.9469932816236120e+04, 5.9478856274215817e+04, 5.9487779865413388e+04, 5.9496703589811092e+04, 5.9505627447391191e+04, 5.9514551438135946e+04, 5.9523475562027634e+04, 5.9532399819048522e+04, 5.9541324209180893e+04, 5.9550248732407024e+04, 5.9559173388709212e+04, 5.9568098178069748e+04, 5.9577023100470928e+04, 5.9585948155895043e+04, 5.9594873344324427e+04, 5.9603798665741371e+04, 5.9612724120128201e+04, 5.9621649707467230e+04, 5.9630575427740790e+04, 5.9639501280931210e+04, 5.9648427267020830e+04, 5.9657353385991984e+04, 5.9666279637827021e+04, 5.9675206022508290e+04, 5.9684132540018138e+04, 5.9693059190338929e+04, 5.9701985973453033e+04, 5.9710912889342813e+04, 5.9719839937990648e+04, 5.9728767119378899e+04, 5.9737694433489953e+04, 5.9746621880306200e+04, 5.9755549459810041e+04, 5.9764477171983868e+04, 5.9773405016810080e+04, 5.9782332994271084e+04, 5.9791261104349280e+04, 5.9800189347027095e+04, 5.9809117722286945e+04, 5.9818046230111257e+04, 5.9826974870482460e+04, 5.9835903643382990e+04, 5.9844832548795275e+04, 5.9853761586701759e+04, 5.9862690757084907e+04, 5.9871620059927147e+04, 5.9880549495210951e+04, 5.9889479062918777e+04, 5.9898408763033083e+04, 5.9907338595536356e+04, 5.9916268560411059e+04, 5.9925198657639674e+04, 5.9934128887204701e+04, 5.9943059249088605e+04, 5.9951989743273902e+04, 5.9960920369743071e+04, 5.9969851128478636e+04, 5.9978782019463077e+04, 5.9987713042678937e+04, 5.9996644198108719e+04, 6.0005575485734946e+04, 6.0014506905540133e+04, 6.0023438457506818e+04, 6.0032370141617554e+04, 6.0041301957854863e+04, 6.0050233906201291e+04, 6.0059165986639389e+04, 6.0068098199151711e+04, 6.0077030543720830e+04, 6.0085963020329291e+04, 6.0094895628959668e+04, 6.0103828369594528e+04, 6.0112761242216468e+04, 6.0121694246808052e+04, 6.0130627383351864e+04, 6.0139560651830492e+04, 6.0148494052226546e+04, 6.0157427584522629e+04, 6.0166361248701331e+04, 6.0175295044745260e+04, 6.0184228972637036e+04, 6.0193163032359291e+04, 6.0202097223894627e+04, 6.0211031547225677e+04, 6.0219966002335081e+04, 6.0228900589205477e+04, 6.0237835307819492e+04, 6.0246770158159779e+04, 6.0255705140208993e+04, 6.0264640253949787e+04, 6.0273575499364808e+04, 6.0282510876436740e+04, 6.0291446385148243e+04, 6.0300382025481995e+04, 6.0309317797420663e+04, 6.0318253700946938e+04, 6.0327189736043503e+04, 6.0336125902693057e+04, 6.0345062200878281e+04, 6.0353998630581889e+04, 6.0362935191786579e+04, 6.0371871884475069e+04, 6.0380808708630073e+04, 6.0389745664234295e+04, 6.0398682751270469e+04, 6.0407619969721323e+04, 6.0416557319569591e+04, 6.0425494800798013e+04, 6.0434432413389317e+04, 6.0443370157326266e+04, 6.0452308032591594e+04, 6.0461246039168058e+04, 6.0470184177038427e+04, 6.0479122446185465e+04, 6.0488060846591943e+04, 6.0496999378240631e+04, 6.0505938041114292e+04, 6.0514876835195726e+04, 6.0523815760467711e+04, 6.0532754816913039e+04, 6.0541694004514517e+04, 6.0550633323254937e+04, 6.0559572773117099e+04, 6.0568512354083825e+04, 6.0577452066137914e+04, 6.0586391909262195e+04, 6.0595331883439489e+04, 6.0604271988652617e+04, 6.0613212224884417e+04, 6.0622152592117731e+04, 6.0631093090335387e+04, 6.0640033719520230e+04, 6.0648974479655117e+04, 6.0657915370722905e+04, 6.0666856392706446e+04, 6.0675797545588604e+04, 6.0684738829352245e+04, 6.0693680243980249e+04, 6.0702621789455501e+04, 6.0711563465760861e+04, 6.0720505272879222e+04, 6.0729447210793485e+04, 6.0738389279486531e+04, 6.0747331478941262e+04, 6.0756273809140584e+04, 6.0765216270067402e+04, 6.0774158861704636e+04, 6.0783101584035197e+04, 6.0792044437042008e+04, 6.0800987420707999e+04, 6.0809930535016087e+04, 6.0818873779949223e+04, 6.0827817155490338e+04, 6.0836760661622378e+04, 6.0845704298328288e+04, 6.0854648065591020e+04, 6.0863591963393534e+04, 6.0872535991718789e+04, 6.0881480150549760e+04, 6.0890424439869414e+04, 6.0899368859660724e+04, 6.0908313409906666e+04, 6.0917258090590221e+04, 6.0926202901694392e+04, 6.0935147843202154e+04, 6.0944092915096509e+04, 6.0953038117360476e+04, 6.0961983449977030e+04, 6.0970928912929216e+04, 6.0979874506200023e+04, 6.0988820229772486e+04, 6.0997766083629613e+04, 6.1006712067754444e+04, 6.1015658182130006e+04, 6.1024604426739330e+04, 6.1033550801565478e+04, 6.1042497306591467e+04, 6.1051443941800375e+04, 6.1060390707175247e+04, 6.1069337602699132e+04, 6.1078284628355104e+04, 6.1087231784126227e+04, 6.1096179069995582e+04, 6.1105126485946224e+04, 6.1114074031961260e+04, 6.1123021708023756e+04, 6.1131969514116812e+04, 6.1140917450223526e+04, 6.1149865516326987e+04, 6.1158813712410301e+04, 6.1167762038456582e+04, 6.1176710494448926e+04, 6.1185659080370460e+04, 6.1194607796204298e+04, 6.1203556641933581e+04, 6.1212505617541414e+04, 6.1221454723010953e+04, 6.1230403958325332e+04, 6.1239353323467694e+04, 6.1248302818421173e+04, 6.1257252443168923e+04, 6.1266202197694111e+04, 6.1275152081979890e+04, 6.1284102096009417e+04, 6.1293052239765879e+04, 6.1302002513232437e+04, 6.1310952916392271e+04, 6.1319903449228565e+04, 6.1328854111724497e+04, 6.1337804903863260e+04, 6.1346755825628061e+04, 6.1355706877002085e+04, 6.1364658057968547e+04, 6.1373609368510646e+04, 6.1382560808611597e+04, 6.1391512378254614e+04, 6.1400464077422926e+04, 6.1409415906099748e+04, 6.1418367864268323e+04, 6.1427319951911864e+04, 6.1436272169013631e+04, 6.1445224515556853e+04, 6.1454176991524793e+04, 6.1463129596900682e+04, 6.1472082331667792e+04, 6.1481035195809367e+04, 6.1489988189308693e+04, 6.1498941312149022e+04, 6.1507894564313625e+04, 6.1516847945785790e+04, 6.1525801456548805e+04, 6.1534755096585926e+04, 6.1543708865880479e+04, 6.1552662764415734e+04, 6.1561616792175009e+04, 6.1570570949141591e+04, 6.1579525235298788e+04, 6.1588479650629924e+04, 6.1597434195118316e+04, 6.1606388868747272e+04, 6.1615343671500115e+04, 6.1624298603360192e+04, 6.1633253664310825e+04, 6.1642208854335353e+04, 6.1651164173417121e+04, 6.1660119621539474e+04, 6.1669075198685750e+04, 6.1678030904839325e+04, 6.1686986739983535e+04, 6.1695942704101777e+04, 6.1704898797177375e+04, 6.1713855019193739e+04, 6.1722811370134230e+04, 6.1731767849982229e+04, 6.1740724458721124e+04, 6.1749681196334306e+04, 6.1758638062805156e+04, 6.1767595058117084e+04, 6.1776552182253487e+04, 6.1785509435197775e+04, 6.1794466816933353e+04, 6.1803424327443652e+04, 6.1812381966712077e+04, 6.1821339734722053e+04, 6.1830297631456997e+04, 6.1839255656900357e+04, 6.1848213811035559e+04, 6.1857172093846071e+04, 6.1866130505315297e+04, 6.1875089045426714e+04, 6.1884047714163768e+04, 6.1893006511509906e+04, 6.1901965437448605e+04, 6.1910924491963327e+04, 6.1919883675037534e+04, 6.1928842986654709e+04, 6.1937802426798313e+04, 6.1946761995451860e+04, 6.1955721692598818e+04, 6.1964681518222686e+04, 6.1973641472306954e+04, 6.1982601554835121e+04, 6.1991561765790691e+04, 6.2000522105157186e+04, 6.2009482572918110e+04, 6.2018443169056969e+04, 6.2027403893557297e+04, 6.2036364746402622e+04, 6.2045325727576470e+04, 6.2054286837062362e+04, 6.2063248074843861e+04, 6.2072209440904488e+04, 6.2081170935227790e+04, 6.2090132557797333e+04, 6.2099094308596665e+04, 6.2108056187609342e+04, 6.2117018194818927e+04, 6.2125980330208993e+04, 6.2134942593763109e+04, 6.2143904985464847e+04, 6.2152867505297792e+04, 6.2161830153245537e+04, 6.2170792929291660e+04, 6.2179755833419746e+04, 6.2188718865613409e+04, 6.2197682025856251e+04, 6.2206645314131849e+04, 6.2215608730423846e+04, 6.2224572274715843e+04, 6.2233535946991462e+04, 6.2242499747234317e+04, 6.2251463675428036e+04, 6.2260427731556250e+04, 6.2269391915602610e+04, 6.2278356227550736e+04, 6.2287320667384280e+04, 6.2296285235086885e+04, 6.2305249930642211e+04, 6.2314214754033885e+04, 6.2323179705245602e+04, 6.2332144784261014e+04, 6.2341109991063786e+04, 6.2350075325637597e+04, 6.2359040787966114e+04, 6.2368006378033024e+04, 6.2376972095822006e+04, 6.2385937941316770e+04, 6.2394903914500981e+04, 6.2403870015358349e+04, 6.2412836243872574e+04, 6.2421802600027368e+04, 6.2430769083806430e+04, 6.2439735695193485e+04, 6.2448702434172243e+04, 6.2457669300726426e+04, 6.2466636294839773e+04, 6.2475603416496007e+04, 6.2484570665678868e+04, 6.2493538042372078e+04, 6.2502505546559398e+04, 6.2511473178224565e+04, 6.2520440937351334e+04, 6.2529408823923462e+04, 6.2538376837924712e+04, 6.2547344979338835e+04, 6.2556313248149600e+04, 6.2565281644340801e+04, 6.2574250167896200e+04, 6.2583218818799578e+04, 6.2592187597034710e+04, 6.2601156502585407e+04, 6.2610125535435436e+04, 6.2619094695568616e+04, 6.2628063982968742e+04, 6.2637033397619605e+04, 6.2646002939505037e+04, 6.2654972608608834e+04, 6.2663942404914815e+04, 6.2672912328406805e+04, 6.2681882379068637e+04, 6.2690852556884121e+04, 6.2699822861837121e+04, 6.2708793293911447e+04, 6.2717763853090953e+04, 6.2726734539359488e+04, 6.2735705352700905e+04, 6.2744676293099044e+04, 6.2753647360537776e+04, 6.2762618555000961e+04, 6.2771589876472462e+04, 6.2780561324936156e+04, 6.2789532900375918e+04, 6.2798504602775625e+04, 6.2807476432119154e+04, 6.2816448388390403e+04, 6.2825420471573256e+04, 6.2834392681651610e+04, 6.2843365018609358e+04, 6.2852337482430412e+04, 6.2861310073098692e+04, 6.2870282790598081e+04, 6.2879255634912523e+04, 6.2888228606025921e+04, 6.2897201703922212e+04, 6.2906174928585315e+04, 6.2915148279999150e+04, 6.2924121758147667e+04, 6.2933095363014814e+04, 6.2942069094584527e+04, 6.2951042952840740e+04, 6.2960016937767439e+04, 6.2968991049348551e+04, 6.2977965287568048e+04, 6.2986939652409892e+04, 6.2995914143858055e+04, 6.3004888761896509e+04, 6.3013863506509231e+04, 6.3022838377680200e+04, 6.3031813375393409e+04, 6.3040788499632836e+04, 6.3049763750382481e+04, 6.3058739127626330e+04, 6.3067714631348397e+04, 6.3076690261532684e+04, 6.3085666018163203e+04, 6.3094641901223964e+04, 6.3103617910698980e+04, 6.3112594046572282e+04, 6.3121570308827897e+04, 6.3130546697449849e+04, 6.3139523212422158e+04, 6.3148499853728885e+04, 6.3157476621354057e+04, 6.3166453515281719e+04, 6.3175430535495922e+04, 6.3184407681980738e+04, 6.3193384954720197e+04, 6.3202362353698387e+04, 6.3211339878899351e+04, 6.3220317530307169e+04, 6.3229295307905908e+04, 6.3238273211679669e+04, 6.3247251241612503e+04, 6.3256229397688512e+04, 6.3265207679891791e+04, 6.3274186088206414e+04, 6.3283164622616503e+04, 6.3292143283106132e+04, 6.3301122069659439e+04, 6.3310100982260512e+04, 6.3319080020893474e+04, 6.3328059185542435e+04, 6.3337038476191527e+04, 6.3346017892824872e+04, 6.3354997435426601e+04, 6.3363977103980847e+04, 6.3372956898471741e+04, 6.3381936818883427e+04, 6.3390916865200059e+04, 6.3399897037405775e+04, 6.3408877335484744e+04, 6.3417857759421124e+04, 6.3426838309199062e+04, 6.3435818984802740e+04, 6.3444799786216317e+04, 6.3453780713423963e+04, 6.3462761766409858e+04, 6.3471742945158199e+04, 6.3480724249653147e+04, 6.3489705679878913e+04, 6.3498687235819685e+04, 6.3507668917459661e+04, 6.3516650724783038e+04, 6.3525632657774018e+04, 6.3534614716416821e+04, 6.3543596900695658e+04, 6.3552579210594748e+04, 6.3561561646098300e+04, 6.3570544207190556e+04, 6.3579526893855749e+04, 6.3588509706078090e+04, 6.3597492643841826e+04, 6.3606475707131212e+04, 6.3615458895930475e+04, 6.3624442210223882e+04, 6.3633425649995661e+04, 6.3642409215230087e+04, 6.3651392905911423e+04, 6.3660376722023917e+04, 6.3669360663551859e+04, 6.3678344730479505e+04, 6.3687328922791152e+04, 6.3696313240471063e+04, 6.3705297683503530e+04, 6.3714282251872835e+04, 6.3723266945563279e+04, 6.3732251764559158e+04, 6.3741236708844757e+04, 6.3750221778404404e+04, 6.3759206973222397e+04, 6.3768192293283049e+04, 6.3777177738570674e+04, 6.3786163309069583e+04, 6.3795149004764113e+04, 6.3804134825638590e+04, 6.3813120771677350e+04, 6.3822106842864720e+04, 6.3831093039185042e+04, 6.3840079360622673e+04, 6.3849065807161940e+04, 6.3858052378787215e+04, 6.3867039075482833e+04, 6.3876025897233165e+04, 6.3885012844022574e+04, 6.3893999915835419e+04, 6.3902987112656083e+04, 6.3911974434468939e+04, 6.3920961881258350e+04, 6.3929949453008725e+04, 6.3938937149704427e+04, 6.3947924971329856e+04, 6.3956912917869406e+04, 6.3965900989307484e+04, 6.3974889185628475e+04, 6.3983877506816796e+04, 6.3992865952856861e+04, 6.4001854523733076e+04, 6.4010843219429866e+04, 6.4019832039931644e+04, 6.4028820985222832e+04, 6.4037810055287875e+04, 6.4046799250111209e+04, 6.4055788569677257e+04, 6.4064778013970463e+04, 6.4073767582975262e+04, 6.4082757276676122e+04, 6.4091747095057493e+04, 6.4100737038103820e+04, 6.4109727105799575e+04, 6.4118717298129210e+04, 6.4127707615077212e+04, 6.4136698056628033e+04, 6.4145688622766167e+04, 6.4154679313476088e+04, 6.4163670128742277e+04, 6.4172661068549220e+04, 6.4181652132881405e+04, 6.4190643321723343e+04, 6.4199634635059512e+04, 6.4208626072874431e+04, 6.4217617635152601e+04, 6.4226609321878532e+04, 6.4235601133036740e+04, 6.4244593068611750e+04, 6.4253585128588078e+04, 6.4262577312950241e+04, 6.4271569621682800e+04, 6.4280562054770257e+04, 6.4289554612197164e+04, 6.4298547293948046e+04, 6.4307540100007471e+04, 6.4316533030359984e+04, 6.4325526084990124e+04, 6.4334519263882466e+04, 6.4343512567021557e+04, 6.4352505994391970e+04, 6.4361499545978266e+04, 6.4370493221765028e+04, 6.4379487021736822e+04, 6.4388480945878240e+04, 6.4397474994173856e+04, 6.4406469166608244e+04, 6.4415463463166023e+04, 6.4424457883831768e+04, 6.4433452428590099e+04, 6.4442447097425596e+04, 6.4451441890322872e+04, 6.4460436807266560e+04, 6.4469431848241242e+04, 6.4478427013231551e+04, 6.4487422302222119e+04, 6.4496417715197553e+04, 6.4505413252142484e+04, 6.4514408913041567e+04, 6.4523404697879414e+04, 6.4532400606640680e+04, 6.4541396639310005e+04, 6.4550392795872038e+04, 6.4559389076311440e+04, 6.4568385480612844e+04, 6.4577382008760935e+04, 6.4586378660740367e+04, 6.4595375436535811e+04, 6.4604372336131928e+04, 6.4613369359513410e+04, 6.4622366506664926e+04, 6.4631363777571161e+04, 6.4640361172216799e+04, 6.4649358690586530e+04, 6.4658356332665055e+04, 6.4667354098437077e+04, 6.4676351987887268e+04, 6.4685350001000370e+04, 6.4694348137761066e+04, 6.4703346398154084e+04, 6.4712344782164131e+04, 6.4721343289775941e+04, 6.4730341920974228e+04, 6.4739340675743726e+04, 6.4748339554069156e+04, 6.4757338555935268e+04, 6.4766337681326797e+04, 6.4775336930228470e+04, 6.4784336302625074e+04, 6.4793335798501321e+04, 6.4802335417841983e+04, 6.4811335160631825e+04, 6.4820335026855588e+04, 6.4829335016498044e+04, 6.4838335129543972e+04, 6.4847335365978150e+04, 6.4856335725785335e+04, 6.4865336208950328e+04, 6.4874336815457893e+04, 6.4883337545292830e+04, 6.4892338398439941e+04, 6.4901339374884010e+04, 6.4910340474609853e+04, 6.4919341697602234e+04, 6.4928343043846005e+04, 6.4937344513325959e+04, 6.4946346106026896e+04, 6.4955347821933661e+04, 6.4964349661031054e+04, 6.4973351623303919e+04, 6.4982353708737079e+04, 6.4991355917315363e+04, 6.5000358249023615e+04, 6.5009360703846665e+04, 6.5018363281769365e+04, 6.5027365982776566e+04, 6.5036368806853112e+04, 6.5045371753983854e+04, 6.5054374824153674e+04, 6.5063378017347422e+04, 6.5072381333549965e+04, 6.5081384772746169e+04, 6.5090388334920914e+04, 6.5099392020059080e+04, 6.5108395828145549e+04, 6.5117399759165208e+04, 6.5126403813102930e+04, 6.5135407989943633e+04, 6.5144412289672196e+04, 6.5153416712273523e+04, 6.5162421257732516e+04, 6.5171425926034091e+04, 6.5180430717163159e+04, 6.5189435631104629e+04, 6.5198440667843410e+04, 6.5207445827364449e+04, 6.5216451109652655e+04, 6.5225456514692960e+04, 6.5234462042470303e+04, 6.5243467692969622e+04, 6.5252473466175856e+04, 6.5261479362073944e+04, 6.5270485380648854e+04, 6.5279491521885517e+04, 6.5288497785768894e+04, 6.5297504172283960e+04, 6.5306510681415646e+04, 6.5315517313148957e+04, 6.5324524067468832e+04, 6.5333530944360260e+04, 6.5342537943808224e+04, 6.5351545065797683e+04, 6.5360552310313651e+04, 6.5369559677341109e+04, 6.5378567166865039e+04, 6.5387574778870432e+04, 6.5396582513342320e+04, 6.5405590370265680e+04, 6.5414598349625521e+04, 6.5423606451406857e+04, 6.5432614675594712e+04, 6.5441623022174092e+04, 6.5450631491130043e+04, 6.5459640082447557e+04, 6.5468648796111687e+04, 6.5477657632107454e+04, 6.5486666590419896e+04, 6.5495675671034056e+04, 6.5504684873934981e+04, 6.5513694199107726e+04, 6.5522703646537317e+04, 6.5531713216208824e+04, 6.5540722908107316e+04, 6.5549732722217843e+04, 6.5558742658525473e+04, 6.5567752717015261e+04, 6.5576762897672306e+04, 6.5585773200481679e+04, 6.5594783625428463e+04, 6.5603794172497728e+04, 6.5612804841674588e+04, 6.5621815632944097e+04, 6.5630826546291370e+04, 6.5639837581701533e+04, 6.5648848739159628e+04, 6.5657860018650812e+04, 6.5666871420160169e+04, 6.5675882943672827e+04, 6.5684894589173884e+04, 6.5693906356648469e+04, 6.5702918246081710e+04, 6.5711930257458764e+04, 6.5720942390764714e+04, 6.5729954645984719e+04, 6.5738967023103905e+04, 6.5747979522107446e+04, 6.5756992142980453e+04, 6.5766004885708084e+04, 6.5775017750275481e+04, 6.5784030736667832e+04, 6.5793043844870277e+04, 6.5802057074867989e+04, 6.5811070426646140e+04, 6.5820083900189871e+04, 6.5829097495484413e+04, 6.5838111212514879e+04, 6.5847125051266499e+04, 6.5856139011724430e+04, 6.5865153093873872e+04, 6.5874167297700013e+04, 6.5883181623188080e+04, 6.5892196070323218e+04, 6.5901210639090656e+04, 6.5910225329475637e+04, 6.5919240141463306e+04, 6.5928255075038964e+04, 6.5937270130187739e+04, 6.5946285306894904e+04, 6.5955300605145676e+04, 6.5964316024925283e+04, 6.5973331566218956e+04, 6.5982347229011939e+04, 6.5991363013289432e+04, 6.6000378919036724e+04, 6.6009394946239059e+04, 6.6018411094881681e+04, 6.6027427364949835e+04, 6.6036443756428780e+04, 6.6045460269303789e+04, 6.6054476903560106e+04, 6.6063493659182990e+04, 6.6072510536157759e+04, 6.6081527534469671e+04, 6.6090544654103986e+04, 6.6099561895045976e+04, 6.6108579257280973e+04, 6.6117596740794237e+04, 6.6126614345571070e+04, 6.6135632071596774e+04, 6.6144649918856652e+04, 6.6153667887335978e+04, 6.6162685977020083e+04, 6.6171704187894269e+04, 6.6180722519943854e+04, 6.6189740973154170e+04, 6.6198759547510519e+04, 6.6207778242998233e+04, 6.6216797059602657e+04, 6.6225815997309095e+04, 6.6234835056102922e+04, 6.6243854235969440e+04, 6.6252873536894011e+04, 6.6261892958861965e+04, 6.6270912501858649e+04, 6.6279932165869453e+04, 6.6288951950879709e+04, 6.6297971856874778e+04, 6.6306991883840048e+04, 6.6316012031760867e+04, 6.6325032300622581e+04, 6.6334052690410594e+04, 6.6343073201110281e+04, 6.6352093832707018e+04, 6.6361114585186224e+04, 6.6370135458533230e+04, 6.6379156452733485e+04, 6.6388177567772364e+04, 6.6397198803635241e+04, 6.6406220160307552e+04, 6.6415241637774670e+04, 6.6424263236022060e+04, 6.6433284955035066e+04, 6.6442306794799180e+04, 6.6451328755299779e+04, 6.6460350836522295e+04, 6.6469373038452162e+04, 6.6478395361074785e+04, 6.6487417804375626e+04, 6.6496440368340118e+04, 6.6505463052953724e+04, 6.6514485858201850e+04, 6.6523508784069956e+04, 6.6532531830543521e+04, 6.6541554997607964e+04, 6.6550578285248776e+04, 6.6559601693451390e+04, 6.6568625222201299e+04, 6.6577648871483965e+04, 6.6586672641284837e+04, 6.6595696531589434e+04, 6.6604720542383220e+04, 6.6613744673651687e+04, 6.6622768925380296e+04, 6.6631793297554555e+04, 6.6640817790159970e+04, 6.6649842403181989e+04, 6.6658867136606175e+04, 6.6667891990418037e+04, 6.6676916964603035e+04, 6.6685942059146706e+04, 6.6694967274034556e+04, 6.6703992609252105e+04, 6.6713018064784876e+04, 6.6722043640618416e+04, 6.6731069336738234e+04, 6.6740095153129863e+04, 6.6749121089778841e+04, 6.6758147146670715e+04, 6.6767173323791023e+04, 6.6776199621125299e+04, 6.6785226038659122e+04, 6.6794252576378014e+04, 6.6803279234267553e+04, 6.6812306012313304e+04, 6.6821332910500816e+04, 6.6830359928815669e+04, 6.6839387067243428e+04, 6.6848414325769671e+04, 6.6857441704379962e+04, 6.6866469203059896e+04, 6.6875496821795066e+04, 6.6884524560571052e+04, 6.6893552419373431e+04, 6.6902580398187813e+04, 6.6911608496999790e+04, 6.6920636715794986e+04, 6.6929665054558980e+04, 6.6938693513277394e+04, 6.6947722091935837e+04, 6.6956750790519916e+04, 6.6965779609015270e+04, 6.6974808547407505e+04, 6.6983837605682274e+04, 6.6992866783825171e+04, 6.7001896081821847e+04, 6.7010925499657926e+04, 6.7019955037319087e+04, 6.7028984694790939e+04, 6.7038014472059120e+04, 6.7047044369109324e+04, 6.7056074385927161e+04, 6.7065104522498325e+04, 6.7074134778808439e+04, 6.7083165154843184e+04, 6.7092195650588241e+04, 6.7101226266029262e+04, 6.7110257001151942e+04, 6.7119287855941948e+04, 6.7128318830384960e+04, 6.7137349924466660e+04, 6.7146381138172728e+04, 6.7155412471488875e+04, 6.7164443924400795e+04, 6.7173475496894171e+04, 6.7182507188954711e+04, 6.7191539000568140e+04, 6.7200570931720154e+04, 6.7209602982396435e+04, 6.7218635152582749e+04, 6.7227667442264777e+04, 6.7236699851428260e+04, 6.7245732380058922e+04, 6.7254765028142487e+04, 6.7263797795664694e+04, 6.7272830682611268e+04, 6.7281863688967962e+04, 6.7290896814720501e+04, 6.7299930059854654e+04, 6.7308963424356160e+04, 6.7317996908210771e+04, 6.7327030511404257e+04, 6.7336064233922341e+04, 6.7345098075750822e+04, 6.7354132036875468e+04, 6.7363166117282017e+04, 6.7372200316956296e+04, 6.7381234635884015e+04, 6.7390269074051015e+04, 6.7399303631443050e+04, 6.7408338308045888e+04, 6.7417373103845355e+04, 6.7426408018827235e+04, 6.7435443052977309e+04, 6.7444478206281390e+04, 6.7453513478725261e+04, 6.7462548870294777e+04, 6.7471584380975706e+04, 6.7480620010753890e+04, 6.7489655759615140e+04, 6.7498691627545268e+04, 6.7507727614530086e+04, 6.7516763720555449e+04, 6.7525799945607170e+04, 6.7534836289671090e+04, 6.7543872752733048e+04, 6.7552909334778888e+04, 6.7561946035794448e+04, 6.7570982855765586e+04, 6.7580019794678155e+04, 6.7589056852517984e+04, 6.7598094029270957e+04, 6.7607131324922928e+04, 6.7616168739459754e+04, 6.7625206272867305e+04, 6.7634243925131464e+04, 6.7643281696238075e+04, 6.7652319586173078e+04, 6.7661357594922272e+04, 6.7670395722471600e+04, 6.7679433968806945e+04, 6.7688472333914178e+04, 6.7697510817779199e+04, 6.7706549420387906e+04, 6.7715588141726243e+04, 6.7724626981780049e+04, 6.7733665940535269e+04, 6.7742705017977816e+04, 6.7751744214093589e+04, 6.7760783528868516e+04, 6.7769822962288526e+04, 6.7778862514339504e+04, 6.7787902185007435e+04, 6.7796941974278205e+04, 6.7805981882137785e+04, 6.7815021908572104e+04, 6.7824062053567090e+04, 6.7833102317108685e+04, 6.7842142699182848e+04, 6.7851183199775536e+04, 6.7860223818872677e+04, 6.7869264556460272e+04, 6.7878305412524249e+04, 6.7887346387050595e+04, 6.7896387480025282e+04, 6.7905428691434237e+04, 6.7914470021263478e+04, 6.7923511469498961e+04, 6.7932553036126686e+04, 6.7941594721132642e+04, 6.7950636524502785e+04, 6.7959678446223130e+04, 6.7968720486279679e+04, 6.7977762644658418e+04, 6.7986804921345349e+04, 6.7995847316326486e+04, 6.8004889829587817e+04, 6.8013932461115357e+04, 6.8022975210895151e+04, 6.8032018078913185e+04, 6.8041061065155489e+04, 6.8050104169608094e+04, 6.8059147392257029e+04, 6.8068190733088297e+04, 6.8077234192087970e+04, 6.8086277769242064e+04, 6.8095321464536639e+04, 6.8104365277957739e+04, 6.8113409209491379e+04, 6.8122453259123620e+04, 6.8131497426840564e+04, 6.8140541712628212e+04, 6.8149586116472652e+04, 6.8158630638359915e+04, 6.8167675278276118e+04, 6.8176720036207305e+04, 6.8185764912139552e+04, 6.8194809906058945e+04, 6.8203855017951559e+04, 6.8212900247803482e+04, 6.8221945595600788e+04, 6.8230991061329565e+04, 6.8240036644975917e+04, 6.8249082346525960e+04, 6.8258128165965769e+04, 6.8267174103281446e+04, 6.8276220158459124e+04, 6.8285266331484876e+04, 6.8294312622344849e+04, 6.8303359031025131e+04, 6.8312405557511855e+04, 6.8321452201791166e+04, 6.8330498963849168e+04, 6.8339545843671993e+04, 6.8348592841245772e+04, 6.8357639956556639e+04, 6.8366687189590753e+04, 6.8375734540334233e+04, 6.8384782008773225e+04, 6.8393829594893890e+04, 6.8402877298682390e+04, 6.8411925120124884e+04, 6.8420973059207492e+04, 6.8430021115916417e+04, 6.8439069290237807e+04, 6.8448117582157822e+04, 6.8457165991662667e+04, 6.8466214518738474e+04, 6.8475263163371448e+04, 6.8484311925547765e+04, 6.8493360805253615e+04, 6.8502409802475173e+04, 6.8511458917198644e+04, 6.8520508149410234e+04, 6.8529557499096103e+04, 6.8538606966242500e+04, 6.8547656550835600e+04, 6.8556706252861622e+04, 6.8565756072306758e+04, 6.8574806009157241e+04, 6.8583856063399289e+04, 6.8592906235019109e+04, 6.8601956524002933e+04, 6.8611006930336982e+04, 6.8620057454007532e+04, 6.8629108095000745e+04, 6.8638158853302899e+04, 6.8647209728900212e+04, 6.8656260721778948e+04, 6.8665311831925355e+04, 6.8674363059325682e+04, 6.8683414403966162e+04, 6.8692465865833074e+04, 6.8701517444912664e+04, 6.8710569141191212e+04, 6.8719620954654951e+04, 6.8728672885290187e+04, 6.8737724933083169e+04, 6.8746777098020175e+04, 6.8755829380087511e+04, 6.8764881779271425e+04, 6.8773934295558196e+04, 6.8782986928934130e+04, 6.8792039679385518e+04, 6.8801092546898668e+04, 6.8810145531459872e+04, 6.8819198633055421e+04, 6.8828251851671608e+04, 6.8837305187294769e+04, 6.8846358639911210e+04, 6.8855412209507238e+04, 6.8864465896069174e+04, 6.8873519699583310e+04, 6.8882573620036026e+04, 6.8891627657413599e+04, 6.8900681811702379e+04, 6.8909736082888703e+04, 6.8918790470958906e+04, 6.8927844975899337e+04, 6.8936899597696320e+04, 6.8945954336336188e+04, 6.8955009191805322e+04, 6.8964064164090058e+04, 6.8973119253176759e+04, 6.8982174459051806e+04, 6.8991229781701506e+04, 6.9000285221112266e+04, 6.9009340777270438e+04, 6.9018396450162414e+04, 6.9027452239774546e+04, 6.9036508146093198e+04, 6.9045564169104793e+04, 6.9054620308795682e+04, 6.9063676565152287e+04, 6.9072732938160960e+04, 6.9081789427808108e+04, 6.9090846034080139e+04, 6.9099902756963449e+04, 6.9108959596444431e+04, 6.9118016552509507e+04, 6.9127073625145087e+04, 6.9136130814337565e+04, 6.9145188120073377e+04, 6.9154245542338933e+04, 6.9163303081120655e+04, 6.9172360736404968e+04, 6.9181418508178293e+04, 6.9190476396427053e+04, 6.9199534401137746e+04, 6.9208592522296734e+04, 6.9217650759890501e+04, 6.9226709113905483e+04, 6.9235767584328118e+04, 6.9244826171144858e+04, 6.9253884874342169e+04, 6.9262943693906505e+04, 6.9272002629824303e+04, 6.9281061682082072e+04, 6.9290120850666251e+04, 6.9299180135563307e+04, 6.9308239536759706e+04, 6.9317299054241957e+04, 6.9326358687996515e+04, 6.9335418438009860e+04, 6.9344478304268472e+04, 6.9353538286758863e+04, 6.9362598385467514e+04, 6.9371658600380906e+04, 6.9380718931485550e+04, 6.9389779378767969e+04, 6.9398839942214618e+04, 6.9407900621812034e+04, 6.9416961417546729e+04, 6.9426022329405227e+04, 6.9435083357374024e+04, 6.9444144501439630e+04, 6.9453205761588586e+04, 6.9462267137807430e+04, 6.9471328630082658e+04, 6.9480390238400825e+04, 6.9489451962748470e+04, 6.9498513803112131e+04, 6.9507575759478350e+04, 6.9516637831833650e+04, 6.9525700020164601e+04, 6.9534762324457755e+04, 6.9543824744699654e+04, 6.9552887280876879e+04, 6.9561949932975942e+04, 6.9571012700983469e+04, 6.9580075584885970e+04, 6.9589138584670058e+04, 6.9598201700322272e+04, 6.9607264931829195e+04, 6.9616328279177440e+04, 6.9625391742353546e+04, 6.9634455321344125e+04, 6.9643519016135760e+04, 6.9652582826715035e+04, 6.9661646753068562e+04, 6.9670710795182924e+04, 6.9679774953044718e+04, 6.9688839226640557e+04, 6.9697903615957053e+04, 6.9706968120980804e+04, 6.9716032741698451e+04, 6.9725097478096563e+04, 6.9734162330161780e+04, 6.9743227297880760e+04, 6.9752292381240070e+04, 6.9761357580226395e+04, 6.9770422894826304e+04, 6.9779488325026497e+04, 6.9788553870813572e+04, 6.9797619532174198e+04, 6.9806685309094973e+04, 6.9815751201562569e+04, 6.9824817209563655e+04, 6.9833883333084872e+04, 6.9842949572112877e+04, 6.9852015926634340e+04, 6.9861082396635902e+04, 6.9870148982104220e+04, 6.9879215683026021e+04, 6.9888282499387904e+04, 6.9897349431176583e+04, 6.9906416478378713e+04, 6.9915483640981009e+04, 6.9924550918970141e+04, 6.9933618312332794e+04, 6.9942685821055667e+04, 6.9951753445125432e+04, 6.9960821184528817e+04, 6.9969889039252521e+04, 6.9978957009283185e+04, 6.9988025094607583e+04, 6.9997093295212399e+04, 7.0006161611084332e+04, 7.0015230042210140e+04, 7.0024298588576479e+04, 7.0033367250170122e+04, 7.0042436026977783e+04, 7.0051504918986189e+04, 7.0060573926182042e+04, 7.0069643048552098e+04, 7.0078712286083100e+04, 7.0087781638761779e+04, 7.0096851106574890e+04, 7.0105920689509148e+04, 7.0114990387551312e+04, 7.0124060200688153e+04, 7.0133130128906414e+04, 7.0142200172192854e+04, 7.0151270330534244e+04, 7.0160340603917313e+04, 7.0169410992328863e+04, 7.0178481495755652e+04, 7.0187552114184451e+04, 7.0196622847602048e+04, 7.0205693695995200e+04, 7.0214764659350709e+04, 7.0223835737655361e+04, 7.0232906930895930e+04, 7.0241978239059201e+04, 7.0251049662132005e+04, 7.0260121200101115e+04, 7.0269192852953318e+04, 7.0278264620675443e+04, 7.0287336503254279e+04, 7.0296408500676611e+04, 7.0305480612929314e+04, 7.0314552839999174e+04, 7.0323625181872994e+04, 7.0332697638537589e+04, 7.0341770209979819e+04, 7.0350842896186499e+04, 7.0359915697144446e+04, 7.0368988612840520e+04, 7.0378061643261535e+04, 7.0387134788394338e+04, 7.0396208048225773e+04, 7.0405281422742657e+04, 7.0414354911931892e+04, 7.0423428515780295e+04, 7.0432502234274725e+04, 7.0441576067402057e+04, 7.0450650015149135e+04, 7.0459724077502804e+04, 7.0468798254449968e+04, 7.0477872545977487e+04, 7.0486946952072220e+04, 7.0496021472721055e+04, 7.0505096107910882e+04, 7.0514170857628545e+04, 7.0523245721860978e+04, 7.0532320700595024e+04, 7.0541395793817588e+04, 7.0550471001515572e+04, 7.0559546323675866e+04, 7.0568621760285387e+04, 7.0577697311331009e+04, 7.0586772976799664e+04, 7.0595848756678242e+04, 7.0604924650953675e+04, 7.0614000659612851e+04, 7.0623076782642718e+04, 7.0632153020030164e+04, 7.0641229371762136e+04, 7.0650305837825566e+04, 7.0659382418207359e+04, 7.0668459112894445e+04, 7.0677535921873787e+04, 7.0686612845132331e+04, 7.0695689882656981e+04, 7.0704767034434713e+04, 7.0713844300452460e+04, 7.0722921680697182e+04, 7.0731999175155841e+04, 7.0741076783815355e+04, 7.0750154506662730e+04, 7.0759232343684896e+04, 7.0768310294868803e+04, 7.0777388360201468e+04, 7.0786466539669840e+04, 7.0795544833260879e+04, 7.0804623240961606e+04, 7.0813701762758967e+04, 7.0822780398639938e+04, 7.0831859148591539e+04, 7.0840938012600731e+04, 7.0850016990654505e+04, 7.0859096082739867e+04, 7.0868175288843806e+04, 7.0877254608953343e+04, 7.0886334043055467e+04, 7.0895413591137185e+04, 7.0904493253185501e+04, 7.0913573029187450e+04, 7.0922652919130051e+04, 7.0931732923000294e+04, 7.0940813040785215e+04, 7.0949893272471832e+04, 7.0958973618047210e+04, 7.0968054077498338e+04, 7.0977134650812251e+04, 7.0986215337975998e+04, 7.0995296138976613e+04, 7.1004377053801145e+04, 7.1013458082436642e+04, 7.1022539224870168e+04, 7.1031620481088743e+04, 7.1040701851079415e+04, 7.1049783334829248e+04, 7.1058864932325348e+04, 7.1067946643554707e+04, 7.1077028468504461e+04, 7.1086110407161628e+04, 7.1095192459513288e+04, 7.1104274625546546e+04, 7.1113356905248467e+04, 7.1122439298606114e+04, 7.1131521805606593e+04, 7.1140604426236969e+04, 7.1149687160484347e+04, 7.1158770008335792e+04, 7.1167852969778425e+04, 7.1176936044799368e+04, 7.1186019233385668e+04, 7.1195102535524478e+04, 7.1204185951202890e+04, 7.1213269480407995e+04, 7.1222353123126930e+04, 7.1231436879346787e+04, 7.1240520749054704e+04, 7.1249604732237800e+04, 7.1258688828883198e+04, 7.1267773038978048e+04, 7.1276857362509443e+04, 7.1285941799464519e+04, 7.1295026349830441e+04, 7.1304111013594331e+04, 7.1313195790743324e+04, 7.1322280681264572e+04, 7.1331365685145240e+04, 7.1340450802372477e+04, 7.1349536032933407e+04, 7.1358621376815208e+04, 7.1367706834005032e+04, 7.1376792404490057e+04, 7.1385878088257436e+04, 7.1394963885294332e+04, 7.1404049795587911e+04, 7.1413135819125382e+04, 7.1422221955893910e+04, 7.1431308205880647e+04, 7.1440394569072800e+04, 7.1449481045457564e+04, 7.1458567635022104e+04, 7.1467654337753629e+04, 7.1476741153639305e+04, 7.1485828082666369e+04, 7.1494915124822001e+04, 7.1504002280093409e+04, 7.1513089548467804e+04, 7.1522176929932379e+04, 7.1531264424474357e+04, 7.1540352032080962e+04, 7.1549439752739374e+04, 7.1558527586436874e+04, 7.1567615533160613e+04, 7.1576703592897888e+04, 7.1585791765635891e+04, 7.1594880051361863e+04, 7.1603968450063025e+04, 7.1613056961726645e+04, 7.1622145586339931e+04, 7.1631234323890138e+04, 7.1640323174364530e+04, 7.1649412137750332e+04, 7.1658501214034841e+04, 7.1667590403205249e+04, 7.1676679705248825e+04, 7.1685769120152865e+04, 7.1694858647904635e+04, 7.1703948288491374e+04, 7.1713038041900363e+04, 7.1722127908118870e+04, 7.1731217887134160e+04, 7.1740307978933561e+04, 7.1749398183504294e+04, 7.1758488500833671e+04, 7.1767578930908960e+04, 7.1776669473717484e+04, 7.1785760129246511e+04, 7.1794850897483353e+04, 7.1803941778415290e+04, 7.1813032772029648e+04, 7.1822123878313680e+04, 7.1831215097254753e+04, 7.1840306428840137e+04, 7.1849397873057169e+04, 7.1858489429893147e+04, 7.1867581099335395e+04, 7.1876672881371240e+04, 7.1885764775988020e+04, 7.1894856783173018e+04, 7.1903948902913602e+04, 7.1913041135197112e+04, 7.1922133480010845e+04, 7.1931225937342155e+04, 7.1940318507178396e+04, 7.1949411189506907e+04, 7.1958503984315030e+04, 7.1967596891590118e+04, 7.1976689911319510e+04, 7.1985783043490577e+04, 7.1994876288090672e+04, 7.2003969645107165e+04, 7.2013063114527409e+04, 7.2022156696338789e+04, 7.2031250390528643e+04, 7.2040344197084356e+04, 7.2049438115993311e+04, 7.2058532147242891e+04, 7.2067626290820466e+04, 7.2076720546713419e+04, 7.2085814914909133e+04, 7.2094909395395021e+04, 7.2104003988158453e+04, 7.2113098693186825e+04, 7.2122193510467550e+04, 7.2131288439987999e+04, 7.2140383481735596e+04, 7.2149478635697771e+04, 7.2158573901861891e+04, 7.2167669280215399e+04, 7.2176764770745707e+04, 7.2185860373440213e+04, 7.2194956088286359e+04, 7.2204051915271542e+04, 7.2213147854383220e+04, 7.2222243905608775e+04, 7.2231340068935693e+04, 7.2240436344351387e+04, 7.2249532731843297e+04, 7.2258629231398852e+04, 7.2267725843005494e+04, 7.2276822566650677e+04, 7.2285919402321859e+04, 7.2295016350006466e+04, 7.2304113409692000e+04, 7.2313210581365856e+04, 7.2322307865015551e+04, 7.2331405260628526e+04, 7.2340502768192222e+04, 7.2349600387694139e+04, 7.2358698119121749e+04, 7.2367795962462493e+04, 7.2376893917703885e+04, 7.2385991984833381e+04, 7.2395090163838468e+04, 7.2404188454706644e+04, 7.2413286857425366e+04, 7.2422385371982164e+04, 7.2431483998364507e+04, 7.2440582736559911e+04, 7.2449681586555846e+04, 7.2458780548339841e+04, 7.2467879621899396e+04, 7.2476978807221996e+04, 7.2486078104295171e+04, 7.2495177513106464e+04, 7.2504277033643331e+04, 7.2513376665893331e+04, 7.2522476409843977e+04, 7.2531576265482770e+04, 7.2540676232797268e+04, 7.2549776311774971e+04, 7.2558876502403451e+04, 7.2567976804670223e+04, 7.2577077218562816e+04, 7.2586177744068787e+04, 7.2595278381175667e+04, 7.2604379129871013e+04, 7.2613479990142368e+04, 7.2622580961977277e+04, 7.2631682045363312e+04, 7.2640783240288030e+04, 7.2649884546738991e+04, 7.2658985964703737e+04, 7.2668087494169857e+04, 7.2677189135124907e+04, 7.2686290887556475e+04, 7.2695392751452120e+04, 7.2704494726799399e+04, 7.2713596813585929e+04, 7.2722699011799283e+04, 7.2731801321427061e+04, 7.2740903742456809e+04, 7.2750006274876141e+04, 7.2759108918672675e+04, 7.2768211673833954e+04, 7.2777314540347637e+04, 7.2786417518201270e+04, 7.2795520607382496e+04, 7.2804623807878917e+04, 7.2813727119678137e+04, 7.2822830542767770e+04, 7.2831934077135418e+04, 7.2841037722768713e+04, 7.2850141479655300e+04, 7.2859245347782751e+04, 7.2868349327138742e+04, 7.2877453417710858e+04, 7.2886557619486775e+04, 7.2895661932454095e+04, 7.2904766356600478e+04, 7.2913870891913553e+04, 7.2922975538380968e+04, 7.2932080295990381e+04, 7.2941185164729395e+04, 7.2950290144585713e+04, 7.2959395235546981e+04, 7.2968500437600829e+04, 7.2977605750734918e+04, 7.2986711174936951e+04, 7.2995816710194558e+04, 7.3004922356495430e+04, 7.3014028113827211e+04, 7.3023133982177576e+04, 7.3032239961534215e+04, 7.3041346051884830e+04, 7.3050452253217052e+04, 7.3059558565518615e+04, 7.3068664988777164e+04, 7.3077771522980416e+04, 7.3086878168116047e+04, 7.3095984924171760e+04, 7.3105091791135259e+04, 7.3114198768994233e+04, 7.3123305857736414e+04, 7.3132413057349462e+04, 7.3141520367821140e+04, 7.3150627789139107e+04, 7.3159735321291111e+04, 7.3168842964264841e+04, 7.3177950718048043e+04, 7.3187058582628437e+04, 7.3196166557993725e+04, 7.3205274644131670e+04, 7.3214382841029961e+04, 7.3223491148676374e+04, 7.3232599567058613e+04, 7.3241708096164439e+04, 7.3250816735981585e+04, 7.3259925486497799e+04, 7.3269034347700843e+04, 7.3278143319578419e+04, 7.3287252402118320e+04, 7.3296361595308277e+04, 7.3305470899136068e+04, 7.3314580313589453e+04, 7.3323689838656195e+04, 7.3332799474324056e+04, 7.3341909220580797e+04, 7.3351019077414196e+04, 7.3360129044812027e+04, 7.3369239122762068e+04, 7.3378349311252110e+04, 7.3387459610269900e+04, 7.3396570019803257e+04, 7.3405680539839945e+04, 7.3414791170367782e+04, 7.3423901911374545e+04, 7.3433012762848011e+04, 7.3442123724776015e+04, 7.3451234797146331e+04, 7.3460345979946767e+04, 7.3469457273165172e+04, 7.3478568676789277e+04, 7.3487680190806946e+04, 7.3496791815206001e+04, 7.3505903549974231e+04, 7.3515015395099457e+04, 7.3524127350569514e+04, 7.3533239416372220e+04, 7.3542351592495412e+04, 7.3551463878926908e+04, 7.3560576275654545e+04, 7.3569688782666184e+04, 7.3578801399949632e+04, 7.3587914127492753e+04, 7.3597026965283367e+04, 7.3606139913309336e+04, 7.3615252971558526e+04, 7.3624366140018741e+04, 7.3633479418677860e+04, 7.3642592807523775e+04, 7.3651706306544293e+04, 7.3660819915727305e+04, 7.3669933635060675e+04, 7.3679047464532283e+04, 7.3688161404129976e+04, 7.3697275453841619e+04, 7.3706389613655134e+04, 7.3715503883558355e+04, 7.3724618263539189e+04, 7.3733732753585515e+04, 7.3742847353685211e+04, 7.3751962063826169e+04, 7.3761076883996284e+04, 7.3770191814183447e+04, 7.3779306854375565e+04, 7.3788422004560533e+04, 7.3797537264726270e+04, 7.3806652634860657e+04, 7.3815768114951614e+04, 7.3824883704987064e+04, 7.3833999404954884e+04, 7.3843115214842997e+04, 7.3852231134639369e+04, 7.3861347164331877e+04, 7.3870463303908444e+04, 7.3879579553357034e+04, 7.3888695912665527e+04, 7.3897812381821888e+04, 7.3906928960814039e+04, 7.3916045649629930e+04, 7.3925162448257499e+04, 7.3934279356684696e+04, 7.3943396374899443e+04, 7.3952513502889691e+04, 7.3961630740643392e+04, 7.3970748088148524e+04, 7.3979865545393026e+04, 7.3988983112364862e+04, 7.3998100789051969e+04, 7.4007218575442326e+04, 7.4016336471523900e+04, 7.4025454477284686e+04, 7.4034572592712619e+04, 7.4043690817795679e+04, 7.4052809152521877e+04, 7.4061927596879148e+04, 7.4071046150855487e+04, 7.4080164814438889e+04, 7.4089283587617349e+04, 7.4098402470378831e+04, 7.4107521462711331e+04, 7.4116640564602887e+04, 7.4125759776041450e+04, 7.4134879097015044e+04, 7.4143998527511649e+04, 7.4153118067519303e+04, 7.4162237717025986e+04, 7.4171357476019723e+04, 7.4180477344488536e+04, 7.4189597322420435e+04, 7.4198717409803430e+04, 7.4207837606625544e+04, 7.4216957912874816e+04, 7.4226078328539268e+04, 7.4235198853606911e+04, 7.4244319488065797e+04, 7.4253440231903980e+04, 7.4262561085109453e+04, 7.4271682047670285e+04, 7.4280803119574499e+04, 7.4289924300810162e+04, 7.4299045591365313e+04, 7.4308166991227990e+04, 7.4317288500386247e+04, 7.4326410118828164e+04, 7.4335531846541795e+04, 7.4344653683515178e+04, 7.4353775629736381e+04, 7.4362897685193500e+04, 7.4372019849874574e+04, 7.4381142123767670e+04, 7.4390264506860884e+04, 7.4399386999142283e+04, 7.4408509600599951e+04, 7.4417632311221954e+04, 7.4426755130996418e+04, 7.4435878059911382e+04, 7.4445001097954955e+04, 7.4454124245115236e+04, 7.4463247501380305e+04, 7.4472370866738289e+04, 7.4481494341177240e+04, 7.4490617924685299e+04, 7.4499741617250562e+04, 7.4508865418861111e+04, 7.4517989329505086e+04, 7.4527113349170599e+04, 7.4536237477845760e+04, 7.4545361715518666e+04, 7.4554486062177471e+04, 7.4563610517810303e+04, 7.4572735082405256e+04, 7.4581859755950471e+04, 7.4590984538434088e+04, 7.4600109429844248e+04, 7.4609234430169061e+04, 7.4618359539396668e+04, 7.4627484757515223e+04, 7.4636610084512882e+04, 7.4645735520377755e+04, 7.4654861065098041e+04, 7.4663986718661850e+04, 7.4673112481057353e+04, 7.4682238352272689e+04, 7.4691364332296056e+04, 7.4700490421115581e+04, 7.4709616618719447e+04, 7.4718742925095823e+04, 7.4727869340232850e+04, 7.4736995864118740e+04, 7.4746122496741649e+04, 7.4755249238089731e+04, 7.4764376088151213e+04, 7.4773503046914237e+04, 7.4782630114367028e+04, 7.4791757290497728e+04, 7.4800884575294564e+04, 7.4810011968745719e+04, 7.4819139470839378e+04, 7.4828267081563768e+04, 7.4837394800907059e+04, 7.4846522628857463e+04, 7.4855650565403179e+04, 7.4864778610532419e+04, 7.4873906764233412e+04, 7.4883035026494355e+04, 7.4892163397303477e+04, 7.4901291876648960e+04, 7.4910420464519077e+04, 7.4919549160902010e+04, 7.4928677965786017e+04, 7.4937806879159296e+04, 7.4946935901010103e+04, 7.4956065031326652e+04, 7.4965194270097199e+04, 7.4974323617309972e+04, 7.4983453072953198e+04, 7.4992582637015163e+04, 7.5001712309484079e+04, 7.5010842090348204e+04, 7.5019971979595779e+04, 7.5029101977215061e+04, 7.5038232083194322e+04, 7.5047362297521831e+04, 7.5056492620185818e+04, 7.5065623051174553e+04, 7.5074753590476321e+04, 7.5083884238079394e+04, 7.5093014993972029e+04, 7.5102145858142510e+04, 7.5111276830579096e+04, 7.5120407911270086e+04, 7.5129539100203750e+04, 7.5138670397368391e+04, 7.5147801802752292e+04, 7.5156933316343726e+04, 7.5166064938130992e+04, 7.5175196668102391e+04, 7.5184328506246209e+04, 7.5193460452550746e+04, 7.5202592507004330e+04, 7.5211724669595249e+04, 7.5220856940311816e+04, 7.5229989319142333e+04, 7.5239121806075098e+04, 7.5248254401098442e+04, 7.5257387104200694e+04, 7.5266519915370154e+04, 7.5275652834595166e+04, 7.5284785861864046e+04, 7.5293918997165107e+04, 7.5303052240486693e+04, 7.5312185591817150e+04, 7.5321319051144790e+04, 7.5330452618457944e+04, 7.5339586293745000e+04, 7.5348720076994228e+04, 7.5357853968194046e+04, 7.5366987967332767e+04, 7.5376122074398736e+04, 7.5385256289380326e+04, 7.5394390612265866e+04, 7.5403525043043745e+04, 7.5412659581702304e+04, 7.5421794228229905e+04, 7.5430928982614932e+04, 7.5440063844845718e+04, 7.5449198814910662e+04, 7.5458333892798139e+04, 7.5467469078496521e+04, 7.5476604371994181e+04, 7.5485739773279507e+04, 7.5494875282340858e+04, 7.5504010899166620e+04, 7.5513146623745211e+04, 7.5522282456065004e+04, 7.5531418396114386e+04, 7.5540554443881789e+04, 7.5549690599355570e+04, 7.5558826862524162e+04, 7.5567963233375936e+04, 7.5577099711899311e+04, 7.5586236298082687e+04, 7.5595372991914497e+04, 7.5604509793383142e+04, 7.5613646702477039e+04, 7.5622783719184590e+04, 7.5631920843494227e+04, 7.5641058075394380e+04, 7.5650195414873480e+04, 7.5659332861919946e+04, 7.5668470416522192e+04, 7.5677608078668665e+04, 7.5686745848347811e+04, 7.5695883725548047e+04, 7.5705021710257832e+04, 7.5714159802465598e+04, 7.5723298002159805e+04, 7.5732436309328885e+04, 7.5741574723961283e+04, 7.5750713246045474e+04, 7.5759851875569890e+04, 7.5768990612523005e+04, 7.5778129456893279e+04, 7.5787268408669173e+04, 7.5796407467839148e+04, 7.5805546634391678e+04, 7.5814685908315223e+04, 7.5823825289598259e+04, 7.5832964778229274e+04, 7.5842104374196715e+04, 7.5851244077489086e+04, 7.5860383888094861e+04, 7.5869523806002544e+04, 7.5878663831200611e+04, 7.5887803963677536e+04, 7.5896944203421852e+04, 7.5906084550421991e+04, 7.5915225004666529e+04, 7.5924365566143882e+04, 7.5933506234842615e+04, 7.5942647010751214e+04, 7.5951787893858185e+04, 7.5960928884152017e+04, 7.5970069981621258e+04, 7.5979211186254397e+04, 7.5988352498039952e+04, 7.5997493916966458e+04, 7.6006635443022431e+04, 7.6015777076196406e+04, 7.6024918816476886e+04, 7.6034060663852419e+04, 7.6043202618311523e+04, 7.6052344679842747e+04, 7.6061486848434623e+04, 7.6070629124075684e+04, 7.6079771506754478e+04, 7.6088913996459538e+04, 7.6098056593179426e+04, 7.6107199296902705e+04, 7.6116342107617878e+04, 7.6125485025313552e+04, 7.6134628049978244e+04, 7.6143771181600518e+04, 7.6152914420168949e+04, 7.6162057765672100e+04, 7.6171201218098533e+04, 7.6180344777436825e+04, 7.6189488443675538e+04, 7.6198632216803235e+04, 7.6207776096808491e+04, 7.6216920083679914e+04, 7.6226064177406064e+04, 7.6235208377975549e+04, 7.6244352685376929e+04, 7.6253497099598768e+04, 7.6262641620629714e+04, 7.6271786248458331e+04, 7.6280930983073209e+04, 7.6290075824462954e+04, 7.6299220772616143e+04, 7.6308365827521411e+04, 7.6317510989167378e+04, 7.6326656257542621e+04, 7.6335801632635732e+04, 7.6344947114435374e+04, 7.6354092702930109e+04, 7.6363238398108588e+04, 7.6372384199959430e+04, 7.6381530108471241e+04, 7.6390676123632671e+04, 7.6399822245432326e+04, 7.6408968473858826e+04, 7.6418114808900835e+04, 7.6427261250546973e+04, 7.6436407798785862e+04, 7.6445554453606150e+04, 7.6454701214996487e+04, 7.6463848082945522e+04, 7.6472995057441862e+04, 7.6482142138474213e+04, 7.6491289326031183e+04, 7.6500436620101435e+04, 7.6509584020673647e+04, 7.6518731527736440e+04, 7.6527879141278507e+04, 7.6537026861288497e+04, 7.6546174687755076e+04, 7.6555322620666921e+04, 7.6564470660012710e+04, 7.6573618805781094e+04, 7.6582767057960751e+04, 7.6591915416540374e+04, 7.6601063881508628e+04, 7.6610212452854204e+04, 7.6619361130565798e+04, 7.6628509914632072e+04, 7.6637658805041734e+04, 7.6646807801783463e+04, 7.6655956904845967e+04, 7.6665106114217953e+04, 7.6674255429888086e+04, 7.6683404851845102e+04, 7.6692554380077694e+04, 7.6701704014574556e+04, 7.6710853755324395e+04, 7.6720003602315948e+04, 7.6729153555537909e+04, 7.6738303614978999e+04, 7.6747453780627926e+04, 7.6756604052473442e+04, 7.6765754430504239e+04, 7.6774904914709055e+04, 7.6784055505076627e+04, 7.6793206201595676e+04, 7.6802357004254911e+04, 7.6811507913043126e+04, 7.6820658927949014e+04, 7.6829810048961328e+04, 7.6838961276068818e+04, 7.6848112609260206e+04, 7.6857264048524274e+04, 7.6866415593849742e+04, 7.6875567245225364e+04, 7.6884719002639904e+04, 7.6893870866082128e+04, 7.6903022835540774e+04, 7.6912174911004608e+04, 7.6921327092462423e+04, 7.6930479379902958e+04, 7.6939631773314992e+04, 7.6948784272687291e+04, 7.6957936878008593e+04, 7.6967089589267736e+04, 7.6976242406453472e+04, 7.6985395329554565e+04, 7.6994548358559827e+04, 7.7003701493458037e+04, 7.7012854734237975e+04, 7.7022008080888438e+04, 7.7031161533398204e+04, 7.7040315091756100e+04, 7.7049468755950875e+04, 7.7058622525971368e+04, 7.7067776401806361e+04, 7.7076930383444662e+04, 7.7086084470875096e+04, 7.7095238664086442e+04, 7.7104392963067541e+04, 7.7113547367807201e+04, 7.7122701878294203e+04, 7.7131856494517415e+04, 7.7141011216465617e+04, 7.7150166044127662e+04, 7.7159320977492374e+04, 7.7168476016548550e+04, 7.7177631161285055e+04, 7.7186786411690729e+04, 7.7195941767754382e+04, 7.7205097229464853e+04, 7.7214252796810979e+04, 7.7223408469781600e+04, 7.7232564248365583e+04, 7.7241720132551782e+04, 7.7250876122329006e+04, 7.7260032217686137e+04, 7.7269188418612001e+04, 7.7278344725095463e+04, 7.7287501137125437e+04, 7.7296657654690716e+04, 7.7305814277780170e+04, 7.7314971006382693e+04, 7.7324127840487141e+04, 7.7333284780082395e+04, 7.7342441825157308e+04, 7.7351598975700777e+04, 7.7360756231701671e+04, 7.7369913593148856e+04, 7.7379071060031230e+04, 7.7388228632337661e+04, 7.7397386310057045e+04, 7.7406544093178280e+04, 7.7415701981690261e+04, 7.7424859975581872e+04, 7.7434018074841995e+04, 7.7443176279459542e+04, 7.7452334589423423e+04, 7.7461493004722535e+04, 7.7470651525345776e+04, 7.7479810151282072e+04, 7.7488968882520319e+04, 7.7498127719049429e+04, 7.7507286660858314e+04, 7.7516445707935898e+04, 7.7525604860271109e+04, 7.7534764117852843e+04, 7.7543923480670055e+04, 7.7553082948711657e+04, 7.7562242521966589e+04, 7.7571402200423763e+04, 7.7580561984072134e+04, 7.7589721872900613e+04, 7.7598881866898155e+04, 7.7608041966053686e+04, 7.7617202170356177e+04, 7.7626362479794538e+04, 7.7635522894357739e+04, 7.7644683414034720e+04, 7.7653844038814452e+04, 7.7663004768685874e+04, 7.7672165603637957e+04, 7.7681326543659641e+04, 7.7690487588739896e+04, 7.7699648738867676e+04, 7.7708809994031952e+04, 7.7717971354221707e+04, 7.7727132819425911e+04, 7.7736294389633520e+04, 7.7745456064833503e+04, 7.7754617845014873e+04, 7.7763779730166585e+04, 7.7772941720277624e+04, 7.7782103815336974e+04, 7.7791266015333633e+04, 7.7800428320256557e+04, 7.7809590730094773e+04, 7.7818753244837266e+04, 7.7827915864473034e+04, 7.7837078588991077e+04, 7.7846241418380378e+04, 7.7855404352629965e+04, 7.7864567391728822e+04, 7.7873730535665978e+04, 7.7882893784430416e+04, 7.7892057138011180e+04, 7.7901220596397266e+04, 7.7910384159577690e+04, 7.7919547827541479e+04, 7.7928711600277631e+04, 7.7937875477775204e+04, 7.7947039460023196e+04, 7.7956203547010664e+04, 7.7965367738726607e+04, 7.7974532035160082e+04, 7.7983696436300132e+04, 7.7992860942135769e+04, 7.8002025552656036e+04, 7.8011190267849975e+04, 7.8020355087706659e+04, 7.8029520012215100e+04, 7.8038685041364370e+04, 7.8047850175143511e+04, 7.8057015413541565e+04, 7.8066180756547605e+04, 7.8075346204150686e+04, 7.8084511756339882e+04, 7.8093677413104233e+04, 7.8102843174432826e+04, 7.8112009040314690e+04, 7.8121175010738923e+04, 7.8130341085694585e+04, 7.8139507265170774e+04, 7.8148673549156549e+04, 7.8157839937641009e+04, 7.8167006430613197e+04, 7.8176173028062229e+04, 7.8185339729977175e+04, 7.8194506536347122e+04, 7.8203673447161171e+04, 7.8212840462408407e+04, 7.8222007582077931e+04, 7.8231174806158844e+04, 7.8240342134640217e+04, 7.8249509567511181e+04, 7.8258677104760820e+04, 7.8267844746378265e+04, 7.8277012492352616e+04, 7.8286180342672975e+04, 7.8295348297328455e+04, 7.8304516356308159e+04, 7.8313684519601229e+04, 7.8322852787196767e+04, 7.8332021159083932e+04, 7.8341189635251809e+04, 7.8350358215689543e+04, 7.8359526900386249e+04, 7.8368695689331071e+04, 7.8377864582513124e+04, 7.8387033579921568e+04, 7.8396202681545517e+04, 7.8405371887374160e+04, 7.8414541197396582e+04, 7.8423710611601957e+04, 7.8432880129979414e+04, 7.8442049752518113e+04, 7.8451219479207197e+04, 7.8460389310035840e+04, 7.8469559244993186e+04, 7.8478729284068395e+04, 7.8487899427250610e+04, 7.8497069674529019e+04, 7.8506240025892796e+04, 7.8515410481331070e+04, 7.8524581040833044e+04, 7.8533751704387847e+04, 7.8542922471984697e+04, 7.8552093343612753e+04, 7.8561264319261216e+04, 7.8570435398919231e+04, 7.8579606582576016e+04, 7.8588777870220729e+04, 7.8597949261842557e+04, 7.8607120757430705e+04, 7.8616292356974373e+04, 7.8625464060462720e+04, 7.8634635867884979e+04, 7.8643807779230337e+04, 7.8652979794487997e+04, 7.8662151913647147e+04, 7.8671324136697018e+04, 7.8680496463626798e+04, 7.8689668894425704e+04, 7.8698841429082924e+04, 7.8708014067587719e+04, 7.8717186809929277e+04, 7.8726359656096829e+04, 7.8735532606079592e+04, 7.8744705659866784e+04, 7.8753878817447621e+04, 7.8763052078811364e+04, 7.8772225443947216e+04, 7.8781398912844408e+04, 7.8790572485492201e+04, 7.8799746161879797e+04, 7.8808919941996472e+04, 7.8818093825831442e+04, 7.8827267813373939e+04, 7.8836441904613253e+04, 7.8845616099538587e+04, 7.8854790398139216e+04, 7.8863964800404385e+04, 7.8873139306323355e+04, 7.8882313915885374e+04, 7.8891488629079700e+04, 7.8900663445895625e+04, 7.8909838366322365e+04, 7.8919013390349210e+04, 7.8928188517965435e+04, 7.8937363749160300e+04, 7.8946539083923068e+04, 7.8955714522243026e+04, 7.8964890064109466e+04, 7.8974065709511633e+04, 7.8983241458438846e+04, 7.8992417310880366e+04, 7.9001593266825483e+04, 7.9010769326263500e+04, 7.9019945489183679e+04, 7.9029121755575310e+04, 7.9038298125427726e+04, 7.9047474598730201e+04, 7.9056651175472012e+04, 7.9065827855642492e+04, 7.9075004639230945e+04, 7.9084181526226661e+04, 7.9093358516618959e+04, 7.9102535610397143e+04, 7.9111712807550517e+04, 7.9120890108068415e+04, 7.9130067511940171e+04, 7.9139245019155045e+04, 7.9148422629702400e+04, 7.9157600343571568e+04, 7.9166778160751841e+04, 7.9175956081232565e+04, 7.9185134105003061e+04, 7.9194312232052675e+04, 7.9203490462370741e+04, 7.9212668795946563e+04, 7.9221847232769534e+04, 7.9231025772828943e+04, 7.9240204416114182e+04, 7.9249383162614569e+04, 7.9258562012319453e+04, 7.9267740965218196e+04, 7.9276920021300146e+04, 7.9286099180554636e+04, 7.9295278442971045e+04, 7.9304457808538704e+04, 7.9313637277247020e+04, 7.9322816849085328e+04, 7.9331996524043003e+04, 7.9341176302109379e+04, 7.9350356183273878e+04, 7.9359536167525832e+04, 7.9368716254854633e+04, 7.9377896445249673e+04, 7.9387076738700285e+04, 7.9396257135195905e+04, 7.9405437634725866e+04, 7.9414618237279574e+04, 7.9423798942846421e+04, 7.9432979751415784e+04, 7.9442160662977083e+04, 7.9451341677519667e+04, 7.9460522795032972e+04, 7.9469704015506388e+04, 7.9478885338929293e+04, 7.9488066765291107e+04, 7.9497248294581223e+04, 7.9506429926789075e+04, 7.9515611661904040e+04, 7.9524793499915540e+04, 7.9533975440813010e+04, 7.9543157484585827e+04, 7.9552339631223440e+04, 7.9561521880715241e+04, 7.9570704233050681e+04, 7.9579886688219165e+04, 7.9589069246210129e+04, 7.9598251907013007e+04, 7.9607434670617222e+04, 7.9616617537012178e+04, 7.9625800506187355e+04, 7.9634983578132174e+04, 7.9644166752836085e+04, 7.9653350030288508e+04, 7.9662533410478893e+04, 7.9671716893396704e+04, 7.9680900479031377e+04, 7.9690084167372363e+04, 7.9699267958409109e+04, 7.9708451852131067e+04, 7.9717635848527701e+04, 7.9726819947588476e+04, 7.9736004149302855e+04, 7.9745188453660274e+04, 7.9754372860650226e+04, 7.9763557370262177e+04, 7.9772741982485590e+04, 7.9781926697309944e+04, 7.9791111514724689e+04, 7.9800296434719319e+04, 7.9809481457283298e+04, 7.9818666582406135e+04, 7.9827851810077278e+04, 7.9837037140286251e+04, 7.9846222573022518e+04, 7.9855408108275587e+04, 7.9864593746034923e+04, 7.9873779486290019e+04, 7.9882965329030398e+04, 7.9892151274245538e+04, 7.9901337321924948e+04, 7.9910523472058121e+04, 7.9919709724634566e+04, 7.9928896079643790e+04, 7.9938082537075301e+04, 7.9947269096918608e+04, 7.9956455759163233e+04, 7.9965642523798684e+04, 7.9974829390814470e+04, 7.9984016360200127e+04, 7.9993203431945163e+04, 8.0002390606039116e+04, 8.0011577882471494e+04, 8.0020765261231834e+04, 8.0029952742309673e+04, 8.0039140325694520e+04, 8.0048328011375925e+04, 8.0057515799343411e+04, 8.0066703689586546e+04, 8.0075891682094865e+04, 8.0085079776857907e+04, 8.0094267973865193e+04, 8.0103456273106291e+04, 8.0112644674570751e+04, 8.0121833178248111e+04, 8.0131021784127937e+04, 8.0140210492199796e+04, 8.0149399302453210e+04, 8.0158588214877775e+04, 8.0167777229463041e+04, 8.0176966346198577e+04, 8.0186155565073932e+04, 8.0195344886078688e+04, 8.0204534309202398e+04, 8.0213723834434655e+04, 8.0222913461765042e+04, 8.0232103191183123e+04, 8.0241293022678467e+04, 8.0250482956240652e+04, 8.0259672991859290e+04, 8.0268863129523947e+04, 8.0278053369224232e+04, 8.0287243710949697e+04, 8.0296434154689952e+04, 8.0305624700434608e+04, 8.0314815348173244e+04, 8.0324006097895457e+04, 8.0333196949590842e+04, 8.0342387903249022e+04, 8.0351578958859580e+04, 8.0360770116412154e+04, 8.0369961375896310e+04, 8.0379152737301672e+04, 8.0388344200617881e+04, 8.0397535765834531e+04, 8.0406727432941232e+04, 8.0415919201927638e+04, 8.0425111072783315e+04, 8.0434303045497931e+04, 8.0443495120061096e+04, 8.0452687296462464e+04, 8.0461879574691629e+04, 8.0471071954738218e+04, 8.0480264436591897e+04, 8.0489457020242276e+04, 8.0498649705679010e+04, 8.0507842492891737e+04, 8.0517035381870126e+04, 8.0526228372603771e+04, 8.0535421465082341e+04, 8.0544614659295476e+04, 8.0553807955232842e+04, 8.0563001352884079e+04, 8.0572194852238841e+04, 8.0581388453286825e+04, 8.0590582156017656e+04, 8.0599775960420986e+04, 8.0608969866486499e+04, 8.0618163874203849e+04, 8.0627357983562732e+04, 8.0636552194552773e+04, 8.0645746507163683e+04, 8.0654940921385103e+04, 8.0664135437206729e+04, 8.0673330054618244e+04, 8.0682524773609301e+04, 8.0691719594169626e+04, 8.0700914516288860e+04, 8.0710109539956728e+04, 8.0719304665162912e+04, 8.0728499891897067e+04, 8.0737695220148933e+04, 8.0746890649908179e+04, 8.0756086181164501e+04, 8.0765281813907597e+04, 8.0774477548127179e+04, 8.0783673383812958e+04, 8.0792869320954633e+04, 8.0802065359541913e+04, 8.0811261499564484e+04, 8.0820457741012084e+04, 8.0829654083874411e+04, 8.0838850528141193e+04, 8.0848047073802139e+04, 8.0857243720846993e+04, 8.0866440469265450e+04, 8.0875637319047222e+04, 8.0884834270182066e+04, 8.0894031322659706e+04, 8.0903228476469885e+04, 8.0912425731602314e+04, 8.0921623088046734e+04, 8.0930820545792871e+04, 8.0940018104830480e+04, 8.0949215765149289e+04, 8.0958413526739067e+04, 8.0967611389589525e+04, 8.0976809353690434e+04, 8.0986007419031535e+04, 8.0995205585602584e+04, 8.1004403853393334e+04, 8.1013602222393529e+04, 8.1022800692592937e+04, 8.1031999263981314e+04, 8.1041197936548429e+04, 8.1050396710284025e+04, 8.1059595585177900e+04, 8.1068794561219809e+04, 8.1077993638399508e+04, 8.1087192816706782e+04, 8.1096392096131400e+04, 8.1105591476663147e+04, 8.1114790958291778e+04, 8.1123990541007093e+04, 8.1133190224798876e+04, 8.1142390009656912e+04, 8.1151589895570971e+04, 8.1160789882530866e+04, 8.1169989970526382e+04, 8.1179190159547303e+04, 8.1188390449583414e+04, 8.1197590840624529e+04, 8.1206791332660461e+04, 8.1215991925680981e+04, 8.1225192619675901e+04, 8.1234393414635022e+04, 8.1243594310548171e+04, 8.1252795307405147e+04, 8.1261996405195765e+04, 8.1271197603909837e+04, 8.1280398903537163e+04, 8.1289600304067571e+04, 8.1298801805490875e+04, 8.1308003407796918e+04, 8.1317205110975512e+04, 8.1326406915016443e+04, 8.1335608819909598e+04, 8.1344810825644789e+04, 8.1354012932211830e+04, 8.1363215139600565e+04, 8.1372417447800850e+04, 8.1381619856802485e+04, 8.1390822366595356e+04, 8.1400024977169262e+04, 8.1409227688514089e+04, 8.1418430500619637e+04, 8.1427633413475793e+04, 8.1436836427072383e+04, 8.1446039541399266e+04, 8.1455242756446300e+04, 8.1464446072203340e+04, 8.1473649488660245e+04, 8.1482853005806872e+04, 8.1492056623633092e+04, 8.1501260342128735e+04, 8.1510464161283729e+04, 8.1519668081087890e+04, 8.1528872101531131e+04, 8.1538076222603282e+04, 8.1547280444294258e+04, 8.1556484766593931e+04, 8.1565689189492143e+04, 8.1574893712978796e+04, 8.1584098337043775e+04, 8.1593303061676968e+04, 8.1602507886868247e+04, 8.1611712812607526e+04, 8.1620917838884678e+04, 8.1630122965689588e+04, 8.1639328193012174e+04, 8.1648533520842335e+04, 8.1657738949169943e+04, 8.1666944477984929e+04, 8.1676150107277179e+04, 8.1685355837036594e+04, 8.1694561667253103e+04, 8.1703767597916594e+04, 8.1712973629016982e+04, 8.1722179760544182e+04, 8.1731385992488111e+04, 8.1740592324838683e+04, 8.1749798757585842e+04, 8.1759005290719477e+04, 8.1768211924229516e+04, 8.1777418658105904e+04, 8.1786625492338557e+04, 8.1795832426917405e+04, 8.1805039461832377e+04, 8.1814246597073405e+04, 8.1823453832630417e+04, 8.1832661168493374e+04, 8.1841868604652205e+04, 8.1851076141096855e+04, 8.1860283777817254e+04, 8.1869491514803347e+04, 8.1878699352045107e+04, 8.1887907289532479e+04, 8.1897115327255393e+04, 8.1906323465203808e+04, 8.1915531703367698e+04, 8.1924740041737008e+04, 8.1933948480301682e+04, 8.1943157019051694e+04, 8.1952365657977032e+04, 8.1961574397067656e+04, 8.1970783236313524e+04, 8.1979992175704581e+04, 8.1989201215230816e+04, 8.1998410354882202e+04, 8.2007619594648742e+04, 8.2016828934520381e+04, 8.2026038374487121e+04, 8.2035247914538922e+04, 8.2044457554665787e+04, 8.2053667294857689e+04, 8.2062877135104645e+04, 8.2072087075396586e+04, 8.2081297115723559e+04, 8.2090507256075522e+04, 8.2099717496442492e+04, 8.2108927836814459e+04, 8.2118138277181439e+04, 8.2127348817533435e+04, 8.2136559457860421e+04, 8.2145770198152430e+04, 8.2154981038399434e+04, 8.2164191978591480e+04, 8.2173403018718585e+04, 8.2182614158770739e+04, 8.2191825398737943e+04, 8.2201036738610244e+04, 8.2210248178377660e+04, 8.2219459718030208e+04, 8.2228671357557905e+04, 8.2237883096950798e+04, 8.2247094936198904e+04, 8.2256306875292226e+04, 8.2265518914220826e+04, 8.2274731052974748e+04, 8.2283943291543997e+04, 8.2293155629918649e+04, 8.2302368068088705e+04, 8.2311580606044226e+04, 8.2320793243775246e+04, 8.2330005981271825e+04, 8.2339218818523994e+04, 8.2348431755521800e+04, 8.2357644792255334e+04, 8.2366857928714613e+04, 8.2376071164889698e+04, 8.2385284500770664e+04, 8.2394497936347543e+04, 8.2403711471610412e+04, 8.2412925106549344e+04, 8.2422138841154389e+04, 8.2431352675415605e+04, 8.2440566609323098e+04, 8.2449780642866914e+04, 8.2458994776037129e+04, 8.2468209008823818e+04, 8.2477423341217072e+04, 8.2486637773206952e+04, 8.2495852304783562e+04, 8.2505066935936949e+04, 8.2514281666657218e+04, 8.2523496496934473e+04, 8.2532711426758789e+04, 8.2541926456120244e+04, 8.2551141585008940e+04, 8.2560356813414997e+04, 8.2569572141328477e+04, 8.2578787568739514e+04, 8.2588003095638182e+04, 8.2597218722014586e+04, 8.2606434447858832e+04, 8.2615650273161038e+04, 8.2624866197911309e+04, 8.2634082222099751e+04, 8.2643298345716481e+04, 8.2652514568751620e+04, 8.2661730891195286e+04, 8.2670947313037570e+04, 8.2680163834268620e+04, 8.2689380454878556e+04, 8.2698597174857496e+04, 8.2707813994195574e+04, 8.2717030912882896e+04, 8.2726247930909623e+04, 8.2735465048265876e+04, 8.2744682264941803e+04, 8.2753899580927522e+04, 8.2763116996213168e+04, 8.2772334510788889e+04, 8.2781552124644833e+04, 8.2790769837771135e+04, 8.2799987650157942e+04, 8.2809205561795432e+04, 8.2818423572673695e+04, 8.2827641682782938e+04, 8.2836859892113309e+04, 8.2846078200654927e+04, 8.2855296608397999e+04, 8.2864515115332644e+04, 8.2873733721449040e+04, 8.2882952426737334e+04, 8.2892171231187735e+04, 8.2901390134790359e+04, 8.2910609137535415e+04, 8.2919828239413066e+04, 8.2929047440413473e+04, 8.2938266740526829e+04, 8.2947486139743312e+04, 8.2956705638053056e+04, 8.2965925235446295e+04, 8.2975144931913208e+04, 8.2984364727443943e+04, 8.2993584622028720e+04, 8.3002804615657718e+04, 8.3012024708321143e+04, 8.3021244900009173e+04, 8.3030465190711999e+04, 8.3039685580419842e+04, 8.3048906069122866e+04, 8.3058126656811306e+04, 8.3067347343475340e+04, 8.3076568129105188e+04, 8.3085789013691043e+04, 8.3095009997223111e+04, 8.3104231079691614e+04, 8.3113452261086772e+04, 8.3122673541398792e+04, 8.3131894920617866e+04, 8.3141116398734244e+04, 8.3150337975738148e+04, 8.3159559651619784e+04, 8.3168781426369373e+04, 8.3178003299977165e+04, 8.3187225272433381e+04, 8.3196447343728229e+04, 8.3205669513851957e+04, 8.3214891782794788e+04, 8.3224114150546986e+04, 8.3233336617098772e+04, 8.3242559182440353e+04, 8.3251781846562037e+04, 8.3261004609454001e+04, 8.3270227471106540e+04, 8.3279450431509875e+04, 8.3288673490654255e+04, 8.3297896648529946e+04, 8.3307119905127198e+04, 8.3316343260436246e+04, 8.3325566714447385e+04, 8.3334790267150820e+04, 8.3344013918536875e+04, 8.3353237668595757e+04, 8.3362461517317744e+04, 8.3371685464693132e+04, 8.3380909510712154e+04, 8.3390133655365105e+04, 8.3399357898642251e+04, 8.3408582240533855e+04, 8.3417806681030226e+04, 8.3427031220121600e+04, 8.3436255857798285e+04, 8.3445480594050532e+04, 8.3454705428868663e+04, 8.3463930362242943e+04, 8.3473155394163667e+04, 8.3482380524621098e+04, 8.3491605753605560e+04, 8.3500831081107361e+04, 8.3510056507116751e+04, 8.3519282031624040e+04, 8.3528507654619549e+04, 8.3537733376093558e+04, 8.3546959196036390e+04, 8.3556185114438340e+04, 8.3565411131289700e+04, 8.3574637246580809e+04, 8.3583863460301931e+04, 8.3593089772443418e+04, 8.3602316182995564e+04, 8.3611542691948707e+04, 8.3620769299293141e+04, 8.3629996005019217e+04, 8.3639222809117229e+04, 8.3648449711577516e+04, 8.3657676712390370e+04, 8.3666903811546159e+04, 8.3676131009035191e+04, 8.3685358304847818e+04, 8.3694585698974348e+04, 8.3703813191405148e+04, 8.3713040782130527e+04, 8.3722268471140836e+04, 8.3731496258426414e+04, 8.3740724143977597e+04, 8.3749952127784753e+04, 8.3759180209838203e+04, 8.3768408390128316e+04, 8.3777636668645428e+04, 8.3786865045379905e+04, 8.3796093520322058e+04, 8.3805322093462310e+04, 8.3814550764790969e+04, 8.3823779534298417e+04, 8.3833008401975007e+04, 8.3842237367811103e+04, 8.3851466431797075e+04, 8.3860695593923301e+04, 8.3869924854180121e+04, 8.3879154212557914e+04, 8.3888383669047093e+04, 8.3897613223637993e+04, 8.3906842876320996e+04, 8.3916072627086483e+04, 8.3925302475924866e+04, 8.3934532422826465e+04, 8.3943762467781722e+04, 8.3952992610780988e+04, 8.3962222851814673e+04, 8.3971453190873159e+04, 8.3980683627946826e+04, 8.3989914163026100e+04, 8.3999144796101347e+04, 8.4008375527162978e+04, 8.4017606356201373e+04, 8.4026837283206958e+04, 8.4036068308170143e+04, 8.4045299431081294e+04, 8.4054530651930851e+04, 8.4063761970709209e+04, 8.4072993387406794e+04, 8.4082224902014001e+04, 8.4091456514521240e+04, 8.4100688224918966e+04, 8.4109920033197559e+04, 8.4119151939347474e+04, 8.4128383943359091e+04, 8.4137616045222851e+04, 8.4146848244929162e+04, 8.4156080542468495e+04, 8.4165312937831259e+04, 8.4174545431007864e+04, 8.4183778021988765e+04, 8.4193010710764414e+04, 8.4202243497325224e+04, 8.4211476381661632e+04, 8.4220709363764094e+04, 8.4229942443623033e+04, 8.4239175621228918e+04, 8.4248408896572175e+04, 8.4257642269643242e+04, 8.4266875740432617e+04, 8.4276109308930696e+04, 8.4285342975127976e+04, 8.4294576739014883e+04, 8.4303810600581914e+04, 8.4313044559819493e+04, 8.4322278616718075e+04, 8.4331512771268142e+04, 8.4340747023460164e+04, 8.4349981373284594e+04, 8.4359215820731915e+04, 8.4368450365792596e+04, 8.4377685008457091e+04, 8.4386919748715896e+04, 8.4396154586559482e+04, 8.4405389521978301e+04, 8.4414624554962866e+04, 8.4423859685503645e+04, 8.4433094913591136e+04, 8.4442330239215822e+04, 8.4451565662368157e+04, 8.4460801183038653e+04, 8.4470036801217822e+04, 8.4479272516896133e+04, 8.4488508330064084e+04, 8.4497744240712171e+04, 8.4506980248830892e+04, 8.4516216354410761e+04, 8.4525452557442244e+04, 8.4534688857915899e+04, 8.4543925255822192e+04, 8.4553161751151652e+04, 8.4562398343894776e+04, 8.4571635034042061e+04, 8.4580871821584064e+04, 8.4590108706511252e+04, 8.4599345688814181e+04, 8.4608582768483335e+04, 8.4617819945509269e+04, 8.4627057219882481e+04, 8.4636294591593512e+04, 8.4645532060632875e+04, 8.4654769626991096e+04, 8.4664007290658730e+04, 8.4673245051626276e+04, 8.4682482909884289e+04, 8.4691720865423296e+04, 8.4700958918233839e+04, 8.4710197068306457e+04, 8.4719435315631665e+04, 8.4728673660200046e+04, 8.4737912102002127e+04, 8.4747150641028464e+04, 8.4756389277269598e+04, 8.4765628010716071e+04, 8.4774866841358453e+04, 8.4784105769187270e+04, 8.4793344794193108e+04, 8.4802583916366508e+04, 8.4811823135698040e+04, 8.4821062452178245e+04, 8.4830301865797694e+04, 8.4839541376546971e+04, 8.4848780984416648e+04, 8.4858020689397250e+04, 8.4867260491479377e+04, 8.4876500390653600e+04, 8.4885740386910489e+04, 8.4894980480240614e+04, 8.4904220670634561e+04, 8.4913460958082898e+04, 8.4922701342576227e+04, 8.4931941824105117e+04, 8.4941182402660139e+04, 8.4950423078231906e+04, 8.4959663850811005e+04, 8.4968904720388018e+04, 8.4978145686953547e+04, 8.4987386750498161e+04, 8.4996627911012460e+04, 8.5005869168487086e+04, 8.5015110522912582e+04, 8.5024351974279591e+04, 8.5033593522578682e+04, 8.5042835167800484e+04, 8.5052076909935597e+04, 8.5061318748974634e+04, 8.5070560684908196e+04, 8.5079802717726925e+04, 8.5089044847421377e+04, 8.5098287073982225e+04, 8.5107529397400052e+04, 8.5116771817665489e+04, 8.5126014334769177e+04, 8.5135256948701703e+04, 8.5144499659453708e+04, 8.5153742467015822e+04, 8.5162985371378658e+04, 8.5172228372532874e+04, 8.5181471470469085e+04, 8.5190714665177933e+04, 8.5199957956650062e+04, 8.5209201344876084e+04, 8.5218444829846645e+04, 8.5227688411552401e+04, 8.5236932089983995e+04, 8.5246175865132056e+04, 8.5255419736987242e+04, 8.5264663705540195e+04, 8.5273907770781560e+04, 8.5283151932702007e+04, 8.5292396191292180e+04, 8.5301640546542752e+04, 8.5310884998444351e+04, 8.5320129546987650e+04, 8.5329374192163319e+04, 8.5338618933962003e+04, 8.5347863772374389e+04, 8.5357108707391133e+04, 8.5366353739002909e+04, 8.5375598867200373e+04, 8.5384844091974199e+04, 8.5394089413315072e+04, 8.5403334831213666e+04, 8.5412580345660666e+04, 8.5421825956646731e+04, 8.5431071664162548e+04, 8.5440317468198802e+04, 8.5449563368746167e+04, 8.5458809365795358e+04, 8.5468055459337047e+04, 8.5477301649361922e+04, 8.5486547935860668e+04, 8.5495794318823988e+04, 8.5505040798242582e+04, 8.5514287374107138e+04, 8.5523534046408357e+04, 8.5532780815136939e+04, 8.5542027680283602e+04, 8.5551274641839031e+04, 8.5560521699793942e+04, 8.5569768854139038e+04, 8.5579016104865033e+04, 8.5588263451962644e+04, 8.5597510895422558e+04, 8.5606758435235533e+04, 8.5616006071392258e+04, 8.5625253803883461e+04, 8.5634501632699845e+04, 8.5643749557832140e+04, 8.5652997579271090e+04, 8.5662245697007413e+04, 8.5671493911031808e+04, 8.5680742221335051e+04, 8.5689990627907828e+04, 8.5699239130740898e+04, 8.5708487729824978e+04, 8.5717736425150841e+04, 8.5726985216709189e+04, 8.5736234104490781e+04, 8.5745483088486362e+04, 8.5754732168686649e+04, 8.5763981345082400e+04, 8.5773230617664391e+04, 8.5782479986423336e+04, 8.5791729451349995e+04, 8.5800979012435113e+04, 8.5810228669669479e+04, 8.5819478423043794e+04, 8.5828728272548862e+04, 8.5837978218175442e+04, 8.5847228259914278e+04, 8.5856478397756131e+04, 8.5865728631691745e+04, 8.5874978961711939e+04, 8.5884229387807442e+04, 8.5893479909969043e+04, 8.5902730528187516e+04, 8.5911981242453636e+04, 8.5921232052758176e+04, 8.5930482959091896e+04, 8.5939733961445585e+04, 8.5948985059810031e+04, 8.5958236254176009e+04, 8.5967487544534320e+04, 8.5976738930875741e+04, 8.5985990413191030e+04, 8.5995241991471019e+04, 8.6004493665706483e+04, 8.6013745435888224e+04, 8.6022997302007017e+04, 8.6032249264053680e+04, 8.6041501322019001e+04, 8.6050753475893784e+04, 8.6060005725668845e+04, 8.6069258071334960e+04, 8.6078510512882945e+04, 8.6087763050303634e+04, 8.6097015683587786e+04, 8.6106268412726247e+04, 8.6115521237709821e+04, 8.6124774158529341e+04, 8.6134027175175594e+04, 8.6143280287639413e+04, 8.6152533495911615e+04, 8.6161786799983034e+04, 8.6171040199844458e+04, 8.6180293695486762e+04, 8.6189547286900735e+04, 8.6198800974077225e+04, 8.6208054757007048e+04, 8.6217308635681038e+04, 8.6226562610090041e+04, 8.6235816680224874e+04, 8.6245070846076400e+04, 8.6254325107635421e+04, 8.6263579464892813e+04, 8.6272833917839424e+04, 8.6282088466466070e+04, 8.6291343110763584e+04, 8.6300597850722872e+04, 8.6309852686334722e+04, 8.6319107617590038e+04, 8.6328362644479639e+04, 8.6337617766994401e+04, 8.6346872985125156e+04, 8.6356128298862764e+04, 8.6365383708198104e+04, 8.6374639213122035e+04, 8.6383894813625418e+04, 8.6393150509699117e+04, 8.6402406301333991e+04, 8.6411662188520917e+04, 8.6420918171250785e+04, 8.6430174249514414e+04, 8.6439430423302736e+04, 8.6448686692606600e+04, 8.6457943057416880e+04, 8.6467199517724468e+04, 8.6476456073520239e+04, 8.6485712724795070e+04, 8.6494969471539851e+04, 8.6504226313745487e+04, 8.6513483251402839e+04, 8.6522740284502812e+04, 8.6531997413036297e+04, 8.6541254636994170e+04, 8.6550511956367351e+04, 8.6559769371146715e+04, 8.6569026881323152e+04, 8.6578284486887584e+04, 8.6587542187830928e+04, 8.6596799984144047e+04, 8.6606057875817889e+04, 8.6615315862843330e+04, 8.6624573945211305e+04, 8.6633832122912703e+04, 8.6643090395938445e+04, 8.6652348764279421e+04, 8.6661607227926579e+04, 8.6670865786870840e+04, 8.6680124441103093e+04, 8.6689383190614273e+04, 8.6698642035395300e+04, 8.6707900975437093e+04, 8.6717160010730615e+04, 8.6726419141266757e+04, 8.6735678367036468e+04, 8.6744937688030666e+04, 8.6754197104240287e+04, 8.6763456615656265e+04, 8.6772716222269548e+04, 8.6781975924071055e+04, 8.6791235721051737e+04, 8.6800495613202555e+04, 8.6809755600514429e+04, 8.6819015682978294e+04, 8.6828275860585112e+04, 8.6837536133325833e+04, 8.6846796501191420e+04, 8.6856056964172793e+04, 8.6865317522260942e+04, 8.6874578175446790e+04, 8.6883838923721327e+04, 8.6893099767075473e+04, 8.6902360705500207e+04, 8.6911621738986505e+04, 8.6920882867525303e+04, 8.6930144091107621e+04, 8.6939405409724364e+04, 8.6948666823366526e+04, 8.6957928332025083e+04, 8.6967189935691000e+04, 8.6976451634355253e+04, 8.6985713428008807e+04, 8.6994975316642667e+04, 8.7004237300247813e+04, 8.7013499378815191e+04, 8.7022761552335796e+04, 8.7032023820800649e+04, 8.7041286184200697e+04, 8.7050548642526905e+04, 8.7059811195770337e+04, 8.7069073843921942e+04, 8.7078336586972684e+04, 8.7087599424913642e+04, 8.7096862357735721e+04, 8.7106125385429972e+04, 8.7115388507987402e+04, 8.7124651725398973e+04, 8.7133915037655723e+04, 8.7143178444748628e+04, 8.7152441946668725e+04, 8.7161705543406992e+04, 8.7170969234954478e+04, 8.7180233021302149e+04, 8.7189496902441053e+04, 8.7198760878362198e+04, 8.7208024949056606e+04, 8.7217289114515268e+04, 8.7226553374729236e+04, 8.7235817729689501e+04, 8.7245082179387129e+04, 8.7254346723813127e+04, 8.7263611362958502e+04, 8.7272876096814318e+04, 8.7282140925371583e+04, 8.7291405848621318e+04, 8.7300670866554603e+04, 8.7309935979162430e+04, 8.7319201186435850e+04, 8.7328466488365899e+04, 8.7337731884943627e+04, 8.7346997376160070e+04, 8.7356262962006265e+04, 8.7365528642473291e+04, 8.7374794417552141e+04, 8.7384060287233922e+04, 8.7393326251509628e+04, 8.7402592310370368e+04, 8.7411858463807162e+04, 8.7421124711811077e+04, 8.7430391054373162e+04, 8.7439657491484497e+04, 8.7448924023136133e+04, 8.7458190649319120e+04, 8.7467457370024524e+04, 8.7476724185243438e+04, 8.7485991094966885e+04, 8.7495258099185972e+04, 8.7504525197891766e+04, 8.7513792391075345e+04, 8.7523059678727761e+04, 8.7532327060840093e+04, 8.7541594537403435e+04, 8.7550862108408866e+04, 8.7560129773847453e+04, 8.7569397533710304e+04, 8.7578665387988483e+04, 8.7587933336673086e+04, 8.7597201379755192e+04, 8.7606469517225894e+04, 8.7615737749076274e+04, 8.7625006075297439e+04, 8.7634274495880498e+04, 8.7643543010816516e+04, 8.7652811620096603e+04, 8.7662080323711882e+04, 8.7671349121653431e+04, 8.7680618013912361e+04, 8.7689887000479765e+04, 8.7699156081346766e+04, 8.7708425256504459e+04, 8.7717694525943967e+04, 8.7726963889656399e+04, 8.7736233347632835e+04, 8.7745502899864456e+04, 8.7754772546342341e+04, 8.7764042287057600e+04, 8.7773312122001371e+04, 8.7782582051164747e+04, 8.7791852074538896e+04, 8.7801122192114897e+04, 8.7810392403883889e+04, 8.7819662709837052e+04, 8.7828933109965437e+04, 8.7838203604260241e+04, 8.7847474192712558e+04, 8.7856744875313540e+04, 8.7866015652054310e+04, 8.7875286522926021e+04, 8.7884557487919825e+04, 8.7893828547026831e+04, 8.7903099700238206e+04, 8.7912370947545103e+04, 8.7921642288938630e+04, 8.7930913724409969e+04, 8.7940185253950258e+04, 8.7949456877550634e+04, 8.7958728595202294e+04, 8.7968000406896361e+04, 8.7977272312623987e+04, 8.7986544312376369e+04, 8.7995816406144615e+04, 8.8005088593919907e+04, 8.8014360875693426e+04, 8.8023633251456326e+04, 8.8032905721199757e+04, 8.8042178284914917e+04, 8.8051450942592943e+04, 8.8060723694225046e+04, 8.8069996539802378e+04, 8.8079269479316106e+04, 8.8088542512757413e+04, 8.8097815640117493e+04, 8.8107088861387514e+04, 8.8116362176558643e+04, 8.8125635585622105e+04, 8.8134909088569038e+04, 8.8144182685390668e+04, 8.8153456376078131e+04, 8.8162730160622668e+04, 8.8172004039015446e+04, 8.8181278011247676e+04, 8.8190552077310538e+04, 8.8199826237195230e+04, 8.8209100490892932e+04, 8.8218374838394884e+04, 8.8227649279692283e+04, 8.8236923814776281e+04, 8.8246198443638146e+04, 8.8255473166269061e+04, 8.8264747982660207e+04, 8.8274022892802837e+04, 8.8283297896688135e+04, 8.8292572994307324e+04, 8.8301848185651630e+04, 8.8311123470712249e+04, 8.8320398849480407e+04, 8.8329674321947328e+04, 8.8338949888104224e+04, 8.8348225547942318e+04, 8.8357501301452852e+04, 8.8366777148627021e+04, 8.8376053089456080e+04, 8.8385329123931253e+04, 8.8394605252043766e+04, 8.8403881473784859e+04, 8.8413157789145771e+04, 8.8422434198117713e+04, 8.8431710700691940e+04, 8.8440987296859690e+04, 8.8450263986612190e+04, 8.8459540769940737e+04, 8.8468817646836513e+04, 8.8478094617290786e+04, 8.8487371681294797e+04, 8.8496648838839814e+04, 8.8505926089917091e+04, 8.8515203434517840e+04, 8.8524480872633358e+04, 8.8533758404254884e+04, 8.8543036029373659e+04, 8.8552313747980967e+04, 8.8561591560068060e+04, 8.8570869465626209e+04, 8.8580147464646652e+04, 8.8589425557120689e+04, 8.8598703743039572e+04, 8.8607982022394543e+04, 8.8617260395176898e+04, 8.8626538861377921e+04, 8.8635817420988853e+04, 8.8645096074001005e+04, 8.8654374820405646e+04, 8.8663653660194046e+04, 8.8672932593357487e+04, 8.8682211619887254e+04, 8.8691490739774628e+04, 8.8700769953010880e+04, 8.8710049259587322e+04, 8.8719328659495237e+04, 8.8728608152725908e+04, 8.8737887739270620e+04, 8.8747167419120669e+04, 8.8756447192267369e+04, 8.8765727058701988e+04, 8.8775007018415854e+04, 8.8784287071400249e+04, 8.8793567217646458e+04, 8.8802847457145806e+04, 8.8812127789889608e+04, 8.8821408215869145e+04, 8.8830688735075746e+04, 8.8839969347500708e+04, 8.8849250053135343e+04, 8.8858530851970965e+04, 8.8867811743998886e+04, 8.8877092729210417e+04, 8.8886373807596901e+04, 8.8895654979149636e+04, 8.8904936243859949e+04, 8.8914217601719167e+04, 8.8923499052718602e+04, 8.8932780596849581e+04, 8.8942062234103447e+04, 8.8951343964471496e+04, 8.8960625787945100e+04, 8.8969907704515557e+04, 8.8979189714174237e+04, 8.8988471816912424e+04, 8.8997754012721503e+04, 8.9007036301592787e+04, 8.9016318683517631e+04, 8.9025601158487378e+04, 8.9034883726493354e+04, 8.9044166387526901e+04, 8.9053449141579389e+04, 8.9062731988642161e+04, 8.9072014928706543e+04, 8.9081297961763921e+04, 8.9090581087805636e+04, 8.9099864306823016e+04, 8.9109147618807459e+04, 8.9118431023750280e+04, 8.9127714521642876e+04, 8.9136998112476576e+04, 8.9146281796242780e+04, 8.9155565572932828e+04, 8.9164849442538107e+04, 8.9174133405049943e+04, 8.9183417460459750e+04, 8.9192701608758856e+04, 8.9201985849938676e+04, 8.9211270183990564e+04, 8.9220554610905892e+04, 8.9229839130676031e+04, 8.9239123743292381e+04, 8.9248408448746311e+04, 8.9257693247029209e+04, 8.9266978138132457e+04, 8.9276263122047429e+04, 8.9285548198765522e+04, 8.9294833368278123e+04, 8.9304118630576617e+04, 8.9313403985652389e+04, 8.9322689433496867e+04, 8.9331974974101395e+04, 8.9341260607457400e+04, 8.9350546333556282e+04, 8.9359832152389441e+04, 8.9369118063948263e+04, 8.9378404068224161e+04, 8.9387690165208536e+04, 8.9396976354892802e+04, 8.9406262637268359e+04, 8.9415549012326621e+04, 8.9424835480058988e+04, 8.9434122040456874e+04, 8.9443408693511694e+04, 8.9452695439214876e+04, 8.9461982277557807e+04, 8.9471269208531943e+04, 8.9480556232128700e+04, 8.9489843348339477e+04, 8.9499130557155702e+04, 8.9508417858568821e+04, 8.9517705252570246e+04, 8.9526992739151392e+04, 8.9536280318303703e+04, 8.9545567990018608e+04, 8.9554855754287550e+04, 8.9564143611101958e+04, 8.9573431560453260e+04, 8.9582719602332902e+04, 8.9592007736732339e+04, 8.9601295963642959e+04, 8.9610584283056261e+04, 8.9619872694963648e+04, 8.9629161199356604e+04, 8.9638449796226545e+04, 8.9647738485564943e+04, 8.9657027267363228e+04, 8.9666316141612842e+04, 8.9675605108305259e+04, 8.9684894167431950e+04, 8.9694183318984360e+04, 8.9703472562953932e+04, 8.9712761899332138e+04, 8.9722051328110436e+04, 8.9731340849280285e+04, 8.9740630462833156e+04, 8.9749920168760509e+04, 8.9759209967053815e+04, 8.9768499857704548e+04, 8.9777789840704179e+04, 8.9787079916044182e+04, 8.9796370083716029e+04, 8.9805660343711163e+04, 8.9814950696021115e+04, 8.9824241140637328e+04, 8.9833531677551306e+04, 8.9842822306754504e+04, 8.9852113028238426e+04, 8.9861403841994557e+04, 8.9870694748014357e+04, 8.9879985746289341e+04, 8.9889276836811012e+04, 8.9898568019570826e+04, 8.9907859294560301e+04, 8.9917150661770909e+04, 8.9926442121194166e+04, 8.9935733672821574e+04, 8.9945025316644620e+04, 8.9954317052654791e+04, 8.9963608880843618e+04, 8.9972900801202602e+04, 8.9982192813723232e+04, 8.9991484918397007e+04, 9.0000777115215446e+04, 9.0010069404170063e+04, 9.0019361785252375e+04, 9.0028654258453898e+04, 9.0037946823766135e+04, 9.0047239481180601e+04, 9.0056532230688812e+04, 9.0065825072282314e+04, 9.0075118005952594e+04, 9.0084411031691183e+04, 9.0093704149489611e+04, 9.0102997359339410e+04, 9.0112290661232109e+04, 9.0121584055159212e+04, 9.0130877541112292e+04, 9.0140171119082821e+04, 9.0149464789062404e+04, 9.0158758551042512e+04, 9.0168052405014707e+04, 9.0177346350970533e+04, 9.0186640388901535e+04, 9.0195934518799244e+04, 9.0205228740655206e+04, 9.0214523054460951e+04, 9.0223817460208054e+04, 9.0233111957888046e+04, 9.0242406547492472e+04, 9.0251701229012891e+04, 9.0260996002440836e+04, 9.0270290867767893e+04, 9.0279585824985596e+04, 9.0288880874085517e+04, 9.0298176015059187e+04, 9.0307471247898182e+04, 9.0316766572594046e+04, 9.0326061989138369e+04, 9.0335357497522724e+04, 9.0344653097738657e+04, 9.0353948789777729e+04, 9.0363244573631528e+04, 9.0372540449291613e+04, 9.0381836416749546e+04, 9.0391132475996914e+04, 9.0400428627025292e+04, 9.0409724869826270e+04, 9.0419021204391407e+04, 9.0428317630712307e+04, 9.0437614148780529e+04, 9.0446910758587648e+04, 9.0456207460125283e+04, 9.0465504253385006e+04, 9.0474801138358380e+04, 9.0484098115037021e+04, 9.0493395183412504e+04, 9.0502692343476418e+04, 9.0511989595220395e+04, 9.0521286938635996e+04, 9.0530584373714810e+04, 9.0539881900448454e+04, 9.0549179518828532e+04, 9.0558477228846648e+04, 9.0567775030494391e+04, 9.0577072923763379e+04, 9.0586370908645200e+04, 9.0595668985131488e+04, 9.0604967153213831e+04, 9.0614265412883833e+04, 9.0623563764133127e+04, 9.0632862206953330e+04, 9.0642160741336033e+04, 9.0651459367272881e+04, 9.0660758084755464e+04, 9.0670056893775429e+04, 9.0679355794324394e+04, 9.0688654786393963e+04, 9.0697953869975769e+04, 9.0707253045061429e+04, 9.0716552311642590e+04, 9.0725851669710901e+04, 9.0735151119257949e+04, 9.0744450660275383e+04, 9.0753750292754834e+04, 9.0763050016687965e+04, 9.0772349832066364e+04, 9.0781649738881708e+04, 9.0790949737125629e+04, 9.0800249826789732e+04, 9.0809550007865706e+04, 9.0818850280345170e+04, 9.0828150644219801e+04, 9.0837451099481215e+04, 9.0846751646121076e+04, 9.0856052284131030e+04, 9.0865353013502739e+04, 9.0874653834227851e+04, 9.0883954746298012e+04, 9.0893255749704898e+04, 9.0902556844440158e+04, 9.0911858030495423e+04, 9.0921159307862399e+04, 9.0930460676532719e+04, 9.0939762136498059e+04, 9.0949063687750095e+04, 9.0958365330280474e+04, 9.0967667064080888e+04, 9.0976968889142983e+04, 9.0986270805458451e+04, 9.0995572813018967e+04, 9.1004874911816194e+04, 9.1014177101841808e+04, 9.1023479383087499e+04, 9.1032781755544900e+04, 9.1042084219205761e+04, 9.1051386774061728e+04, 9.1060689420104478e+04, 9.1069992157325716e+04, 9.1079294985717148e+04, 9.1088597905270421e+04, 9.1097900915977239e+04, 9.1107204017829295e+04, 9.1116507210818279e+04, 9.1125810494935897e+04, 9.1135113870173838e+04, 9.1144417336523809e+04, 9.1153720893977486e+04, 9.1163024542526604e+04, 9.1172328282162838e+04, 9.1181632112877895e+04, 9.1190936034663493e+04, 9.1200240047511325e+04, 9.1209544151413109e+04, 9.1218848346360566e+04, 9.1228152632345387e+04, 9.1237457009359292e+04, 9.1246761477393986e+04, 9.1256066036441189e+04, 9.1265370686492635e+04, 9.1274675427540031e+04, 9.1283980259575095e+04, 9.1293285182589563e+04, 9.1302590196575125e+04, 9.1311895301523546e+04, 9.1321200497426544e+04, 9.1330505784275825e+04, 9.1339811162063139e+04, 9.1349116630780220e+04, 9.1358422190418787e+04, 9.1367727840970561e+04, 9.1377033582427292e+04, 9.1386339414780727e+04, 9.1395645338022601e+04, 9.1404951352144650e+04, 9.1414257457138607e+04, 9.1423563652996236e+04, 9.1432869939709257e+04, 9.1442176317269434e+04, 9.1451482785668515e+04, 9.1460789344898250e+04, 9.1470095994950345e+04, 9.1479402735816620e+04, 9.1488709567488782e+04, 9.1498016489958594e+04, 9.1507323503217820e+04, 9.1516630607258237e+04, 9.1525937802071581e+04, 9.1535245087649600e+04, 9.1544552463984088e+04, 9.1553859931066792e+04, 9.1563167488889478e+04, 9.1572475137443922e+04, 9.1581782876721903e+04, 9.1591090706715142e+04, 9.1600398627415445e+04, 9.1609706638814605e+04, 9.1619014740904357e+04, 9.1628322933676507e+04, 9.1637631217122806e+04, 9.1646939591235045e+04, 9.1656248056005032e+04, 9.1665556611424501e+04, 9.1674865257485275e+04, 9.1684173994179102e+04, 9.1693482821497790e+04, 9.1702791739433131e+04, 9.1712100747976903e+04, 9.1721409847120900e+04, 9.1730719036856928e+04, 9.1740028317176751e+04, 9.1749337688072192e+04, 9.1758647149535042e+04, 9.1767956701557079e+04, 9.1777266344130127e+04, 9.1786576077245976e+04, 9.1795885900896435e+04, 9.1805195815073312e+04, 9.1814505819768412e+04, 9.1823815914973529e+04, 9.1833126100680471e+04, 9.1842436376881058e+04, 9.1851746743567099e+04, 9.1861057200730414e+04, 9.1870367748362798e+04, 9.1879678386456086e+04, 9.1888989115002085e+04, 9.1898299933992632e+04, 9.1907610843419534e+04, 9.1916921843274598e+04, 9.1926232933549676e+04, 9.1935544114236574e+04, 9.1944855385327130e+04, 9.1954166746813164e+04, 9.1963478198686513e+04, 9.1972789740938999e+04, 9.1982101373562473e+04, 9.1991413096548742e+04, 9.2000724909889657e+04, 9.2010036813577040e+04, 9.2019348807602742e+04, 9.2028660891958600e+04, 9.2037973066636463e+04, 9.2047285331628169e+04, 9.2056597686925554e+04, 9.2065910132520468e+04, 9.2075222668404749e+04, 9.2084535294570291e+04, 9.2093848011008886e+04, 9.2103160817712414e+04, 9.2112473714672713e+04, 9.2121786701881632e+04, 9.2131099779331067e+04, 9.2140412947012839e+04, 9.2149726204918828e+04, 9.2159039553040886e+04, 9.2168352991370863e+04, 9.2177666519900624e+04, 9.2186980138622050e+04, 9.2196293847527006e+04, 9.2205607646607357e+04, 9.2214921535854955e+04, 9.2224235515261680e+04, 9.2233549584819426e+04, 9.2242863744520044e+04, 9.2252177994355414e+04, 9.2261492334317401e+04, 9.2270806764397901e+04, 9.2280121284588793e+04, 9.2289435894881928e+04, 9.2298750595269230e+04, 9.2308065385742564e+04, 9.2317380266293825e+04, 9.2326695236914864e+04, 9.2336010297597619e+04, 9.2345325448333970e+04, 9.2354640689115768e+04, 9.2363956019934922e+04, 9.2373271440783341e+04, 9.2382586951652935e+04, 9.2391902552535554e+04, 9.2401218243423151e+04, 9.2410534024307577e+04, 9.2419849895180756e+04, 9.2429165856034582e+04, 9.2438481906860994e+04, 9.2447798047651857e+04, 9.2457114278399094e+04, 9.2466430599094601e+04, 9.2475747009730316e+04, 9.2485063510298103e+04, 9.2494380100789931e+04, 9.2503696781197679e+04, 9.2513013551513286e+04, 9.2522330411728646e+04, 9.2531647361835683e+04, 9.2540964401826306e+04, 9.2550281531692468e+04, 9.2559598751426078e+04, 9.2568916061019059e+04, 9.2578233460463336e+04, 9.2587550949750817e+04, 9.2596868528873470e+04, 9.2606186197823205e+04, 9.2615503956591958e+04, 9.2624821805171639e+04, 9.2634139743554217e+04, 9.2643457771731613e+04, 9.2652775889695768e+04, 9.2662094097438603e+04, 9.2671412394952087e+04, 9.2680730782228158e+04, 9.2690049259258740e+04, 9.2699367826035785e+04, 9.2708686482551260e+04, 9.2718005228797061e+04, 9.2727324064765184e+04, 9.2736642990447595e+04, 9.2745962005836191e+04, 9.2755281110922951e+04, 9.2764600305699845e+04, 9.2773919590158825e+04, 9.2783238964291813e+04, 9.2792558428090822e+04, 9.2801877981547746e+04, 9.2811197624654640e+04, 9.2820517357403383e+04, 9.2829837179785987e+04, 9.2839157091794390e+04, 9.2848477093420574e+04, 9.2857797184656505e+04, 9.2867117365494152e+04, 9.2876437635925497e+04, 9.2885757995942506e+04, 9.2895078445537161e+04, 9.2904398984701445e+04, 9.2913719613427311e+04, 9.2923040331706754e+04, 9.2932361139531757e+04, 9.2941682036894301e+04, 9.2951003023786368e+04, 9.2960324100199941e+04, 9.2969645266127001e+04, 9.2978966521559531e+04, 9.2988287866489554e+04, 9.2997609300909040e+04, 9.3006930824809970e+04, 9.3016252438184340e+04, 9.3025574141024175e+04, 9.3034895933321444e+04, 9.3044217815068158e+04, 9.3053539786256311e+04, 9.3062861846877888e+04, 9.3072183996924912e+04, 9.3081506236389396e+04, 9.3090828565263320e+04, 9.3100150983538697e+04, 9.3109473491207536e+04, 9.3118796088261850e+04, 9.3128118774693663e+04, 9.3137441550494972e+04, 9.3146764415657788e+04, 9.3156087370174137e+04, 9.3165410414036029e+04, 9.3174733547235490e+04, 9.3184056769764517e+04, 9.3193380081615149e+04, 9.3202703482779398e+04, 9.3212026973249303e+04, 9.3221350553016906e+04, 9.3230674222074187e+04, 9.3239997980413187e+04, 9.3249321828025961e+04, 9.3258645764904504e+04, 9.3267969791040887e+04, 9.3277293906427120e+04, 9.3286618111055242e+04, 9.3295942404917310e+04, 9.3305266788005334e+04, 9.3314591260311354e+04, 9.3323915821827424e+04, 9.3333240472545585e+04, 9.3342565212457877e+04, 9.3351890041556340e+04, 9.3361214959833014e+04, 9.3370539967279969e+04, 9.3379865063889229e+04, 9.3389190249652878e+04, 9.3398515524562943e+04, 9.3407840888611492e+04, 9.3417166341790566e+04, 9.3426491884092233e+04, 9.3435817515508548e+04, 9.3445143236031552e+04, 9.3454469045653328e+04, 9.3463794944365945e+04, 9.3473120932161459e+04, 9.3482447009031908e+04, 9.3491773174969378e+04, 9.3501099429965936e+04, 9.3510425774013667e+04, 9.3519752207104641e+04, 9.3529078729230896e+04, 9.3538405340384546e+04, 9.3547732040557632e+04, 9.3557058829742251e+04, 9.3566385707930473e+04, 9.3575712675114395e+04, 9.3585039731286073e+04, 9.3594366876437591e+04, 9.3603694110561046e+04, 9.3613021433648508e+04, 9.3622348845692075e+04, 9.3631676346683817e+04, 9.3641003936615845e+04, 9.3650331615480231e+04, 9.3659659383269071e+04, 9.3668987239974478e+04, 9.3678315185588537e+04, 9.3687643220103331e+04, 9.3696971343510973e+04, 9.3706299555803533e+04, 9.3715627856973151e+04, 9.3724956247011898e+04, 9.3734284725911872e+04, 9.3743613293665199e+04, 9.3752941950263994e+04, 9.3762270695700354e+04, 9.3771599529966377e+04, 9.3780928453054163e+04, 9.3790257464955866e+04, 9.3799586565663572e+04, 9.3808915755169393e+04, 9.3818245033465428e+04, 9.3827574400543846e+04, 9.3836903856396704e+04, 9.3846233401016158e+04, 9.3855563034394319e+04, 9.3864892756523317e+04, 9.3874222567395263e+04, 9.3883552467002301e+04, 9.3892882455336527e+04, 9.3902212532390113e+04, 9.3911542698155157e+04, 9.3920872952623802e+04, 9.3930203295788175e+04, 9.3939533727640403e+04, 9.3948864248172627e+04, 9.3958194857376991e+04, 9.3967525555245636e+04, 9.3976856341770690e+04, 9.3986187216944294e+04, 9.3995518180758605e+04, 9.4004849233205721e+04, 9.4014180374277857e+04, 9.4023511603967112e+04, 9.4032842922265641e+04, 9.4042174329165602e+04, 9.4051505824659151e+04, 9.4060837408738400e+04, 9.4070169081395550e+04, 9.4079500842622743e+04, 9.4088832692412121e+04, 9.4098164630755840e+04, 9.4107496657646072e+04, 9.4116828773074973e+04, 9.4126160977034699e+04, 9.4135493269517421e+04, 9.4144825650515297e+04, 9.4154158120020511e+04, 9.4163490678025220e+04, 9.4172823324521582e+04, 9.4182156059501751e+04, 9.4191488882957958e+04, 9.4200821794882315e+04, 9.4210154795267037e+04, 9.4219487884104266e+04, 9.4228821061386188e+04, 9.4238154327105003e+04, 9.4247487681252882e+04, 9.4256821123822010e+04, 9.4266154654804530e+04, 9.4275488274192670e+04, 9.4284821981978603e+04, 9.4294155778154498e+04, 9.4303489662712585e+04, 9.4312823635645007e+04, 9.4322157696943963e+04, 9.4331491846601682e+04, 9.4340826084610293e+04, 9.4350160410962053e+04, 9.4359494825649119e+04, 9.4368829328663734e+04, 9.4378163919998042e+04, 9.4387498599644270e+04, 9.4396833367594620e+04, 9.4406168223841290e+04, 9.4415503168376483e+04, 9.4424838201192411e+04, 9.4434173322281291e+04, 9.4443508531635307e+04, 9.4452843829246689e+04, 9.4462179215107637e+04, 9.4471514689210366e+04, 9.4480850251547090e+04, 9.4490185902110024e+04, 9.4499521640891398e+04, 9.4508857467883427e+04, 9.4518193383078309e+04, 9.4527529386468304e+04, 9.4536865478045569e+04, 9.4546201657802376e+04, 9.4555537925730940e+04, 9.4564874281823490e+04, 9.4574210726072270e+04, 9.4583547258469480e+04, 9.4592883879007350e+04, 9.4602220587678137e+04, 9.4611557384474072e+04, 9.4620894269387354e+04, 9.4630231242410257e+04, 9.4639568303535023e+04, 9.4648905452753839e+04, 9.4658242690058993e+04, 9.4667580015442712e+04, 9.4676917428897243e+04, 9.4686254930414842e+04, 9.4695592519987738e+04, 9.4704930197608177e+04, 9.4714267963268416e+04, 9.4723605816960684e+04, 9.4732943758677255e+04, 9.4742281788410372e+04, 9.4751619906152293e+04, 9.4760958111895277e+04, 9.4770296405631554e+04, 9.4779634787353411e+04, 9.4788973257053105e+04, 9.4798311814722867e+04, 9.4807650460355013e+04, 9.4816989193941743e+04, 9.4826328015475359e+04, 9.4835666924948149e+04, 9.4845005922352328e+04, 9.4854345007680196e+04, 9.4863684180924014e+04, 9.4873023442076053e+04, 9.4882362791128588e+04, 9.4891702228073889e+04, 9.4901041752904261e+04, 9.4910381365611946e+04, 9.4919721066189231e+04, 9.4929060854628406e+04, 9.4938400730921712e+04, 9.4947740695061497e+04, 9.4957080747039989e+04, 9.4966420886849504e+04, 9.4975761114482317e+04, 9.4985101429930699e+04, 9.4994441833186982e+04, 9.5003782324243439e+04, 9.5013122903092342e+04, 9.5022463569725980e+04, 9.5031804324136683e+04, 9.5041145166316695e+04, 9.5050486096258377e+04, 9.5059827113954001e+04, 9.5069168219395855e+04, 9.5078509412576255e+04, 9.5087850693487489e+04, 9.5097192062121860e+04, 9.5106533518471697e+04, 9.5115875062529289e+04, 9.5125216694286937e+04, 9.5134558413736959e+04, 9.5143900220871670e+04, 9.5153242115683403e+04, 9.5162584098164429e+04, 9.5171926168307080e+04, 9.5181268326103673e+04, 9.5190610571546538e+04, 9.5199952904627964e+04, 9.5209295325340310e+04, 9.5218637833675850e+04, 9.5227980429626943e+04, 9.5237323113185907e+04, 9.5246665884345071e+04, 9.5256008743096740e+04, 9.5265351689433257e+04, 9.5274694723346955e+04, 9.5284037844830164e+04, 9.5293381053875230e+04, 9.5302724350474455e+04, 9.5312067734620185e+04, 9.5321411206304765e+04, 9.5330754765520542e+04, 9.5340098412259817e+04, 9.5349442146514964e+04, 9.5358785968278316e+04, 9.5368129877542218e+04, 9.5377473874299016e+04, 9.5386817958541025e+04, 9.5396162130260622e+04, 9.5405506389450165e+04, 9.5414850736101987e+04, 9.5424195170208433e+04, 9.5433539691761864e+04, 9.5442884300754653e+04, 9.5452228997179118e+04, 9.5461573781027633e+04, 9.5470918652292559e+04, 9.5480263610966271e+04, 9.5489608657041099e+04, 9.5498953790509418e+04, 9.5508299011363575e+04, 9.5517644319595958e+04, 9.5526989715198943e+04, 9.5536335198164859e+04, 9.5545680768486083e+04, 9.5555026426155004e+04, 9.5564372171163996e+04, 9.5573718003505404e+04, 9.5583063923171620e+04, 9.5592409930155030e+04, 9.5601756024447968e+04, 9.5611102206042866e+04, 9.5620448474932069e+04, 9.5629794831107967e+04, 9.5639141274562920e+04, 9.5648487805289347e+04, 9.5657834423279623e+04, 9.5667181128526121e+04, 9.5676527921021232e+04, 9.5685874800757345e+04, 9.5695221767726864e+04, 9.5704568821922148e+04, 9.5713915963335618e+04, 9.5723263191959646e+04, 9.5732610507786652e+04, 9.5741957910808997e+04, 9.5751305401019126e+04, 9.5760652978409402e+04, 9.5770000642972256e+04, 9.5779348394700050e+04, 9.5788696233585186e+04, 9.5798044159620127e+04, 9.5807392172797234e+04, 9.5816740273108910e+04, 9.5826088460547588e+04, 9.5835436735105643e+04, 9.5844785096775522e+04, 9.5854133545549630e+04, 9.5863482081420356e+04, 9.5872830704380132e+04, 9.5882179414421378e+04, 9.5891528211536483e+04, 9.5900877095717908e+04, 9.5910226066958057e+04, 9.5919575125249350e+04, 9.5928924270584204e+04, 9.5938273502955053e+04, 9.5947622822354300e+04, 9.5956972228774393e+04, 9.5966321722207780e+04, 9.5975671302646835e+04, 9.5985020970084035e+04, 9.5994370724511798e+04, 9.6003720565922544e+04, 9.6013070494308733e+04, 9.6022420509662799e+04, 9.6031770611977161e+04, 9.6041120801244251e+04, 9.6050471077456517e+04, 9.6059821440606407e+04, 9.6069171890686368e+04, 9.6078522427688848e+04, 9.6087873051606264e+04, 9.6097223762431080e+04, 9.6106574560155757e+04, 9.6115925444772714e+04, 9.6125276416274428e+04, 9.6134627474653360e+04, 9.6143978619901914e+04, 9.6153329852012597e+04, 9.6162681170977827e+04, 9.6172032576790079e+04, 9.6181384069441789e+04, 9.6190735648925474e+04, 9.6200087315233541e+04, 9.6209439068358479e+04, 9.6218790908292736e+04, 9.6228142835028761e+04, 9.6237494848559058e+04, 9.6246846948876075e+04, 9.6256199135972260e+04, 9.6265551409840133e+04, 9.6274903770472127e+04, 9.6284256217860733e+04, 9.6293608751998414e+04, 9.6302961372877646e+04, 9.6312314080490905e+04, 9.6321666874830684e+04, 9.6331019755889443e+04, 9.6340372723659661e+04, 9.6349725778133827e+04, 9.6359078919304433e+04, 9.6368432147163956e+04, 9.6377785461704872e+04, 9.6387138862919703e+04, 9.6396492350800880e+04, 9.6405845925340953e+04, 9.6415199586532384e+04, 9.6424553334367665e+04, 9.6433907168839272e+04, 9.6443261089939740e+04, 9.6452615097661517e+04, 9.6461969191997137e+04, 9.6471323372939092e+04, 9.6480677640479873e+04, 9.6490031994611985e+04, 9.6499386435327935e+04, 9.6508740962620228e+04, 9.6518095576481355e+04, 9.6527450276903837e+04, 9.6536805063880180e+04, 9.6546159937402903e+04, 9.6555514897464498e+04, 9.6564869944057486e+04, 9.6574225077174371e+04, 9.6583580296807675e+04, 9.6592935602949889e+04, 9.6602290995593576e+04, 9.6611646474731213e+04, 9.6621002040355350e+04, 9.6630357692458478e+04, 9.6639713431033146e+04, 9.6649069256071874e+04, 9.6658425167567169e+04, 9.6667781165511580e+04, 9.6677137249897612e+04, 9.6686493420717787e+04, 9.6695849677964681e+04, 9.6705206021630787e+04, 9.6714562451708625e+04, 9.6723918968190745e+04, 9.6733275571069695e+04, 9.6742632260338010e+04, 9.6751989035988212e+04, 9.6761345898012834e+04, 9.6770702846404456e+04, 9.6780059881155583e+04, 9.6789417002258764e+04, 9.6798774209706549e+04, 9.6808131503491473e+04, 9.6817488883606085e+04, 9.6826846350042964e+04, 9.6836203902794616e+04, 9.6845561541853618e+04, 9.6854919267212492e+04, 9.6864277078863830e+04, 9.6873634976800153e+04, 9.6882992961014039e+04, 9.6892351031498038e+04, 9.6901709188244713e+04, 9.6911067431246614e+04, 9.6920425760496306e+04, 9.6929784175986337e+04, 9.6939142677709300e+04, 9.6948501265657731e+04, 9.6957859939824237e+04, 9.6967218700201338e+04, 9.6976577546781613e+04, 9.6985936479557655e+04, 9.6995295498521999e+04, 9.7004654603667252e+04, 9.7014013794985964e+04, 9.7023373072470713e+04, 9.7032732436114078e+04, 9.7042091885908652e+04, 9.7051451421846999e+04, 9.7060811043921683e+04, 9.7070170752125327e+04, 9.7079530546450464e+04, 9.7088890426889717e+04, 9.7098250393435665e+04, 9.7107610446080871e+04, 9.7116970584817929e+04, 9.7126330809639447e+04, 9.7135691120538017e+04, 9.7145051517506203e+04, 9.7154412000536613e+04, 9.7163772569621840e+04, 9.7173133224754478e+04, 9.7182493965927119e+04, 9.7191854793132356e+04, 9.7201215706362811e+04, 9.7210576705611078e+04, 9.7219937790869735e+04, 9.7229298962131419e+04, 9.7238660219388708e+04, 9.7248021562634225e+04, 9.7257382991860562e+04, 9.7266744507060328e+04, 9.7276106108226144e+04, 9.7285467795350618e+04, 9.7294829568426358e+04, 9.7304191427445970e+04, 9.7313553372402079e+04, 9.7322915403287290e+04, 9.7332277520094241e+04, 9.7341639722815526e+04, 9.7351002011443779e+04, 9.7360364385971610e+04, 9.7369726846391655e+04, 9.7379089392696522e+04, 9.7388452024878847e+04, 9.7397814742931238e+04, 9.7407177546846360e+04, 9.7416540436616793e+04, 9.7425903412235188e+04, 9.7435266473694195e+04, 9.7444629620986409e+04, 9.7453992854104494e+04, 9.7463356173041073e+04, 9.7472719577788768e+04, 9.7482083068340246e+04, 9.7491446644688127e+04, 9.7500810306825049e+04, 9.7510174054743664e+04, 9.7519537888436593e+04, 9.7528901807896487e+04, 9.7538265813116013e+04, 9.7547629904087793e+04, 9.7556994080804478e+04, 9.7566358343258733e+04, 9.7575722691443181e+04, 9.7585087125350503e+04, 9.7594451644973320e+04, 9.7603816250304313e+04, 9.7613180941336119e+04, 9.7622545718061388e+04, 9.7631910580472802e+04, 9.7641275528563012e+04, 9.7650640562324654e+04, 9.7660005681750423e+04, 9.7669370886832941e+04, 9.7678736177564919e+04, 9.7688101553939006e+04, 9.7697467015947841e+04, 9.7706832563584117e+04, 9.7716198196840516e+04, 9.7725563915709659e+04, 9.7734929720184271e+04, 9.7744295610256973e+04, 9.7753661585920476e+04, 9.7763027647167459e+04, 9.7772393793990559e+04, 9.7781760026382486e+04, 9.7791126344335920e+04, 9.7800492747843527e+04, 9.7809859236897988e+04, 9.7819225811491997e+04, 9.7828592471618234e+04, 9.7837959217269366e+04, 9.7847326048438117e+04, 9.7856692965117138e+04, 9.7866059967299123e+04, 9.7875427054976768e+04, 9.7884794228142782e+04, 9.7894161486789846e+04, 9.7903528830910625e+04, 9.7912896260497844e+04, 9.7922263775544197e+04, 9.7931631376042365e+04, 9.7940999061985058e+04, 9.7950366833364984e+04, 9.7959734690174839e+04, 9.7969102632407317e+04, 9.7978470660055129e+04, 9.7987838773110969e+04, 9.7997206971567561e+04, 9.8006575255417600e+04, 9.8015943624653795e+04, 9.8025312079268871e+04, 9.8034680619255523e+04, 9.8044049244606460e+04, 9.8053417955314406e+04, 9.8062786751372070e+04, 9.8072155632772177e+04, 9.8081524599507437e+04, 9.8090893651570557e+04, 9.8100262788954278e+04, 9.8109632011651309e+04, 9.8119001319654388e+04, 9.8128370712956224e+04, 9.8137740191549543e+04, 9.8147109755427067e+04, 9.8156479404581522e+04, 9.8165849139005644e+04, 9.8175218958692174e+04, 9.8184588863633820e+04, 9.8193958853823322e+04, 9.8203328929253417e+04, 9.8212699089916845e+04, 9.8222069335806329e+04, 9.8231439666914594e+04, 9.8240810083234421e+04, 9.8250180584758520e+04, 9.8259551171479630e+04, 9.8268921843390504e+04, 9.8278292600483896e+04, 9.8287663442752528e+04, 9.8297034370189154e+04, 9.8306405382786514e+04, 9.8315776480537359e+04, 9.8325147663434458e+04, 9.8334518931470549e+04, 9.8343890284638357e+04, 9.8353261722930692e+04, 9.8362633246340250e+04, 9.8372004854859828e+04, 9.8381376548482163e+04, 9.8390748327200010e+04, 9.8400120191006150e+04, 9.8409492139893307e+04, 9.8418864173854294e+04, 9.8428236292881833e+04, 9.8437608496968707e+04, 9.8446980786107670e+04, 9.8456353160291503e+04, 9.8465725619512959e+04, 9.8475098163764807e+04, 9.8484470793039844e+04, 9.8493843507330806e+04, 9.8503216306630478e+04, 9.8512589190931641e+04, 9.8521962160227049e+04, 9.8531335214509498e+04, 9.8540708353771784e+04, 9.8550081578006677e+04, 9.8559454887206928e+04, 9.8568828281365350e+04, 9.8578201760474709e+04, 9.8587575324527774e+04, 9.8596948973517370e+04, 9.8606322707436237e+04, 9.8615696526277228e+04, 9.8625070430033054e+04, 9.8634444418696570e+04, 9.8643818492260514e+04, 9.8653192650717727e+04, 9.8662566894060976e+04, 9.8671941222283058e+04, 9.8681315635376799e+04, 9.8690690133334938e+04, 9.8700064716150300e+04, 9.8709439383815683e+04, 9.8718814136323927e+04, 9.8728188973667769e+04, 9.8737563895840081e+04, 9.8746938902833615e+04, 9.8756313994641197e+04, 9.8765689171255624e+04, 9.8775064432669707e+04, 9.8784439778876273e+04, 9.8793815209868117e+04, 9.8803190725638065e+04, 9.8812566326178901e+04, 9.8821942011483465e+04, 9.8831317781544567e+04, 9.8840693636355034e+04, 9.8850069575907648e+04, 9.8859445600195279e+04, 9.8868821709210708e+04, 9.8878197902946777e+04, 9.8887574181396296e+04, 9.8896950544552092e+04, 9.8906326992406990e+04, 9.8915703524953831e+04, 9.8925080142185427e+04, 9.8934456844094631e+04, 9.8943833630674242e+04, 9.8953210501917099e+04, 9.8962587457816044e+04, 9.8971964498363915e+04, 9.8981341623553526e+04, 9.8990718833377745e+04, 9.9000096127829369e+04, 9.9009473506901268e+04, 9.9018850970586282e+04, 9.9028228518877237e+04, 9.9037606151766988e+04, 9.9046983869248375e+04, 9.9056361671314226e+04, 9.9065739557957422e+04, 9.9075117529170777e+04, 9.9084495584947159e+04, 9.9093873725279424e+04, 9.9103251950160382e+04, 9.9112630259582933e+04, 9.9122008653539902e+04, 9.9131387132024160e+04, 9.9140765695028560e+04, 9.9150144342545944e+04, 9.9159523074569181e+04, 9.9168901891091140e+04, 9.9178280792104662e+04, 9.9187659777602617e+04, 9.9197038847577860e+04, 9.9206418002023274e+04, 9.9215797240931715e+04, 9.9225176564296053e+04, 9.9234555972109141e+04, 9.9243935464363865e+04, 9.9253315041053080e+04, 9.9262694702169669e+04, 9.9272074447706502e+04, 9.9281454277656463e+04, 9.9290834192012408e+04, 9.9300214190767234e+04, 9.9309594273913797e+04, 9.9318974441444981e+04, 9.9328354693353685e+04, 9.9337735029632764e+04, 9.9347115450275101e+04, 9.9356495955273611e+04, 9.9365876544621133e+04, 9.9375257218310580e+04, 9.9384637976334838e+04, 9.9394018818686804e+04, 9.9403399745359333e+04, 9.9412780756345368e+04, 9.9422161851637778e+04, 9.9431543031229434e+04, 9.9440924295113262e+04, 9.9450305643282147e+04, 9.9459687075728958e+04, 9.9469068592446638e+04, 9.9478450193428042e+04, 9.9487831878666111e+04, 9.9497213648153702e+04, 9.9506595501883770e+04, 9.9515977439849172e+04, 9.9525359462042834e+04, 9.9534741568457655e+04, 9.9544123759086564e+04, 9.9553506033922458e+04, 9.9562888392958223e+04, 9.9572270836186799e+04, 9.9581653363601086e+04, 9.9591035975193998e+04, 9.9600418670958461e+04, 9.9609801450887389e+04, 9.9619184314973681e+04, 9.9628567263210265e+04, 9.9637950295590053e+04, 9.9647333412105989e+04, 9.9656716612750984e+04, 9.9666099897517939e+04, 9.9675483266399795e+04, 9.9684866719389480e+04, 9.9694250256479922e+04, 9.9703633877664048e+04, 9.9713017582934772e+04, 9.9722401372285050e+04, 9.9731785245707768e+04, 9.9741169203195910e+04, 9.9750553244742390e+04, 9.9759937370340136e+04, 9.9769321579982105e+04, 9.9778705873661209e+04, 9.9788090251370391e+04, 9.9797474713102609e+04, 9.9806859258850774e+04, 9.9816243888607845e+04, 9.9825628602366764e+04, 9.9835013400120471e+04, 9.9844398281861926e+04, 9.9853783247584055e+04, 9.9863168297279830e+04, 9.9872553430942178e+04, 9.9881938648564043e+04, 9.9891323950138394e+04, 9.9900709335658205e+04, 9.9910094805116401e+04, 9.9919480358505942e+04, 9.9928865995819768e+04, 9.9938251717050851e+04, 9.9947637522192163e+04, 9.9957023411236631e+04, 9.9966409384177241e+04, 9.9975795441006965e+04, 9.9985181581718731e+04, 9.9994567806305538e+04, 1.0000395411476036e+05, 1.0001334050707611e+05, 1.0002272698324580e+05, 1.0003211354326237e+05, 1.0004150018711881e+05, 1.0005088691480808e+05, 1.0006027372632317e+05, 1.0006966062165705e+05, 1.0007904760080267e+05, 1.0008843466375304e+05, 1.0009782181050112e+05, 1.0010720904103990e+05, 1.0011659635536233e+05, 1.0012598375346142e+05, 1.0013537123533015e+05, 1.0014475880096148e+05, 1.0015414645034845e+05, 1.0016353418348398e+05, 1.0017292200036108e+05, 1.0018230990097274e+05, 1.0019169788531195e+05, 1.0020108595337172e+05, 1.0021047410514503e+05, 1.0021986234062485e+05, 1.0022925065980418e+05, 1.0023863906267603e+05, 1.0024802754923339e+05, 1.0025741611946926e+05, 1.0026680477337661e+05, 1.0027619351094851e+05, 1.0028558233217790e+05, 1.0029497123705779e+05, 1.0030436022558119e+05, 1.0031374929774110e+05, 1.0032313845353054e+05, 1.0033252769294252e+05, 1.0034191701597002e+05, 1.0035130642260607e+05, 1.0036069591284367e+05, 1.0037008548667585e+05, 1.0037947514409560e+05, 1.0038886488509594e+05, 1.0039825470966988e+05, 1.0040764461781047e+05, 1.0041703460951067e+05, 1.0042642468476354e+05, 1.0043581484356208e+05, 1.0044520508589933e+05, 1.0045459541176831e+05, 1.0046398582116203e+05, 1.0047337631407351e+05, 1.0048276689049578e+05, 1.0049215755042188e+05, 1.0050154829384483e+05, 1.0051093912075764e+05, 1.0052033003115338e+05, 1.0052972102502505e+05, 1.0053911210236569e+05, 1.0054850326316834e+05, 1.0055789450742603e+05, 1.0056728583513177e+05, 1.0057667724627865e+05, 1.0058606874085964e+05, 1.0059546031886786e+05, 1.0060485198029628e+05, 1.0061424372513800e+05, 1.0062363555338600e+05, 1.0063302746503337e+05, 1.0064241946007316e+05, 1.0065181153849840e+05, 1.0066120370030213e+05, 1.0067059594547740e+05, 1.0067998827401729e+05, 1.0068938068591479e+05, 1.0069877318116301e+05, 1.0070816575975496e+05, 1.0071755842168373e+05, 1.0072695116694235e+05, 1.0073634399552390e+05, 1.0074573690742144e+05, 1.0075512990262800e+05, 1.0076452298113666e+05, 1.0077391614294048e+05, 1.0078330938803250e+05, 1.0079270271640579e+05, 1.0080209612805345e+05, 1.0081148962296853e+05, 1.0082088320114407e+05, 1.0083027686257317e+05, 1.0083967060724889e+05, 1.0084906443516430e+05, 1.0085845834631249e+05, 1.0086785234068648e+05, 1.0087724641827941e+05, 1.0088664057908430e+05, 1.0089603482309423e+05, 1.0090542915030234e+05, 1.0091482356070164e+05, 1.0092421805428524e+05, 1.0093361263104621e+05, 1.0094300729097764e+05, 1.0095240203407261e+05, 1.0096179686032420e+05, 1.0097119176972551e+05, 1.0098058676226962e+05, 1.0098998183794961e+05, 1.0099937699675858e+05, 1.0100877223868961e+05, 1.0101816756373580e+05, 1.0102756297189023e+05, 1.0103695846314602e+05, 1.0104635403749623e+05, 1.0105574969493398e+05, 1.0106514543545237e+05, 1.0107454125904449e+05, 1.0108393716570344e+05, 1.0109333315542230e+05, 1.0110272922819419e+05, 1.0111212538401222e+05, 1.0112152162286948e+05, 1.0113091794475909e+05, 1.0114031434967414e+05, 1.0114971083760775e+05, 1.0115910740855301e+05, 1.0116850406250305e+05, 1.0117790079945097e+05, 1.0118729761938989e+05, 1.0119669452231292e+05, 1.0120609150821317e+05, 1.0121548857708374e+05, 1.0122488572891781e+05, 1.0123428296370844e+05, 1.0124368028144875e+05, 1.0125307768213186e+05, 1.0126247516575090e+05, 1.0127187273229900e+05, 1.0128127038176928e+05, 1.0129066811415486e+05, 1.0130006592944886e+05, 1.0130946382764445e+05, 1.0131886180873470e+05, 1.0132825987271276e+05, 1.0133765801957175e+05, 1.0134705624930482e+05, 1.0135645456190511e+05, 1.0136585295736573e+05, 1.0137525143567982e+05, 1.0138464999684053e+05, 1.0139404864084098e+05, 1.0140344736767431e+05, 1.0141284617733369e+05, 1.0142224506981220e+05, 1.0143164404510304e+05, 1.0144104310319931e+05, 1.0145044224409419e+05, 1.0145984146778079e+05, 1.0146924077425229e+05, 1.0147864016350180e+05, 1.0148803963552249e+05, 1.0149743919030752e+05, 1.0150683882785004e+05, 1.0151623854814317e+05, 1.0152563835118010e+05, 1.0153503823695396e+05, 1.0154443820545792e+05, 1.0155383825668512e+05, 1.0156323839062874e+05, 1.0157263860728192e+05, 1.0158203890663780e+05, 1.0159143928868958e+05, 1.0160083975343042e+05, 1.0161024030085346e+05, 1.0161964093095188e+05, 1.0162904164371885e+05, 1.0163844243914752e+05, 1.0164784331723106e+05, 1.0165724427796266e+05, 1.0166664532133545e+05, 1.0167604644734263e+05, 1.0168544765597736e+05, 1.0169484894723284e+05, 1.0170425032110223e+05, 1.0171365177757868e+05, 1.0172305331665539e+05, 1.0173245493832554e+05, 1.0174185664258229e+05, 1.0175125842941883e+05, 1.0176066029882834e+05, 1.0177006225080404e+05, 1.0177946428533904e+05, 1.0178886640242657e+05, 1.0179826860205980e+05, 1.0180767088423195e+05, 1.0181707324893618e+05, 1.0182647569616568e+05, 1.0183587822591365e+05, 1.0184528083817328e+05, 1.0185468353293774e+05, 1.0186408631020025e+05, 1.0187348916995400e+05, 1.0188289211219217e+05, 1.0189229513690798e+05, 1.0190169824409462e+05, 1.0191110143374528e+05, 1.0192050470585315e+05, 1.0192990806041146e+05, 1.0193931149741338e+05, 1.0194871501685216e+05, 1.0195811861872097e+05, 1.0196752230301301e+05, 1.0197692606972151e+05, 1.0198632991883966e+05, 1.0199573385036069e+05, 1.0200513786427779e+05, 1.0201454196058419e+05, 1.0202394613927309e+05, 1.0203335040033770e+05, 1.0204275474377125e+05, 1.0205215916956693e+05, 1.0206156367771798e+05, 1.0207096826821764e+05, 1.0208037294105906e+05, 1.0208977769623553e+05, 1.0209918253374023e+05, 1.0210858745356640e+05, 1.0211799245570725e+05, 1.0212739754015603e+05, 1.0213680270690595e+05, 1.0214620795595022e+05, 1.0215561328728209e+05, 1.0216501870089478e+05, 1.0217442419678152e+05, 1.0218382977493556e+05, 1.0219323543535011e+05, 1.0220264117801841e+05, 1.0221204700293372e+05, 1.0222145291008921e+05, 1.0223085889947820e+05, 1.0224026497109390e+05, 1.0224967112492952e+05, 1.0225907736097831e+05, 1.0226848367923353e+05, 1.0227789007968840e+05, 1.0228729656233618e+05, 1.0229670312717011e+05, 1.0230610977418347e+05, 1.0231551650336944e+05, 1.0232492331472132e+05, 1.0233433020823234e+05, 1.0234373718389575e+05, 1.0235314424170481e+05, 1.0236255138165278e+05, 1.0237195860373287e+05, 1.0238136590793838e+05, 1.0239077329426256e+05, 1.0240018076269864e+05, 1.0240958831323989e+05, 1.0241899594587959e+05, 1.0242840366061099e+05, 1.0243781145742732e+05, 1.0244721933632190e+05, 1.0245662729728797e+05, 1.0246603534031876e+05, 1.0247544346540759e+05, 1.0248485167254769e+05, 1.0249425996173236e+05, 1.0250366833295483e+05, 1.0251307678620839e+05, 1.0252248532148630e+05, 1.0253189393878187e+05, 1.0254130263808832e+05, 1.0255071141939897e+05, 1.0256012028270707e+05, 1.0256952922800591e+05, 1.0257893825528875e+05, 1.0258834736454891e+05, 1.0259775655577962e+05, 1.0260716582897420e+05, 1.0261657518412590e+05, 1.0262598462122804e+05, 1.0263539414027387e+05, 1.0264480374125669e+05, 1.0265421342416979e+05, 1.0266362318900647e+05, 1.0267303303576000e+05, 1.0268244296442367e+05, 1.0269185297499078e+05, 1.0270126306745461e+05, 1.0271067324180847e+05, 1.0272008349804563e+05, 1.0272949383615941e+05, 1.0273890425614311e+05, 1.0274831475799003e+05, 1.0275772534169343e+05, 1.0276713600724663e+05, 1.0277654675464294e+05, 1.0278595758387567e+05, 1.0279536849493810e+05, 1.0280477948782356e+05, 1.0281419056252533e+05, 1.0282360171903673e+05, 1.0283301295735108e+05, 1.0284242427746167e+05, 1.0285183567936181e+05, 1.0286124716304483e+05, 1.0287065872850403e+05, 1.0288007037573270e+05, 1.0288948210472420e+05, 1.0289889391547181e+05, 1.0290830580796886e+05, 1.0291771778220867e+05, 1.0292712983818454e+05, 1.0293654197588981e+05, 1.0294595419531780e+05, 1.0295536649646182e+05, 1.0296477887931518e+05, 1.0297419134387125e+05, 1.0298360389012333e+05, 1.0299301651806473e+05, 1.0300242922768880e+05, 1.0301184201898886e+05, 1.0302125489195826e+05, 1.0303066784659030e+05, 1.0304008088287832e+05, 1.0304949400081566e+05, 1.0305890720039564e+05, 1.0306832048161160e+05, 1.0307773384445692e+05, 1.0308714728892487e+05, 1.0309656081500881e+05, 1.0310597442270209e+05, 1.0311538811199807e+05, 1.0312480188289005e+05, 1.0313421573537140e+05, 1.0314362966943545e+05, 1.0315304368507552e+05, 1.0316245778228501e+05, 1.0317187196105724e+05, 1.0318128622138557e+05, 1.0319070056326334e+05, 1.0320011498668390e+05, 1.0320952949164060e+05, 1.0321894407812678e+05, 1.0322835874613581e+05, 1.0323777349566107e+05, 1.0324718832669585e+05, 1.0325660323923355e+05, 1.0326601823326752e+05, 1.0327543330879112e+05, 1.0328484846579771e+05, 1.0329426370428066e+05, 1.0330367902423332e+05, 1.0331309442564906e+05, 1.0332250990852124e+05, 1.0333192547284323e+05, 1.0334134111860838e+05, 1.0335075684581009e+05, 1.0336017265444169e+05, 1.0336958854449658e+05, 1.0337900451596813e+05, 1.0338842056884969e+05, 1.0339783670313464e+05, 1.0340725291881638e+05, 1.0341666921588825e+05, 1.0342608559434363e+05, 1.0343550205417591e+05, 1.0344491859537848e+05, 1.0345433521794470e+05, 1.0346375192186795e+05, 1.0347316870714162e+05, 1.0348258557375909e+05, 1.0349200252171376e+05, 1.0350141955099898e+05, 1.0351083666160816e+05, 1.0352025385353467e+05, 1.0352967112677191e+05, 1.0353908848131326e+05, 1.0354850591715213e+05, 1.0355792343428188e+05, 1.0356734103269594e+05, 1.0357675871238769e+05, 1.0358617647335051e+05, 1.0359559431557781e+05, 1.0360501223906298e+05, 1.0361443024379942e+05, 1.0362384832978054e+05, 1.0363326649699973e+05, 1.0364268474545039e+05, 1.0365210307512589e+05, 1.0366152148601970e+05, 1.0367093997812517e+05, 1.0368035855143571e+05, 1.0368977720594476e+05, 1.0369919594164571e+05, 1.0370861475853196e+05, 1.0371803365659693e+05, 1.0372745263583402e+05, 1.0373687169623665e+05, 1.0374629083779821e+05, 1.0375571006051215e+05, 1.0376512936437188e+05, 1.0377454874937078e+05, 1.0378396821550230e+05, 1.0379338776275984e+05, 1.0380280739113686e+05, 1.0381222710062673e+05, 1.0382164689122289e+05, 1.0383106676291875e+05, 1.0384048671570774e+05, 1.0384990674958332e+05, 1.0385932686453887e+05, 1.0386874706056784e+05, 1.0387816733766363e+05, 1.0388758769581970e+05, 1.0389700813502946e+05, 1.0390642865528636e+05, 1.0391584925658382e+05, 1.0392526993891528e+05, 1.0393469070227417e+05, 1.0394411154665390e+05, 1.0395353247204794e+05, 1.0396295347844971e+05, 1.0397237456585266e+05, 1.0398179573425023e+05, 1.0399121698363587e+05, 1.0400063831400296e+05, 1.0401005972534501e+05, 1.0401948121765546e+05, 1.0402890279092772e+05, 1.0403832444515525e+05, 1.0404774618033151e+05, 1.0405716799644993e+05, 1.0406658989350397e+05, 1.0407601187148709e+05, 1.0408543393039271e+05, 1.0409485607021430e+05, 1.0410427829094531e+05, 1.0411370059257920e+05, 1.0412312297510942e+05, 1.0413254543852941e+05, 1.0414196798283268e+05, 1.0415139060801263e+05, 1.0416081331406275e+05, 1.0417023610097650e+05, 1.0417965896874733e+05, 1.0418908191736871e+05, 1.0419850494683409e+05, 1.0420792805713695e+05, 1.0421735124827073e+05, 1.0422677452022894e+05, 1.0423619787300505e+05, 1.0424562130659247e+05, 1.0425504482098474e+05, 1.0426446841617528e+05, 1.0427389209215758e+05, 1.0428331584892511e+05, 1.0429273968647134e+05, 1.0430216360478978e+05, 1.0431158760387385e+05, 1.0432101168371706e+05, 1.0433043584431290e+05, 1.0433986008565483e+05, 1.0434928440773633e+05, 1.0435870881055089e+05, 1.0436813329409198e+05, 1.0437755785835310e+05, 1.0438698250332772e+05, 1.0439640722900935e+05, 1.0440583203539143e+05, 1.0441525692246750e+05, 1.0442468189023102e+05, 1.0443410693867550e+05, 1.0444353206779440e+05, 1.0445295727758126e+05, 1.0446238256802950e+05, 1.0447180793913269e+05, 1.0448123339088428e+05, 1.0449065892327778e+05, 1.0450008453630669e+05, 1.0450951022996449e+05, 1.0451893600424471e+05, 1.0452836185914083e+05, 1.0453778779464636e+05, 1.0454721381075481e+05, 1.0455663990745966e+05, 1.0456606608475443e+05, 1.0457549234263264e+05, 1.0458491868108777e+05, 1.0459434510011334e+05, 1.0460377159970283e+05, 1.0461319817984980e+05, 1.0462262484054775e+05, 1.0463205158179019e+05, 1.0464147840357061e+05, 1.0465090530588254e+05, 1.0466033228871950e+05, 1.0466975935207499e+05, 1.0467918649594254e+05, 1.0468861372031565e+05, 1.0469804102518788e+05, 1.0470746841055271e+05, 1.0471689587640367e+05, 1.0472632342273429e+05, 1.0473575104953808e+05, 1.0474517875680861e+05, 1.0475460654453936e+05, 1.0476403441272386e+05, 1.0477346236135566e+05, 1.0478289039042826e+05, 1.0479231849993521e+05, 1.0480174668987004e+05, 1.0481117496022627e+05, 1.0482060331099745e+05, 1.0483003174217712e+05, 1.0483946025375879e+05, 1.0484888884573600e+05, 1.0485831751810228e+05, 1.0486774627085120e+05, 1.0487717510397627e+05, 1.0488660401747105e+05, 1.0489603301132907e+05, 1.0490546208554387e+05, 1.0491489124010899e+05, 1.0492432047501799e+05, 1.0493374979026440e+05, 1.0494317918584179e+05, 1.0495260866174367e+05, 1.0496203821796362e+05, 1.0497146785449519e+05, 1.0498089757133191e+05, 1.0499032736846733e+05, 1.0499975724589503e+05, 1.0500918720360856e+05, 1.0501861724160145e+05, 1.0502804735986727e+05, 1.0503747755839957e+05, 1.0504690783719190e+05, 1.0505633819623785e+05, 1.0506576863553096e+05, 1.0507519915506477e+05, 1.0508462975483289e+05, 1.0509406043482885e+05, 1.0510349119504621e+05, 1.0511292203547855e+05, 1.0512235295611943e+05, 1.0513178395696242e+05, 1.0514121503800110e+05, 1.0515064619922901e+05, 1.0516007744063975e+05, 1.0516950876222688e+05, 1.0517894016398396e+05, 1.0518837164590460e+05, 1.0519780320798232e+05, 1.0520723485021073e+05, 1.0521666657258340e+05, 1.0522609837509390e+05, 1.0523553025773582e+05, 1.0524496222050276e+05, 1.0525439426338824e+05, 1.0526382638638590e+05, 1.0527325858948928e+05, 1.0528269087269199e+05, 1.0529212323598762e+05, 1.0530155567936972e+05, 1.0531098820283193e+05, 1.0532042080636778e+05, 1.0532985348997089e+05, 1.0533928625363485e+05, 1.0534871909735326e+05, 1.0535815202111969e+05, 1.0536758502492774e+05, 1.0537701810877102e+05, 1.0538645127264310e+05, 1.0539588451653760e+05, 1.0540531784044810e+05, 1.0541475124436819e+05, 1.0542418472829147e+05, 1.0543361829221157e+05, 1.0544305193612205e+05, 1.0545248566001654e+05, 1.0546191946388864e+05, 1.0547135334773194e+05, 1.0548078731154009e+05, 1.0549022135530662e+05, 1.0549965547902521e+05, 1.0550908968268943e+05, 1.0551852396629289e+05, 1.0552795832982921e+05, 1.0553739277329201e+05, 1.0554682729667488e+05, 1.0555626189997143e+05, 1.0556569658317530e+05, 1.0557513134628009e+05, 1.0558456618927942e+05, 1.0559400111216689e+05, 1.0560343611493615e+05, 1.0561287119758080e+05, 1.0562230636009446e+05, 1.0563174160247078e+05, 1.0564117692470334e+05, 1.0565061232678578e+05, 1.0566004780871174e+05, 1.0566948337047482e+05, 1.0567891901206865e+05, 1.0568835473348688e+05, 1.0569779053472311e+05, 1.0570722641577102e+05, 1.0571666237662418e+05, 1.0572609841727624e+05, 1.0573553453772086e+05, 1.0574497073795165e+05, 1.0575440701796225e+05, 1.0576384337774629e+05, 1.0577327981729741e+05, 1.0578271633660923e+05, 1.0579215293567543e+05, 1.0580158961448963e+05, 1.0581102637304545e+05, 1.0582046321133657e+05, 1.0582990012935660e+05, 1.0583933712709919e+05, 1.0584877420455799e+05, 1.0585821136172666e+05, 1.0586764859859883e+05, 1.0587708591516814e+05, 1.0588652331142825e+05, 1.0589596078737281e+05, 1.0590539834299548e+05, 1.0591483597828988e+05, 1.0592427369324971e+05, 1.0593371148786858e+05, 1.0594314936214018e+05, 1.0595258731605812e+05, 1.0596202534961612e+05, 1.0597146346280778e+05, 1.0598090165562680e+05, 1.0599033992806682e+05, 1.0599977828012152e+05, 1.0600921671178454e+05, 1.0601865522304953e+05, 1.0602809381391019e+05, 1.0603753248436017e+05, 1.0604697123439315e+05, 1.0605641006400275e+05, 1.0606584897318268e+05, 1.0607528796192660e+05, 1.0608472703022818e+05, 1.0609416617808110e+05, 1.0610360540547901e+05, 1.0611304471241561e+05, 1.0612248409888455e+05, 1.0613192356487951e+05, 1.0614136311039419e+05, 1.0615080273542223e+05, 1.0616024243995734e+05, 1.0616968222399318e+05, 1.0617912208752343e+05, 1.0618856203054178e+05, 1.0619800205304190e+05, 1.0620744215501749e+05, 1.0621688233646222e+05, 1.0622632259736977e+05, 1.0623576293773386e+05, 1.0624520335754816e+05, 1.0625464385680633e+05, 1.0626408443550211e+05, 1.0627352509362914e+05, 1.0628296583118115e+05, 1.0629240664815179e+05, 1.0630184754453480e+05, 1.0631128852032385e+05, 1.0632072957551264e+05, 1.0633017071009486e+05, 1.0633961192406420e+05, 1.0634905321741437e+05, 1.0635849459013909e+05, 1.0636793604223202e+05, 1.0637737757368690e+05, 1.0638681918449741e+05, 1.0639626087465724e+05, 1.0640570264416012e+05, 1.0641514449299975e+05, 1.0642458642116982e+05, 1.0643402842866407e+05, 1.0644347051547616e+05, 1.0645291268159983e+05, 1.0646235492702879e+05, 1.0647179725175675e+05, 1.0648123965577740e+05, 1.0649068213908449e+05, 1.0650012470167171e+05, 1.0650956734353278e+05, 1.0651901006466142e+05, 1.0652845286505135e+05, 1.0653789574469629e+05, 1.0654733870358994e+05, 1.0655678174172602e+05, 1.0656622485909826e+05, 1.0657566805570041e+05, 1.0658511133152616e+05, 1.0659455468656926e+05, 1.0660399812082339e+05, 1.0661344163428232e+05, 1.0662288522693976e+05, 1.0663232889878945e+05, 1.0664177264982510e+05, 1.0665121648004044e+05, 1.0666066038942922e+05, 1.0667010437798516e+05, 1.0667954844570199e+05, 1.0668899259257344e+05, 1.0669843681859327e+05, 1.0670788112375519e+05, 1.0671732550805295e+05, 1.0672676997148027e+05, 1.0673621451403093e+05, 1.0674565913569863e+05, 1.0675510383647712e+05, 1.0676454861636015e+05, 1.0677399347534146e+05, 1.0678343841341478e+05, 1.0679288343057389e+05, 1.0680232852681250e+05, 1.0681177370212438e+05, 1.0682121895650326e+05, 1.0683066428994293e+05, 1.0684010970243708e+05, 1.0684955519397948e+05, 1.0685900076456391e+05, 1.0686844641418410e+05, 1.0687789214283381e+05, 1.0688733795050679e+05, 1.0689678383719678e+05, 1.0690622980289758e+05, 1.0691567584760291e+05, 1.0692512197130655e+05, 1.0693456817400224e+05, 1.0694401445568376e+05, 1.0695346081634487e+05, 1.0696290725597933e+05, 1.0697235377458089e+05, 1.0698180037214333e+05, 1.0699124704866040e+05, 1.0700069380412588e+05, 1.0701014063853353e+05, 1.0701958755187714e+05, 1.0702903454415046e+05, 1.0703848161534726e+05, 1.0704792876546131e+05, 1.0705737599448640e+05, 1.0706682330241628e+05, 1.0707627068924473e+05, 1.0708571815496555e+05, 1.0709516569957249e+05, 1.0710461332305934e+05, 1.0711406102541985e+05, 1.0712350880664786e+05, 1.0713295666673710e+05, 1.0714240460568135e+05, 1.0715185262347442e+05, 1.0716130072011007e+05, 1.0717074889558210e+05, 1.0718019714988430e+05, 1.0718964548301045e+05, 1.0719909389495432e+05, 1.0720854238570972e+05, 1.0721799095527042e+05, 1.0722743960363024e+05, 1.0723688833078292e+05, 1.0724633713672232e+05, 1.0725578602144217e+05, 1.0726523498493631e+05, 1.0727468402719851e+05, 1.0728413314822258e+05, 1.0729358234800231e+05, 1.0730303162653151e+05, 1.0731248098380395e+05, 1.0732193041981345e+05, 1.0733137993455380e+05, 1.0734082952801882e+05, 1.0735027920020232e+05, 1.0735972895109808e+05, 1.0736917878069991e+05, 1.0737862868900161e+05, 1.0738807867599699e+05, 1.0739752874167988e+05, 1.0740697888604406e+05, 1.0741642910908334e+05, 1.0742587941079157e+05, 1.0743532979116253e+05, 1.0744478025019003e+05, 1.0745423078786790e+05, 1.0746368140418993e+05, 1.0747313209914995e+05, 1.0748258287274178e+05, 1.0749203372495924e+05, 1.0750148465579613e+05, 1.0751093566524629e+05, 1.0752038675330354e+05, 1.0752983791996166e+05, 1.0753928916521452e+05, 1.0754874048905594e+05, 1.0755819189147973e+05, 1.0756764337247972e+05, 1.0757709493204973e+05, 1.0758654657018361e+05, 1.0759599828687515e+05, 1.0760545008211822e+05, 1.0761490195590661e+05, 1.0762435390823419e+05, 1.0763380593909476e+05, 1.0764325804848218e+05, 1.0765271023639025e+05, 1.0766216250281285e+05, 1.0767161484774378e+05, 1.0768106727117690e+05, 1.0769051977310602e+05, 1.0769997235352501e+05, 1.0770942501242767e+05, 1.0771887774980789e+05, 1.0772833056565946e+05, 1.0773778345997629e+05, 1.0774723643275216e+05, 1.0775668948398094e+05, 1.0776614261365647e+05, 1.0777559582177261e+05, 1.0778504910832319e+05, 1.0779450247330207e+05, 1.0780395591670311e+05, 1.0781340943852015e+05, 1.0782286303874702e+05, 1.0783231671737757e+05, 1.0784177047440571e+05, 1.0785122430982525e+05, 1.0786067822363005e+05, 1.0787013221581395e+05, 1.0787958628637085e+05, 1.0788904043529458e+05, 1.0789849466257900e+05, 1.0790794896821798e+05, 1.0791740335220538e+05, 1.0792685781453506e+05, 1.0793631235520088e+05, 1.0794576697419670e+05, 1.0795522167151638e+05, 1.0796467644715380e+05, 1.0797413130110283e+05, 1.0798358623335732e+05, 1.0799304124391115e+05, 1.0800249633275819e+05, 1.0801195149989234e+05, 1.0802140674530741e+05, 1.0803086206899730e+05, 1.0804031747095588e+05, 1.0804977295117704e+05, 1.0805922850965465e+05, 1.0806868414638260e+05, 1.0807813986135474e+05, 1.0808759565456495e+05, 1.0809705152600713e+05, 1.0810650747567513e+05, 1.0811596350356285e+05, 1.0812541960966418e+05, 1.0813487579397298e+05, 1.0814433205648318e+05, 1.0815378839718860e+05, 1.0816324481608318e+05, 1.0817270131316078e+05, 1.0818215788841530e+05, 1.0819161454184062e+05, 1.0820107127343064e+05, 1.0821052808317923e+05, 1.0821998497108030e+05, 1.0822944193712773e+05, 1.0823889898131545e+05, 1.0824835610363730e+05, 1.0825781330408722e+05, 1.0826727058265907e+05, 1.0827672793934676e+05, 1.0828618537414422e+05, 1.0829564288704532e+05, 1.0830510047804394e+05, 1.0831455814713403e+05, 1.0832401589430944e+05, 1.0833347371956412e+05, 1.0834293162289196e+05, 1.0835238960428686e+05, 1.0836184766374271e+05, 1.0837130580125345e+05, 1.0838076401681296e+05, 1.0839022231041515e+05, 1.0839968068205395e+05, 1.0840913913172326e+05, 1.0841859765941699e+05, 1.0842805626512905e+05, 1.0843751494885336e+05, 1.0844697371058383e+05, 1.0845643255031439e+05, 1.0846589146803894e+05, 1.0847535046375140e+05, 1.0848480953744569e+05, 1.0849426868911572e+05, 1.0850372791875544e+05, 1.0851318722635873e+05, 1.0852264661191954e+05, 1.0853210607543179e+05, 1.0854156561688939e+05, 1.0855102523628628e+05, 1.0856048493361639e+05, 1.0856994470887365e+05, 1.0857940456205196e+05, 1.0858886449314527e+05, 1.0859832450214749e+05, 1.0860778458905258e+05, 1.0861724475385444e+05, 1.0862670499654702e+05, 1.0863616531712427e+05, 1.0864562571558010e+05, 1.0865508619190846e+05, 1.0866454674610325e+05, 1.0867400737815847e+05, 1.0868346808806801e+05, 1.0869292887582582e+05, 1.0870238974142585e+05, 1.0871185068486205e+05, 1.0872131170612831e+05, 1.0873077280521863e+05, 1.0874023398212694e+05, 1.0874969523684718e+05, 1.0875915656937328e+05, 1.0876861797969922e+05, 1.0877807946781891e+05, 1.0878754103372635e+05, 1.0879700267741544e+05, 1.0880646439888015e+05, 1.0881592619811442e+05, 1.0882538807511221e+05, 1.0883485002986748e+05, 1.0884431206237417e+05, 1.0885377417262625e+05, 1.0886323636061768e+05, 1.0887269862634239e+05, 1.0888216096979435e+05, 1.0889162339096754e+05, 1.0890108588985589e+05, 1.0891054846645339e+05, 1.0892001112075396e+05, 1.0892947385275162e+05, 1.0893893666244029e+05, 1.0894839954981394e+05, 1.0895786251486656e+05, 1.0896732555759208e+05, 1.0897678867798451e+05, 1.0898625187603777e+05, 1.0899571515174588e+05, 1.0900517850510277e+05, 1.0901464193610243e+05, 1.0902410544473885e+05, 1.0903356903100594e+05, 1.0904303269489775e+05, 1.0905249643640820e+05, 1.0906196025553127e+05, 1.0907142415226098e+05, 1.0908088812659126e+05, 1.0909035217851613e+05, 1.0909981630802953e+05, 1.0910928051512547e+05, 1.0911874479979790e+05, 1.0912820916204086e+05, 1.0913767360184828e+05, 1.0914713811921416e+05, 1.0915660271413247e+05, 1.0916606738659722e+05, 1.0917553213660239e+05, 1.0918499696414196e+05, 1.0919446186920992e+05, 1.0920392685180026e+05, 1.0921339191190698e+05, 1.0922285704952406e+05, 1.0923232226464553e+05, 1.0924178755726533e+05, 1.0925125292737749e+05, 1.0926071837497597e+05, 1.0927018390005481e+05, 1.0927964950260799e+05, 1.0928911518262950e+05, 1.0929858094011334e+05, 1.0930804677505353e+05, 1.0931751268744403e+05, 1.0932697867727886e+05, 1.0933644474455206e+05, 1.0934591088925760e+05, 1.0935537711138948e+05, 1.0936484341094171e+05, 1.0937430978790832e+05, 1.0938377624228327e+05, 1.0939324277406061e+05, 1.0940270938323432e+05, 1.0941217606979844e+05, 1.0942164283374697e+05, 1.0943110967507392e+05, 1.0944057659377328e+05, 1.0945004358983912e+05, 1.0945951066326541e+05, 1.0946897781404617e+05, 1.0947844504217543e+05, 1.0948791234764720e+05, 1.0949737973045549e+05, 1.0950684719059433e+05, 1.0951631472805775e+05, 1.0952578234283975e+05, 1.0953525003493436e+05, 1.0954471780433564e+05, 1.0955418565103755e+05, 1.0956365357503416e+05, 1.0957312157631945e+05, 1.0958258965488750e+05, 1.0959205781073232e+05, 1.0960152604384793e+05, 1.0961099435422837e+05, 1.0962046274186767e+05, 1.0962993120675985e+05, 1.0963939974889894e+05, 1.0964886836827900e+05, 1.0965833706489405e+05, 1.0966780583873813e+05, 1.0967727468980527e+05, 1.0968674361808949e+05, 1.0969621262358486e+05, 1.0970568170628537e+05, 1.0971515086618513e+05, 1.0972462010327812e+05, 1.0973408941755841e+05, 1.0974355880902005e+05, 1.0975302827765705e+05, 1.0976249782346349e+05, 1.0977196744643341e+05, 1.0978143714656083e+05, 1.0979090692383981e+05, 1.0980037677826440e+05, 1.0980984670982866e+05, 1.0981931671852662e+05, 1.0982878680435236e+05, 1.0983825696729989e+05, 1.0984772720736331e+05, 1.0985719752453665e+05, 1.0986666791881395e+05, 1.0987613839018930e+05, 1.0988560893865672e+05, 1.0989507956421030e+05, 1.0990455026684407e+05, 1.0991402104655211e+05, 1.0992349190332847e+05, 1.0993296283716722e+05, 1.0994243384806238e+05, 1.0995190493600807e+05, 1.0996137610099834e+05, 1.0997084734302723e+05, 1.0998031866208882e+05, 1.0998979005817718e+05, 1.0999926153128636e+05, 1.1000873308141045e+05, 1.1001820470854352e+05, 1.1002767641267962e+05, 1.1003714819381286e+05, 1.1004662005193725e+05, 1.1005609198704691e+05, 1.1006556399913591e+05, 1.1007503608819831e+05, 1.1008450825422820e+05, 1.1009398049721964e+05, 1.1010345281716672e+05, 1.1011292521406351e+05, 1.1012239768790410e+05, 1.1013187023868256e+05, 1.1014134286639297e+05, 1.1015081557102942e+05, 1.1016028835258598e+05, 1.1016976121105674e+05, 1.1017923414643579e+05, 1.1018870715871721e+05, 1.1019818024789510e+05, 1.1020765341396352e+05, 1.1021712665691657e+05, 1.1022659997674836e+05, 1.1023607337345295e+05, 1.1024554684702445e+05, 1.1025502039745696e+05, 1.1026449402474455e+05, 1.1027396772888133e+05, 1.1028344150986137e+05, 1.1029291536767880e+05, 1.1030238930232771e+05, 1.1031186331380217e+05, 1.1032133740209631e+05, 1.1033081156720422e+05, 1.1034028580911997e+05, 1.1034976012783771e+05, 1.1035923452335151e+05, 1.1036870899565547e+05, 1.1037818354474372e+05, 1.1038765817061033e+05, 1.1039713287324943e+05, 1.1040660765265512e+05, 1.1041608250882152e+05, 1.1042555744174273e+05, 1.1043503245141284e+05, 1.1044450753782599e+05, 1.1045398270097628e+05, 1.1046345794085781e+05, 1.1047293325746471e+05, 1.1048240865079107e+05, 1.1049188412083103e+05, 1.1050135966757870e+05, 1.1051083529102818e+05, 1.1052031099117362e+05, 1.1052978676800911e+05, 1.1053926262152876e+05, 1.1054873855172671e+05, 1.1055821455859707e+05, 1.1056769064213398e+05, 1.1057716680233156e+05, 1.1058664303918391e+05, 1.1059611935268519e+05, 1.1060559574282948e+05, 1.1061507220961095e+05, 1.1062454875302372e+05, 1.1063402537306187e+05, 1.1064350206971960e+05, 1.1065297884299098e+05, 1.1066245569287018e+05, 1.1067193261935131e+05, 1.1068140962242852e+05, 1.1069088670209593e+05, 1.1070036385834769e+05, 1.1070984109117791e+05, 1.1071931840058074e+05, 1.1072879578655031e+05, 1.1073827324908078e+05, 1.1074775078816626e+05, 1.1075722840380089e+05, 1.1076670609597885e+05, 1.1077618386469423e+05, 1.1078566170994123e+05, 1.1079513963171396e+05, 1.1080461763000654e+05, 1.1081409570481317e+05, 1.1082357385612796e+05, 1.1083305208394506e+05, 1.1084253038825863e+05, 1.1085200876906281e+05, 1.1086148722635175e+05, 1.1087096576011958e+05, 1.1088044437036049e+05, 1.1088992305706862e+05, 1.1089940182023810e+05, 1.1090888065986311e+05, 1.1091835957593779e+05, 1.1092783856845633e+05, 1.1093731763741282e+05, 1.1094679678280148e+05, 1.1095627600461645e+05, 1.1096575530285187e+05, 1.1097523467750193e+05, 1.1098471412856076e+05, 1.1099419365602257e+05, 1.1100367325988147e+05, 1.1101315294013165e+05, 1.1102263269676728e+05, 1.1103211252978251e+05, 1.1104159243917151e+05, 1.1105107242492847e+05, 1.1106055248704751e+05, 1.1107003262552284e+05, 1.1107951284034862e+05, 1.1108899313151902e+05, 1.1109847349902822e+05, 1.1110795394287037e+05, 1.1111743446303967e+05, 1.1112691505953028e+05, 1.1113639573233636e+05, 1.1114587648145213e+05, 1.1115535730687172e+05, 1.1116483820858932e+05, 1.1117431918659914e+05, 1.1118380024089533e+05, 1.1119328137147210e+05, 1.1120276257832359e+05, 1.1121224386144402e+05, 1.1122172522082755e+05, 1.1123120665646836e+05, 1.1124068816836066e+05, 1.1125016975649861e+05, 1.1125965142087643e+05, 1.1126913316148828e+05, 1.1127861497832836e+05, 1.1128809687139087e+05, 1.1129757884066997e+05, 1.1130706088615989e+05, 1.1131654300785478e+05, 1.1132602520574885e+05, 1.1133550747983632e+05, 1.1134498983011134e+05, 1.1135447225656816e+05, 1.1136395475920092e+05, 1.1137343733800385e+05, 1.1138291999297116e+05, 1.1139240272409702e+05, 1.1140188553137564e+05, 1.1141136841480123e+05, 1.1142085137436797e+05, 1.1143033441007008e+05, 1.1143981752190177e+05, 1.1144930070985726e+05, 1.1145878397393071e+05, 1.1146826731411637e+05, 1.1147775073040841e+05, 1.1148723422280104e+05, 1.1149671779128851e+05, 1.1150620143586499e+05, 1.1151568515652471e+05, 1.1152516895326188e+05, 1.1153465282607070e+05, 1.1154413677494541e+05, 1.1155362079988020e+05, 1.1156310490086928e+05, 1.1157258907790687e+05, 1.1158207333098722e+05, 1.1159155766010450e+05, 1.1160104206525297e+05, 1.1161052654642683e+05, 1.1162001110362030e+05, 1.1162949573682761e+05, 1.1163898044604297e+05, 1.1164846523126062e+05, 1.1165795009247477e+05, 1.1166743502967965e+05, 1.1167692004286949e+05, 1.1168640513203850e+05, 1.1169589029718093e+05, 1.1170537553829099e+05, 1.1171486085536292e+05, 1.1172434624839094e+05, 1.1173383171736929e+05, 1.1174331726229221e+05, 1.1175280288315393e+05, 1.1176228857994865e+05, 1.1177177435267065e+05, 1.1178126020131414e+05, 1.1179074612587337e+05, 1.1180023212634257e+05, 1.1180971820271599e+05, 1.1181920435498783e+05, 1.1182869058315239e+05, 1.1183817688720388e+05, 1.1184766326713652e+05, 1.1185714972294458e+05, 1.1186663625462230e+05, 1.1187612286216392e+05, 1.1188560954556367e+05, 1.1189509630481582e+05, 1.1190458313991460e+05, 1.1191407005085426e+05, 1.1192355703762907e+05, 1.1193304410023325e+05, 1.1194253123866106e+05, 1.1195201845290678e+05, 1.1196150574296464e+05, 1.1197099310882889e+05, 1.1198048055049377e+05, 1.1198996806795355e+05, 1.1199945566120249e+05, 1.1200894333023485e+05, 1.1201843107504486e+05, 1.1202791889562679e+05, 1.1203740679197492e+05, 1.1204689476408350e+05, 1.1205638281194679e+05, 1.1206587093555904e+05, 1.1207535913491451e+05, 1.1208484741000748e+05, 1.1209433576083223e+05, 1.1210382418738298e+05, 1.1211331268965403e+05, 1.1212280126763962e+05, 1.1213228992133406e+05, 1.1214177865073156e+05, 1.1215126745582644e+05, 1.1216075633661296e+05, 1.1217024529308539e+05, 1.1217973432523798e+05, 1.1218922343306501e+05, 1.1219871261656078e+05, 1.1220820187571953e+05, 1.1221769121053556e+05, 1.1222718062100312e+05, 1.1223667010711653e+05, 1.1224615966887004e+05, 1.1225564930625792e+05, 1.1226513901927449e+05, 1.1227462880791399e+05, 1.1228411867217071e+05, 1.1229360861203894e+05, 1.1230309862751298e+05, 1.1231258871858708e+05, 1.1232207888525556e+05, 1.1233156912751266e+05, 1.1234105944535270e+05, 1.1235054983876995e+05, 1.1236004030775873e+05, 1.1236953085231330e+05, 1.1237902147242796e+05, 1.1238851216809700e+05, 1.1239800293931470e+05, 1.1240749378607537e+05, 1.1241698470837332e+05, 1.1242647570620282e+05, 1.1243596677955816e+05, 1.1244545792843365e+05, 1.1245494915282358e+05, 1.1246444045272225e+05, 1.1247393182812394e+05, 1.1248342327902299e+05, 1.1249291480541367e+05, 1.1250240640729028e+05, 1.1251189808464715e+05, 1.1252138983747858e+05, 1.1253088166577884e+05, 1.1254037356954227e+05, 1.1254986554876315e+05, 1.1255935760343580e+05, 1.1256884973355451e+05, 1.1257834193911363e+05, 1.1258783422010743e+05, 1.1259732657653023e+05, 1.1260681900837636e+05, 1.1261631151564009e+05, 1.1262580409831578e+05, 1.1263529675639770e+05, 1.1264478948988020e+05, 1.1265428229875758e+05, 1.1266377518302415e+05, 1.1267326814267423e+05, 1.1268276117770215e+05, 1.1269225428810222e+05, 1.1270174747386877e+05, 1.1271124073499610e+05, 1.1272073407147854e+05, 1.1273022748331040e+05, 1.1273972097048604e+05, 1.1274921453299974e+05, 1.1275870817084586e+05, 1.1276820188401868e+05, 1.1277769567251256e+05, 1.1278718953632185e+05, 1.1279668347544083e+05, 1.1280617748986388e+05, 1.1281567157958528e+05, 1.1282516574459938e+05, 1.1283465998490050e+05, 1.1284415430048300e+05, 1.1285364869134118e+05, 1.1286314315746940e+05, 1.1287263769886199e+05, 1.1288213231551328e+05, 1.1289162700741763e+05, 1.1290112177456931e+05, 1.1291061661696274e+05, 1.1292011153459222e+05, 1.1292960652745208e+05, 1.1293910159553669e+05, 1.1294859673884037e+05, 1.1295809195735746e+05, 1.1296758725108234e+05, 1.1297708262000930e+05, 1.1298657806413270e+05, 1.1299607358344691e+05, 1.1300556917794626e+05, 1.1301506484762509e+05, 1.1302456059247779e+05, 1.1303405641249867e+05, 1.1304355230768207e+05, 1.1305304827802238e+05, 1.1306254432351392e+05, 1.1307204044415106e+05, 1.1308153663992813e+05, 1.1309103291083954e+05, 1.1310052925687958e+05, 1.1311002567804266e+05, 1.1311952217432308e+05, 1.1312901874571525e+05, 1.1313851539221351e+05, 1.1314801211381219e+05, 1.1315750891050571e+05, 1.1316700578228838e+05, 1.1317650272915458e+05, 1.1318599975109867e+05, 1.1319549684811503e+05, 1.1320499402019801e+05, 1.1321449126734196e+05, 1.1322398858954127e+05, 1.1323348598679031e+05, 1.1324298345908342e+05, 1.1325248100641501e+05, 1.1326197862877941e+05, 1.1327147632617102e+05, 1.1328097409858418e+05, 1.1329047194601329e+05, 1.1329996986845272e+05, 1.1330946786589683e+05, 1.1331896593834001e+05, 1.1332846408577664e+05, 1.1333796230820107e+05, 1.1334746060560769e+05, 1.1335695897799087e+05, 1.1336645742534501e+05, 1.1337595594766449e+05, 1.1338545454494368e+05, 1.1339495321717694e+05, 1.1340445196435868e+05, 1.1341395078648326e+05, 1.1342344968354510e+05, 1.1343294865553857e+05, 1.1344244770245806e+05, 1.1345194682429792e+05, 1.1346144602105257e+05, 1.1347094529271639e+05, 1.1348044463928379e+05, 1.1348994406074913e+05, 1.1349944355710682e+05, 1.1350894312835124e+05, 1.1351844277447677e+05, 1.1352794249547784e+05, 1.1353744229134881e+05, 1.1354694216208410e+05, 1.1355644210767809e+05, 1.1356594212812519e+05, 1.1357544222341977e+05, 1.1358494239355628e+05, 1.1359444263852906e+05, 1.1360394295833255e+05, 1.1361344335296113e+05, 1.1362294382240919e+05, 1.1363244436667116e+05, 1.1364194498574144e+05, 1.1365144567961442e+05, 1.1366094644828454e+05, 1.1367044729174614e+05, 1.1367994820999367e+05, 1.1368944920302156e+05, 1.1369895027082418e+05, 1.1370845141339593e+05, 1.1371795263073125e+05, 1.1372745392282454e+05, 1.1373695528967022e+05, 1.1374645673126269e+05, 1.1375595824759638e+05, 1.1376545983866569e+05, 1.1377496150446503e+05, 1.1378446324498882e+05, 1.1379396506023149e+05, 1.1380346695018745e+05, 1.1381296891485110e+05, 1.1382247095421689e+05, 1.1383197306827921e+05, 1.1384147525703251e+05, 1.1385097752047121e+05, 1.1386047985858971e+05, 1.1386998227138244e+05, 1.1387948475884383e+05, 1.1388898732096830e+05, 1.1389848995775028e+05, 1.1390799266918418e+05, 1.1391749545526446e+05, 1.1392699831598555e+05, 1.1393650125134183e+05, 1.1394600426132778e+05, 1.1395550734593780e+05, 1.1396501050516634e+05, 1.1397451373900782e+05, 1.1398401704745668e+05, 1.1399352043050734e+05, 1.1400302388815426e+05, 1.1401252742039185e+05, 1.1402203102721459e+05, 1.1403153470861685e+05, 1.1404103846459312e+05, 1.1405054229513783e+05, 1.1406004620024543e+05, 1.1406955017991034e+05, 1.1407905423412699e+05, 1.1408855836288986e+05, 1.1409806256619336e+05, 1.1410756684403194e+05, 1.1411707119640007e+05, 1.1412657562329216e+05, 1.1413608012470270e+05, 1.1414558470062610e+05, 1.1415508935105681e+05, 1.1416459407598931e+05, 1.1417409887541801e+05, 1.1418360374933740e+05, 1.1419310869774190e+05, 1.1420261372062597e+05, 1.1421211881798408e+05, 1.1422162398981067e+05, 1.1423112923610018e+05, 1.1424063455684710e+05, 1.1425013995204587e+05, 1.1425964542169093e+05, 1.1426915096577675e+05, 1.1427865658429782e+05, 1.1428816227724854e+05, 1.1429766804462342e+05, 1.1430717388641690e+05, 1.1431667980262345e+05, 1.1432618579323753e+05, 1.1433569185825360e+05, 1.1434519799766612e+05, 1.1435470421146958e+05, 1.1436421049965844e+05, 1.1437371686222714e+05, 1.1438322329917017e+05, 1.1439272981048199e+05, 1.1440223639615708e+05, 1.1441174305618991e+05, 1.1442124979057493e+05, 1.1443075659930665e+05, 1.1444026348237951e+05, 1.1444977043978799e+05, 1.1445927747152656e+05, 1.1446878457758971e+05, 1.1447829175797194e+05, 1.1448779901266769e+05, 1.1449730634167143e+05, 1.1450681374497766e+05, 1.1451632122258085e+05, 1.1452582877447549e+05, 1.1453533640065607e+05, 1.1454484410111707e+05, 1.1455435187585292e+05, 1.1456385972485818e+05, 1.1457336764812727e+05, 1.1458287564565473e+05, 1.1459238371743502e+05, 1.1460189186346262e+05, 1.1461140008373202e+05, 1.1462090837823771e+05, 1.1463041674697418e+05, 1.1463992518993595e+05, 1.1464943370711747e+05, 1.1465894229851326e+05, 1.1466845096411777e+05, 1.1467795970392555e+05, 1.1468746851793103e+05, 1.1469697740612878e+05, 1.1470648636851324e+05, 1.1471599540507894e+05, 1.1472550451582034e+05, 1.1473501370073197e+05, 1.1474452295980831e+05, 1.1475403229304389e+05, 1.1476354170043318e+05, 1.1477305118197070e+05, 1.1478256073765094e+05, 1.1479207036746840e+05, 1.1480158007141763e+05, 1.1481108984949306e+05, 1.1482059970168925e+05, 1.1483010962800069e+05, 1.1483961962842186e+05, 1.1484912970294732e+05, 1.1485863985157153e+05, 1.1486815007428905e+05, 1.1487766037109436e+05, 1.1488717074198197e+05, 1.1489668118694640e+05, 1.1490619170598216e+05, 1.1491570229908376e+05, 1.1492521296624573e+05, 1.1493472370746257e+05, 1.1494423452272879e+05, 1.1495374541203893e+05, 1.1496325637538749e+05, 1.1497276741276900e+05, 1.1498227852417796e+05, 1.1499178970960893e+05, 1.1500130096905639e+05, 1.1501081230251487e+05, 1.1502032370997893e+05, 1.1502983519144302e+05, 1.1503934674690173e+05, 1.1504885837634956e+05, 1.1505837007978105e+05, 1.1506788185719070e+05, 1.1507739370857306e+05, 1.1508690563392267e+05, 1.1509641763323401e+05, 1.1510592970650166e+05, 1.1511544185372014e+05, 1.1512495407488396e+05, 1.1513446636998767e+05, 1.1514397873902581e+05, 1.1515349118199291e+05, 1.1516300369888349e+05, 1.1517251628969210e+05, 1.1518202895441325e+05, 1.1519154169304152e+05, 1.1520105450557143e+05, 1.1521056739199749e+05, 1.1522008035231428e+05, 1.1522959338651634e+05, 1.1523910649459816e+05, 1.1524861967655434e+05, 1.1525813293237941e+05, 1.1526764626206791e+05, 1.1527715966561438e+05, 1.1528667314301335e+05, 1.1529618669425941e+05, 1.1530570031934707e+05, 1.1531521401827088e+05, 1.1532472779102541e+05, 1.1533424163760518e+05, 1.1534375555800476e+05, 1.1535326955221870e+05, 1.1536278362024153e+05, 1.1537229776206784e+05, 1.1538181197769217e+05, 1.1539132626710905e+05, 1.1540084063031308e+05, 1.1541035506729878e+05, 1.1541986957806071e+05, 1.1542938416259344e+05, 1.1543889882089151e+05, 1.1544841355294951e+05, 1.1545792835876197e+05, 1.1546744323832345e+05, 1.1547695819162854e+05, 1.1548647321867177e+05, 1.1549598831944774e+05, 1.1550550349395099e+05, 1.1551501874217609e+05, 1.1552453406411759e+05, 1.1553404945977007e+05, 1.1554356492912810e+05, 1.1555308047218625e+05, 1.1556259608893908e+05, 1.1557211177938117e+05, 1.1558162754350707e+05, 1.1559114338131137e+05, 1.1560065929278864e+05, 1.1561017527793342e+05, 1.1561969133674033e+05, 1.1562920746920395e+05, 1.1563872367531880e+05, 1.1564823995507949e+05, 1.1565775630848062e+05, 1.1566727273551672e+05, 1.1567678923618239e+05, 1.1568630581047223e+05, 1.1569582245838078e+05, 1.1570533917990263e+05, 1.1571485597503238e+05, 1.1572437284376462e+05, 1.1573388978609387e+05, 1.1574340680201480e+05, 1.1575292389152195e+05, 1.1576244105460990e+05, 1.1577195829127326e+05, 1.1578147560150659e+05, 1.1579099298530450e+05, 1.1580051044266157e+05, 1.1581002797357239e+05, 1.1581954557803154e+05, 1.1582906325603362e+05, 1.1583858100757323e+05, 1.1584809883264496e+05, 1.1585761673124341e+05, 1.1586713470336313e+05, 1.1587665274899876e+05, 1.1588617086814489e+05, 1.1589568906079612e+05, 1.1590520732694703e+05, 1.1591472566659222e+05, 1.1592424407972630e+05, 1.1593376256634385e+05, 1.1594328112643950e+05, 1.1595279976000784e+05, 1.1596231846704346e+05, 1.1597183724754097e+05, 1.1598135610149497e+05, 1.1599087502890007e+05, 1.1600039402975090e+05, 1.1600991310404202e+05, 1.1601943225176807e+05, 1.1602895147292363e+05, 1.1603847076750334e+05, 1.1604799013550178e+05, 1.1605750957691358e+05, 1.1606702909173333e+05, 1.1607654867995568e+05, 1.1608606834157521e+05, 1.1609558807658654e+05, 1.1610510788498429e+05, 1.1611462776676308e+05, 1.1612414772191750e+05, 1.1613366775044220e+05, 1.1614318785233177e+05, 1.1615270802758084e+05, 1.1616222827618402e+05, 1.1617174859813595e+05, 1.1618126899343122e+05, 1.1619078946206447e+05, 1.1620031000403033e+05, 1.1620983061932340e+05, 1.1621935130793833e+05, 1.1622887206986970e+05, 1.1623839290511220e+05, 1.1624791381366040e+05, 1.1625743479550896e+05, 1.1626695585065248e+05, 1.1627647697908564e+05, 1.1628599818080301e+05, 1.1629551945579925e+05, 1.1630504080406899e+05, 1.1631456222560685e+05, 1.1632408372040746e+05, 1.1633360528846546e+05, 1.1634312692977549e+05, 1.1635264864433216e+05, 1.1636217043213015e+05, 1.1637169229316406e+05, 1.1638121422742853e+05, 1.1639073623491822e+05, 1.1640025831562774e+05, 1.1640978046955173e+05, 1.1641930269668486e+05, 1.1642882499702174e+05, 1.1643834737055702e+05, 1.1644786981728536e+05, 1.1645739233720140e+05, 1.1646691493029975e+05, 1.1647643759657510e+05, 1.1648596033602206e+05, 1.1649548314863529e+05, 1.1650500603440945e+05, 1.1651452899333916e+05, 1.1652405202541909e+05, 1.1653357513064389e+05, 1.1654309830900820e+05, 1.1655262156050667e+05, 1.1656214488513395e+05, 1.1657166828288474e+05, 1.1658119175375362e+05, 1.1659071529773527e+05, 1.1660023891482435e+05, 1.1660976260501552e+05, 1.1661928636830342e+05, 1.1662881020468274e+05, 1.1663833411414812e+05, 1.1664785809669422e+05, 1.1665738215231569e+05, 1.1666690628100722e+05, 1.1667643048276345e+05, 1.1668595475757902e+05, 1.1669547910544861e+05, 1.1670500352636691e+05, 1.1671452802032854e+05, 1.1672405258732819e+05, 1.1673357722736053e+05, 1.1674310194042021e+05, 1.1675262672650190e+05, 1.1676215158560028e+05, 1.1677167651771000e+05, 1.1678120152282577e+05, 1.1679072660094222e+05, 1.1680025175205403e+05, 1.1680977697615589e+05, 1.1681930227324244e+05, 1.1682882764330840e+05, 1.1683835308634839e+05, 1.1684787860235713e+05, 1.1685740419132927e+05, 1.1686692985325951e+05, 1.1687645558814250e+05, 1.1688598139597294e+05, 1.1689550727674550e+05, 1.1690503323045486e+05, 1.1691455925709568e+05, 1.1692408535666268e+05, 1.1693361152915053e+05, 1.1694313777455389e+05, 1.1695266409286746e+05, 1.1696219048408596e+05, 1.1697171694820402e+05, 1.1698124348521633e+05, 1.1699077009511762e+05, 1.1700029677790253e+05, 1.1700982353356580e+05, 1.1701935036210206e+05, 1.1702887726350604e+05, 1.1703840423777244e+05, 1.1704793128489589e+05, 1.1705745840487115e+05, 1.1706698559769288e+05, 1.1707651286335578e+05, 1.1708604020185454e+05, 1.1709556761318387e+05, 1.1710509509733845e+05, 1.1711462265431297e+05, 1.1712415028410214e+05, 1.1713367798670067e+05, 1.1714320576210326e+05, 1.1715273361030458e+05, 1.1716226153129934e+05, 1.1717178952508228e+05, 1.1718131759164807e+05, 1.1719084573099139e+05, 1.1720037394310700e+05, 1.1720990222798954e+05, 1.1721943058563377e+05, 1.1722895901603438e+05, 1.1723848751918605e+05, 1.1724801609508351e+05, 1.1725754474372149e+05, 1.1726707346509467e+05, 1.1727660225919777e+05, 1.1728613112602550e+05, 1.1729566006557256e+05, 1.1730518907783368e+05, 1.1731471816280356e+05, 1.1732424732047692e+05, 1.1733377655084847e+05, 1.1734330585391294e+05, 1.1735283522966501e+05, 1.1736236467809945e+05, 1.1737189419921093e+05, 1.1738142379299420e+05, 1.1739095345944396e+05, 1.1740048319855494e+05, 1.1741001301032187e+05, 1.1741954289473944e+05, 1.1742907285180240e+05, 1.1743860288150546e+05, 1.1744813298384334e+05, 1.1745766315881078e+05, 1.1746719340640251e+05, 1.1747672372661324e+05, 1.1748625411943770e+05, 1.1749578458487061e+05, 1.1750531512290672e+05, 1.1751484573354074e+05, 1.1752437641676741e+05, 1.1753390717258147e+05, 1.1754343800097762e+05, 1.1755296890195062e+05, 1.1756249987549518e+05, 1.1757203092160606e+05, 1.1758156204027798e+05, 1.1759109323150567e+05, 1.1760062449528389e+05, 1.1761015583160734e+05, 1.1761968724047080e+05, 1.1762921872186897e+05, 1.1763875027579660e+05, 1.1764828190224843e+05, 1.1765781360121923e+05, 1.1766734537270370e+05, 1.1767687721669659e+05, 1.1768640913319266e+05, 1.1769594112218666e+05, 1.1770547318367331e+05, 1.1771500531764737e+05, 1.1772453752410358e+05, 1.1773406980303669e+05, 1.1774360215444144e+05, 1.1775313457831259e+05, 1.1776266707464487e+05, 1.1777219964343306e+05, 1.1778173228467189e+05, 1.1779126499835610e+05, 1.1780079778448046e+05, 1.1781033064303972e+05, 1.1781986357402863e+05, 1.1782939657744196e+05, 1.1783892965327445e+05, 1.1784846280152087e+05, 1.1785799602217594e+05, 1.1786752931523447e+05, 1.1787706268069118e+05, 1.1788659611854081e+05, 1.1789612962877817e+05, 1.1790566321139800e+05, 1.1791519686639505e+05, 1.1792473059376411e+05, 1.1793426439349991e+05, 1.1794379826559722e+05, 1.1795333221005082e+05, 1.1796286622685548e+05, 1.1797240031600594e+05, 1.1798193447749698e+05, 1.1799146871132337e+05, 1.1800100301747987e+05, 1.1801053739596125e+05, 1.1802007184676229e+05, 1.1802960636987776e+05, 1.1803914096530242e+05, 1.1804867563303105e+05, 1.1805821037305841e+05, 1.1806774518537927e+05, 1.1807728006998844e+05, 1.1808681502688066e+05, 1.1809635005605072e+05, 1.1810588515749338e+05, 1.1811542033120342e+05, 1.1812495557717564e+05, 1.1813449089540483e+05, 1.1814402628588573e+05, 1.1815356174861314e+05, 1.1816309728358182e+05, 1.1817263289078658e+05, 1.1818216857022219e+05, 1.1819170432188343e+05, 1.1820124014576510e+05, 1.1821077604186196e+05, 1.1822031201016880e+05, 1.1822984805068040e+05, 1.1823938416339159e+05, 1.1824892034829711e+05, 1.1825845660539177e+05, 1.1826799293467034e+05, 1.1827752933612764e+05, 1.1828706580975844e+05, 1.1829660235555752e+05, 1.1830613897351970e+05, 1.1831567566363978e+05, 1.1832521242591251e+05, 1.1833474926033271e+05, 1.1834428616689518e+05, 1.1835382314559471e+05, 1.1836336019642612e+05, 1.1837289731938417e+05, 1.1838243451446368e+05, 1.1839197178165942e+05, 1.1840150912096622e+05, 1.1841104653237888e+05, 1.1842058401589219e+05, 1.1843012157150095e+05, 1.1843965919919996e+05, 1.1844919689898402e+05, 1.1845873467084797e+05, 1.1846827251478657e+05, 1.1847781043079466e+05, 1.1848734841886703e+05, 1.1849688647899848e+05, 1.1850642461118383e+05, 1.1851596281541791e+05, 1.1852550109169546e+05, 1.1853503944001137e+05, 1.1854457786036041e+05, 1.1855411635273740e+05, 1.1856365491713713e+05, 1.1857319355355443e+05, 1.1858273226198413e+05, 1.1859227104242102e+05, 1.1860180989485994e+05, 1.1861134881929569e+05, 1.1862088781572309e+05, 1.1863042688413696e+05, 1.1863996602453213e+05, 1.1864950523690338e+05, 1.1865904452124557e+05, 1.1866858387755348e+05, 1.1867812330582195e+05, 1.1868766280604583e+05, 1.1869720237821992e+05, 1.1870674202233905e+05, 1.1871628173839803e+05, 1.1872582152639168e+05, 1.1873536138631486e+05, 1.1874490131816237e+05, 1.1875444132192904e+05, 1.1876398139760971e+05, 1.1877352154519920e+05, 1.1878306176469232e+05, 1.1879260205608392e+05, 1.1880214241936886e+05, 1.1881168285454191e+05, 1.1882122336159795e+05, 1.1883076394053178e+05, 1.1884030459133827e+05, 1.1884984531401221e+05, 1.1885938610854850e+05, 1.1886892697494190e+05, 1.1887846791318730e+05, 1.1888800892327951e+05, 1.1889755000521342e+05, 1.1890709115898379e+05, 1.1891663238458554e+05, 1.1892617368201345e+05, 1.1893571505126238e+05, 1.1894525649232717e+05, 1.1895479800520268e+05, 1.1896433958988373e+05, 1.1897388124636518e+05, 1.1898342297464189e+05, 1.1899296477470866e+05, 1.1900250664656036e+05, 1.1901204859019184e+05, 1.1902159060559796e+05, 1.1903113269277356e+05, 1.1904067485171345e+05, 1.1905021708241256e+05, 1.1905975938486568e+05, 1.1906930175906769e+05, 1.1907884420501342e+05, 1.1908838672269773e+05, 1.1909792931211549e+05, 1.1910747197326153e+05, 1.1911701470613074e+05, 1.1912655751071793e+05, 1.1913610038701798e+05, 1.1914564333502577e+05, 1.1915518635473614e+05, 1.1916472944614395e+05, 1.1917427260924404e+05, 1.1918381584403128e+05, 1.1919335915050055e+05, 1.1920290252864671e+05, 1.1921244597846459e+05, 1.1922198949994911e+05, 1.1923153309309506e+05, 1.1924107675789736e+05, 1.1925062049435086e+05, 1.1926016430245044e+05, 1.1926970818219095e+05, 1.1927925213356728e+05, 1.1928879615657426e+05, 1.1929834025120678e+05, 1.1930788441745973e+05, 1.1931742865532797e+05, 1.1932697296480635e+05, 1.1933651734588976e+05, 1.1934606179857308e+05, 1.1935560632285116e+05, 1.1936515091871889e+05, 1.1937469558617113e+05, 1.1938424032520279e+05, 1.1939378513580871e+05, 1.1940333001798380e+05, 1.1941287497172291e+05, 1.1942241999702093e+05, 1.1943196509387276e+05, 1.1944151026227324e+05, 1.1945105550221728e+05, 1.1946060081369976e+05, 1.1947014619671555e+05, 1.1947969165125953e+05, 1.1948923717732659e+05, 1.1949878277491164e+05, 1.1950832844400953e+05, 1.1951787418461518e+05, 1.1952741999672345e+05, 1.1953696588032924e+05, 1.1954651183542740e+05, 1.1955605786201289e+05, 1.1956560396008054e+05, 1.1957515012962528e+05, 1.1958469637064198e+05, 1.1959424268312553e+05, 1.1960378906707083e+05, 1.1961333552247277e+05, 1.1962288204932626e+05, 1.1963242864762621e+05, 1.1964197531736745e+05, 1.1965152205854491e+05, 1.1966106887115352e+05, 1.1967061575518813e+05, 1.1968016271064366e+05, 1.1968970973751502e+05, 1.1969925683579710e+05, 1.1970880400548480e+05, 1.1971835124657302e+05, 1.1972789855905666e+05, 1.1973744594293063e+05, 1.1974699339818982e+05, 1.1975654092482917e+05, 1.1976608852284354e+05, 1.1977563619222786e+05, 1.1978518393297704e+05, 1.1979473174508597e+05, 1.1980427962854957e+05, 1.1981382758336276e+05, 1.1982337560952042e+05, 1.1983292370701750e+05, 1.1984247187584886e+05, 1.1985202011600946e+05, 1.1986156842749417e+05, 1.1987111681029793e+05, 1.1988066526441567e+05, 1.1989021378984227e+05, 1.1989976238657265e+05, 1.1990931105460174e+05, 1.1991885979392446e+05, 1.1992840860453573e+05, 1.1993795748643044e+05, 1.1994750643960355e+05, 1.1995705546404993e+05, 1.1996660455976456e+05, 1.1997615372674231e+05, 1.1998570296497812e+05, 1.1999525227446693e+05, 1.2000480165520363e+05, 1.2001435110718315e+05, 1.2002390063040044e+05, 1.2003345022485043e+05, 1.2004299989052799e+05, 1.2005254962742812e+05, 1.2006209943554568e+05, 1.2007164931487564e+05, 1.2008119926541293e+05, 1.2009074928715249e+05, 1.2010029938008920e+05, 1.2010984954421804e+05, 1.2011939977953392e+05, 1.2012895008603178e+05, 1.2013850046370654e+05, 1.2014805091255315e+05, 1.2015760143256655e+05, 1.2016715202374164e+05, 1.2017670268607339e+05, 1.2018625341955673e+05, 1.2019580422418659e+05, 1.2020535509995793e+05, 1.2021490604686567e+05, 1.2022445706490475e+05, 1.2023400815407012e+05, 1.2024355931435672e+05, 1.2025311054575948e+05, 1.2026266184827335e+05, 1.2027221322189327e+05, 1.2028176466661418e+05, 1.2029131618243105e+05, 1.2030086776933882e+05, 1.2031041942733240e+05, 1.2031997115640678e+05, 1.2032952295655689e+05, 1.2033907482777769e+05, 1.2034862677006410e+05, 1.2035817878341109e+05, 1.2036773086781362e+05, 1.2037728302326665e+05, 1.2038683524976508e+05, 1.2039638754730391e+05, 1.2040593991587807e+05, 1.2041549235548254e+05, 1.2042504486611225e+05, 1.2043459744776218e+05, 1.2044415010042724e+05, 1.2045370282410244e+05, 1.2046325561878271e+05, 1.2047280848446304e+05, 1.2048236142113832e+05, 1.2049191442880360e+05, 1.2050146750745377e+05, 1.2051102065708382e+05, 1.2052057387768870e+05, 1.2053012716926339e+05, 1.2053968053180285e+05, 1.2054923396530205e+05, 1.2055878746975595e+05, 1.2056834104515950e+05, 1.2057789469150768e+05, 1.2058744840879543e+05, 1.2059700219701776e+05, 1.2060655605616963e+05, 1.2061610998624600e+05, 1.2062566398724183e+05, 1.2063521805915210e+05, 1.2064477220197181e+05, 1.2065432641569589e+05, 1.2066388070031932e+05, 1.2067343505583708e+05, 1.2068298948224416e+05, 1.2069254397953551e+05, 1.2070209854770612e+05, 1.2071165318675098e+05, 1.2072120789666502e+05, 1.2073076267744326e+05, 1.2074031752908068e+05, 1.2074987245157223e+05, 1.2075942744491292e+05, 1.2076898250909771e+05, 1.2077853764412161e+05, 1.2078809284997956e+05, 1.2079764812666658e+05, 1.2080720347417762e+05, 1.2081675889250770e+05, 1.2082631438165177e+05, 1.2083586994160485e+05, 1.2084542557236190e+05, 1.2085498127391792e+05, 1.2086453704626788e+05, 1.2087409288940679e+05, 1.2088364880332963e+05, 1.2089320478803141e+05, 1.2090276084350709e+05, 1.2091231696975167e+05, 1.2092187316676015e+05, 1.2093142943452753e+05, 1.2094098577304877e+05, 1.2095054218231891e+05, 1.2096009866233289e+05, 1.2096965521308576e+05, 1.2097921183457249e+05, 1.2098876852678809e+05, 1.2099832528972754e+05, 1.2100788212338585e+05, 1.2101743902775800e+05, 1.2102699600283902e+05, 1.2103655304862390e+05, 1.2104611016510763e+05, 1.2105566735228522e+05, 1.2106522461015168e+05, 1.2107478193870203e+05, 1.2108433933793122e+05, 1.2109389680783430e+05, 1.2110345434840625e+05, 1.2111301195964208e+05, 1.2112256964153681e+05, 1.2113212739408544e+05, 1.2114168521728300e+05, 1.2115124311112447e+05, 1.2116080107560486e+05, 1.2117035911071919e+05, 1.2117991721646248e+05, 1.2118947539282971e+05, 1.2119903363981594e+05, 1.2120859195741615e+05, 1.2121815034562537e+05, 1.2122770880443860e+05, 1.2123726733385086e+05, 1.2124682593385718e+05, 1.2125638460445254e+05, 1.2126594334563200e+05, 1.2127550215739057e+05, 1.2128506103972327e+05, 1.2129461999262508e+05, 1.2130417901609108e+05, 1.2131373811011625e+05, 1.2132329727469562e+05, 1.2133285650982421e+05, 1.2134241581549705e+05, 1.2135197519170918e+05, 1.2136153463845560e+05, 1.2137109415573133e+05, 1.2138065374353142e+05, 1.2139021340185090e+05, 1.2139977313068476e+05, 1.2140933293002806e+05, 1.2141889279987580e+05, 1.2142845274022306e+05, 1.2143801275106482e+05, 1.2144757283239612e+05, 1.2145713298421203e+05, 1.2146669320650752e+05, 1.2147625349927769e+05, 1.2148581386251752e+05, 1.2149537429622206e+05, 1.2150493480038634e+05, 1.2151449537500543e+05, 1.2152405602007433e+05, 1.2153361673558810e+05, 1.2154317752154177e+05, 1.2155273837793038e+05, 1.2156229930474896e+05, 1.2157186030199255e+05, 1.2158142136965618e+05, 1.2159098250773492e+05, 1.2160054371622381e+05, 1.2161010499511787e+05, 1.2161966634441215e+05, 1.2162922776410170e+05, 1.2163878925418157e+05, 1.2164835081464681e+05, 1.2165791244549246e+05, 1.2166747414671355e+05, 1.2167703591830513e+05, 1.2168659776026227e+05, 1.2169615967258000e+05, 1.2170572165525339e+05, 1.2171528370827748e+05, 1.2172484583164731e+05, 1.2173440802535796e+05, 1.2174397028940446e+05, 1.2175353262378187e+05, 1.2176309502848522e+05, 1.2177265750350959e+05, 1.2178222004885005e+05, 1.2179178266450162e+05, 1.2180134535045938e+05, 1.2181090810671837e+05, 1.2182047093327368e+05, 1.2183003383012033e+05, 1.2183959679725341e+05, 1.2184915983466795e+05, 1.2185872294235903e+05, 1.2186828612032170e+05, 1.2187784936855105e+05, 1.2188741268704213e+05, 1.2189697607578999e+05, 1.2190653953478968e+05, 1.2191610306403633e+05, 1.2192566666352493e+05, 1.2193523033325060e+05, 1.2194479407320838e+05, 1.2195435788339334e+05, 1.2196392176380055e+05, 1.2197348571442509e+05, 1.2198304973526202e+05, 1.2199261382630642e+05, 1.2200217798755335e+05, 1.2201174221899788e+05, 1.2202130652063510e+05, 1.2203087089246005e+05, 1.2204043533446784e+05, 1.2204999984665355e+05, 1.2205956442901223e+05, 1.2206912908153895e+05, 1.2207869380422881e+05, 1.2208825859707686e+05, 1.2209782346007822e+05, 1.2210738839322793e+05, 1.2211695339652109e+05, 1.2212651846995276e+05, 1.2213608361351804e+05, 1.2214564882721200e+05, 1.2215521411102974e+05, 1.2216477946496633e+05, 1.2217434488901685e+05, 1.2218391038317639e+05, 1.2219347594744006e+05, 1.2220304158180289e+05, 1.2221260728626001e+05, 1.2222217306080648e+05, 1.2223173890543742e+05, 1.2224130482014790e+05, 1.2225087080493299e+05, 1.2226043685978781e+05, 1.2227000298470743e+05, 1.2227956917968697e+05, 1.2228913544472148e+05, 1.2229870177980610e+05, 1.2230826818493589e+05, 1.2231783466010595e+05, 1.2232740120531138e+05, 1.2233696782054727e+05, 1.2234653450580873e+05, 1.2235610126109084e+05, 1.2236566808638872e+05, 1.2237523498169742e+05, 1.2238480194701209e+05, 1.2239436898232781e+05, 1.2240393608763968e+05, 1.2241350326294280e+05, 1.2242307050823227e+05, 1.2243263782350319e+05, 1.2244220520875067e+05, 1.2245177266396981e+05, 1.2246134018915572e+05, 1.2247090778430350e+05, 1.2248047544940827e+05, 1.2249004318446512e+05, 1.2249961098946915e+05, 1.2250917886441549e+05, 1.2251874680929922e+05, 1.2252831482411547e+05, 1.2253788290885933e+05, 1.2254745106352594e+05, 1.2255701928811039e+05, 1.2256658758260781e+05, 1.2257615594701329e+05, 1.2258572438132195e+05, 1.2259529288552892e+05, 1.2260486145962930e+05, 1.2261443010361822e+05, 1.2262399881749078e+05, 1.2263356760124209e+05, 1.2264313645486727e+05, 1.2265270537836145e+05, 1.2266227437171974e+05, 1.2267184343493728e+05, 1.2268141256800918e+05, 1.2269098177093054e+05, 1.2270055104369650e+05, 1.2271012038630217e+05, 1.2271968979874266e+05, 1.2272925928101313e+05, 1.2273882883310870e+05, 1.2274839845502448e+05, 1.2275796814675559e+05, 1.2276753790829719e+05, 1.2277710773964437e+05, 1.2278667764079226e+05, 1.2279624761173601e+05, 1.2280581765247072e+05, 1.2281538776299154e+05, 1.2282495794329360e+05, 1.2283452819337201e+05, 1.2284409851322194e+05, 1.2285366890283849e+05, 1.2286323936221682e+05, 1.2287280989135202e+05, 1.2288238049023927e+05, 1.2289195115887368e+05, 1.2290152189725039e+05, 1.2291109270536456e+05, 1.2292066358321127e+05, 1.2293023453078570e+05, 1.2293980554808301e+05, 1.2294937663509829e+05, 1.2295894779182669e+05, 1.2296851901826338e+05, 1.2297809031440347e+05, 1.2298766168024213e+05, 1.2299723311577448e+05, 1.2300680462099565e+05, 1.2301637619590083e+05, 1.2302594784048513e+05, 1.2303551955474369e+05, 1.2304509133867169e+05, 1.2305466319226424e+05, 1.2306423511551651e+05, 1.2307380710842364e+05, 1.2308337917098077e+05, 1.2309295130318307e+05, 1.2310252350502567e+05, 1.2311209577650373e+05, 1.2312166811761241e+05, 1.2313124052834685e+05, 1.2314081300870221e+05, 1.2315038555867362e+05, 1.2315995817825626e+05, 1.2316953086744527e+05, 1.2317910362623581e+05, 1.2318867645462303e+05, 1.2319824935260210e+05, 1.2320782232016818e+05, 1.2321739535731640e+05, 1.2322696846404196e+05, 1.2323654164033999e+05, 1.2324611488620566e+05, 1.2325568820163413e+05, 1.2326526158662055e+05, 1.2327483504116008e+05, 1.2328440856524791e+05, 1.2329398215887918e+05, 1.2330355582204908e+05, 1.2331312955475276e+05, 1.2332270335698535e+05, 1.2333227722874205e+05, 1.2334185117001805e+05, 1.2335142518080848e+05, 1.2336099926110853e+05, 1.2337057341091335e+05, 1.2338014763021811e+05, 1.2338972191901799e+05, 1.2339929627730817e+05, 1.2340887070508380e+05, 1.2341844520234005e+05, 1.2342801976907213e+05, 1.2343759440527517e+05, 1.2344716911094436e+05, 1.2345674388607490e+05, 1.2346631873066192e+05, 1.2347589364470064e+05, 1.2348546862818621e+05, 1.2349504368111381e+05, 1.2350461880347862e+05, 1.2351419399527581e+05, 1.2352376925650057e+05, 1.2353334458714808e+05, 1.2354291998721355e+05, 1.2355249545669211e+05, 1.2356207099557896e+05, 1.2357164660386929e+05, 1.2358122228155829e+05, 1.2359079802864113e+05, 1.2360037384511300e+05, 1.2360994973096908e+05, 1.2361952568620455e+05, 1.2362910171081463e+05, 1.2363867780479448e+05, 1.2364825396813930e+05, 1.2365783020084426e+05, 1.2366740650290456e+05, 1.2367698287431539e+05, 1.2368655931507195e+05, 1.2369613582516943e+05, 1.2370571240460302e+05, 1.2371528905336792e+05, 1.2372486577145927e+05, 1.2373444255887234e+05, 1.2374401941560229e+05, 1.2375359634164430e+05, 1.2376317333699360e+05, 1.2377275040164536e+05, 1.2378232753559480e+05, 1.2379190473883711e+05, 1.2380148201136747e+05, 1.2381105935318109e+05, 1.2382063676427318e+05, 1.2383021424463895e+05, 1.2383979179427358e+05, 1.2384936941317227e+05, 1.2385894710133025e+05, 1.2386852485874269e+05, 1.2387810268540483e+05, 1.2388768058131184e+05, 1.2389725854645894e+05, 1.2390683658084135e+05, 1.2391641468445424e+05, 1.2392599285729286e+05, 1.2393557109935240e+05, 1.2394514941062807e+05, 1.2395472779111507e+05, 1.2396430624080861e+05, 1.2397388475970391e+05, 1.2398346334779618e+05, 1.2399304200508063e+05, 1.2400262073155248e+05, 1.2401219952720693e+05, 1.2402177839203921e+05, 1.2403135732604453e+05, 1.2404093632921808e+05, 1.2405051540155512e+05, 1.2406009454305084e+05, 1.2406967375370045e+05, 1.2407925303349918e+05, 1.2408883238244226e+05, 1.2409841180052489e+05, 1.2410799128774230e+05, 1.2411757084408970e+05, 1.2412715046956233e+05, 1.2413673016415541e+05, 1.2414630992786413e+05, 1.2415588976068376e+05, 1.2416546966260950e+05, 1.2417504963363658e+05, 1.2418462967376021e+05, 1.2419420978297562e+05, 1.2420378996127807e+05, 1.2421337020866275e+05, 1.2422295052512489e+05, 1.2423253091065976e+05, 1.2424211136526254e+05, 1.2425169188892849e+05, 1.2426127248165281e+05, 1.2427085314343075e+05, 1.2428043387425755e+05, 1.2429001467412844e+05, 1.2429959554303864e+05, 1.2430917648098340e+05, 1.2431875748795795e+05, 1.2432833856395750e+05, 1.2433791970897731e+05, 1.2434750092301263e+05, 1.2435708220605868e+05, 1.2436666355811071e+05, 1.2437624497916394e+05, 1.2438582646921360e+05, 1.2439540802825497e+05, 1.2440498965628324e+05, 1.2441457135329369e+05, 1.2442415311928155e+05, 1.2443373495424206e+05, 1.2444331685817047e+05, 1.2445289883106202e+05, 1.2446248087291194e+05, 1.2447206298371550e+05, 1.2448164516346794e+05, 1.2449122741216447e+05, 1.2450080972980040e+05, 1.2451039211637093e+05, 1.2451997457187135e+05, 1.2452955709629686e+05, 1.2453913968964275e+05, 1.2454872235190425e+05, 1.2455830508307660e+05, 1.2456788788315505e+05, 1.2457747075213489e+05, 1.2458705369001134e+05, 1.2459663669677966e+05, 1.2460621977243514e+05, 1.2461580291697297e+05, 1.2462538613038848e+05, 1.2463496941267683e+05, 1.2464455276383336e+05, 1.2465413618385330e+05, 1.2466371967273191e+05, 1.2467330323046444e+05, 1.2468288685704616e+05, 1.2469247055247232e+05, 1.2470205431673821e+05, 1.2471163814983905e+05, 1.2472122205177012e+05, 1.2473080602252670e+05, 1.2474039006210402e+05, 1.2474997417049736e+05, 1.2475955834770197e+05, 1.2476914259371314e+05, 1.2477872690852614e+05, 1.2478831129213620e+05, 1.2479789574453863e+05, 1.2480748026572866e+05, 1.2481706485570160e+05, 1.2482664951445267e+05, 1.2483623424197719e+05, 1.2484581903827039e+05, 1.2485540390332758e+05, 1.2486498883714399e+05, 1.2487457383971491e+05, 1.2488415891103564e+05, 1.2489374405110140e+05, 1.2490332925990752e+05, 1.2491291453744924e+05, 1.2492249988372182e+05, 1.2493208529872059e+05, 1.2494167078244079e+05, 1.2495125633487770e+05, 1.2496084195602662e+05, 1.2497042764588280e+05, 1.2498001340444153e+05, 1.2498959923169810e+05, 1.2499918512764778e+05, 1.2500877109228585e+05, 1.2501835712560761e+05, 1.2502794322760832e+05, 1.2503752939828328e+05, 1.2504711563762776e+05, 1.2505670194563705e+05, 1.2506628832230646e+05, 1.2507587476763123e+05, 1.2508546128160668e+05, 1.2509504786422808e+05, 1.2510463451549073e+05, 1.2511422123538992e+05, 1.2512380802392095e+05, 1.2513339488107906e+05, 1.2514298180685961e+05, 1.2515256880125785e+05, 1.2516215586426907e+05, 1.2517174299588856e+05, 1.2518133019611162e+05, 1.2519091746493356e+05, 1.2520050480234968e+05, 1.2521009220835523e+05, 1.2521967968294554e+05, 1.2522926722611592e+05, 1.2523885483786163e+05, 1.2524844251817798e+05, 1.2525803026706028e+05, 1.2526761808450382e+05, 1.2527720597050391e+05, 1.2528679392505584e+05, 1.2529638194815491e+05, 1.2530597003979643e+05, 1.2531555819997570e+05, 1.2532514642868801e+05, 1.2533473472592866e+05, 1.2534432309169299e+05, 1.2535391152597629e+05, 1.2536350002877385e+05, 1.2537308860008097e+05, 1.2538267723989298e+05, 1.2539226594820518e+05, 1.2540185472501289e+05, 1.2541144357031137e+05, 1.2542103248409599e+05, 1.2543062146636203e+05, 1.2544021051710480e+05, 1.2544979963631963e+05, 1.2545938882400181e+05, 1.2546897808014665e+05, 1.2547856740474948e+05, 1.2548815679780561e+05, 1.2549774625931036e+05, 1.2550733578925903e+05, 1.2551692538764696e+05, 1.2552651505446943e+05, 1.2553610478972178e+05, 1.2554569459339933e+05, 1.2555528446549739e+05, 1.2556487440601129e+05, 1.2557446441493632e+05, 1.2558405449226784e+05, 1.2559364463800115e+05, 1.2560323485213156e+05, 1.2561282513465441e+05, 1.2562241548556503e+05, 1.2563200590485871e+05, 1.2564159639253080e+05, 1.2565118694857664e+05, 1.2566077757299154e+05, 1.2567036826577081e+05, 1.2567995902690978e+05, 1.2568954985640380e+05, 1.2569914075424818e+05, 1.2570873172043824e+05, 1.2571832275496933e+05, 1.2572791385783677e+05, 1.2573750502903591e+05, 1.2574709626856205e+05, 1.2575668757641054e+05, 1.2576627895257670e+05, 1.2577587039705589e+05, 1.2578546190984339e+05, 1.2579505349093459e+05, 1.2580464514032479e+05, 1.2581423685800936e+05, 1.2582382864398358e+05, 1.2583342049824285e+05, 1.2584301242078247e+05, 1.2585260441159779e+05, 1.2586219647068415e+05, 1.2587178859803689e+05, 1.2588138079365135e+05, 1.2589097305752283e+05, 1.2590056538964674e+05, 1.2591015779001838e+05, 1.2591975025863310e+05, 1.2592934279548623e+05, 1.2593893540057316e+05, 1.2594852807388919e+05, 1.2595812081542968e+05, 1.2596771362518998e+05, 1.2597730650316543e+05, 1.2598689944935135e+05, 1.2599649246374314e+05, 1.2600608554633611e+05, 1.2601567869712564e+05, 1.2602527191610706e+05, 1.2603486520327571e+05, 1.2604445855862695e+05, 1.2605405198215615e+05, 1.2606364547385863e+05, 1.2607323903372977e+05, 1.2608283266176490e+05, 1.2609242635795940e+05, 1.2610202012230860e+05, 1.2611161395480788e+05, 1.2612120785545258e+05, 1.2613080182423805e+05, 1.2614039586115966e+05, 1.2614998996621276e+05, 1.2615958413939271e+05, 1.2616917838069487e+05, 1.2617877269011461e+05, 1.2618836706764728e+05, 1.2619796151328825e+05, 1.2620755602703287e+05, 1.2621715060887650e+05, 1.2622674525881452e+05, 1.2623633997684227e+05, 1.2624593476295512e+05, 1.2625552961714847e+05, 1.2626512453941762e+05, 1.2627471952975799e+05, 1.2628431458816493e+05, 1.2629390971463382e+05, 1.2630350490915998e+05, 1.2631310017173883e+05, 1.2632269550236572e+05, 1.2633229090103602e+05, 1.2634188636774510e+05, 1.2635148190248833e+05, 1.2636107750526110e+05, 1.2637067317605877e+05, 1.2638026891487671e+05, 1.2638986472171028e+05, 1.2639946059655488e+05, 1.2640905653940586e+05, 1.2641865255025863e+05, 1.2642824862910854e+05, 1.2643784477595096e+05, 1.2644744099078128e+05, 1.2645703727359489e+05, 1.2646663362438715e+05, 1.2647623004315344e+05, 1.2648582652988912e+05, 1.2649542308458962e+05, 1.2650501970725029e+05, 1.2651461639786651e+05, 1.2652421315643369e+05, 1.2653380998294718e+05, 1.2654340687740236e+05, 1.2655300383979465e+05, 1.2656260087011941e+05, 1.2657219796837204e+05, 1.2658179513454789e+05, 1.2659139236864241e+05, 1.2660098967065093e+05, 1.2661058704056886e+05, 1.2662018447839159e+05, 1.2662978198411450e+05, 1.2663937955773299e+05, 1.2664897719924244e+05, 1.2665857490863826e+05, 1.2666817268591581e+05, 1.2667777053107052e+05, 1.2668736844409777e+05, 1.2669696642499293e+05, 1.2670656447375144e+05, 1.2671616259036864e+05, 1.2672576077483998e+05, 1.2673535902716081e+05, 1.2674495734732655e+05, 1.2675455573533260e+05, 1.2676415419117436e+05, 1.2677375271484719e+05, 1.2678335130634654e+05, 1.2679294996566778e+05, 1.2680254869280632e+05, 1.2681214748775755e+05, 1.2682174635051691e+05, 1.2683134528107978e+05, 1.2684094427944154e+05, 1.2685054334559760e+05, 1.2686014247954340e+05, 1.2686974168127432e+05, 1.2687934095078577e+05, 1.2688894028807315e+05, 1.2689853969313187e+05, 1.2690813916595733e+05, 1.2691773870654493e+05, 1.2692733831489012e+05, 1.2693693799098827e+05, 1.2694653773483481e+05, 1.2695613754642513e+05, 1.2696573742575466e+05, 1.2697533737281882e+05, 1.2698493738761300e+05, 1.2699453747013262e+05, 1.2700413762037308e+05, 1.2701373783832982e+05, 1.2702333812399827e+05, 1.2703293847737378e+05, 1.2704253889845182e+05, 1.2705213938722778e+05, 1.2706173994369710e+05, 1.2707134056785519e+05, 1.2708094125969744e+05, 1.2709054201921933e+05, 1.2710014284641625e+05, 1.2710974374128360e+05, 1.2711934470381682e+05, 1.2712894573401133e+05, 1.2713854683186255e+05, 1.2714814799736590e+05, 1.2715774923051680e+05, 1.2716735053131067e+05, 1.2717695189974295e+05, 1.2718655333580908e+05, 1.2719615483950445e+05, 1.2720575641082450e+05, 1.2721535804976468e+05, 1.2722495975632040e+05, 1.2723456153048706e+05, 1.2724416337226014e+05, 1.2725376528163503e+05, 1.2726336725860719e+05, 1.2727296930317202e+05, 1.2728257141532499e+05, 1.2729217359506147e+05, 1.2730177584237696e+05, 1.2731137815726687e+05, 1.2732098053972662e+05, 1.2733058298975165e+05, 1.2734018550733740e+05, 1.2734978809247933e+05, 1.2735939074517283e+05, 1.2736899346541337e+05, 1.2737859625319639e+05, 1.2738819910851728e+05, 1.2739780203137154e+05, 1.2740740502175457e+05, 1.2741700807966183e+05, 1.2742661120508876e+05, 1.2743621439803079e+05, 1.2744581765848339e+05, 1.2745542098644195e+05, 1.2746502438190197e+05, 1.2747462784485884e+05, 1.2748423137530807e+05, 1.2749383497324503e+05, 1.2750343863866522e+05, 1.2751304237156408e+05, 1.2752264617193703e+05, 1.2753225003977954e+05, 1.2754185397508706e+05, 1.2755145797785500e+05, 1.2756106204807888e+05, 1.2757066618575407e+05, 1.2758027039087610e+05, 1.2758987466344038e+05, 1.2759947900344236e+05, 1.2760908341087747e+05, 1.2761868788574122e+05, 1.2762829242802902e+05, 1.2763789703773634e+05, 1.2764750171485863e+05, 1.2765710645939136e+05, 1.2766671127132997e+05, 1.2767631615066991e+05, 1.2768592109740667e+05, 1.2769552611153566e+05, 1.2770513119305238e+05, 1.2771473634195227e+05, 1.2772434155823079e+05, 1.2773394684188340e+05, 1.2774355219290558e+05, 1.2775315761129279e+05, 1.2776276309704046e+05, 1.2777236865014408e+05, 1.2778197427059911e+05, 1.2779157995840101e+05, 1.2780118571354523e+05, 1.2781079153602727e+05, 1.2782039742584257e+05, 1.2783000338298660e+05, 1.2783960940745482e+05, 1.2784921549924271e+05, 1.2785882165834574e+05, 1.2786842788475938e+05, 1.2787803417847908e+05, 1.2788764053950034e+05, 1.2789724696781862e+05, 1.2790685346342938e+05, 1.2791646002632810e+05, 1.2792606665651024e+05, 1.2793567335397130e+05, 1.2794528011870671e+05, 1.2795488695071200e+05, 1.2796449384998262e+05, 1.2797410081651404e+05, 1.2798370785030173e+05, 1.2799331495134118e+05, 1.2800292211962787e+05, 1.2801252935515727e+05, 1.2802213665792487e+05, 1.2803174402792612e+05, 1.2804135146515653e+05, 1.2805095896961159e+05, 1.2806056654128672e+05, 1.2807017418017748e+05, 1.2807978188627931e+05, 1.2808938965958769e+05, 1.2809899750009811e+05, 1.2810860540780606e+05, 1.2811821338270702e+05, 1.2812782142479648e+05, 1.2813742953406993e+05, 1.2814703771052285e+05, 1.2815664595415072e+05, 1.2816625426494905e+05, 1.2817586264291331e+05, 1.2818547108803898e+05, 1.2819507960032158e+05, 1.2820468817975657e+05, 1.2821429682633947e+05, 1.2822390554006574e+05, 1.2823351432093089e+05, 1.2824312316893041e+05, 1.2825273208405980e+05, 1.2826234106631453e+05, 1.2827195011569012e+05, 1.2828155923218206e+05, 1.2829116841578584e+05, 1.2830077766649696e+05, 1.2831038698431090e+05, 1.2831999636922318e+05, 1.2832960582122930e+05, 1.2833921534032475e+05, 1.2834882492650501e+05, 1.2835843457976561e+05, 1.2836804430010203e+05, 1.2837765408750979e+05, 1.2838726394198438e+05, 1.2839687386352131e+05, 1.2840648385211609e+05, 1.2841609390776418e+05, 1.2842570403046113e+05, 1.2843531422020242e+05, 1.2844492447698356e+05, 1.2845453480080007e+05, 1.2846414519164746e+05, 1.2847375564952122e+05, 1.2848336617441686e+05, 1.2849297676632990e+05, 1.2850258742525583e+05, 1.2851219815119018e+05, 1.2852180894412844e+05, 1.2853141980406614e+05, 1.2854103073099878e+05, 1.2855064172492188e+05, 1.2856025278583093e+05, 1.2856986391372146e+05, 1.2857947510858899e+05, 1.2858908637042902e+05, 1.2859869769923708e+05, 1.2860830909500868e+05, 1.2861792055773933e+05, 1.2862753208742455e+05, 1.2863714368405985e+05, 1.2864675534764078e+05, 1.2865636707816283e+05, 1.2866597887562151e+05, 1.2867559074001237e+05, 1.2868520267133090e+05, 1.2869481466957263e+05, 1.2870442673473312e+05, 1.2871403886680781e+05, 1.2872365106579231e+05, 1.2873326333168207e+05, 1.2874287566447270e+05, 1.2875248806415965e+05, 1.2876210053073846e+05, 1.2877171306420467e+05, 1.2878132566455382e+05, 1.2879093833178141e+05, 1.2880055106588299e+05, 1.2881016386685404e+05, 1.2881977673469012e+05, 1.2882938966938679e+05, 1.2883900267093955e+05, 1.2884861573934391e+05, 1.2885822887459544e+05, 1.2886784207668965e+05, 1.2887745534562207e+05, 1.2888706868138826e+05, 1.2889668208398373e+05, 1.2890629555340404e+05, 1.2891590908964467e+05, 1.2892552269270120e+05, 1.2893513636256916e+05, 1.2894475009924408e+05, 1.2895436390272150e+05, 1.2896397777299696e+05, 1.2897359171006597e+05, 1.2898320571392412e+05, 1.2899281978456692e+05, 1.2900243392198990e+05, 1.2901204812618862e+05, 1.2902166239715861e+05, 1.2903127673489542e+05, 1.2904089113939459e+05, 1.2905050561065166e+05, 1.2906012014866216e+05, 1.2906973475342168e+05, 1.2907934942492572e+05, 1.2908896416316985e+05, 1.2909857896814961e+05, 1.2910819383986051e+05, 1.2911780877829816e+05, 1.2912742378345807e+05, 1.2913703885533579e+05, 1.2914665399392690e+05, 1.2915626919922691e+05, 1.2916588447123137e+05, 1.2917549980993586e+05, 1.2918511521533593e+05, 1.2919473068742710e+05, 1.2920434622620496e+05, 1.2921396183166504e+05, 1.2922357750380292e+05, 1.2923319324261411e+05, 1.2924280904809420e+05, 1.2925242492023873e+05, 1.2926204085904326e+05, 1.2927165686450335e+05, 1.2928127293661455e+05, 1.2929088907537241e+05, 1.2930050528077252e+05, 1.2931012155281041e+05, 1.2931973789148165e+05, 1.2932935429678179e+05, 1.2933897076870641e+05, 1.2934858730725109e+05, 1.2935820391241134e+05, 1.2936782058418276e+05, 1.2937743732256089e+05, 1.2938705412754131e+05, 1.2939667099911957e+05, 1.2940628793729126e+05, 1.2941590494205193e+05, 1.2942552201339712e+05, 1.2943513915132244e+05, 1.2944475635582344e+05, 1.2945437362689567e+05, 1.2946399096453472e+05, 1.2947360836873618e+05, 1.2948322583949559e+05, 1.2949284337680852e+05, 1.2950246098067054e+05, 1.2951207865107723e+05, 1.2952169638802417e+05, 1.2953131419150691e+05, 1.2954093206152106e+05, 1.2955054999806215e+05, 1.2956016800112577e+05, 1.2956978607070752e+05, 1.2957940420680292e+05, 1.2958902240940763e+05, 1.2959864067851714e+05, 1.2960825901412708e+05, 1.2961787741623302e+05, 1.2962749588483053e+05, 1.2963711441991517e+05, 1.2964673302148255e+05, 1.2965635168952825e+05, 1.2966597042404783e+05, 1.2967558922503689e+05, 1.2968520809249101e+05, 1.2969482702640576e+05, 1.2970444602677674e+05, 1.2971406509359954e+05, 1.2972368422686972e+05, 1.2973330342658288e+05, 1.2974292269273459e+05, 1.2975254202532045e+05, 1.2976216142433604e+05, 1.2977178088977694e+05, 1.2978140042163875e+05, 1.2979102001991709e+05, 1.2980063968460749e+05, 1.2981025941570559e+05, 1.2981987921320694e+05, 1.2982949907710715e+05, 1.2983911900740182e+05, 1.2984873900408651e+05, 1.2985835906715685e+05, 1.2986797919660842e+05, 1.2987759939243681e+05, 1.2988721965463761e+05, 1.2989683998320643e+05, 1.2990646037813884e+05, 1.2991608083943046e+05, 1.2992570136707688e+05, 1.2993532196107370e+05, 1.2994494262141650e+05, 1.2995456334810091e+05, 1.2996418414112252e+05, 1.2997380500047690e+05, 1.2998342592615969e+05, 1.2999304691816645e+05, 1.3000266797649281e+05, 1.3001228910113437e+05, 1.3002191029208673e+05, 1.3003153154934550e+05, 1.3004115287290626e+05, 1.3005077426276464e+05, 1.3006039571891623e+05, 1.3007001724135665e+05, 1.3007963883008149e+05, 1.3008926048508637e+05, 1.3009888220636689e+05, 1.3010850399391865e+05, 1.3011812584773728e+05, 1.3012774776781836e+05, 1.3013736975415751e+05, 1.3014699180675036e+05, 1.3015661392559249e+05, 1.3016623611067954e+05, 1.3017585836200711e+05, 1.3018548067957081e+05, 1.3019510306336627e+05, 1.3020472551338907e+05, 1.3021434802963486e+05, 1.3022397061209922e+05, 1.3023359326077778e+05, 1.3024321597566616e+05, 1.3025283875675997e+05, 1.3026246160405484e+05, 1.3027208451754639e+05, 1.3028170749723021e+05, 1.3029133054310195e+05, 1.3030095365515721e+05, 1.3031057683339162e+05, 1.3032020007780081e+05, 1.3032982338838038e+05, 1.3033944676512596e+05, 1.3034907020803317e+05, 1.3035869371709765e+05, 1.3036831729231500e+05, 1.3037794093368083e+05, 1.3038756464119082e+05, 1.3039718841484057e+05, 1.3040681225462568e+05, 1.3041643616054179e+05, 1.3042606013258452e+05, 1.3043568417074955e+05, 1.3044530827503245e+05, 1.3045493244542887e+05, 1.3046455668193444e+05, 1.3047418098454477e+05, 1.3048380535325552e+05, 1.3049342978806231e+05, 1.3050305428896076e+05, 1.3051267885594651e+05, 1.3052230348901519e+05, 1.3053192818816246e+05, 1.3054155295338391e+05, 1.3055117778467522e+05, 1.3056080268203198e+05, 1.3057042764544985e+05, 1.3058005267492447e+05, 1.3058967777045145e+05, 1.3059930293202646e+05, 1.3060892815964510e+05, 1.3061855345330307e+05, 1.3062817881299598e+05, 1.3063780423871944e+05, 1.3064742973046911e+05, 1.3065705528824063e+05, 1.3066668091202965e+05, 1.3067630660183181e+05, 1.3068593235764276e+05, 1.3069555817945811e+05, 1.3070518406727353e+05, 1.3071481002108468e+05, 1.3072443604088716e+05, 1.3073406212667664e+05, 1.3074368827844878e+05, 1.3075331449619921e+05, 1.3076294077992358e+05, 1.3077256712961754e+05, 1.3078219354527675e+05, 1.3079182002689684e+05, 1.3080144657447346e+05, 1.3081107318800226e+05, 1.3082069986747891e+05, 1.3083032661289902e+05, 1.3083995342425829e+05, 1.3084958030155234e+05, 1.3085920724477683e+05, 1.3086883425392743e+05, 1.3087846132899978e+05, 1.3088808846998953e+05, 1.3089771567689233e+05, 1.3090734294970385e+05, 1.3091697028841973e+05, 1.3092659769303566e+05, 1.3093622516354728e+05, 1.3094585269995025e+05, 1.3095548030224020e+05, 1.3096510797041282e+05, 1.3097473570446376e+05, 1.3098436350438870e+05, 1.3099399137018328e+05, 1.3100361930184317e+05, 1.3101324729936401e+05, 1.3102287536274150e+05, 1.3103250349197126e+05, 1.3104213168704900e+05, 1.3105175994797035e+05, 1.3106138827473098e+05, 1.3107101666732658e+05, 1.3108064512575284e+05, 1.3109027365000531e+05, 1.3109990224007980e+05, 1.3110953089597184e+05, 1.3111915961767721e+05, 1.3112878840519153e+05, 1.3113841725851048e+05, 1.3114804617762970e+05, 1.3115767516254491e+05, 1.3116730421325174e+05, 1.3117693332974587e+05, 1.3118656251202300e+05, 1.3119619176007877e+05, 1.3120582107390885e+05, 1.3121545045350899e+05, 1.3122507989887477e+05, 1.3123470941000193e+05, 1.3124433898688608e+05, 1.3125396862952298e+05, 1.3126359833790819e+05, 1.3127322811203750e+05, 1.3128285795190651e+05, 1.3129248785751098e+05, 1.3130211782884656e+05, 1.3131174786590892e+05, 1.3132137796869371e+05, 1.3133100813719659e+05, 1.3134063837141334e+05, 1.3135026867133955e+05, 1.3135989903697098e+05, 1.3136952946830323e+05, 1.3137915996533202e+05, 1.3138879052805310e+05, 1.3139842115646205e+05, 1.3140805185055465e+05, 1.3141768261032648e+05, 1.3142731343577334e+05, 1.3143694432689084e+05, 1.3144657528367470e+05, 1.3145620630612058e+05, 1.3146583739422419e+05, 1.3147546854798126e+05, 1.3148509976738741e+05, 1.3149473105243832e+05, 1.3150436240312972e+05, 1.3151399381945733e+05, 1.3152362530141679e+05, 1.3153325684900381e+05, 1.3154288846221409e+05, 1.3155252014104329e+05, 1.3156215188548720e+05, 1.3157178369554141e+05, 1.3158141557120168e+05, 1.3159104751246367e+05, 1.3160067951932311e+05, 1.3161031159177565e+05, 1.3161994372981702e+05, 1.3162957593344292e+05, 1.3163920820264900e+05, 1.3164884053743104e+05, 1.3165847293778471e+05, 1.3166810540370567e+05, 1.3167773793518968e+05, 1.3168737053223242e+05, 1.3169700319482957e+05, 1.3170663592297686e+05, 1.3171626871666996e+05, 1.3172590157590463e+05, 1.3173553450067656e+05, 1.3174516749098140e+05, 1.3175480054681489e+05, 1.3176443366817274e+05, 1.3177406685505062e+05, 1.3178370010744434e+05, 1.3179333342534953e+05, 1.3180296680876188e+05, 1.3181260025767714e+05, 1.3182223377209101e+05, 1.3183186735199924e+05, 1.3184150099739741e+05, 1.3185113470828134e+05, 1.3186076848464677e+05, 1.3187040232648936e+05, 1.3188003623380480e+05, 1.3188967020658884e+05, 1.3189930424483720e+05, 1.3190893834854555e+05, 1.3191857251770963e+05, 1.3192820675232520e+05, 1.3193784105238793e+05, 1.3194747541789355e+05, 1.3195710984883775e+05, 1.3196674434521626e+05, 1.3197637890702483e+05, 1.3198601353425917e+05, 1.3199564822691496e+05, 1.3200528298498798e+05, 1.3201491780847390e+05, 1.3202455269736846e+05, 1.3203418765166737e+05, 1.3204382267136639e+05, 1.3205345775646120e+05, 1.3206309290694754e+05, 1.3207272812282117e+05, 1.3208236340407777e+05, 1.3209199875071310e+05, 1.3210163416272280e+05, 1.3211126964010269e+05, 1.3212090518284845e+05, 1.3213054079095583e+05, 1.3214017646442054e+05, 1.3214981220323834e+05, 1.3215944800740489e+05, 1.3216908387691603e+05, 1.3217871981176743e+05, 1.3218835581195480e+05, 1.3219799187747386e+05, 1.3220762800832040e+05, 1.3221726420449011e+05, 1.3222690046597878e+05, 1.3223653679278208e+05, 1.3224617318489574e+05, 1.3225580964231558e+05, 1.3226544616503723e+05, 1.3227508275305649e+05, 1.3228471940636908e+05, 1.3229435612497077e+05, 1.3230399290885724e+05, 1.3231362975802427e+05, 1.3232326667246758e+05, 1.3233290365218290e+05, 1.3234254069716597e+05, 1.3235217780741255e+05, 1.3236181498291838e+05, 1.3237145222367920e+05, 1.3238108952969074e+05, 1.3239072690094874e+05, 1.3240036433744893e+05, 1.3241000183918708e+05, 1.3241963940615894e+05, 1.3242927703836025e+05, 1.3243891473578673e+05, 1.3244855249843415e+05, 1.3245819032629827e+05, 1.3246782821937479e+05, 1.3247746617765949e+05, 1.3248710420114815e+05, 1.3249674228983646e+05, 1.3250638044372023e+05, 1.3251601866279519e+05, 1.3252565694705700e+05, 1.3253529529650154e+05, 1.3254493371112450e+05, 1.3255457219092161e+05, 1.3256421073588866e+05, 1.3257384934602139e+05, 1.3258348802131554e+05, 1.3259312676176691e+05, 1.3260276556737121e+05, 1.3261240443812424e+05, 1.3262204337402171e+05, 1.3263168237505940e+05, 1.3264132144123307e+05, 1.3265096057253843e+05, 1.3266059976897133e+05, 1.3267023903052742e+05, 1.3267987835720257e+05, 1.3268951774899245e+05, 1.3269915720589287e+05, 1.3270879672789961e+05, 1.3271843631500835e+05, 1.3272807596721489e+05, 1.3273771568451502e+05, 1.3274735546690450e+05, 1.3275699531437905e+05, 1.3276663522693451e+05, 1.3277627520456654e+05, 1.3278591524727098e+05, 1.3279555535504359e+05, 1.3280519552788010e+05, 1.3281483576577631e+05, 1.3282447606872802e+05, 1.3283411643673092e+05, 1.3284375686978080e+05, 1.3285339736787349e+05, 1.3286303793100471e+05, 1.3287267855917022e+05, 1.3288231925236582e+05, 1.3289196001058727e+05, 1.3290160083383031e+05, 1.3291124172209075e+05, 1.3292088267536438e+05, 1.3293052369364697e+05, 1.3294016477693422e+05, 1.3294980592522200e+05, 1.3295944713850602e+05, 1.3296908841678203e+05, 1.3297872976004591e+05, 1.3298837116829338e+05, 1.3299801264152024e+05, 1.3300765417972222e+05, 1.3301729578289515e+05, 1.3302693745103476e+05, 1.3303657918413682e+05, 1.3304622098219718e+05, 1.3305586284521158e+05, 1.3306550477317581e+05, 1.3307514676608564e+05, 1.3308478882393689e+05, 1.3309443094672528e+05, 1.3310407313444666e+05, 1.3311371538709677e+05, 1.3312335770467139e+05, 1.3313300008716632e+05, 1.3314264253457738e+05, 1.3315228504690027e+05, 1.3316192762413085e+05, 1.3317157026626490e+05, 1.3318121297329816e+05, 1.3319085574522646e+05, 1.3320049858204561e+05, 1.3321014148375133e+05, 1.3321978445033950e+05, 1.3322942748180585e+05, 1.3323907057814617e+05, 1.3324871373935626e+05, 1.3325835696543189e+05, 1.3326800025636892e+05, 1.3327764361216308e+05, 1.3328728703281019e+05, 1.3329693051830604e+05, 1.3330657406864644e+05, 1.3331621768382716e+05, 1.3332586136384401e+05, 1.3333550510869277e+05, 1.3334514891836926e+05, 1.3335479279286930e+05, 1.3336443673218862e+05, 1.3337408073632309e+05, 1.3338372480526843e+05, 1.3339336893902050e+05, 1.3340301313757512e+05, 1.3341265740092803e+05, 1.3342230172907506e+05, 1.3343194612201198e+05, 1.3344159057973465e+05, 1.3345123510223883e+05, 1.3346087968952034e+05, 1.3347052434157496e+05, 1.3348016905839855e+05, 1.3348981383998686e+05, 1.3349945868633571e+05, 1.3350910359744093e+05, 1.3351874857329833e+05, 1.3352839361390367e+05, 1.3353803871925280e+05, 1.3354768388934151e+05, 1.3355732912416564e+05, 1.3356697442372094e+05, 1.3357661978800324e+05, 1.3358626521700842e+05, 1.3359591071073219e+05, 1.3360555626917043e+05, 1.3361520189231893e+05, 1.3362484758017349e+05, 1.3363449333272994e+05, 1.3364413914998408e+05, 1.3365378503193174e+05, 1.3366343097856871e+05, 1.3367307698989086e+05, 1.3368272306589395e+05, 1.3369236920657379e+05, 1.3370201541192626e+05, 1.3371166168194718e+05, 1.3372130801663228e+05, 1.3373095441597747e+05, 1.3374060087997848e+05, 1.3375024740863120e+05, 1.3375989400193145e+05, 1.3376954065987503e+05, 1.3377918738245778e+05, 1.3378883416967548e+05, 1.3379848102152397e+05, 1.3380812793799906e+05, 1.3381777491909661e+05, 1.3382742196481244e+05, 1.3383706907514238e+05, 1.3384671625008219e+05, 1.3385636348962778e+05, 1.3386601079377491e+05, 1.3387565816251945e+05, 1.3388530559585721e+05, 1.3389495309378402e+05, 1.3390460065629569e+05, 1.3391424828338812e+05, 1.3392389597505706e+05, 1.3393354373129836e+05, 1.3394319155210786e+05, 1.3395283943748142e+05, 1.3396248738741482e+05, 1.3397213540190394e+05, 1.3398178348094458e+05, 1.3399143162453256e+05, 1.3400107983266376e+05, 1.3401072810533398e+05, 1.3402037644253910e+05, 1.3403002484427489e+05, 1.3403967331053721e+05, 1.3404932184132191e+05, 1.3405897043662480e+05, 1.3406861909644175e+05, 1.3407826782076858e+05, 1.3408791660960115e+05, 1.3409756546293528e+05, 1.3410721438076679e+05, 1.3411686336309157e+05, 1.3412651240990544e+05, 1.3413616152120422e+05, 1.3414581069698377e+05, 1.3415545993723991e+05, 1.3416510924196852e+05, 1.3417475861116545e+05, 1.3418440804482650e+05, 1.3419405754294756e+05, 1.3420370710552446e+05, 1.3421335673255299e+05, 1.3422300642402904e+05, 1.3423265617994848e+05, 1.3424230600030714e+05, 1.3425195588510091e+05, 1.3426160583432557e+05, 1.3427125584797698e+05, 1.3428090592605097e+05, 1.3429055606854346e+05, 1.3430020627545024e+05, 1.3430985654676723e+05, 1.3431950688249021e+05, 1.3432915728261508e+05, 1.3433880774713765e+05, 1.3434845827605380e+05, 1.3435810886935939e+05, 1.3436775952705022e+05, 1.3437741024912221e+05, 1.3438706103557118e+05, 1.3439671188639302e+05, 1.3440636280158354e+05, 1.3441601378113864e+05, 1.3442566482505415e+05, 1.3443531593332591e+05, 1.3444496710594985e+05, 1.3445461834292175e+05, 1.3446426964423750e+05, 1.3447392100989301e+05, 1.3448357243988407e+05, 1.3449322393420653e+05, 1.3450287549285631e+05, 1.3451252711582926e+05, 1.3452217880312123e+05, 1.3453183055472808e+05, 1.3454148237064568e+05, 1.3455113425086992e+05, 1.3456078619539662e+05, 1.3457043820422169e+05, 1.3458009027734093e+05, 1.3458974241475028e+05, 1.3459939461644556e+05, 1.3460904688242261e+05, 1.3461869921267743e+05, 1.3462835160720575e+05, 1.3463800406600349e+05, 1.3464765658906649e+05, 1.3465730917639070e+05, 1.3466696182797191e+05, 1.3467661454380600e+05, 1.3468626732388890e+05, 1.3469592016821646e+05, 1.3470557307678447e+05, 1.3471522604958891e+05, 1.3472487908662559e+05, 1.3473453218789044e+05, 1.3474418535337929e+05, 1.3475383858308804e+05, 1.3476349187701254e+05, 1.3477314523514867e+05, 1.3478279865749233e+05, 1.3479245214403939e+05, 1.3480210569478571e+05, 1.3481175930972720e+05, 1.3482141298885969e+05, 1.3483106673217914e+05, 1.3484072053968138e+05, 1.3485037441136225e+05, 1.3486002834721771e+05, 1.3486968234724362e+05, 1.3487933641143583e+05, 1.3488899053979022e+05, 1.3489864473230272e+05, 1.3490829898896918e+05, 1.3491795330978549e+05, 1.3492760769474754e+05, 1.3493726214385123e+05, 1.3494691665709240e+05, 1.3495657123446697e+05, 1.3496622587597082e+05, 1.3497588058159989e+05, 1.3498553535134997e+05, 1.3499519018521704e+05, 1.3500484508319691e+05, 1.3501450004528553e+05, 1.3502415507147880e+05, 1.3503381016177253e+05, 1.3504346531616268e+05, 1.3505312053464510e+05, 1.3506277581721576e+05, 1.3507243116387047e+05, 1.3508208657460514e+05, 1.3509174204941571e+05, 1.3510139758829796e+05, 1.3511105319124792e+05, 1.3512070885826147e+05, 1.3513036458933441e+05, 1.3514002038446270e+05, 1.3514967624364223e+05, 1.3515933216686893e+05, 1.3516898815413861e+05, 1.3517864420544726e+05, 1.3518830032079076e+05, 1.3519795650016496e+05, 1.3520761274356581e+05, 1.3521726905098918e+05, 1.3522692542243103e+05, 1.3523658185788718e+05, 1.3524623835735361e+05, 1.3525589492082616e+05, 1.3526555154830078e+05, 1.3527520823977335e+05, 1.3528486499523977e+05, 1.3529452181469597e+05, 1.3530417869813781e+05, 1.3531383564556122e+05, 1.3532349265696210e+05, 1.3533314973233640e+05, 1.3534280687167999e+05, 1.3535246407498879e+05, 1.3536212134225867e+05, 1.3537177867348556e+05, 1.3538143606866541e+05, 1.3539109352779409e+05, 1.3540075105086755e+05, 1.3541040863788166e+05, 1.3542006628883234e+05, 1.3542972400371550e+05, 1.3543938178252708e+05, 1.3544903962526296e+05, 1.3545869753191905e+05, 1.3546835550249126e+05, 1.3547801353697557e+05, 1.3548767163536785e+05, 1.3549732979766399e+05, 1.3550698802385997e+05, 1.3551664631395167e+05, 1.3552630466793501e+05, 1.3553596308580591e+05, 1.3554562156756027e+05, 1.3555528011319402e+05, 1.3556493872270311e+05, 1.3557459739608341e+05, 1.3558425613333087e+05, 1.3559391493444145e+05, 1.3560357379941098e+05, 1.3561323272823545e+05, 1.3562289172091073e+05, 1.3563255077743283e+05, 1.3564220989779761e+05, 1.3565186908200101e+05, 1.3566152833003894e+05, 1.3567118764190734e+05, 1.3568084701760212e+05, 1.3569050645711922e+05, 1.3570016596045456e+05, 1.3570982552760406e+05, 1.3571948515856368e+05, 1.3572914485332934e+05, 1.3573880461189695e+05, 1.3574846443426248e+05, 1.3575812432042183e+05, 1.3576778427037087e+05, 1.3577744428410564e+05, 1.3578710436162198e+05, 1.3579676450291587e+05, 1.3580642470798324e+05, 1.3581608497682004e+05, 1.3582574530942220e+05, 1.3583540570578558e+05, 1.3584506616590620e+05, 1.3585472668977999e+05, 1.3586438727740283e+05, 1.3587404792877071e+05, 1.3588370864387954e+05, 1.3589336942272526e+05, 1.3590303026530382e+05, 1.3591269117161116e+05, 1.3592235214164318e+05, 1.3593201317539590e+05, 1.3594167427286517e+05, 1.3595133543404695e+05, 1.3596099665893725e+05, 1.3597065794753193e+05, 1.3598031929982692e+05, 1.3598998071581827e+05, 1.3599964219550183e+05, 1.3600930373887360e+05, 1.3601896534592952e+05, 1.3602862701666544e+05, 1.3603828875107740e+05, 1.3604795054916132e+05, 1.3605761241091319e+05, 1.3606727433632888e+05, 1.3607693632540439e+05, 1.3608659837813562e+05, 1.3609626049451856e+05, 1.3610592267454916e+05, 1.3611558491822335e+05, 1.3612524722553711e+05, 1.3613490959648634e+05, 1.3614457203106701e+05, 1.3615423452927507e+05, 1.3616389709110648e+05, 1.3617355971655718e+05, 1.3618322240562318e+05, 1.3619288515830037e+05, 1.3620254797458471e+05, 1.3621221085447219e+05, 1.3622187379795875e+05, 1.3623153680504029e+05, 1.3624119987571283e+05, 1.3625086300997232e+05, 1.3626052620781472e+05, 1.3627018946923595e+05, 1.3627985279423202e+05, 1.3628951618279886e+05, 1.3629917963493240e+05, 1.3630884315062867e+05, 1.3631850672988355e+05, 1.3632817037269307e+05, 1.3633783407905314e+05, 1.3634749784895973e+05, 1.3635716168240883e+05, 1.3636682557939636e+05, 1.3637648953991832e+05, 1.3638615356397067e+05, 1.3639581765154938e+05, 1.3640548180265041e+05, 1.3641514601726973e+05, 1.3642481029540327e+05, 1.3643447463704701e+05, 1.3644413904219694e+05, 1.3645380351084899e+05, 1.3646346804299916e+05, 1.3647313263864338e+05, 1.3648279729777767e+05, 1.3649246202039800e+05, 1.3650212680650031e+05, 1.3651179165608058e+05, 1.3652145656913475e+05, 1.3653112154565885e+05, 1.3654078658564881e+05, 1.3655045168910059e+05, 1.3656011685601019e+05, 1.3656978208637363e+05, 1.3657944738018682e+05, 1.3658911273744571e+05, 1.3659877815814634e+05, 1.3660844364228466e+05, 1.3661810918985662e+05, 1.3662777480085826e+05, 1.3663744047528549e+05, 1.3664710621313434e+05, 1.3665677201440072e+05, 1.3666643787908068e+05, 1.3667610380717018e+05, 1.3668576979866516e+05, 1.3669543585356168e+05, 1.3670510197185562e+05, 1.3671476815354303e+05, 1.3672443439861989e+05, 1.3673410070708214e+05, 1.3674376707892577e+05, 1.3675343351414680e+05, 1.3676310001274123e+05, 1.3677276657470496e+05, 1.3678243320003405e+05, 1.3679209988872445e+05, 1.3680176664077214e+05, 1.3681143345617314e+05, 1.3682110033492342e+05, 1.3683076727701898e+05, 1.3684043428245577e+05, 1.3685010135122979e+05, 1.3685976848333707e+05, 1.3686943567877350e+05, 1.3687910293753521e+05, 1.3688877025961806e+05, 1.3689843764501813e+05, 1.3690810509373137e+05, 1.3691777260575379e+05, 1.3692744018108136e+05, 1.3693710781971007e+05, 1.3694677552163595e+05, 1.3695644328685495e+05, 1.3696611111536311e+05, 1.3697577900715638e+05, 1.3698544696223075e+05, 1.3699511498058232e+05, 1.3700478306220696e+05, 1.3701445120710070e+05, 1.3702411941525957e+05, 1.3703378768667957e+05, 1.3704345602135666e+05, 1.3705312441928685e+05, 1.3706279288046615e+05, 1.3707246140489055e+05, 1.3708212999255606e+05, 1.3709179864345872e+05, 1.3710146735759446e+05, 1.3711113613495929e+05, 1.3712080497554922e+05, 1.3713047387936033e+05, 1.3714014284638851e+05, 1.3714981187662986e+05, 1.3715948097008030e+05, 1.3716915012673588e+05, 1.3717881934659262e+05, 1.3718848862964645e+05, 1.3719815797589347e+05, 1.3720782738532967e+05, 1.3721749685795102e+05, 1.3722716639375352e+05, 1.3723683599273319e+05, 1.3724650565488604e+05, 1.3725617538020815e+05, 1.3726584516869544e+05, 1.3727551502034391e+05, 1.3728518493514968e+05, 1.3729485491310866e+05, 1.3730452495421690e+05, 1.3731419505847039e+05, 1.3732386522586516e+05, 1.3733353545639722e+05, 1.3734320575006257e+05, 1.3735287610685726e+05, 1.3736254652677730e+05, 1.3737221700981871e+05, 1.3738188755597745e+05, 1.3739155816524956e+05, 1.3740122883763112e+05, 1.3741089957311805e+05, 1.3742057037170645e+05, 1.3743024123339230e+05, 1.3743991215817159e+05, 1.3744958314604039e+05, 1.3745925419699471e+05, 1.3746892531103056e+05, 1.3747859648814399e+05, 1.3748826772833095e+05, 1.3749793903158754e+05, 1.3750761039790971e+05, 1.3751728182729357e+05, 1.3752695331973504e+05, 1.3753662487523022e+05, 1.3754629649377512e+05, 1.3755596817536579e+05, 1.3756563991999818e+05, 1.3757531172766842e+05, 1.3758498359837240e+05, 1.3759465553210626e+05, 1.3760432752886601e+05, 1.3761399958864765e+05, 1.3762367171144721e+05, 1.3763334389726075e+05, 1.3764301614608426e+05, 1.3765268845791381e+05, 1.3766236083274538e+05, 1.3767203327057505e+05, 1.3768170577139879e+05, 1.3769137833521268e+05, 1.3770105096201276e+05, 1.3771072365179504e+05, 1.3772039640455556e+05, 1.3773006922029037e+05, 1.3773974209899546e+05, 1.3774941504066691e+05, 1.3775908804530077e+05, 1.3776876111289300e+05, 1.3777843424343970e+05, 1.3778810743693690e+05, 1.3779778069338063e+05, 1.3780745401276692e+05, 1.3781712739509181e+05, 1.3782680084035135e+05, 1.3783647434854155e+05, 1.3784614791965848e+05, 1.3785582155369816e+05, 1.3786549525065668e+05, 1.3787516901053002e+05, 1.3788484283331424e+05, 1.3789451671900539e+05, 1.3790419066759953e+05, 1.3791386467909266e+05, 1.3792353875348085e+05, 1.3793321289076019e+05, 1.3794288709092664e+05, 1.3795256135397631e+05, 1.3796223567990522e+05, 1.3797191006870940e+05, 1.3798158452038490e+05, 1.3799125903492776e+05, 1.3800093361233408e+05, 1.3801060825259989e+05, 1.3802028295572122e+05, 1.3802995772169411e+05, 1.3803963255051462e+05, 1.3804930744217883e+05, 1.3805898239668275e+05, 1.3806865741402251e+05, 1.3807833249419407e+05, 1.3808800763719351e+05, 1.3809768284301690e+05, 1.3810735811166023e+05, 1.3811703344311961e+05, 1.3812670883739114e+05, 1.3813638429447077e+05, 1.3814605981435464e+05, 1.3815573539703875e+05, 1.3816541104251920e+05, 1.3817508675079202e+05, 1.3818476252185326e+05, 1.3819443835569904e+05, 1.3820411425232529e+05, 1.3821379021172819e+05, 1.3822346623390377e+05, 1.3823314231884805e+05, 1.3824281846655716e+05, 1.3825249467702711e+05, 1.3826217095025393e+05, 1.3827184728623377e+05, 1.3828152368496262e+05, 1.3829120014643660e+05, 1.3830087667065175e+05, 1.3831055325760410e+05, 1.3832022990728973e+05, 1.3832990661970471e+05, 1.3833958339484513e+05, 1.3834926023270702e+05, 1.3835893713328644e+05, 1.3836861409657952e+05, 1.3837829112258228e+05, 1.3838796821129078e+05, 1.3839764536270112e+05, 1.3840732257680935e+05, 1.3841699985361152e+05, 1.3842667719310374e+05, 1.3843635459528206e+05, 1.3844603206014255e+05, 1.3845570958768125e+05, 1.3846538717789430e+05, 1.3847506483077767e+05, 1.3848474254632753e+05, 1.3849442032453994e+05, 1.3850409816541092e+05, 1.3851377606893660e+05, 1.3852345403511307e+05, 1.3853313206393633e+05, 1.3854281015540252e+05, 1.3855248830950767e+05, 1.3856216652624786e+05, 1.3857184480561921e+05, 1.3858152314761776e+05, 1.3859120155223963e+05, 1.3860088001948086e+05, 1.3861055854933753e+05, 1.3862023714180573e+05, 1.3862991579688151e+05, 1.3863959451456097e+05, 1.3864927329484021e+05, 1.3865895213771533e+05, 1.3866863104318237e+05, 1.3867831001123742e+05, 1.3868798904187657e+05, 1.3869766813509588e+05, 1.3870734729089146e+05, 1.3871702650925939e+05, 1.3872670579019579e+05, 1.3873638513369669e+05, 1.3874606453975820e+05, 1.3875574400837635e+05, 1.3876542353954731e+05, 1.3877510313326717e+05, 1.3878478278953195e+05, 1.3879446250833778e+05, 1.3880414228968075e+05, 1.3881382213355691e+05, 1.3882350203996239e+05, 1.3883318200889335e+05, 1.3884286204034573e+05, 1.3885254213431571e+05, 1.3886222229079934e+05, 1.3887190250979277e+05, 1.3888158279129202e+05, 1.3889126313529324e+05, 1.3890094354179251e+05, 1.3891062401078589e+05, 1.3892030454226956e+05, 1.3892998513623953e+05, 1.3893966579269193e+05, 1.3894934651162283e+05, 1.3895902729302837e+05, 1.3896870813690466e+05, 1.3897838904324771e+05, 1.3898807001205368e+05, 1.3899775104331868e+05, 1.3900743213703882e+05, 1.3901711329321016e+05, 1.3902679451182880e+05, 1.3903647579289082e+05, 1.3904615713639237e+05, 1.3905583854232953e+05, 1.3906552001069844e+05, 1.3907520154149510e+05, 1.3908488313471572e+05, 1.3909456479035638e+05, 1.3910424650841314e+05, 1.3911392828888213e+05, 1.3912361013175949e+05, 1.3913329203704125e+05, 1.3914297400472363e+05, 1.3915265603480258e+05, 1.3916233812727433e+05, 1.3917202028213497e+05, 1.3918170249938057e+05, 1.3919138477900723e+05, 1.3920106712101115e+05, 1.3921074952538835e+05, 1.3922043199213492e+05, 1.3923011452124704e+05, 1.3923979711272081e+05, 1.3924947976655231e+05, 1.3925916248273765e+05, 1.3926884526127303e+05, 1.3927852810215447e+05, 1.3928821100537808e+05, 1.3929789397094000e+05, 1.3930757699883636e+05, 1.3931726008906323e+05, 1.3932694324161680e+05, 1.3933662645649308e+05, 1.3934630973368828e+05, 1.3935599307319848e+05, 1.3936567647501981e+05, 1.3937535993914833e+05, 1.3938504346558027e+05, 1.3939472705431166e+05, 1.3940441070533864e+05, 1.3941409441865733e+05, 1.3942377819426387e+05, 1.3943346203215432e+05, 1.3944314593232487e+05, 1.3945282989477160e+05, 1.3946251391949068e+05, 1.3947219800647817e+05, 1.3948188215573021e+05, 1.3949156636724295e+05, 1.3950125064101254e+05, 1.3951093497703504e+05, 1.3952061937530659e+05, 1.3953030383582332e+05, 1.3953998835858138e+05, 1.3954967294357685e+05, 1.3955935759080588e+05, 1.3956904230026461e+05, 1.3957872707194916e+05, 1.3958841190585567e+05, 1.3959809680198022e+05, 1.3960778176031899e+05, 1.3961746678086810e+05, 1.3962715186362367e+05, 1.3963683700858185e+05, 1.3964652221573878e+05, 1.3965620748509050e+05, 1.3966589281663322e+05, 1.3967557821036308e+05, 1.3968526366627618e+05, 1.3969494918436871e+05, 1.3970463476463672e+05, 1.3971432040707639e+05, 1.3972400611168388e+05, 1.3973369187845528e+05, 1.3974337770738671e+05, 1.3975306359847437e+05, 1.3976274955171437e+05, 1.3977243556710283e+05, 1.3978212164463589e+05, 1.3979180778430973e+05, 1.3980149398612045e+05, 1.3981118025006418e+05, 1.3982086657613708e+05, 1.3983055296433528e+05, 1.3984023941465496e+05, 1.3984992592709220e+05, 1.3985961250164316e+05, 1.3986929913830402e+05, 1.3987898583707088e+05, 1.3988867259793991e+05, 1.3989835942090722e+05, 1.3990804630596898e+05, 1.3991773325312135e+05, 1.3992742026236042e+05, 1.3993710733368239e+05, 1.3994679446708338e+05, 1.3995648166255958e+05, 1.3996616892010710e+05, 1.3997585623972202e+05, 1.3998554362140060e+05, 1.3999523106513891e+05, 1.4000491857093317e+05, 1.4001460613877946e+05, 1.4002429376867396e+05, 1.4003398146061285e+05, 1.4004366921459226e+05, 1.4005335703060831e+05, 1.4006304490865718e+05, 1.4007273284873506e+05, 1.4008242085083801e+05, 1.4009210891496224e+05, 1.4010179704110391e+05, 1.4011148522925915e+05, 1.4012117347942412e+05, 1.4013086179159500e+05, 1.4014055016576793e+05, 1.4015023860193905e+05, 1.4015992710010448e+05, 1.4016961566026049e+05, 1.4017930428240317e+05, 1.4018899296652866e+05, 1.4019868171263314e+05, 1.4020837052071278e+05, 1.4021805939076372e+05, 1.4022774832278214e+05, 1.4023743731676415e+05, 1.4024712637270600e+05, 1.4025681549060377e+05, 1.4026650467045367e+05, 1.4027619391225182e+05, 1.4028588321599443e+05, 1.4029557258167761e+05, 1.4030526200929758e+05, 1.4031495149885051e+05, 1.4032464105033249e+05, 1.4033433066373973e+05, 1.4034402033906835e+05, 1.4035371007631457e+05, 1.4036339987547460e+05, 1.4037308973654450e+05, 1.4038277965952051e+05, 1.4039246964439875e+05, 1.4040215969117545e+05, 1.4041184979984671e+05, 1.4042153997040875e+05, 1.4043123020285773e+05, 1.4044092049718980e+05, 1.4045061085340116e+05, 1.4046030127148793e+05, 1.4046999175144633e+05, 1.4047968229327252e+05, 1.4048937289696268e+05, 1.4049906356251295e+05, 1.4050875428991957e+05, 1.4051844507917864e+05, 1.4052813593028637e+05, 1.4053782684323893e+05, 1.4054751781803250e+05, 1.4055720885466327e+05, 1.4056689995312743e+05, 1.4057659111342108e+05, 1.4058628233554045e+05, 1.4059597361948169e+05, 1.4060566496524104e+05, 1.4061535637281463e+05, 1.4062504784219863e+05, 1.4063473937338925e+05, 1.4064443096638267e+05, 1.4065412262117502e+05, 1.4066381433776251e+05, 1.4067350611614133e+05, 1.4068319795630773e+05, 1.4069288985825781e+05, 1.4070258182198773e+05, 1.4071227384749372e+05, 1.4072196593477196e+05, 1.4073165808381862e+05, 1.4074135029462990e+05, 1.4075104256720198e+05, 1.4076073490153105e+05, 1.4077042729761329e+05, 1.4078011975544487e+05, 1.4078981227502204e+05, 1.4079950485634091e+05, 1.4080919749939771e+05, 1.4081889020418865e+05, 1.4082858297070986e+05, 1.4083827579895756e+05, 1.4084796868892794e+05, 1.4085766164061721e+05, 1.4086735465402150e+05, 1.4087704772913706e+05, 1.4088674086596011e+05, 1.4089643406448676e+05, 1.4090612732471325e+05, 1.4091582064663575e+05, 1.4092551403025049e+05, 1.4093520747555365e+05, 1.4094490098254138e+05, 1.4095459455120994e+05, 1.4096428818155549e+05, 1.4097398187357423e+05, 1.4098367562726236e+05, 1.4099336944261607e+05, 1.4100306331963153e+05, 1.4101275725830504e+05, 1.4102245125863273e+05, 1.4103214532061078e+05, 1.4104183944423540e+05, 1.4105153362950281e+05, 1.4106122787640919e+05, 1.4107092218495076e+05, 1.4108061655512368e+05, 1.4109031098692422e+05, 1.4110000548034854e+05, 1.4110970003539286e+05, 1.4111939465205336e+05, 1.4112908933032627e+05, 1.4113878407020777e+05, 1.4114847887169410e+05, 1.4115817373478139e+05, 1.4116786865946592e+05, 1.4117756364574388e+05, 1.4118725869361145e+05, 1.4119695380306485e+05, 1.4120664897410030e+05, 1.4121634420671401e+05, 1.4122603950090217e+05, 1.4123573485666103e+05, 1.4124543027398674e+05, 1.4125512575287552e+05, 1.4126482129332362e+05, 1.4127451689532725e+05, 1.4128421255888257e+05, 1.4129390828398583e+05, 1.4130360407063324e+05, 1.4131329991882099e+05, 1.4132299582854533e+05, 1.4133269179980241e+05, 1.4134238783258852e+05, 1.4135208392689982e+05, 1.4136178008273256e+05, 1.4137147630008296e+05, 1.4138117257894718e+05, 1.4139086891932148e+05, 1.4140056532120213e+05, 1.4141026178458522e+05, 1.4141995830946704e+05, 1.4142965489584382e+05, 1.4143935154371176e+05, 1.4144904825306707e+05, 1.4145874502390597e+05, 1.4146844185622470e+05, 1.4147813875001948e+05, 1.4148783570528650e+05, 1.4149753272202201e+05, 1.4150722980022224e+05, 1.4151692693988336e+05, 1.4152662414100167e+05, 1.4153632140357330e+05, 1.4154601872759455e+05, 1.4155571611306159e+05, 1.4156541355997071e+05, 1.4157511106831810e+05, 1.4158480863809996e+05, 1.4159450626931258e+05, 1.4160420396195215e+05, 1.4161390171601487e+05, 1.4162359953149702e+05, 1.4163329740839478e+05, 1.4164299534670438e+05, 1.4165269334642211e+05, 1.4166239140754414e+05, 1.4167208953006673e+05, 1.4168178771398606e+05, 1.4169148595929841e+05, 1.4170118426600003e+05, 1.4171088263408709e+05, 1.4172058106355587e+05, 1.4173027955440257e+05, 1.4173997810662346e+05, 1.4174967672021472e+05, 1.4175937539517262e+05, 1.4176907413149340e+05, 1.4177877292917328e+05, 1.4178847178820847e+05, 1.4179817070859531e+05, 1.4180786969032991e+05, 1.4181756873340855e+05, 1.4182726783782750e+05, 1.4183696700358295e+05, 1.4184666623067117e+05, 1.4185636551908840e+05, 1.4186606486883084e+05, 1.4187576427989479e+05, 1.4188546375227644e+05, 1.4189516328597208e+05, 1.4190486288097789e+05, 1.4191456253729016e+05, 1.4192426225490507e+05, 1.4193396203381894e+05, 1.4194366187402798e+05, 1.4195336177552841e+05, 1.4196306173831649e+05, 1.4197276176238849e+05, 1.4198246184774057e+05, 1.4199216199436909e+05, 1.4200186220227025e+05, 1.4201156247144027e+05, 1.4202126280187545e+05, 1.4203096319357195e+05, 1.4204066364652611e+05, 1.4205036416073408e+05, 1.4206006473619223e+05, 1.4206976537289668e+05, 1.4207946607084380e+05, 1.4208916683002978e+05, 1.4209886765045082e+05, 1.4210856853210327e+05, 1.4211826947498333e+05, 1.4212797047908721e+05, 1.4213767154441122e+05, 1.4214737267095162e+05, 1.4215707385870465e+05, 1.4216677510766656e+05, 1.4217647641783359e+05, 1.4218617778920199e+05, 1.4219587922176806e+05, 1.4220558071552799e+05, 1.4221528227047808e+05, 1.4222498388661456e+05, 1.4223468556393369e+05, 1.4224438730243174e+05, 1.4225408910210503e+05, 1.4226379096294972e+05, 1.4227349288496209e+05, 1.4228319486813841e+05, 1.4229289691247494e+05, 1.4230259901796796e+05, 1.4231230118461369e+05, 1.4232200341240840e+05, 1.4233170570134840e+05, 1.4234140805142987e+05, 1.4235111046264911e+05, 1.4236081293500241e+05, 1.4237051546848597e+05, 1.4238021806309614e+05, 1.4238992071882915e+05, 1.4239962343568122e+05, 1.4240932621364866e+05, 1.4241902905272771e+05, 1.4242873195291462e+05, 1.4243843491420572e+05, 1.4244813793659720e+05, 1.4245784102008541e+05, 1.4246754416466656e+05, 1.4247724737033693e+05, 1.4248695063709278e+05, 1.4249665396493039e+05, 1.4250635735384602e+05, 1.4251606080383595e+05, 1.4252576431489646e+05, 1.4253546788702381e+05, 1.4254517152021427e+05, 1.4255487521446409e+05, 1.4256457896976959e+05, 1.4257428278612698e+05, 1.4258398666353259e+05, 1.4259369060198267e+05, 1.4260339460147350e+05, 1.4261309866200131e+05, 1.4262280278356242e+05, 1.4263250696615313e+05, 1.4264221120976968e+05, 1.4265191551440835e+05, 1.4266161988006541e+05, 1.4267132430673717e+05, 1.4268102879441984e+05, 1.4269073334310975e+05, 1.4270043795280319e+05, 1.4271014262349642e+05, 1.4271984735518569e+05, 1.4272955214786730e+05, 1.4273925700153754e+05, 1.4274896191619273e+05, 1.4275866689182905e+05, 1.4276837192844288e+05, 1.4277807702603043e+05, 1.4278778218458800e+05, 1.4279748740411192e+05, 1.4280719268459841e+05, 1.4281689802604381e+05, 1.4282660342844436e+05, 1.4283630889179639e+05, 1.4284601441609615e+05, 1.4285572000133991e+05, 1.4286542564752401e+05, 1.4287513135464469e+05, 1.4288483712269823e+05, 1.4289454295168098e+05, 1.4290424884158914e+05, 1.4291395479241910e+05, 1.4292366080416707e+05, 1.4293336687682936e+05, 1.4294307301040224e+05, 1.4295277920488207e+05, 1.4296248546026507e+05, 1.4297219177654761e+05, 1.4298189815372587e+05, 1.4299160459179617e+05, 1.4300131109075490e+05, 1.4301101765059825e+05, 1.4302072427132254e+05, 1.4303043095292407e+05, 1.4304013769539914e+05, 1.4304984449874406e+05, 1.4305955136295510e+05, 1.4306925828802859e+05, 1.4307896527396076e+05, 1.4308867232074792e+05, 1.4309837942838640e+05, 1.4310808659687248e+05, 1.4311779382620248e+05, 1.4312750111637267e+05, 1.4313720846737942e+05, 1.4314691587921893e+05, 1.4315662335188757e+05, 1.4316633088538158e+05, 1.4317603847969731e+05, 1.4318574613483102e+05, 1.4319545385077901e+05, 1.4320516162753766e+05, 1.4321486946510317e+05, 1.4322457736347194e+05, 1.4323428532264018e+05, 1.4324399334260428e+05, 1.4325370142336047e+05, 1.4326340956490507e+05, 1.4327311776723442e+05, 1.4328282603034482e+05, 1.4329253435423254e+05, 1.4330224273889392e+05, 1.4331195118432527e+05, 1.4332165969052288e+05, 1.4333136825748303e+05, 1.4334107688520208e+05, 1.4335078557367632e+05, 1.4336049432290200e+05, 1.4337020313287553e+05, 1.4337991200359317e+05, 1.4338962093505121e+05, 1.4339932992724603e+05, 1.4340903898017388e+05, 1.4341874809383109e+05, 1.4342845726821400e+05, 1.4343816650331888e+05, 1.4344787579914203e+05, 1.4345758515567979e+05, 1.4346729457292848e+05, 1.4347700405088442e+05, 1.4348671358954394e+05, 1.4349642318890331e+05, 1.4350613284895883e+05, 1.4351584256970690e+05, 1.4352555235114376e+05, 1.4353526219326578e+05, 1.4354497209606922e+05, 1.4355468205955048e+05, 1.4356439208370578e+05, 1.4357410216853151e+05, 1.4358381231402393e+05, 1.4359352252017945e+05, 1.4360323278699431e+05, 1.4361294311446484e+05, 1.4362265350258743e+05, 1.4363236395135833e+05, 1.4364207446077390e+05, 1.4365178503083042e+05, 1.4366149566152427e+05, 1.4367120635285170e+05, 1.4368091710480908e+05, 1.4369062791739276e+05, 1.4370033879059899e+05, 1.4371004972442417e+05, 1.4371976071886462e+05, 1.4372947177391662e+05, 1.4373918288957648e+05, 1.4374889406584061e+05, 1.4375860530270528e+05, 1.4376831660016681e+05, 1.4377802795822156e+05, 1.4378773937686588e+05, 1.4379745085609602e+05, 1.4380716239590835e+05, 1.4381687399629920e+05, 1.4382658565726492e+05, 1.4383629737880189e+05, 1.4384600916090631e+05, 1.4385572100357458e+05, 1.4386543290680303e+05, 1.4387514487058800e+05, 1.4388485689492582e+05, 1.4389456897981284e+05, 1.4390428112524535e+05, 1.4391399333121971e+05, 1.4392370559773227e+05, 1.4393341792477935e+05, 1.4394313031235730e+05, 1.4395284276046243e+05, 1.4396255526909107e+05, 1.4397226783823955e+05, 1.4398198046790430e+05, 1.4399169315808156e+05, 1.4400140590876772e+05, 1.4401111871995905e+05, 1.4402083159165198e+05, 1.4403054452384281e+05, 1.4404025751652784e+05, 1.4404997056970350e+05, 1.4405968368336605e+05, 1.4406939685751186e+05, 1.4407911009213727e+05, 1.4408882338723863e+05, 1.4409853674281231e+05, 1.4410825015885459e+05, 1.4411796363536181e+05, 1.4412767717233041e+05, 1.4413739076975666e+05, 1.4414710442763689e+05, 1.4415681814596750e+05, 1.4416653192474484e+05, 1.4417624576396515e+05, 1.4418595966362490e+05, 1.4419567362372042e+05, 1.4420538764424800e+05, 1.4421510172520400e+05, 1.4422481586658480e+05, 1.4423453006838678e+05, 1.4424424433060622e+05, 1.4425395865323950e+05, 1.4426367303628294e+05, 1.4427338747973295e+05, 1.4428310198358580e+05, 1.4429281654783789e+05, 1.4430253117248559e+05, 1.4431224585752521e+05, 1.4432196060295316e+05, 1.4433167540876573e+05, 1.4434139027495933e+05, 1.4435110520153027e+05, 1.4436082018847496e+05, 1.4437053523578969e+05, 1.4438025034347083e+05, 1.4438996551151475e+05, 1.4439968073991782e+05, 1.4440939602867636e+05, 1.4441911137778679e+05, 1.4442882678724540e+05, 1.4443854225704860e+05, 1.4444825778719271e+05, 1.4445797337767412e+05, 1.4446768902848920e+05, 1.4447740473963425e+05, 1.4448712051110569e+05, 1.4449683634289986e+05, 1.4450655223501311e+05, 1.4451626818744180e+05, 1.4452598420018231e+05, 1.4453570027323096e+05, 1.4454541640658421e+05, 1.4455513260023831e+05, 1.4456484885418968e+05, 1.4457456516843470e+05, 1.4458428154296969e+05, 1.4459399797779106e+05, 1.4460371447289517e+05, 1.4461343102827834e+05, 1.4462314764393700e+05, 1.4463286431986743e+05, 1.4464258105606612e+05, 1.4465229785252930e+05, 1.4466201470925350e+05, 1.4467173162623492e+05, 1.4468144860347002e+05, 1.4469116564095515e+05, 1.4470088273868669e+05, 1.4471059989666106e+05, 1.4472031711487458e+05, 1.4473003439332359e+05, 1.4473975173200452e+05, 1.4474946913091367e+05, 1.4475918659004749e+05, 1.4476890410940233e+05, 1.4477862168897453e+05, 1.4478833932876048e+05, 1.4479805702875656e+05, 1.4480777478895918e+05, 1.4481749260936465e+05, 1.4482721048996941e+05, 1.4483692843076980e+05, 1.4484664643176220e+05, 1.4485636449294299e+05, 1.4486608261430854e+05, 1.4487580079585523e+05, 1.4488551903757948e+05, 1.4489523733947758e+05, 1.4490495570154596e+05, 1.4491467412378101e+05, 1.4492439260617911e+05, 1.4493411114873664e+05, 1.4494382975144999e+05, 1.4495354841431550e+05, 1.4496326713732962e+05, 1.4497298592048869e+05, 1.4498270476378905e+05, 1.4499242366722715e+05, 1.4500214263079932e+05, 1.4501186165450199e+05, 1.4502158073833151e+05, 1.4503129988228428e+05, 1.4504101908635671e+05, 1.4505073835054514e+05, 1.4506045767484599e+05, 1.4507017705925563e+05, 1.4507989650377046e+05, 1.4508961600838689e+05, 1.4509933557310124e+05, 1.4510905519790997e+05, 1.4511877488280940e+05, 1.4512849462779603e+05, 1.4513821443286608e+05, 1.4514793429801610e+05, 1.4515765422324243e+05, 1.4516737420854141e+05, 1.4517709425390948e+05, 1.4518681435934306e+05, 1.4519653452483847e+05, 1.4520625475039214e+05, 1.4521597503600045e+05, 1.4522569538165981e+05, 1.4523541578736663e+05, 1.4524513625311726e+05, 1.4525485677890811e+05, 1.4526457736473560e+05, 1.4527429801059610e+05, 1.4528401871648600e+05, 1.4529373948240172e+05, 1.4530346030833962e+05, 1.4531318119429614e+05, 1.4532290214026766e+05, 1.4533262314625061e+05, 1.4534234421224133e+05, 1.4535206533823628e+05, 1.4536178652423187e+05, 1.4537150777022436e+05, 1.4538122907621035e+05, 1.4539095044218606e+05, 1.4540067186814800e+05, 1.4541039335409258e+05, 1.4542011490001617e+05, 1.4542983650591510e+05, 1.4543955817178587e+05, 1.4544927989762489e+05, 1.4545900168342851e+05, 1.4546872352919317e+05, 1.4547844543491525e+05, 1.4548816740059116e+05, 1.4549788942621733e+05, 1.4550761151179014e+05, 1.4551733365730601e+05, 1.4552705586276134e+05, 1.4553677812815251e+05, 1.4554650045347598e+05, 1.4555622283872811e+05, 1.4556594528390534e+05, 1.4557566778900407e+05, 1.4558539035402075e+05, 1.4559511297895171e+05, 1.4560483566379343e+05, 1.4561455840854227e+05, 1.4562428121319468e+05, 1.4563400407774706e+05, 1.4564372700219581e+05, 1.4565344998653737e+05, 1.4566317303076814e+05, 1.4567289613488450e+05, 1.4568261929888290e+05, 1.4569234252275978e+05, 1.4570206580651150e+05, 1.4571178915013449e+05, 1.4572151255362516e+05, 1.4573123601697996e+05, 1.4574095954019530e+05, 1.4575068312326755e+05, 1.4576040676619316e+05, 1.4577013046896853e+05, 1.4577985423159014e+05, 1.4578957805405435e+05, 1.4579930193635757e+05, 1.4580902587849626e+05, 1.4581874988046681e+05, 1.4582847394226567e+05, 1.4583819806388923e+05, 1.4584792224533396e+05, 1.4585764648659620e+05, 1.4586737078767244e+05, 1.4587709514855908e+05, 1.4588681956925252e+05, 1.4589654404974924e+05, 1.4590626859004566e+05, 1.4591599319013814e+05, 1.4592571785002312e+05, 1.4593544256969704e+05, 1.4594516734915637e+05, 1.4595489218839753e+05, 1.4596461708741685e+05, 1.4597434204621080e+05, 1.4598406706477585e+05, 1.4599379214310841e+05, 1.4600351728120490e+05, 1.4601324247906177e+05, 1.4602296773667540e+05, 1.4603269305404229e+05, 1.4604241843115879e+05, 1.4605214386802138e+05, 1.4606186936462647e+05, 1.4607159492097053e+05, 1.4608132053704993e+05, 1.4609104621286117e+05, 1.4610077194840065e+05, 1.4611049774366475e+05, 1.4612022359864996e+05, 1.4612994951335274e+05, 1.4613967548776948e+05, 1.4614940152189659e+05, 1.4615912761573054e+05, 1.4616885376926782e+05, 1.4617857998250477e+05, 1.4618830625543787e+05, 1.4619803258806359e+05, 1.4620775898037828e+05, 1.4621748543237845e+05, 1.4622721194406055e+05, 1.4623693851542094e+05, 1.4624666514645613e+05, 1.4625639183716255e+05, 1.4626611858753656e+05, 1.4627584539757468e+05, 1.4628557226727335e+05, 1.4629529919662903e+05, 1.4630502618563807e+05, 1.4631475323429695e+05, 1.4632448034260215e+05, 1.4633420751055007e+05, 1.4634393473813718e+05, 1.4635366202535992e+05, 1.4636338937221473e+05, 1.4637311677869808e+05, 1.4638284424480636e+05, 1.4639257177053604e+05, 1.4640229935588359e+05, 1.4641202700084544e+05, 1.4642175470541802e+05, 1.4643148246959780e+05, 1.4644121029338118e+05, 1.4645093817676467e+05, 1.4646066611974471e+05, 1.4647039412231772e+05, 1.4648012218448016e+05, 1.4648985030622847e+05, 1.4649957848755908e+05, 1.4650930672846848e+05, 1.4651903502895311e+05, 1.4652876338900943e+05, 1.4653849180863390e+05, 1.4654822028782294e+05, 1.4655794882657300e+05, 1.4656767742488056e+05, 1.4657740608274206e+05, 1.4658713480015393e+05, 1.4659686357711264e+05, 1.4660659241361468e+05, 1.4661632130965646e+05, 1.4662605026523449e+05, 1.4663577928034513e+05, 1.4664550835498495e+05, 1.4665523748915037e+05, 1.4666496668283781e+05, 1.4667469593604375e+05, 1.4668442524876460e+05, 1.4669415462099691e+05, 1.4670388405273706e+05, 1.4671361354398151e+05, 1.4672334309472679e+05, 1.4673307270496929e+05, 1.4674280237470556e+05, 1.4675253210393197e+05, 1.4676226189264498e+05, 1.4677199174084110e+05, 1.4678172164851677e+05, 1.4679145161566845e+05, 1.4680118164229265e+05, 1.4681091172838578e+05, 1.4682064187394429e+05, 1.4683037207896469e+05, 1.4684010234344340e+05, 1.4684983266737693e+05, 1.4685956305076173e+05, 1.4686929349359428e+05, 1.4687902399587099e+05, 1.4688875455758837e+05, 1.4689848517874288e+05, 1.4690821585933099e+05, 1.4691794659934918e+05, 1.4692767739879389e+05, 1.4693740825766156e+05, 1.4694713917594875e+05, 1.4695687015365186e+05, 1.4696660119076740e+05, 1.4697633228729182e+05, 1.4698606344322159e+05, 1.4699579465855318e+05, 1.4700552593328306e+05, 1.4701525726740769e+05, 1.4702498866092358e+05, 1.4703472011382721e+05, 1.4704445162611501e+05, 1.4705418319778345e+05, 1.4706391482882906e+05, 1.4707364651924826e+05, 1.4708337826903752e+05, 1.4709311007819333e+05, 1.4710284194671223e+05, 1.4711257387459060e+05, 1.4712230586182501e+05, 1.4713203790841182e+05, 1.4714177001434760e+05, 1.4715150217962882e+05, 1.4716123440425191e+05, 1.4717096668821338e+05, 1.4718069903150966e+05, 1.4719043143413731e+05, 1.4720016389609282e+05, 1.4720989641737260e+05, 1.4721962899797314e+05, 1.4722936163789095e+05, 1.4723909433712251e+05, 1.4724882709566431e+05, 1.4725855991351278e+05, 1.4726829279066445e+05, 1.4727802572711580e+05, 1.4728775872286328e+05, 1.4729749177790343e+05, 1.4730722489223268e+05, 1.4731695806584755e+05, 1.4732669129874447e+05, 1.4733642459092004e+05, 1.4734615794237063e+05, 1.4735589135309280e+05, 1.4736562482308300e+05, 1.4737535835233773e+05, 1.4738509194085348e+05, 1.4739482558862676e+05, 1.4740455929565398e+05, 1.4741429306193170e+05, 1.4742402688745639e+05, 1.4743376077222454e+05, 1.4744349471623264e+05, 1.4745322871947716e+05, 1.4746296278195464e+05, 1.4747269690366159e+05, 1.4748243108459437e+05, 1.4749216532474963e+05, 1.4750189962412376e+05, 1.4751163398271328e+05, 1.4752136840051471e+05, 1.4753110287752451e+05, 1.4754083741373921e+05, 1.4755057200915526e+05, 1.4756030666376915e+05, 1.4757004137757744e+05, 1.4757977615057660e+05, 1.4758951098276311e+05, 1.4759924587413346e+05, 1.4760898082468417e+05, 1.4761871583441176e+05, 1.4762845090331268e+05, 1.4763818603138343e+05, 1.4764792121862053e+05, 1.4765765646502047e+05, 1.4766739177057976e+05, 1.4767712713529490e+05, 1.4768686255916234e+05, 1.4769659804217870e+05, 1.4770633358434038e+05, 1.4771606918564395e+05, 1.4772580484608584e+05, 1.4773554056566261e+05, 1.4774527634437074e+05, 1.4775501218220670e+05, 1.4776474807916704e+05, 1.4777448403524823e+05, 1.4778422005044683e+05, 1.4779395612475934e+05, 1.4780369225818219e+05, 1.4781342845071197e+05, 1.4782316470234515e+05, 1.4783290101307817e+05, 1.4784263738290765e+05, 1.4785237381183007e+05, 1.4786211029984188e+05, 1.4787184684693968e+05, 1.4788158345311994e+05, 1.4789132011837914e+05, 1.4790105684271379e+05, 1.4791079362612046e+05, 1.4792053046859556e+05, 1.4793026737013567e+05, 1.4794000433073731e+05, 1.4794974135039700e+05, 1.4795947842911122e+05, 1.4796921556687646e+05, 1.4797895276368930e+05, 1.4798869001954616e+05, 1.4799842733444364e+05, 1.4800816470837823e+05, 1.4801790214134639e+05, 1.4802763963334475e+05, 1.4803737718436975e+05, 1.4804711479441790e+05, 1.4805685246348576e+05, 1.4806659019156979e+05, 1.4807632797866658e+05, 1.4808606582477258e+05, 1.4809580372988433e+05, 1.4810554169399830e+05, 1.4811527971711109e+05, 1.4812501779921923e+05, 1.4813475594031916e+05, 1.4814449414040745e+05, 1.4815423239948065e+05, 1.4816397071753521e+05, 1.4817370909456766e+05, 1.4818344753057457e+05, 1.4819318602555242e+05, 1.4820292457949781e+05, 1.4821266319240714e+05, 1.4822240186427702e+05, 1.4823214059510393e+05, 1.4824187938488441e+05, 1.4825161823361501e+05, 1.4826135714129222e+05, 1.4827109610791257e+05, 1.4828083513347260e+05, 1.4829057421796888e+05, 1.4830031336139786e+05, 1.4831005256375609e+05, 1.4831979182504010e+05, 1.4832953114524644e+05, 1.4833927052437159e+05, 1.4834900996241215e+05, 1.4835874945936454e+05, 1.4836848901522538e+05, 1.4837822862999118e+05, 1.4838796830365848e+05, 1.4839770803622378e+05, 1.4840744782768365e+05, 1.4841718767803459e+05, 1.4842692758727315e+05, 1.4843666755539586e+05, 1.4844640758239923e+05, 1.4845614766827985e+05, 1.4846588781303418e+05, 1.4847562801665877e+05, 1.4848536827915022e+05, 1.4849510860050499e+05, 1.4850484898071966e+05, 1.4851458941979072e+05, 1.4852432991771479e+05, 1.4853407047448828e+05, 1.4854381109010783e+05, 1.4855355176456997e+05, 1.4856329249787121e+05, 1.4857303329000808e+05, 1.4858277414097713e+05, 1.4859251505077491e+05, 1.4860225601939793e+05, 1.4861199704684276e+05, 1.4862173813310594e+05, 1.4863147927818404e+05, 1.4864122048207349e+05, 1.4865096174477096e+05, 1.4866070306627292e+05, 1.4867044444657594e+05, 1.4868018588567653e+05, 1.4868992738357125e+05, 1.4869966894025661e+05, 1.4870941055572924e+05, 1.4871915222998560e+05, 1.4872889396302230e+05, 1.4873863575483582e+05, 1.4874837760542278e+05, 1.4875811951477965e+05, 1.4876786148290304e+05, 1.4877760350978942e+05, 1.4878734559543541e+05, 1.4879708773983756e+05, 1.4880682994299236e+05, 1.4881657220489639e+05, 1.4882631452554619e+05, 1.4883605690493830e+05, 1.4884579934306935e+05, 1.4885554183993576e+05, 1.4886528439553417e+05, 1.4887502700986111e+05, 1.4888476968291315e+05, 1.4889451241468682e+05, 1.4890425520517866e+05, 1.4891399805438524e+05, 1.4892374096230313e+05, 1.4893348392892879e+05, 1.4894322695425892e+05, 1.4895297003828996e+05, 1.4896271318101848e+05, 1.4897245638244110e+05, 1.4898219964255430e+05, 1.4899194296135471e+05, 1.4900168633883880e+05, 1.4901142977500323e+05, 1.4902117326984444e+05, 1.4903091682335906e+05, 1.4904066043554363e+05, 1.4905040410639471e+05, 1.4906014783590887e+05, 1.4906989162408264e+05, 1.4907963547091259e+05, 1.4908937937639526e+05, 1.4909912334052727e+05, 1.4910886736330515e+05, 1.4911861144472545e+05, 1.4912835558478476e+05, 1.4913809978347961e+05, 1.4914784404080658e+05, 1.4915758835676219e+05, 1.4916733273134308e+05, 1.4917707716454574e+05, 1.4918682165636681e+05, 1.4919656620680276e+05, 1.4920631081585024e+05, 1.4921605548350580e+05, 1.4922580020976593e+05, 1.4923554499462730e+05, 1.4924528983808641e+05, 1.4925503474013985e+05, 1.4926477970078416e+05, 1.4927452472001594e+05, 1.4928426979783172e+05, 1.4929401493422809e+05, 1.4930376012920163e+05, 1.4931350538274887e+05, 1.4932325069486647e+05, 1.4933299606555089e+05, 1.4934274149479877e+05, 1.4935248698260664e+05, 1.4936223252897113e+05, 1.4937197813388874e+05, 1.4938172379735607e+05, 1.4939146951936971e+05, 1.4940121529992623e+05, 1.4941096113902217e+05, 1.4942070703665409e+05, 1.4943045299281864e+05, 1.4944019900751236e+05, 1.4944994508073179e+05, 1.4945969121247358e+05, 1.4946943740273421e+05, 1.4947918365151034e+05, 1.4948892995879852e+05, 1.4949867632459529e+05, 1.4950842274889725e+05, 1.4951816923170097e+05, 1.4952791577300304e+05, 1.4953766237280003e+05, 1.4954740903108855e+05, 1.4955715574786515e+05, 1.4956690252312639e+05, 1.4957664935686887e+05, 1.4958639624908919e+05, 1.4959614319978390e+05, 1.4960589020894960e+05, 1.4961563727658291e+05, 1.4962538440268030e+05, 1.4963513158723843e+05, 1.4964487883025387e+05, 1.4965462613172320e+05, 1.4966437349164303e+05, 1.4967412091000992e+05, 1.4968386838682048e+05, 1.4969361592207124e+05, 1.4970336351575883e+05, 1.4971311116787980e+05, 1.4972285887843076e+05, 1.4973260664740831e+05, 1.4974235447480899e+05, 1.4975210236062942e+05, 1.4976185030486618e+05, 1.4977159830751587e+05, 1.4978134636857509e+05, 1.4979109448804037e+05, 1.4980084266590836e+05, 1.4981059090217561e+05, 1.4982033919683876e+05, 1.4983008754989435e+05, 1.4983983596133898e+05, 1.4984958443116923e+05, 1.4985933295938175e+05, 1.4986908154597308e+05, 1.4987883019093980e+05, 1.4988857889427853e+05, 1.4989832765598586e+05, 1.4990807647605840e+05, 1.4991782535449267e+05, 1.4992757429128539e+05, 1.4993732328643303e+05, 1.4994707233993229e+05, 1.4995682145177966e+05, 1.4996657062197183e+05, 1.4997631985050533e+05, 1.4998606913737679e+05, 1.4999581848258281e+05, 1.5000556788611994e+05, 1.5001531734798485e+05, 1.5002506686817409e+05, 1.5003481644668430e+05, 1.5004456608351201e+05, 1.5005431577865384e+05, 1.5006406553210644e+05, 1.5007381534386639e+05, 1.5008356521393027e+05, 1.5009331514229465e+05, 1.5010306512895622e+05, 1.5011281517391151e+05, 1.5012256527715712e+05, 1.5013231543868972e+05, 1.5014206565850586e+05, 1.5015181593660210e+05, 1.5016156627297518e+05, 1.5017131666762158e+05, 1.5018106712053795e+05, 1.5019081763172089e+05, 1.5020056820116699e+05, 1.5021031882887287e+05, 1.5022006951483514e+05, 1.5022982025905041e+05, 1.5023957106151528e+05, 1.5024932192222634e+05, 1.5025907284118026e+05, 1.5026882381837355e+05, 1.5027857485380289e+05, 1.5028832594746488e+05, 1.5029807709935610e+05, 1.5030782830947320e+05, 1.5031757957781278e+05, 1.5032733090437145e+05, 1.5033708228914576e+05, 1.5034683373213239e+05, 1.5035658523332796e+05, 1.5036633679272904e+05, 1.5037608841033222e+05, 1.5038584008613415e+05, 1.5039559182013146e+05, 1.5040534361232075e+05, 1.5041509546269866e+05, 1.5042484737126171e+05, 1.5043459933800661e+05, 1.5044435136292994e+05, 1.5045410344602831e+05, 1.5046385558729834e+05, 1.5047360778673668e+05, 1.5048336004433988e+05, 1.5049311236010460e+05, 1.5050286473402745e+05, 1.5051261716610508e+05, 1.5052236965633405e+05, 1.5053212220471102e+05, 1.5054187481123261e+05, 1.5055162747589540e+05, 1.5056138019869605e+05, 1.5057113297963116e+05, 1.5058088581869737e+05, 1.5059063871589125e+05, 1.5060039167120951e+05, 1.5061014468464869e+05, 1.5061989775620546e+05, 1.5062965088587638e+05, 1.5063940407365817e+05, 1.5064915731954738e+05, 1.5065891062354061e+05, 1.5066866398563457e+05, 1.5067841740582586e+05, 1.5068817088411105e+05, 1.5069792442048682e+05, 1.5070767801494978e+05, 1.5071743166749657e+05, 1.5072718537812380e+05, 1.5073693914682808e+05, 1.5074669297360608e+05, 1.5075644685845441e+05, 1.5076620080136962e+05, 1.5077595480234845e+05, 1.5078570886138751e+05, 1.5079546297848341e+05, 1.5080521715363275e+05, 1.5081497138683219e+05, 1.5082472567807836e+05, 1.5083448002736786e+05, 1.5084423443469740e+05, 1.5085398890006356e+05, 1.5086374342346293e+05, 1.5087349800489223e+05, 1.5088325264434805e+05, 1.5089300734182697e+05, 1.5090276209732570e+05, 1.5091251691084084e+05, 1.5092227178236903e+05, 1.5093202671190695e+05, 1.5094178169945115e+05, 1.5095153674499830e+05, 1.5096129184854508e+05, 1.5097104701008808e+05, 1.5098080222962395e+05, 1.5099055750714932e+05, 1.5100031284266082e+05, 1.5101006823615509e+05, 1.5101982368762881e+05, 1.5102957919707859e+05, 1.5103933476450102e+05, 1.5104909038989287e+05, 1.5105884607325064e+05, 1.5106860181457101e+05, 1.5107835761385065e+05, 1.5108811347108620e+05, 1.5109786938627425e+05, 1.5110762535941153e+05, 1.5111738139049459e+05, 1.5112713747952011e+05, 1.5113689362648479e+05, 1.5114664983138515e+05, 1.5115640609421791e+05, 1.5116616241497974e+05, 1.5117591879366722e+05, 1.5118567523027703e+05, 1.5119543172480582e+05, 1.5120518827725024e+05, 1.5121494488760692e+05, 1.5122470155587248e+05, 1.5123445828204363e+05, 1.5124421506611697e+05, 1.5125397190808918e+05, 1.5126372880795688e+05, 1.5127348576571670e+05, 1.5128324278136532e+05, 1.5129299985489942e+05, 1.5130275698631559e+05, 1.5131251417561053e+05, 1.5132227142278085e+05, 1.5133202872782323e+05, 1.5134178609073427e+05, 1.5135154351151068e+05, 1.5136130099014906e+05, 1.5137105852664611e+05, 1.5138081612099844e+05, 1.5139057377320275e+05, 1.5140033148325569e+05, 1.5141008925115393e+05, 1.5141984707689405e+05, 1.5142960496047273e+05, 1.5143936290188663e+05, 1.5144912090113241e+05, 1.5145887895820674e+05, 1.5146863707310622e+05, 1.5147839524582762e+05, 1.5148815347636747e+05, 1.5149791176472252e+05, 1.5150767011088939e+05, 1.5151742851486476e+05, 1.5152718697664523e+05, 1.5153694549622750e+05, 1.5154670407360824e+05, 1.5155646270878409e+05, 1.5156622140175174e+05, 1.5157598015250781e+05, 1.5158573896104895e+05, 1.5159549782737190e+05, 1.5160525675147321e+05, 1.5161501573334966e+05, 1.5162477477299783e+05, 1.5163453387041442e+05, 1.5164429302559607e+05, 1.5165405223853944e+05, 1.5166381150924126e+05, 1.5167357083769806e+05, 1.5168333022390661e+05, 1.5169308966786353e+05, 1.5170284916956557e+05, 1.5171260872900928e+05, 1.5172236834619139e+05, 1.5173212802110860e+05, 1.5174188775375750e+05, 1.5175164754413476e+05, 1.5176140739223710e+05, 1.5177116729806116e+05, 1.5178092726160362e+05, 1.5179068728286112e+05, 1.5180044736183036e+05, 1.5181020749850801e+05, 1.5181996769289073e+05, 1.5182972794497517e+05, 1.5183948825475803e+05, 1.5184924862223599e+05, 1.5185900904740565e+05, 1.5186876953026379e+05, 1.5187853007080700e+05, 1.5188829066903202e+05, 1.5189805132493543e+05, 1.5190781203851398e+05, 1.5191757280976427e+05, 1.5192733363868308e+05, 1.5193709452526705e+05, 1.5194685546951278e+05, 1.5195661647141702e+05, 1.5196637753097643e+05, 1.5197613864818765e+05, 1.5198589982304742e+05, 1.5199566105555234e+05, 1.5200542234569919e+05, 1.5201518369348458e+05, 1.5202494509890518e+05, 1.5203470656195772e+05, 1.5204446808263881e+05, 1.5205422966094516e+05, 1.5206399129687346e+05, 1.5207375299042035e+05, 1.5208351474158259e+05, 1.5209327655035679e+05, 1.5210303841673967e+05, 1.5211280034072787e+05, 1.5212256232231812e+05, 1.5213232436150708e+05, 1.5214208645829139e+05, 1.5215184861266782e+05, 1.5216161082463295e+05, 1.5217137309418357e+05, 1.5218113542131628e+05, 1.5219089780602785e+05, 1.5220066024831490e+05, 1.5221042274817414e+05, 1.5222018530560224e+05, 1.5222994792059588e+05, 1.5223971059315177e+05, 1.5224947332326657e+05, 1.5225923611093702e+05, 1.5226899895615978e+05, 1.5227876185893148e+05, 1.5228852481924885e+05, 1.5229828783710863e+05, 1.5230805091250746e+05, 1.5231781404544204e+05, 1.5232757723590903e+05, 1.5233734048390514e+05, 1.5234710378942706e+05, 1.5235686715247150e+05, 1.5236663057303516e+05, 1.5237639405111465e+05, 1.5238615758670674e+05, 1.5239592117980818e+05, 1.5240568483041550e+05, 1.5241544853852550e+05, 1.5242521230413485e+05, 1.5243497612724028e+05, 1.5244474000783841e+05, 1.5245450394592600e+05, 1.5246426794149971e+05, 1.5247403199455625e+05, 1.5248379610509230e+05, 1.5249356027310461e+05, 1.5250332449858982e+05, 1.5251308878154465e+05, 1.5252285312196577e+05, 1.5253261751984994e+05, 1.5254238197519377e+05, 1.5255214648799400e+05, 1.5256191105824735e+05, 1.5257167568595053e+05, 1.5258144037110018e+05, 1.5259120511369305e+05, 1.5260096991372583e+05, 1.5261073477119519e+05, 1.5262049968609787e+05, 1.5263026465843056e+05, 1.5264002968818991e+05, 1.5264979477537272e+05, 1.5265955991997564e+05, 1.5266932512199538e+05, 1.5267909038142866e+05, 1.5268885569827215e+05, 1.5269862107252257e+05, 1.5270838650417663e+05, 1.5271815199323100e+05, 1.5272791753968241e+05, 1.5273768314352757e+05, 1.5274744880476317e+05, 1.5275721452338598e+05, 1.5276698029939266e+05, 1.5277674613277990e+05, 1.5278651202354443e+05, 1.5279627797168295e+05, 1.5280604397719217e+05, 1.5281581004006878e+05, 1.5282557616030955e+05, 1.5283534233791113e+05, 1.5284510857287023e+05, 1.5285487486518361e+05, 1.5286464121484794e+05, 1.5287440762185992e+05, 1.5288417408621628e+05, 1.5289394060791374e+05, 1.5290370718694903e+05, 1.5291347382331884e+05, 1.5292324051701985e+05, 1.5293300726804879e+05, 1.5294277407640242e+05, 1.5295254094207741e+05, 1.5296230786507050e+05, 1.5297207484537837e+05, 1.5298184188299775e+05, 1.5299160897792538e+05, 1.5300137613015797e+05, 1.5301114333969224e+05, 1.5302091060652485e+05, 1.5303067793065257e+05, 1.5304044531207209e+05, 1.5305021275078016e+05, 1.5305998024677348e+05, 1.5306974780004876e+05, 1.5307951541060273e+05, 1.5308928307843211e+05, 1.5309905080353364e+05, 1.5310881858590402e+05, 1.5311858642553995e+05, 1.5312835432243819e+05, 1.5313812227659542e+05, 1.5314789028800841e+05, 1.5315765835667387e+05, 1.5316742648258849e+05, 1.5317719466574900e+05, 1.5318696290615216e+05, 1.5319673120379462e+05, 1.5320649955867315e+05, 1.5321626797078451e+05, 1.5322603644012535e+05, 1.5323580496669246e+05, 1.5324557355048251e+05, 1.5325534219149229e+05, 1.5326511088971849e+05, 1.5327487964515784e+05, 1.5328464845780705e+05, 1.5329441732766284e+05, 1.5330418625472198e+05, 1.5331395523898120e+05, 1.5332372428043719e+05, 1.5333349337908669e+05, 1.5334326253492641e+05, 1.5335303174795315e+05, 1.5336280101816359e+05, 1.5337257034555441e+05, 1.5338233973012242e+05, 1.5339210917186434e+05, 1.5340187867077687e+05, 1.5341164822685675e+05, 1.5342141784010074e+05, 1.5343118751050552e+05, 1.5344095723806787e+05, 1.5345072702278453e+05, 1.5346049686465220e+05, 1.5347026676366763e+05, 1.5348003671982756e+05, 1.5348980673312867e+05, 1.5349957680356776e+05, 1.5350934693114154e+05, 1.5351911711584678e+05, 1.5352888735768019e+05, 1.5353865765663845e+05, 1.5354842801271839e+05, 1.5355819842591669e+05, 1.5356796889623010e+05, 1.5357773942365541e+05, 1.5358751000818930e+05, 1.5359728064982849e+05, 1.5360705134856977e+05, 1.5361682210440986e+05, 1.5362659291734549e+05, 1.5363636378737341e+05, 1.5364613471449038e+05, 1.5365590569869307e+05, 1.5366567673997831e+05, 1.5367544783834281e+05, 1.5368521899378329e+05, 1.5369499020629653e+05, 1.5370476147587920e+05, 1.5371453280252812e+05, 1.5372430418624004e+05, 1.5373407562701165e+05, 1.5374384712483973e+05, 1.5375361867972097e+05, 1.5376339029165218e+05, 1.5377316196063010e+05, 1.5378293368665143e+05, 1.5379270546971296e+05, 1.5380247730981145e+05, 1.5381224920694358e+05, 1.5382202116110615e+05, 1.5383179317229587e+05, 1.5384156524050955e+05, 1.5385133736574388e+05, 1.5386110954799561e+05, 1.5387088178726155e+05, 1.5388065408353839e+05, 1.5389042643682289e+05, 1.5390019884711181e+05, 1.5390997131440192e+05, 1.5391974383868993e+05, 1.5392951641997261e+05, 1.5393928905824674e+05, 1.5394906175350904e+05, 1.5395883450575624e+05, 1.5396860731498513e+05, 1.5397838018119248e+05, 1.5398815310437497e+05, 1.5399792608452943e+05, 1.5400769912165261e+05, 1.5401747221574123e+05, 1.5402724536679199e+05, 1.5403701857480177e+05, 1.5404679183976728e+05, 1.5405656516168523e+05, 1.5406633854055239e+05, 1.5407611197636553e+05, 1.5408588546912142e+05, 1.5409565901881683e+05, 1.5410543262544848e+05, 1.5411520628901318e+05, 1.5412498000950762e+05, 1.5413475378692866e+05, 1.5414452762127298e+05, 1.5415430151253738e+05, 1.5416407546071857e+05, 1.5417384946581334e+05, 1.5418362352781845e+05, 1.5419339764673065e+05, 1.5420317182254672e+05, 1.5421294605526337e+05, 1.5422272034487742e+05, 1.5423249469138565e+05, 1.5424226909478477e+05, 1.5425204355507161e+05, 1.5426181807224287e+05, 1.5427159264629529e+05, 1.5428136727722571e+05, 1.5429114196503084e+05, 1.5430091670970747e+05, 1.5431069151125237e+05, 1.5432046636966229e+05, 1.5433024128493405e+05, 1.5434001625706433e+05, 1.5434979128604996e+05, 1.5435956637188772e+05, 1.5436934151457431e+05, 1.5437911671410652e+05, 1.5438889197048117e+05, 1.5439866728369496e+05, 1.5440844265374474e+05, 1.5441821808062721e+05, 1.5442799356433915e+05, 1.5443776910487737e+05, 1.5444754470223861e+05, 1.5445732035641960e+05, 1.5446709606741721e+05, 1.5447687183522811e+05, 1.5448664765984917e+05, 1.5449642354127715e+05, 1.5450619947950874e+05, 1.5451597547454078e+05, 1.5452575152636998e+05, 1.5453552763499317e+05, 1.5454530380040716e+05, 1.5455508002260863e+05, 1.5456485630159444e+05, 1.5457463263736130e+05, 1.5458440902990600e+05, 1.5459418547922536e+05, 1.5460396198531613e+05, 1.5461373854817508e+05, 1.5462351516779902e+05, 1.5463329184418469e+05, 1.5464306857732890e+05, 1.5465284536722838e+05, 1.5466262221387995e+05, 1.5467239911728041e+05, 1.5468217607742644e+05, 1.5469195309431496e+05, 1.5470173016794259e+05, 1.5471150729830624e+05, 1.5472128448540266e+05, 1.5473106172922862e+05, 1.5474083902978091e+05, 1.5475061638705633e+05, 1.5476039380105160e+05, 1.5477017127176357e+05, 1.5477994879918898e+05, 1.5478972638332457e+05, 1.5479950402416725e+05, 1.5480928172171372e+05, 1.5481905947596079e+05, 1.5482883728690521e+05, 1.5483861515454386e+05, 1.5484839307887346e+05, 1.5485817105989077e+05, 1.5486794909759256e+05, 1.5487772719197566e+05, 1.5488750534303690e+05, 1.5489728355077299e+05, 1.5490706181518079e+05, 1.5491684013625703e+05, 1.5492661851399852e+05, 1.5493639694840205e+05, 1.5494617543946445e+05, 1.5495595398718238e+05, 1.5496573259155281e+05, 1.5497551125257241e+05, 1.5498528997023802e+05, 1.5499506874454638e+05, 1.5500484757549432e+05, 1.5501462646307860e+05, 1.5502440540729609e+05, 1.5503418440814348e+05, 1.5504396346561762e+05, 1.5505374257971533e+05, 1.5506352175043334e+05, 1.5507330097776849e+05, 1.5508308026171755e+05, 1.5509285960227734e+05, 1.5510263899944467e+05, 1.5511241845321629e+05, 1.5512219796358902e+05, 1.5513197753055961e+05, 1.5514175715412493e+05, 1.5515153683428172e+05, 1.5516131657102684e+05, 1.5517109636435704e+05, 1.5518087621426911e+05, 1.5519065612075987e+05, 1.5520043608382615e+05, 1.5521021610346469e+05, 1.5521999617967228e+05, 1.5522977631244581e+05, 1.5523955650178200e+05, 1.5524933674767768e+05, 1.5525911705012963e+05, 1.5526889740913472e+05, 1.5527867782468966e+05, 1.5528845829679130e+05, 1.5529823882543642e+05, 1.5530801941062184e+05, 1.5531780005234439e+05, 1.5532758075060084e+05, 1.5533736150538796e+05, 1.5534714231670261e+05, 1.5535692318454161e+05, 1.5536670410890173e+05, 1.5537648508977977e+05, 1.5538626612717254e+05, 1.5539604722107685e+05, 1.5540582837148954e+05, 1.5541560957840737e+05, 1.5542539084182712e+05, 1.5543517216174569e+05, 1.5544495353815984e+05, 1.5545473497106638e+05, 1.5546451646046209e+05, 1.5547429800634383e+05, 1.5548407960870839e+05, 1.5549386126755259e+05, 1.5550364298287316e+05, 1.5551342475466701e+05, 1.5552320658293093e+05, 1.5553298846766172e+05, 1.5554277040885619e+05, 1.5555255240651115e+05, 1.5556233446062339e+05, 1.5557211657118975e+05, 1.5558189873820703e+05, 1.5559168096167210e+05, 1.5560146324158169e+05, 1.5561124557793268e+05, 1.5562102797072183e+05, 1.5563081041994598e+05, 1.5564059292560193e+05, 1.5565037548768651e+05, 1.5566015810619658e+05, 1.5566994078112891e+05, 1.5567972351248030e+05, 1.5568950630024759e+05, 1.5569928914442763e+05, 1.5570907204501718e+05, 1.5571885500201309e+05, 1.5572863801541214e+05, 1.5573842108521124e+05, 1.5574820421140711e+05, 1.5575798739399656e+05, 1.5576777063297652e+05, 1.5577755392834372e+05, 1.5578733728009503e+05, 1.5579712068822721e+05, 1.5580690415273712e+05, 1.5581668767362161e+05, 1.5582647125087745e+05, 1.5583625488450154e+05, 1.5584603857449064e+05, 1.5585582232084154e+05, 1.5586560612355114e+05, 1.5587538998261627e+05, 1.5588517389803365e+05, 1.5589495786980022e+05, 1.5590474189791270e+05, 1.5591452598236801e+05, 1.5592431012316293e+05, 1.5593409432029427e+05, 1.5594387857375891e+05, 1.5595366288355360e+05, 1.5596344724967523e+05, 1.5597323167212057e+05, 1.5598301615088657e+05, 1.5599280068596991e+05, 1.5600258527736753e+05, 1.5601236992507617e+05, 1.5602215462909275e+05, 1.5603193938941401e+05, 1.5604172420603683e+05, 1.5605150907895801e+05, 1.5606129400817445e+05, 1.5607107899368290e+05, 1.5608086403548025e+05, 1.5609064913356327e+05, 1.5610043428792886e+05, 1.5611021949857383e+05, 1.5612000476549496e+05, 1.5612979008868913e+05, 1.5613957546815323e+05, 1.5614936090388399e+05, 1.5615914639587828e+05, 1.5616893194413296e+05, 1.5617871754864484e+05, 1.5618850320941076e+05, 1.5619828892642757e+05, 1.5620807469969208e+05, 1.5621786052920113e+05, 1.5622764641495157e+05, 1.5623743235694023e+05, 1.5624721835516396e+05, 1.5625700440961961e+05, 1.5626679052030400e+05, 1.5627657668721396e+05, 1.5628636291034633e+05, 1.5629614918969793e+05, 1.5630593552526564e+05, 1.5631572191704629e+05, 1.5632550836503672e+05, 1.5633529486923377e+05, 1.5634508142963427e+05, 1.5635486804623509e+05, 1.5636465471903305e+05, 1.5637444144802494e+05, 1.5638422823320772e+05, 1.5639401507457811e+05, 1.5640380197213305e+05, 1.5641358892586932e+05, 1.5642337593578378e+05, 1.5643316300187330e+05, 1.5644295012413469e+05, 1.5645273730256481e+05, 1.5646252453716053e+05, 1.5647231182791866e+05, 1.5648209917483607e+05, 1.5649188657790961e+05, 1.5650167403713608e+05, 1.5651146155251237e+05, 1.5652124912403530e+05, 1.5653103675170173e+05, 1.5654082443550855e+05, 1.5655061217545252e+05, 1.5656039997153057e+05, 1.5657018782373946e+05, 1.5657997573207616e+05, 1.5658976369653744e+05, 1.5659955171712016e+05, 1.5660933979382119e+05, 1.5661912792663736e+05, 1.5662891611556555e+05, 1.5663870436060260e+05, 1.5664849266174529e+05, 1.5665828101899056e+05, 1.5666806943233527e+05, 1.5667785790177621e+05, 1.5668764642731028e+05, 1.5669743500893432e+05, 1.5670722364664517e+05, 1.5671701234043972e+05, 1.5672680109031475e+05, 1.5673658989626722e+05, 1.5674637875829390e+05, 1.5675616767639166e+05, 1.5676595665055743e+05, 1.5677574568078801e+05, 1.5678553476708021e+05, 1.5679532390943097e+05, 1.5680511310783710e+05, 1.5681490236229546e+05, 1.5682469167280296e+05, 1.5683448103935638e+05, 1.5684427046195263e+05, 1.5685405994058854e+05, 1.5686384947526097e+05, 1.5687363906596683e+05, 1.5688342871270291e+05, 1.5689321841546614e+05, 1.5690300817425334e+05, 1.5691279798906136e+05, 1.5692258785988711e+05, 1.5693237778672742e+05, 1.5694216776957916e+05, 1.5695195780843918e+05, 1.5696174790330438e+05, 1.5697153805417157e+05, 1.5698132826103765e+05, 1.5699111852389947e+05, 1.5700090884275388e+05, 1.5701069921759781e+05, 1.5702048964842805e+05, 1.5703028013524151e+05, 1.5704007067803503e+05, 1.5704986127680549e+05, 1.5705965193154977e+05, 1.5706944264226471e+05, 1.5707923340894718e+05, 1.5708902423159408e+05, 1.5709881511020221e+05, 1.5710860604476850e+05, 1.5711839703528982e+05, 1.5712818808176302e+05, 1.5713797918418495e+05, 1.5714777034255251e+05, 1.5715756155686255e+05, 1.5716735282711199e+05, 1.5717714415329765e+05, 1.5718693553541641e+05, 1.5719672697346512e+05, 1.5720651846744068e+05, 1.5721631001733997e+05, 1.5722610162315986e+05, 1.5723589328489723e+05, 1.5724568500254894e+05, 1.5725547677611184e+05, 1.5726526860558282e+05, 1.5727506049095880e+05, 1.5728485243223660e+05, 1.5729464442941308e+05, 1.5730443648248521e+05, 1.5731422859144976e+05, 1.5732402075630362e+05, 1.5733381297704374e+05, 1.5734360525366693e+05, 1.5735339758617012e+05, 1.5736318997455013e+05, 1.5737298241880388e+05, 1.5738277491892819e+05, 1.5739256747492004e+05, 1.5740236008677626e+05, 1.5741215275449370e+05, 1.5742194547806925e+05, 1.5743173825749979e+05, 1.5744153109278224e+05, 1.5745132398391343e+05, 1.5746111693089028e+05, 1.5747090993370963e+05, 1.5748070299236837e+05, 1.5749049610686346e+05, 1.5750028927719168e+05, 1.5751008250334996e+05, 1.5751987578533517e+05, 1.5752966912314421e+05, 1.5753946251677396e+05, 1.5754925596622130e+05, 1.5755904947148310e+05, 1.5756884303255624e+05, 1.5757863664943763e+05, 1.5758843032212416e+05, 1.5759822405061268e+05, 1.5760801783490012e+05, 1.5761781167498333e+05, 1.5762760557085925e+05, 1.5763739952252471e+05, 1.5764719352997662e+05, 1.5765698759321184e+05, 1.5766678171222733e+05, 1.5767657588701992e+05, 1.5768637011758657e+05, 1.5769616440392405e+05, 1.5770595874602933e+05, 1.5771575314389929e+05, 1.5772554759753076e+05, 1.5773534210692075e+05, 1.5774513667206606e+05, 1.5775493129296359e+05, 1.5776472596961030e+05, 1.5777452070200298e+05, 1.5778431549013860e+05, 1.5779411033401403e+05, 1.5780390523362617e+05, 1.5781370018897188e+05, 1.5782349520004808e+05, 1.5783329026685169e+05, 1.5784308538937953e+05, 1.5785288056762854e+05, 1.5786267580159564e+05, 1.5787247109127769e+05, 1.5788226643667158e+05, 1.5789206183777424e+05, 1.5790185729458259e+05, 1.5791165280709346e+05, 1.5792144837530379e+05, 1.5793124399921045e+05, 1.5794103967881037e+05, 1.5795083541410041e+05, 1.5796063120507749e+05, 1.5797042705173852e+05, 1.5798022295408035e+05, 1.5799001891209994e+05, 1.5799981492579414e+05, 1.5800961099515992e+05, 1.5801940712019414e+05, 1.5802920330089369e+05, 1.5803899953725547e+05, 1.5804879582927638e+05, 1.5805859217695336e+05, 1.5806838858028330e+05, 1.5807818503926310e+05, 1.5808798155388964e+05, 1.5809777812415981e+05, 1.5810757475007055e+05, 1.5811737143161878e+05, 1.5812716816880135e+05, 1.5813696496161522e+05, 1.5814676181005727e+05, 1.5815655871412440e+05, 1.5816635567381349e+05, 1.5817615268912152e+05, 1.5818594976004539e+05, 1.5819574688658194e+05, 1.5820554406872811e+05, 1.5821534130648081e+05, 1.5822513859983691e+05, 1.5823493594879340e+05, 1.5824473335334714e+05, 1.5825453081349505e+05, 1.5826432832923398e+05, 1.5827412590056090e+05, 1.5828392352747271e+05, 1.5829372120996634e+05, 1.5830351894803872e+05, 1.5831331674168669e+05, 1.5832311459090721e+05, 1.5833291249569718e+05, 1.5834271045605349e+05, 1.5835250847197309e+05, 1.5836230654345284e+05, 1.5837210467048970e+05, 1.5838190285308060e+05, 1.5839170109122241e+05, 1.5840149938491211e+05, 1.5841129773414653e+05, 1.5842109613892261e+05, 1.5843089459923728e+05, 1.5844069311508743e+05, 1.5845049168647002e+05, 1.5846029031338196e+05, 1.5847008899582012e+05, 1.5847988773378145e+05, 1.5848968652726285e+05, 1.5849948537626129e+05, 1.5850928428077363e+05, 1.5851908324079681e+05, 1.5852888225632775e+05, 1.5853868132736339e+05, 1.5854848045390059e+05, 1.5855827963593631e+05, 1.5856807887346749e+05, 1.5857787816649099e+05, 1.5858767751500377e+05, 1.5859747691900277e+05, 1.5860727637848491e+05, 1.5861707589344704e+05, 1.5862687546388616e+05, 1.5863667508979916e+05, 1.5864647477118298e+05, 1.5865627450803452e+05, 1.5866607430035071e+05, 1.5867587414812847e+05, 1.5868567405136474e+05, 1.5869547401005647e+05, 1.5870527402420051e+05, 1.5871507409379387e+05, 1.5872487421883340e+05, 1.5873467439931611e+05, 1.5874447463523885e+05, 1.5875427492659853e+05, 1.5876407527339214e+05, 1.5877387567561661e+05, 1.5878367613326886e+05, 1.5879347664634575e+05, 1.5880327721484427e+05, 1.5881307783876135e+05, 1.5882287851809390e+05, 1.5883267925283886e+05, 1.5884248004299318e+05, 1.5885228088855371e+05, 1.5886208178951746e+05, 1.5887188274588130e+05, 1.5888168375764228e+05, 1.5889148482479720e+05, 1.5890128594734307e+05, 1.5891108712527677e+05, 1.5892088835859529e+05, 1.5893068964729551e+05, 1.5894049099137439e+05, 1.5895029239082886e+05, 1.5896009384565579e+05, 1.5896989535585223e+05, 1.5897969692141507e+05, 1.5898949854234120e+05, 1.5899930021862758e+05, 1.5900910195027117e+05, 1.5901890373726891e+05, 1.5902870557961770e+05, 1.5903850747731450e+05, 1.5904830943035620e+05, 1.5905811143873975e+05, 1.5906791350246215e+05, 1.5907771562152030e+05, 1.5908751779591115e+05, 1.5909732002563163e+05, 1.5910712231067865e+05, 1.5911692465104916e+05, 1.5912672704674015e+05, 1.5913652949774847e+05, 1.5914633200407116e+05, 1.5915613456570511e+05, 1.5916593718264729e+05, 1.5917573985489458e+05, 1.5918554258244400e+05, 1.5919534536529239e+05, 1.5920514820343678e+05, 1.5921495109687408e+05, 1.5922475404560124e+05, 1.5923455704961519e+05, 1.5924436010891286e+05, 1.5925416322349125e+05, 1.5926396639334728e+05, 1.5927376961847788e+05, 1.5928357289888000e+05, 1.5929337623455058e+05, 1.5930317962548658e+05, 1.5931298307168495e+05, 1.5932278657314260e+05, 1.5933259012985651e+05, 1.5934239374182362e+05, 1.5935219740904088e+05, 1.5936200113150521e+05, 1.5937180490921362e+05, 1.5938160874216302e+05, 1.5939141263035030e+05, 1.5940121657377252e+05, 1.5941102057242655e+05, 1.5942082462630939e+05, 1.5943062873541794e+05, 1.5944043289974920e+05, 1.5945023711930006e+05, 1.5946004139406752e+05, 1.5946984572404853e+05, 1.5947965010924006e+05, 1.5948945454963896e+05, 1.5949925904524227e+05, 1.5950906359604697e+05, 1.5951886820204998e+05, 1.5952867286324821e+05, 1.5953847757963865e+05, 1.5954828235121822e+05, 1.5955808717798392e+05, 1.5956789205993272e+05, 1.5957769699706152e+05, 1.5958750198936730e+05, 1.5959730703684702e+05, 1.5960711213949765e+05, 1.5961691729731610e+05, 1.5962672251029938e+05, 1.5963652777844443e+05, 1.5964633310174817e+05, 1.5965613848020759e+05, 1.5966594391381965e+05, 1.5967574940258131e+05, 1.5968555494648952e+05, 1.5969536054554125e+05, 1.5970516619973342e+05, 1.5971497190906302e+05, 1.5972477767352705e+05, 1.5973458349312239e+05, 1.5974438936784607e+05, 1.5975419529769503e+05, 1.5976400128266617e+05, 1.5977380732275650e+05, 1.5978361341796303e+05, 1.5979341956828267e+05, 1.5980322577371239e+05, 1.5981303203424916e+05, 1.5982283834988996e+05, 1.5983264472063171e+05, 1.5984245114647140e+05, 1.5985225762740601e+05, 1.5986206416343246e+05, 1.5987187075454777e+05, 1.5988167740074886e+05, 1.5989148410203270e+05, 1.5990129085839627e+05, 1.5991109766983654e+05, 1.5992090453635048e+05, 1.5993071145793507e+05, 1.5994051843458723e+05, 1.5995032546630394e+05, 1.5996013255308219e+05, 1.5996993969491895e+05, 1.5997974689181117e+05, 1.5998955414375581e+05, 1.5999936145074989e+05, 1.6000916881279030e+05, 1.6001897622987410e+05, 1.6002878370199821e+05, 1.6003859122915959e+05, 1.6004839881135523e+05, 1.6005820644858212e+05, 1.6006801414083719e+05, 1.6007782188811747e+05, 1.6008762969041988e+05, 1.6009743754774140e+05, 1.6010724546007902e+05, 1.6011705342742970e+05, 1.6012686144979042e+05, 1.6013666952715817e+05, 1.6014647765952989e+05, 1.6015628584690258e+05, 1.6016609408927319e+05, 1.6017590238663874e+05, 1.6018571073899616e+05, 1.6019551914634244e+05, 1.6020532760867459e+05, 1.6021513612598952e+05, 1.6022494469828426e+05, 1.6023475332555579e+05, 1.6024456200780108e+05, 1.6025437074501708e+05, 1.6026417953720075e+05, 1.6027398838434916e+05, 1.6028379728645922e+05, 1.6029360624352787e+05, 1.6030341525555222e+05, 1.6031322432252913e+05, 1.6032303344445562e+05, 1.6033284262132869e+05, 1.6034265185314530e+05, 1.6035246113990244e+05, 1.6036227048159711e+05, 1.6037207987822624e+05, 1.6038188932978685e+05, 1.6039169883627590e+05, 1.6040150839769040e+05, 1.6041131801402738e+05, 1.6042112768528372e+05, 1.6043093741145643e+05, 1.6044074719254253e+05, 1.6045055702853898e+05, 1.6046036691944278e+05, 1.6047017686525089e+05, 1.6047998686596032e+05, 1.6048979692156805e+05, 1.6049960703207107e+05, 1.6050941719746639e+05, 1.6051922741775098e+05, 1.6052903769292176e+05, 1.6053884802297581e+05, 1.6054865840791009e+05, 1.6055846884772155e+05, 1.6056827934240727e+05, 1.6057808989196413e+05, 1.6058790049638919e+05, 1.6059771115567943e+05, 1.6060752186983181e+05, 1.6061733263884336e+05, 1.6062714346271104e+05, 1.6063695434143182e+05, 1.6064676527500278e+05, 1.6065657626342084e+05, 1.6066638730668300e+05, 1.6067619840478624e+05, 1.6068600955772755e+05, 1.6069582076550394e+05, 1.6070563202811245e+05, 1.6071544334555001e+05, 1.6072525471781363e+05, 1.6073506614490031e+05, 1.6074487762680708e+05, 1.6075468916353086e+05, 1.6076450075506870e+05, 1.6077431240141756e+05, 1.6078412410257448e+05, 1.6079393585853640e+05, 1.6080374766930038e+05, 1.6081355953486337e+05, 1.6082337145522234e+05, 1.6083318343037437e+05, 1.6084299546031645e+05, 1.6085280754504551e+05, 1.6086261968455857e+05, 1.6087243187885266e+05, 1.6088224412792476e+05, 1.6089205643177187e+05, 1.6090186879039099e+05, 1.6091168120377912e+05, 1.6092149367193328e+05, 1.6093130619485045e+05, 1.6094111877252761e+05, 1.6095093140496183e+05, 1.6096074409215004e+05, 1.6097055683408925e+05, 1.6098036963077652e+05, 1.6099018248220879e+05, 1.6099999538838310e+05, 1.6100980834929648e+05, 1.6101962136494584e+05, 1.6102943443532821e+05, 1.6103924756044068e+05, 1.6104906074028020e+05, 1.6105887397484374e+05, 1.6106868726412838e+05, 1.6107850060813106e+05, 1.6108831400684881e+05, 1.6109812746027863e+05, 1.6110794096841753e+05, 1.6111775453126256e+05, 1.6112756814881068e+05, 1.6113738182105887e+05, 1.6114719554800421e+05, 1.6115700932964362e+05, 1.6116682316597423e+05, 1.6117663705699294e+05, 1.6118645100269682e+05, 1.6119626500308287e+05, 1.6120607905814805e+05, 1.6121589316788944e+05, 1.6122570733230401e+05, 1.6123552155138878e+05, 1.6124533582514076e+05, 1.6125515015355698e+05, 1.6126496453663445e+05, 1.6127477897437016e+05, 1.6128459346676114e+05, 1.6129440801380438e+05, 1.6130422261549687e+05, 1.6131403727183572e+05, 1.6132385198281787e+05, 1.6133366674844036e+05, 1.6134348156870017e+05, 1.6135329644359433e+05, 1.6136311137311987e+05, 1.6137292635727386e+05, 1.6138274139605317e+05, 1.6139255648945496e+05, 1.6140237163747617e+05, 1.6141218684011383e+05, 1.6142200209736498e+05, 1.6143181740922661e+05, 1.6144163277569573e+05, 1.6145144819676943e+05, 1.6146126367244468e+05, 1.6147107920271842e+05, 1.6148089478758781e+05, 1.6149071042704975e+05, 1.6150052612110131e+05, 1.6151034186973958e+05, 1.6152015767296148e+05, 1.6152997353076405e+05, 1.6153978944314434e+05, 1.6154960541009935e+05, 1.6155942143162608e+05, 1.6156923750772161e+05, 1.6157905363838293e+05, 1.6158886982360709e+05, 1.6159868606339104e+05, 1.6160850235773189e+05, 1.6161831870662660e+05, 1.6162813511007224e+05, 1.6163795156806579e+05, 1.6164776808060432e+05, 1.6165758464768482e+05, 1.6166740126930433e+05, 1.6167721794545988e+05, 1.6168703467614850e+05, 1.6169685146136722e+05, 1.6170666830111304e+05, 1.6171648519538299e+05, 1.6172630214417409e+05, 1.6173611914748338e+05, 1.6174593620530789e+05, 1.6175575331764467e+05, 1.6176557048449077e+05, 1.6177538770584317e+05, 1.6178520498169886e+05, 1.6179502231205491e+05, 1.6180483969690840e+05, 1.6181465713625631e+05, 1.6182447463009568e+05, 1.6183429217842358e+05, 1.6184410978123694e+05, 1.6185392743853285e+05, 1.6186374515030833e+05, 1.6187356291656048e+05, 1.6188338073728626e+05, 1.6189319861248272e+05, 1.6190301654214691e+05, 1.6191283452627584e+05, 1.6192265256486656e+05, 1.6193247065791607e+05, 1.6194228880542144e+05, 1.6195210700737967e+05, 1.6196192526378785e+05, 1.6197174357464301e+05, 1.6198156193994213e+05, 1.6199138035968228e+05, 1.6200119883386049e+05, 1.6201101736247379e+05, 1.6202083594551927e+05, 1.6203065458299386e+05, 1.6204047327489470e+05, 1.6205029202121880e+05, 1.6206011082196317e+05, 1.6206992967712483e+05, 1.6207974858670091e+05, 1.6208956755068837e+05, 1.6209938656908425e+05, 1.6210920564188567e+05, 1.6211902476908956e+05, 1.6212884395069306e+05, 1.6213866318669316e+05, 1.6214848247708686e+05, 1.6215830182187128e+05, 1.6216812122104340e+05, 1.6217794067460034e+05, 1.6218776018253906e+05, 1.6219757974485666e+05, 1.6220739936155017e+05, 1.6221721903261662e+05, 1.6222703875805304e+05, 1.6223685853785652e+05, 1.6224667837202406e+05, 1.6225649826055270e+05, 1.6226631820343956e+05, 1.6227613820068157e+05, 1.6228595825227589e+05, 1.6229577835821948e+05, 1.6230559851850942e+05, 1.6231541873314275e+05, 1.6232523900211652e+05, 1.6233505932542781e+05, 1.6234487970307365e+05, 1.6235470013505104e+05, 1.6236452062135708e+05, 1.6237434116198879e+05, 1.6238416175694321e+05, 1.6239398240621746e+05, 1.6240380310980853e+05, 1.6241362386771350e+05, 1.6242344467992935e+05, 1.6243326554645321e+05, 1.6244308646728206e+05, 1.6245290744241304e+05, 1.6246272847184312e+05, 1.6247254955556942e+05, 1.6248237069358892e+05, 1.6249219188589870e+05, 1.6250201313249586e+05, 1.6251183443337740e+05, 1.6252165578854040e+05, 1.6253147719798191e+05, 1.6254129866169897e+05, 1.6255112017968862e+05, 1.6256094175194798e+05, 1.6257076337847402e+05, 1.6258058505926386e+05, 1.6259040679431453e+05, 1.6260022858362307e+05, 1.6261005042718659e+05, 1.6261987232500207e+05, 1.6262969427706659e+05, 1.6263951628337728e+05, 1.6264933834393110e+05, 1.6265916045872515e+05, 1.6266898262775652e+05, 1.6267880485102220e+05, 1.6268862712851932e+05, 1.6269844946024488e+05, 1.6270827184619597e+05, 1.6271809428636968e+05, 1.6272791678076301e+05, 1.6273773932937306e+05, 1.6274756193219684e+05, 1.6275738458923146e+05, 1.6276720730047394e+05, 1.6277703006592143e+05, 1.6278685288557093e+05, 1.6279667575941945e+05, 1.6280649868746413e+05, 1.6281632166970201e+05, 1.6282614470613017e+05, 1.6283596779674562e+05, 1.6284579094154548e+05, 1.6285561414052677e+05, 1.6286543739368659e+05, 1.6287526070102199e+05, 1.6288508406253005e+05, 1.6289490747820781e+05, 1.6290473094805234e+05, 1.6291455447206076e+05, 1.6292437805023004e+05, 1.6293420168255732e+05, 1.6294402536903965e+05, 1.6295384910967408e+05, 1.6296367290445772e+05, 1.6297349675338759e+05, 1.6298332065646077e+05, 1.6299314461367435e+05, 1.6300296862502539e+05, 1.6301279269051095e+05, 1.6302261681012812e+05, 1.6303244098387394e+05, 1.6304226521174549e+05, 1.6305208949373980e+05, 1.6306191382985402e+05, 1.6307173822008522e+05, 1.6308156266443041e+05, 1.6309138716288671e+05, 1.6310121171545115e+05, 1.6311103632212084e+05, 1.6312086098289277e+05, 1.6313068569776413e+05, 1.6314051046673194e+05, 1.6315033528979329e+05, 1.6316016016694525e+05, 1.6316998509818487e+05, 1.6317981008350925e+05, 1.6318963512291547e+05, 1.6319946021640056e+05, 1.6320928536396165e+05, 1.6321911056559579e+05, 1.6322893582130005e+05, 1.6323876113107154e+05, 1.6324858649490730e+05, 1.6325841191280438e+05, 1.6326823738475997e+05, 1.6327806291077103e+05, 1.6328788849083468e+05, 1.6329771412494802e+05, 1.6330753981310807e+05, 1.6331736555531199e+05, 1.6332719135155680e+05, 1.6333701720183963e+05, 1.6334684310615752e+05, 1.6335666906450753e+05, 1.6336649507688681e+05, 1.6337632114329239e+05, 1.6338614726372133e+05, 1.6339597343817077e+05, 1.6340579966663770e+05, 1.6341562594911933e+05, 1.6342545228561270e+05, 1.6343527867611486e+05, 1.6344510512062290e+05, 1.6345493161913389e+05, 1.6346475817164496e+05, 1.6347458477815316e+05, 1.6348441143865557e+05, 1.6349423815314932e+05, 1.6350406492163145e+05, 1.6351389174409906e+05, 1.6352371862054925e+05, 1.6353354555097909e+05, 1.6354337253538563e+05, 1.6355319957376597e+05, 1.6356302666611725e+05, 1.6357285381243654e+05, 1.6358268101272089e+05, 1.6359250826696743e+05, 1.6360233557517320e+05, 1.6361216293733535e+05, 1.6362199035345094e+05, 1.6363181782351705e+05, 1.6364164534753078e+05, 1.6365147292548919e+05, 1.6366130055738942e+05, 1.6367112824322854e+05, 1.6368095598300363e+05, 1.6369078377671179e+05, 1.6370061162435013e+05, 1.6371043952591572e+05, 1.6372026748140564e+05, 1.6373009549081701e+05, 1.6373992355414690e+05, 1.6374975167139241e+05, 1.6375957984255065e+05, 1.6376940806761867e+05, 1.6377923634659359e+05, 1.6378906467947256e+05, 1.6379889306625258e+05, 1.6380872150693080e+05, 1.6381855000150425e+05, 1.6382837854997013e+05, 1.6383820715232549e+05, 1.6384803580856739e+05, 1.6385786451869298e+05, 1.6386769328269933e+05, 1.6387752210058353e+05, 1.6388735097234268e+05, 1.6389717989797393e+05, 1.6390700887747429e+05, 1.6391683791084090e+05, 1.6392666699807087e+05, 1.6393649613916129e+05, 1.6394632533410925e+05, 1.6395615458291181e+05, 1.6396598388556618e+05, 1.6397581324206936e+05, 1.6398564265241849e+05, 1.6399547211661070e+05, 1.6400530163464305e+05, 1.6401513120651260e+05, 1.6402496083221654e+05, 1.6403479051175193e+05, 1.6404462024511589e+05, 1.6405445003230550e+05, 1.6406427987331786e+05, 1.6407410976815008e+05, 1.6408393971679930e+05, 1.6409376971926258e+05, 1.6410359977553700e+05, 1.6411342988561976e+05, 1.6412326004950787e+05, 1.6413309026719845e+05, 1.6414292053868863e+05, 1.6415275086397555e+05, 1.6416258124305628e+05, 1.6417241167592790e+05, 1.6418224216258756e+05, 1.6419207270303232e+05, 1.6420190329725933e+05, 1.6421173394526567e+05, 1.6422156464704848e+05, 1.6423139540260486e+05, 1.6424122621193188e+05, 1.6425105707502671e+05, 1.6426088799188641e+05, 1.6427071896250808e+05, 1.6428054998688889e+05, 1.6429038106502593e+05, 1.6430021219691628e+05, 1.6431004338255705e+05, 1.6431987462194540e+05, 1.6432970591507838e+05, 1.6433953726195314e+05, 1.6434936866256676e+05, 1.6435920011691639e+05, 1.6436903162499913e+05, 1.6437886318681206e+05, 1.6438869480235240e+05, 1.6439852647161714e+05, 1.6440835819460344e+05, 1.6441818997130843e+05, 1.6442802180172919e+05, 1.6443785368586288e+05, 1.6444768562370661e+05, 1.6445751761525744e+05, 1.6446734966051253e+05, 1.6447718175946898e+05, 1.6448701391212389e+05, 1.6449684611847441e+05, 1.6450667837851765e+05, 1.6451651069225074e+05, 1.6452634305967079e+05, 1.6453617548077489e+05, 1.6454600795556017e+05, 1.6455584048402376e+05, 1.6456567306616277e+05, 1.6457550570197433e+05, 1.6458533839145556e+05, 1.6459517113460356e+05, 1.6460500393141544e+05, 1.6461483678188839e+05, 1.6462466968601948e+05, 1.6463450264380584e+05, 1.6464433565524459e+05, 1.6465416872033285e+05, 1.6466400183906770e+05, 1.6467383501144632e+05, 1.6468366823746581e+05, 1.6469350151712328e+05, 1.6470333485041591e+05, 1.6471316823734075e+05, 1.6472300167789499e+05, 1.6473283517207572e+05, 1.6474266871988005e+05, 1.6475250232130513e+05, 1.6476233597634808e+05, 1.6477216968500597e+05, 1.6478200344727602e+05, 1.6479183726315532e+05, 1.6480167113264097e+05, 1.6481150505573014e+05, 1.6482133903241990e+05, 1.6483117306270738e+05, 1.6484100714658981e+05, 1.6485084128406420e+05, 1.6486067547512773e+05, 1.6487050971977753e+05, 1.6488034401801068e+05, 1.6489017836982440e+05, 1.6490001277521573e+05, 1.6490984723418186e+05, 1.6491968174671993e+05, 1.6492951631282701e+05, 1.6493935093250024e+05, 1.6494918560573680e+05, 1.6495902033253375e+05, 1.6496885511288827e+05, 1.6497868994679747e+05, 1.6498852483425854e+05, 1.6499835977526853e+05, 1.6500819476982462e+05, 1.6501802981792390e+05, 1.6502786491956358e+05, 1.6503770007474071e+05, 1.6504753528345251e+05, 1.6505737054569603e+05, 1.6506720586146848e+05, 1.6507704123076695e+05, 1.6508687665358858e+05, 1.6509671212993050e+05, 1.6510654765978982e+05, 1.6511638324316376e+05, 1.6512621888004939e+05, 1.6513605457044384e+05, 1.6514589031434432e+05, 1.6515572611174788e+05, 1.6516556196265173e+05, 1.6517539786705290e+05, 1.6518523382494863e+05, 1.6519506983633604e+05, 1.6520490590121227e+05, 1.6521474201957442e+05, 1.6522457819141969e+05, 1.6523441441674516e+05, 1.6524425069554802e+05, 1.6525408702782539e+05, 1.6526392341357440e+05, 1.6527375985279217e+05, 1.6528359634547587e+05, 1.6529343289162268e+05, 1.6530326949122967e+05, 1.6531310614429403e+05, 1.6532294285081286e+05, 1.6533277961078336e+05, 1.6534261642420260e+05, 1.6535245329106780e+05, 1.6536229021137609e+05, 1.6537212718512455e+05, 1.6538196421231038e+05, 1.6539180129293073e+05, 1.6540163842698274e+05, 1.6541147561446353e+05, 1.6542131285537023e+05, 1.6543115014970006e+05, 1.6544098749745006e+05, 1.6545082489861746e+05, 1.6546066235319938e+05, 1.6547049986119298e+05, 1.6548033742259542e+05, 1.6549017503740379e+05, 1.6550001270561531e+05, 1.6550985042722712e+05, 1.6551968820223628e+05, 1.6552952603064003e+05, 1.6553936391243545e+05, 1.6554920184761976e+05, 1.6555903983619006e+05, 1.6556887787814354e+05, 1.6557871597347732e+05, 1.6558855412218857e+05, 1.6559839232427441e+05, 1.6560823057973202e+05, 1.6561806888855854e+05, 1.6562790725075113e+05, 1.6563774566630693e+05, 1.6564758413522309e+05, 1.6565742265749679e+05, 1.6566726123312517e+05, 1.6567709986210536e+05, 1.6568693854443455e+05, 1.6569677728010985e+05, 1.6570661606912845e+05, 1.6571645491148750e+05, 1.6572629380718415e+05, 1.6573613275621555e+05, 1.6574597175857885e+05, 1.6575581081427124e+05, 1.6576564992328983e+05, 1.6577548908563182e+05, 1.6578532830129433e+05, 1.6579516757027456e+05, 1.6580500689256963e+05, 1.6581484626817668e+05, 1.6582468569709288e+05, 1.6583452517931547e+05, 1.6584436471484153e+05, 1.6585420430366820e+05, 1.6586404394579268e+05, 1.6587388364121213e+05, 1.6588372338992369e+05, 1.6589356319192454e+05, 1.6590340304721182e+05, 1.6591324295578271e+05, 1.6592308291763437e+05, 1.6593292293276396e+05, 1.6594276300116864e+05, 1.6595260312284555e+05, 1.6596244329779188e+05, 1.6597228352600476e+05, 1.6598212380748137e+05, 1.6599196414221893e+05, 1.6600180453021452e+05, 1.6601164497146534e+05, 1.6602148546596852e+05, 1.6603132601372129e+05, 1.6604116661472077e+05, 1.6605100726896414e+05, 1.6606084797644857e+05, 1.6607068873717121e+05, 1.6608052955112918e+05, 1.6609037041831974e+05, 1.6610021133874002e+05, 1.6611005231238715e+05, 1.6611989333925839e+05, 1.6612973441935077e+05, 1.6613957555266155e+05, 1.6614941673918793e+05, 1.6615925797892699e+05, 1.6616909927187592e+05, 1.6617894061803192e+05, 1.6618878201739214e+05, 1.6619862346995375e+05, 1.6620846497571393e+05, 1.6621830653466986e+05, 1.6622814814681871e+05, 1.6623798981215761e+05, 1.6624783153068376e+05, 1.6625767330239431e+05, 1.6626751512728646e+05, 1.6627735700535736e+05, 1.6628719893660417e+05, 1.6629704092102411e+05, 1.6630688295861435e+05, 1.6631672504937201e+05, 1.6632656719329432e+05, 1.6633640939037842e+05, 1.6634625164062149e+05, 1.6635609394402069e+05, 1.6636593630057320e+05, 1.6637577871027624e+05, 1.6638562117312691e+05, 1.6639546368912244e+05, 1.6640530625825998e+05, 1.6641514888053676e+05, 1.6642499155594991e+05, 1.6643483428449658e+05, 1.6644467706617402e+05, 1.6645451990097930e+05, 1.6646436278890973e+05, 1.6647420572996241e+05, 1.6648404872413451e+05, 1.6649389177142325e+05, 1.6650373487182579e+05, 1.6651357802533926e+05, 1.6652342123196091e+05, 1.6653326449168791e+05, 1.6654310780451741e+05, 1.6655295117044658e+05, 1.6656279458947264e+05, 1.6657263806159279e+05, 1.6658248158680415e+05, 1.6659232516510395e+05, 1.6660216879648934e+05, 1.6661201248095746e+05, 1.6662185621850559e+05, 1.6663170000913084e+05, 1.6664154385283045e+05, 1.6665138774960156e+05, 1.6666123169944138e+05, 1.6667107570234712e+05, 1.6668091975831587e+05, 1.6669076386734488e+05, 1.6670060802943134e+05, 1.6671045224457243e+05, 1.6672029651276532e+05, 1.6673014083400718e+05, 1.6673998520829520e+05, 1.6674982963562658e+05, 1.6675967411599855e+05, 1.6676951864940824e+05, 1.6677936323585286e+05, 1.6678920787532956e+05, 1.6679905256783558e+05, 1.6680889731336810e+05, 1.6681874211192428e+05, 1.6682858696350132e+05, 1.6683843186809641e+05, 1.6684827682570677e+05, 1.6685812183632952e+05, 1.6686796689996193e+05, 1.6687781201660115e+05, 1.6688765718624435e+05, 1.6689750240888877e+05, 1.6690734768453156e+05, 1.6691719301316995e+05, 1.6692703839480106e+05, 1.6693688382942218e+05, 1.6694672931703040e+05, 1.6695657485762297e+05, 1.6696642045119710e+05, 1.6697626609774993e+05, 1.6698611179727866e+05, 1.6699595754978055e+05, 1.6700580335525275e+05, 1.6701564921369246e+05, 1.6702549512509687e+05, 1.6703534108946315e+05, 1.6704518710678854e+05, 1.6705503317707020e+05, 1.6706487930030533e+05, 1.6707472547649118e+05, 1.6708457170562487e+05, 1.6709441798770364e+05, 1.6710426432272463e+05, 1.6711411071068514e+05, 1.6712395715158232e+05, 1.6713380364541331e+05, 1.6714365019217538e+05, 1.6715349679186571e+05, 1.6716334344448149e+05, 1.6717319015001995e+05, 1.6718303690847824e+05, 1.6719288371985359e+05, 1.6720273058414317e+05, 1.6721257750134423e+05, 1.6722242447145394e+05, 1.6723227149446949e+05, 1.6724211857038812e+05, 1.6725196569920695e+05, 1.6726181288092324e+05, 1.6727166011553424e+05, 1.6728150740303707e+05, 1.6729135474342899e+05, 1.6730120213670717e+05, 1.6731104958286884e+05, 1.6732089708191116e+05, 1.6733074463383137e+05, 1.6734059223862668e+05, 1.6735043989629421e+05, 1.6736028760683126e+05, 1.6737013537023505e+05, 1.6737998318650268e+05, 1.6738983105563148e+05, 1.6739967897761858e+05, 1.6740952695246116e+05, 1.6741937498015648e+05, 1.6742922306070177e+05, 1.6743907119409417e+05, 1.6744891938033092e+05, 1.6745876761940922e+05, 1.6746861591132628e+05, 1.6747846425607931e+05, 1.6748831265366552e+05, 1.6749816110408213e+05, 1.6750800960732636e+05, 1.6751785816339534e+05, 1.6752770677228639e+05, 1.6753755543399663e+05, 1.6754740414852332e+05, 1.6755725291586365e+05, 1.6756710173601483e+05, 1.6757695060897409e+05, 1.6758679953473862e+05, 1.6759664851330564e+05, 1.6760649754467237e+05, 1.6761634662883604e+05, 1.6762619576579382e+05, 1.6763604495554292e+05, 1.6764589419808061e+05, 1.6765574349340406e+05, 1.6766559284151052e+05, 1.6767544224239714e+05, 1.6768529169606118e+05, 1.6769514120249986e+05, 1.6770499076171039e+05, 1.6771484037368995e+05, 1.6772469003843580e+05, 1.6773453975594515e+05, 1.6774438952621521e+05, 1.6775423934924314e+05, 1.6776408922502626e+05, 1.6777393915356169e+05, 1.6778378913484671e+05, 1.6779363916887849e+05, 1.6780348925565431e+05, 1.6781333939517135e+05, 1.6782318958742684e+05, 1.6783303983241797e+05, 1.6784289013014198e+05, 1.6785274048059611e+05, 1.6786259088377754e+05, 1.6787244133968352e+05, 1.6788229184831126e+05, 1.6789214240965800e+05, 1.6790199302372089e+05, 1.6791184369049728e+05, 1.6792169440998425e+05, 1.6793154518217914e+05, 1.6794139600707908e+05, 1.6795124688468134e+05, 1.6796109781498316e+05, 1.6797094879798172e+05, 1.6798079983367425e+05, 1.6799065092205800e+05, 1.6800050206313017e+05, 1.6801035325688799e+05, 1.6802020450332868e+05, 1.6803005580244950e+05, 1.6803990715424763e+05, 1.6804975855872029e+05, 1.6805961001586475e+05, 1.6806946152567820e+05, 1.6807931308815785e+05, 1.6808916470330098e+05, 1.6809901637110481e+05, 1.6810886809156652e+05, 1.6811871986468337e+05, 1.6812857169045258e+05, 1.6813842356887137e+05, 1.6814827549993698e+05, 1.6815812748364665e+05, 1.6816797951999761e+05, 1.6817783160898703e+05, 1.6818768375061222e+05, 1.6819753594487035e+05, 1.6820738819175871e+05, 1.6821724049127448e+05, 1.6822709284341490e+05, 1.6823694524817722e+05, 1.6824679770555865e+05, 1.6825665021555641e+05, 1.6826650277816775e+05, 1.6827635539338988e+05, 1.6828620806122010e+05, 1.6829606078165557e+05, 1.6830591355469354e+05, 1.6831576638033130e+05, 1.6832561925856603e+05, 1.6833547218939493e+05, 1.6834532517281530e+05, 1.6835517820882433e+05, 1.6836503129741928e+05, 1.6837488443859736e+05, 1.6838473763235586e+05, 1.6839459087869193e+05, 1.6840444417760288e+05, 1.6841429752908592e+05, 1.6842415093313830e+05, 1.6843400438975720e+05, 1.6844385789893995e+05, 1.6845371146068373e+05, 1.6846356507498576e+05, 1.6847341874184331e+05, 1.6848327246125363e+05, 1.6849312623321393e+05, 1.6850298005772143e+05, 1.6851283393477343e+05, 1.6852268786436712e+05, 1.6853254184649975e+05, 1.6854239588116857e+05, 1.6855224996837080e+05, 1.6856210410810373e+05, 1.6857195830036452e+05, 1.6858181254515049e+05, 1.6859166684245883e+05, 1.6860152119228680e+05, 1.6861137559463162e+05, 1.6862123004949064e+05, 1.6863108455686094e+05, 1.6864093911673987e+05, 1.6865079372912462e+05, 1.6866064839401245e+05, 1.6867050311140061e+05, 1.6868035788128636e+05, 1.6869021270366691e+05, 1.6870006757853954e+05, 1.6870992250590146e+05, 1.6871977748574995e+05, 1.6872963251808222e+05, 1.6873948760289553e+05, 1.6874934274018713e+05, 1.6875919792995424e+05, 1.6876905317219420e+05, 1.6877890846690413e+05, 1.6878876381408135e+05, 1.6879861921372308e+05, 1.6880847466582659e+05, 1.6881833017038909e+05, 1.6882818572740789e+05, 1.6883804133688018e+05, 1.6884789699880325e+05, 1.6885775271317433e+05, 1.6886760847999062e+05, 1.6887746429924943e+05, 1.6888732017094805e+05, 1.6889717609508368e+05, 1.6890703207165352e+05, 1.6891688810065491e+05, 1.6892674418208504e+05, 1.6893660031594118e+05, 1.6894645650222059e+05, 1.6895631274092052e+05, 1.6896616903203825e+05, 1.6897602537557099e+05, 1.6898588177151600e+05, 1.6899573821987052e+05, 1.6900559472063181e+05, 1.6901545127379714e+05, 1.6902530787936377e+05, 1.6903516453732896e+05, 1.6904502124768993e+05, 1.6905487801044396e+05, 1.6906473482558833e+05, 1.6907459169312022e+05, 1.6908444861303695e+05, 1.6909430558533577e+05, 1.6910416261001388e+05, 1.6911401968706862e+05, 1.6912387681649716e+05, 1.6913373399829681e+05, 1.6914359123246482e+05, 1.6915344851899846e+05, 1.6916330585789500e+05, 1.6917316324915164e+05, 1.6918302069276568e+05, 1.6919287818873438e+05, 1.6920273573705499e+05, 1.6921259333772474e+05, 1.6922245099074094e+05, 1.6923230869610081e+05, 1.6924216645380162e+05, 1.6925202426384063e+05, 1.6926188212621515e+05, 1.6927174004092239e+05, 1.6928159800795963e+05, 1.6929145602732411e+05, 1.6930131409901314e+05, 1.6931117222302392e+05, 1.6932103039935377e+05, 1.6933088862799990e+05, 1.6934074690895958e+05, 1.6935060524223014e+05, 1.6936046362780876e+05, 1.6937032206569274e+05, 1.6938018055587934e+05, 1.6939003909836584e+05, 1.6939989769314951e+05, 1.6940975634022758e+05, 1.6941961503959735e+05, 1.6942947379125605e+05, 1.6943933259520095e+05, 1.6944919145142933e+05, 1.6945905035993850e+05, 1.6946890932072565e+05, 1.6947876833378809e+05, 1.6948862739912307e+05, 1.6949848651672786e+05, 1.6950834568659976e+05, 1.6951820490873599e+05, 1.6952806418313386e+05, 1.6953792350979062e+05, 1.6954778288870354e+05, 1.6955764231986989e+05, 1.6956750180328693e+05, 1.6957736133895192e+05, 1.6958722092686218e+05, 1.6959708056701493e+05, 1.6960694025940748e+05, 1.6961680000403704e+05, 1.6962665980090096e+05, 1.6963651964999651e+05, 1.6964637955132089e+05, 1.6965623950487142e+05, 1.6966609951064532e+05, 1.6967595956863993e+05, 1.6968581967885248e+05, 1.6969567984128025e+05, 1.6970554005592054e+05, 1.6971540032277061e+05, 1.6972526064182777e+05, 1.6973512101308923e+05, 1.6974498143655228e+05, 1.6975484191221424e+05, 1.6976470244007232e+05, 1.6977456302012387e+05, 1.6978442365236609e+05, 1.6979428433679629e+05, 1.6980414507341178e+05, 1.6981400586220980e+05, 1.6982386670318758e+05, 1.6983372759634248e+05, 1.6984358854167175e+05, 1.6985344953917270e+05, 1.6986331058884252e+05, 1.6987317169067854e+05, 1.6988303284467809e+05, 1.6989289405083834e+05, 1.6990275530915667e+05, 1.6991261661963034e+05, 1.6992247798225659e+05, 1.6993233939703272e+05, 1.6994220086395598e+05, 1.6995206238302367e+05, 1.6996192395423312e+05, 1.6997178557758156e+05, 1.6998164725306627e+05, 1.6999150898068456e+05, 1.7000137076043367e+05, 1.7001123259231093e+05, 1.7002109447631365e+05, 1.7003095641243903e+05, 1.7004081840068442e+05, 1.7005068044104706e+05, 1.7006054253352422e+05, 1.7007040467811321e+05, 1.7008026687481129e+05, 1.7009012912361580e+05, 1.7009999142452399e+05, 1.7010985377753316e+05, 1.7011971618264061e+05, 1.7012957863984359e+05, 1.7013944114913940e+05, 1.7014930371052533e+05, 1.7015916632399862e+05, 1.7016902898955665e+05, 1.7017889170719666e+05, 1.7018875447691590e+05, 1.7019861729871170e+05, 1.7020848017258136e+05, 1.7021834309852216e+05, 1.7022820607653135e+05, 1.7023806910660624e+05, 1.7024793218874413e+05, 1.7025779532294231e+05, 1.7026765850919808e+05, 1.7027752174750870e+05, 1.7028738503787149e+05, 1.7029724838028377e+05, 1.7030711177474272e+05, 1.7031697522124572e+05, 1.7032683871979002e+05, 1.7033670227037297e+05, 1.7034656587299181e+05, 1.7035642952764384e+05, 1.7036629323432635e+05, 1.7037615699303665e+05, 1.7038602080377203e+05, 1.7039588466652978e+05, 1.7040574858130721e+05, 1.7041561254810158e+05, 1.7042547656691022e+05, 1.7043534063773041e+05, 1.7044520476055943e+05, 1.7045506893539458e+05, 1.7046493316223315e+05, 1.7047479744107247e+05, 1.7048466177190983e+05, 1.7049452615474249e+05, 1.7050439058956777e+05, 1.7051425507638295e+05, 1.7052411961518537e+05, 1.7053398420597226e+05, 1.7054384884874098e+05, 1.7055371354348879e+05, 1.7056357829021302e+05, 1.7057344308891095e+05, 1.7058330793957986e+05, 1.7059317284221706e+05, 1.7060303779681993e+05, 1.7061290280338563e+05, 1.7062276786191156e+05, 1.7063263297239496e+05, 1.7064249813483315e+05, 1.7065236334922348e+05, 1.7066222861556316e+05, 1.7067209393384962e+05, 1.7068195930408002e+05, 1.7069182472625177e+05, 1.7070169020036212e+05, 1.7071155572640838e+05, 1.7072142130438783e+05, 1.7073128693429779e+05, 1.7074115261613563e+05, 1.7075101834989854e+05, 1.7076088413558391e+05, 1.7077074997318897e+05, 1.7078061586271107e+05, 1.7079048180414757e+05, 1.7080034779749569e+05, 1.7081021384275274e+05, 1.7082007993991603e+05, 1.7082994608898292e+05, 1.7083981228995067e+05, 1.7084967854281660e+05, 1.7085954484757801e+05, 1.7086941120423222e+05, 1.7087927761277652e+05, 1.7088914407320824e+05, 1.7089901058552464e+05, 1.7090887714972306e+05, 1.7091874376580078e+05, 1.7092861043375518e+05, 1.7093847715358352e+05, 1.7094834392528312e+05, 1.7095821074885127e+05, 1.7096807762428530e+05, 1.7097794455158254e+05, 1.7098781153074029e+05, 1.7099767856175581e+05, 1.7100754564462643e+05, 1.7101741277934951e+05, 1.7102727996592232e+05, 1.7103714720434215e+05, 1.7104701449460638e+05, 1.7105688183671224e+05, 1.7106674923065712e+05, 1.7107661667643831e+05, 1.7108648417405313e+05, 1.7109635172349884e+05, 1.7110621932477283e+05, 1.7111608697787236e+05, 1.7112595468279475e+05, 1.7113582243953733e+05, 1.7114569024809741e+05, 1.7115555810847232e+05, 1.7116542602065936e+05, 1.7117529398465581e+05, 1.7118516200045904e+05, 1.7119503006806632e+05, 1.7120489818747502e+05, 1.7121476635868242e+05, 1.7122463458168585e+05, 1.7123450285648264e+05, 1.7124437118307009e+05, 1.7125423956144552e+05, 1.7126410799160626e+05, 1.7127397647354956e+05, 1.7128384500727284e+05, 1.7129371359277336e+05, 1.7130358223004846e+05, 1.7131345091909546e+05, 1.7132331965991168e+05, 1.7133318845249442e+05, 1.7134305729684103e+05, 1.7135292619294880e+05, 1.7136279514081505e+05, 1.7137266414043715e+05, 1.7138253319181237e+05, 1.7139240229493802e+05, 1.7140227144981147e+05, 1.7141214065643001e+05, 1.7142200991479101e+05, 1.7143187922489175e+05, 1.7144174858672958e+05, 1.7145161800030179e+05, 1.7146148746560572e+05, 1.7147135698263868e+05, 1.7148122655139802e+05, 1.7149109617188104e+05, 1.7150096584408509e+05, 1.7151083556800749e+05, 1.7152070534364556e+05, 1.7153057517099660e+05, 1.7154044505005793e+05, 1.7155031498082692e+05, 1.7156018496330091e+05, 1.7157005499747719e+05, 1.7157992508335310e+05, 1.7158979522092594e+05, 1.7159966541019306e+05, 1.7160953565115182e+05, 1.7161940594379950e+05, 1.7162927628813343e+05, 1.7163914668415097e+05, 1.7164901713184942e+05, 1.7165888763122616e+05, 1.7166875818227846e+05, 1.7167862878500365e+05, 1.7168849943939908e+05, 1.7169837014546210e+05, 1.7170824090318999e+05, 1.7171811171258014e+05, 1.7172798257362985e+05, 1.7173785348633648e+05, 1.7174772445069731e+05, 1.7175759546670967e+05, 1.7176746653437096e+05, 1.7177733765367846e+05, 1.7178720882462955e+05, 1.7179708004722148e+05, 1.7180695132145166e+05, 1.7181682264731740e+05, 1.7182669402481601e+05, 1.7183656545394484e+05, 1.7184643693470125e+05, 1.7185630846708256e+05, 1.7186618005108609e+05, 1.7187605168670916e+05, 1.7188592337394916e+05, 1.7189579511280337e+05, 1.7190566690326916e+05, 1.7191553874534383e+05, 1.7192541063902478e+05, 1.7193528258430929e+05, 1.7194515458119474e+05, 1.7195502662967841e+05, 1.7196489872975767e+05, 1.7197477088142990e+05, 1.7198464308469242e+05, 1.7199451533954250e+05, 1.7200438764597758e+05, 1.7201426000399492e+05, 1.7202413241359187e+05, 1.7203400487476579e+05, 1.7204387738751402e+05, 1.7205374995183389e+05, 1.7206362256772275e+05, 1.7207349523517792e+05, 1.7208336795419682e+05, 1.7209324072477670e+05, 1.7210311354691494e+05, 1.7211298642060885e+05, 1.7212285934585580e+05, 1.7213273232265312e+05, 1.7214260535099817e+05, 1.7215247843088830e+05, 1.7216235156232081e+05, 1.7217222474529306e+05, 1.7218209797980246e+05, 1.7219197126584628e+05, 1.7220184460342187e+05, 1.7221171799252660e+05, 1.7222159143315780e+05, 1.7223146492531282e+05, 1.7224133846898898e+05, 1.7225121206418367e+05, 1.7226108571089420e+05, 1.7227095940911790e+05, 1.7228083315885218e+05, 1.7229070696009434e+05, 1.7230058081284177e+05, 1.7231045471709178e+05, 1.7232032867284172e+05, 1.7233020268008890e+05, 1.7234007673883077e+05, 1.7234995084906460e+05, 1.7235982501078775e+05, 1.7236969922399760e+05, 1.7237957348869144e+05, 1.7238944780486671e+05, 1.7239932217252065e+05, 1.7240919659165072e+05, 1.7241907106225420e+05, 1.7242894558432844e+05, 1.7243882015787082e+05, 1.7244869478287868e+05, 1.7245856945934938e+05, 1.7246844418728026e+05, 1.7247831896666862e+05, 1.7248819379751192e+05, 1.7249806867980745e+05, 1.7250794361355255e+05, 1.7251781859874463e+05, 1.7252769363538103e+05, 1.7253756872345903e+05, 1.7254744386297607e+05, 1.7255731905392947e+05, 1.7256719429631656e+05, 1.7257706959013475e+05, 1.7258694493538133e+05, 1.7259682033205370e+05, 1.7260669578014920e+05, 1.7261657127966516e+05, 1.7262644683059899e+05, 1.7263632243294807e+05, 1.7264619808670966e+05, 1.7265607379188118e+05, 1.7266594954845999e+05, 1.7267582535644341e+05, 1.7268570121582886e+05, 1.7269557712661361e+05, 1.7270545308879507e+05, 1.7271532910237060e+05, 1.7272520516733755e+05, 1.7273508128369329e+05, 1.7274495745143518e+05, 1.7275483367056059e+05, 1.7276470994106683e+05, 1.7277458626295128e+05, 1.7278446263621136e+05, 1.7279433906084430e+05, 1.7280421553684759e+05, 1.7281409206421857e+05, 1.7282396864295451e+05, 1.7283384527305290e+05, 1.7284372195451104e+05, 1.7285359868732627e+05, 1.7286347547149597e+05, 1.7287335230701751e+05, 1.7288322919388826e+05, 1.7289310613210555e+05, 1.7290298312166677e+05, 1.7291286016256930e+05, 1.7292273725481049e+05, 1.7293261439838767e+05, 1.7294249159329827e+05, 1.7295236883953959e+05, 1.7296224613710903e+05, 1.7297212348600398e+05, 1.7298200088622177e+05, 1.7299187833775976e+05, 1.7300175584061531e+05, 1.7301163339478584e+05, 1.7302151100026866e+05, 1.7303138865706118e+05, 1.7304126636516070e+05, 1.7305114412456466e+05, 1.7306102193527040e+05, 1.7307089979727528e+05, 1.7308077771057672e+05, 1.7309065567517205e+05, 1.7310053369105858e+05, 1.7311041175823379e+05, 1.7312028987669499e+05, 1.7313016804643956e+05, 1.7314004626746484e+05, 1.7314992453976825e+05, 1.7315980286334714e+05, 1.7316968123819886e+05, 1.7317955966432078e+05, 1.7318943814171033e+05, 1.7319931667036482e+05, 1.7320919525028163e+05, 1.7321907388145817e+05, 1.7322895256389177e+05, 1.7323883129757983e+05, 1.7324871008251977e+05, 1.7325858891870885e+05, 1.7326846780614450e+05, 1.7327834674482411e+05, 1.7328822573474506e+05, 1.7329810477590471e+05, 1.7330798386830045e+05, 1.7331786301192958e+05, 1.7332774220678955e+05, 1.7333762145287771e+05, 1.7334750075019142e+05, 1.7335738009872811e+05, 1.7336725949848513e+05, 1.7337713894945983e+05, 1.7338701845164961e+05, 1.7339689800505189e+05, 1.7340677760966393e+05, 1.7341665726548320e+05, 1.7342653697250708e+05, 1.7343641673073289e+05, 1.7344629654015810e+05, 1.7345617640078001e+05, 1.7346605631259602e+05, 1.7347593627560354e+05, 1.7348581628979987e+05, 1.7349569635518247e+05, 1.7350557647174867e+05, 1.7351545663949591e+05, 1.7352533685842151e+05, 1.7353521712852287e+05, 1.7354509744979738e+05, 1.7355497782224239e+05, 1.7356485824585531e+05, 1.7357473872063356e+05, 1.7358461924657444e+05, 1.7359449982367540e+05, 1.7360438045193380e+05, 1.7361426113134704e+05, 1.7362414186191245e+05, 1.7363402264362745e+05, 1.7364390347648942e+05, 1.7365378436049575e+05, 1.7366366529564382e+05, 1.7367354628193102e+05, 1.7368342731935473e+05, 1.7369330840791229e+05, 1.7370318954760119e+05, 1.7371307073841873e+05, 1.7372295198036230e+05, 1.7373283327342934e+05, 1.7374271461761720e+05, 1.7375259601292326e+05, 1.7376247745934490e+05, 1.7377235895687953e+05, 1.7378224050552456e+05, 1.7379212210527738e+05, 1.7380200375613529e+05, 1.7381188545809578e+05, 1.7382176721115614e+05, 1.7383164901531386e+05, 1.7384153087056629e+05, 1.7385141277691079e+05, 1.7386129473434479e+05, 1.7387117674286562e+05, 1.7388105880247074e+05, 1.7389094091315754e+05, 1.7390082307492333e+05, 1.7391070528776560e+05, 1.7392058755168167e+05, 1.7393046986666898e+05, 1.7394035223272489e+05, 1.7395023464984683e+05, 1.7396011711803215e+05, 1.7396999963727823e+05, 1.7397988220758253e+05, 1.7398976482894240e+05, 1.7399964750135521e+05, 1.7400953022481839e+05, 1.7401941299932933e+05, 1.7402929582488543e+05, 1.7403917870148405e+05, 1.7404906162912262e+05, 1.7405894460779853e+05, 1.7406882763750918e+05, 1.7407871071825194e+05, 1.7408859385002419e+05, 1.7409847703282337e+05, 1.7410836026664689e+05, 1.7411824355149211e+05, 1.7412812688735645e+05, 1.7413801027423731e+05, 1.7414789371213203e+05, 1.7415777720103809e+05, 1.7416766074095279e+05, 1.7417754433187359e+05, 1.7418742797379789e+05, 1.7419731166672311e+05, 1.7420719541064661e+05, 1.7421707920556582e+05, 1.7422696305147812e+05, 1.7423684694838087e+05, 1.7424673089627153e+05, 1.7425661489514750e+05, 1.7426649894500614e+05, 1.7427638304584488e+05, 1.7428626719766110e+05, 1.7429615140045225e+05, 1.7430603565421567e+05, 1.7431591995894877e+05, 1.7432580431464899e+05, 1.7433568872131369e+05, 1.7434557317894034e+05, 1.7435545768752627e+05, 1.7436534224706894e+05, 1.7437522685756569e+05, 1.7438511151901400e+05, 1.7439499623141121e+05, 1.7440488099475476e+05, 1.7441476580904203e+05, 1.7442465067427044e+05, 1.7443453559043742e+05, 1.7444442055754032e+05, 1.7445430557557658e+05, 1.7446419064454359e+05, 1.7447407576443878e+05, 1.7448396093525953e+05, 1.7449384615700328e+05, 1.7450373142966739e+05, 1.7451361675324931e+05, 1.7452350212774644e+05, 1.7453338755315615e+05, 1.7454327302947585e+05, 1.7455315855670298e+05, 1.7456304413483499e+05, 1.7457292976386921e+05, 1.7458281544380309e+05, 1.7459270117463407e+05, 1.7460258695635948e+05, 1.7461247278897677e+05, 1.7462235867248333e+05, 1.7463224460687660e+05, 1.7464213059215399e+05, 1.7465201662831294e+05, 1.7466190271535076e+05, 1.7467178885326497e+05, 1.7468167504205293e+05, 1.7469156128171203e+05, 1.7470144757223973e+05, 1.7471133391363345e+05, 1.7472122030589057e+05, 1.7473110674900853e+05, 1.7474099324298467e+05, 1.7475087978781646e+05, 1.7476076638350132e+05, 1.7477065303003666e+05, 1.7478053972741988e+05, 1.7479042647564842e+05, 1.7480031327471967e+05, 1.7481020012463102e+05, 1.7482008702537997e+05, 1.7482997397696384e+05, 1.7483986097938014e+05, 1.7484974803262620e+05, 1.7485963513669948e+05, 1.7486952229159742e+05, 1.7487940949731736e+05, 1.7488929675385679e+05, 1.7489918406121308e+05, 1.7490907141938369e+05, 1.7491895882836598e+05, 1.7492884628815745e+05, 1.7493873379875545e+05, 1.7494862136015744e+05, 1.7495850897236078e+05, 1.7496839663536297e+05, 1.7497828434916140e+05, 1.7498817211375345e+05, 1.7499805992913659e+05, 1.7500794779530825e+05, 1.7501783571226581e+05, 1.7502772368000669e+05, 1.7503761169852832e+05, 1.7504749976782812e+05, 1.7505738788790352e+05, 1.7506727605875194e+05, 1.7507716428037078e+05, 1.7508705255275749e+05, 1.7509694087590949e+05, 1.7510682924982425e+05, 1.7511671767449912e+05, 1.7512660614993155e+05, 1.7513649467611895e+05, 1.7514638325305877e+05, 1.7515627188074836e+05, 1.7516616055918529e+05, 1.7517604928836680e+05, 1.7518593806829047e+05, 1.7519582689895367e+05, 1.7520571578035384e+05, 1.7521560471248836e+05, 1.7522549369535467e+05, 1.7523538272895027e+05, 1.7524527181327247e+05, 1.7525516094831878e+05, 1.7526505013408660e+05, 1.7527493937057335e+05, 1.7528482865777650e+05, 1.7529471799569341e+05, 1.7530460738432160e+05, 1.7531449682365838e+05, 1.7532438631370128e+05, 1.7533427585444768e+05, 1.7534416544589505e+05, 1.7535405508804080e+05, 1.7536394478088227e+05, 1.7537383452441703e+05, 1.7538372431864246e+05, 1.7539361416355596e+05, 1.7540350405915498e+05, 1.7541339400543695e+05, 1.7542328400239930e+05, 1.7543317405003947e+05, 1.7544306414835487e+05, 1.7545295429734295e+05, 1.7546284449700115e+05, 1.7547273474732690e+05, 1.7548262504831763e+05, 1.7549251539997075e+05, 1.7550240580228376e+05, 1.7551229625525398e+05, 1.7552218675887899e+05, 1.7553207731315613e+05, 1.7554196791808281e+05, 1.7555185857365653e+05, 1.7556174927987473e+05, 1.7557164003673478e+05, 1.7558153084423416e+05, 1.7559142170237031e+05, 1.7560131261114066e+05, 1.7561120357054260e+05, 1.7562109458057364e+05, 1.7563098564123115e+05, 1.7564087675251262e+05, 1.7565076791441548e+05, 1.7566065912693716e+05, 1.7567055039007508e+05, 1.7568044170382668e+05, 1.7569033306818939e+05, 1.7570022448316071e+05, 1.7571011594873801e+05, 1.7572000746491877e+05, 1.7572989903170042e+05, 1.7573979064908039e+05, 1.7574968231705611e+05, 1.7575957403562506e+05, 1.7576946580478465e+05, 1.7577935762453233e+05, 1.7578924949486551e+05, 1.7579914141578169e+05, 1.7580903338727826e+05, 1.7581892540935270e+05, 1.7582881748200240e+05, 1.7583870960522487e+05, 1.7584860177901751e+05, 1.7585849400337777e+05, 1.7586838627830308e+05, 1.7587827860379091e+05, 1.7588817097983870e+05, 1.7589806340644386e+05, 1.7590795588360389e+05, 1.7591784841131620e+05, 1.7592774098957825e+05, 1.7593763361838745e+05, 1.7594752629774131e+05, 1.7595741902763717e+05, 1.7596731180807258e+05, 1.7597720463904491e+05, 1.7598709752055167e+05, 1.7599699045259025e+05, 1.7600688343515815e+05, 1.7601677646825276e+05, 1.7602666955187160e+05, 1.7603656268601207e+05, 1.7604645587067158e+05, 1.7605634910584765e+05, 1.7606624239153773e+05, 1.7607613572773919e+05, 1.7608602911444954e+05, 1.7609592255166621e+05, 1.7610581603938667e+05, 1.7611570957760836e+05, 1.7612560316632868e+05, 1.7613549680554520e+05, 1.7614539049525527e+05, 1.7615528423545632e+05, 1.7616517802614588e+05, 1.7617507186732139e+05, 1.7618496575898022e+05, 1.7619485970111992e+05, 1.7620475369373790e+05, 1.7621464773683160e+05, 1.7622454183039849e+05, 1.7623443597443600e+05, 1.7624433016894161e+05, 1.7625422441391277e+05, 1.7626411870934695e+05, 1.7627401305524155e+05, 1.7628390745159410e+05, 1.7629380189840199e+05, 1.7630369639566270e+05, 1.7631359094337365e+05, 1.7632348554153234e+05, 1.7633338019013620e+05, 1.7634327488918271e+05, 1.7635316963866929e+05, 1.7636306443859343e+05, 1.7637295928895255e+05, 1.7638285418974416e+05, 1.7639274914096564e+05, 1.7640264414261453e+05, 1.7641253919468826e+05, 1.7642243429718429e+05, 1.7643232945010002e+05, 1.7644222465343296e+05, 1.7645211990718057e+05, 1.7646201521134030e+05, 1.7647191056590961e+05, 1.7648180597088597e+05, 1.7649170142626681e+05, 1.7650159693204961e+05, 1.7651149248823180e+05, 1.7652138809481086e+05, 1.7653128375178430e+05, 1.7654117945914951e+05, 1.7655107521690396e+05, 1.7656097102504518e+05, 1.7657086688357050e+05, 1.7658076279247753e+05, 1.7659065875176364e+05, 1.7660055476142632e+05, 1.7661045082146305e+05, 1.7662034693187123e+05, 1.7663024309264836e+05, 1.7664013930379189e+05, 1.7665003556529933e+05, 1.7665993187716813e+05, 1.7666982823939575e+05, 1.7667972465197963e+05, 1.7668962111491724e+05, 1.7669951762820603e+05, 1.7670941419184351e+05, 1.7671931080582709e+05, 1.7672920747015427e+05, 1.7673910418482250e+05, 1.7674900094982929e+05, 1.7675889776517209e+05, 1.7676879463084831e+05, 1.7677869154685549e+05, 1.7678858851319103e+05, 1.7679848552985245e+05, 1.7680838259683715e+05, 1.7681827971414267e+05, 1.7682817688176647e+05, 1.7683807409970602e+05, 1.7684797136795873e+05, 1.7685786868652212e+05, 1.7686776605539367e+05, 1.7687766347457081e+05, 1.7688756094405099e+05, 1.7689745846383178e+05, 1.7690735603391056e+05, 1.7691725365428481e+05, 1.7692715132495202e+05, 1.7693704904590966e+05, 1.7694694681715523e+05, 1.7695684463868610e+05, 1.7696674251049987e+05, 1.7697664043259394e+05, 1.7698653840496577e+05, 1.7699643642761288e+05, 1.7700633450053271e+05, 1.7701623262372272e+05, 1.7702613079718046e+05, 1.7703602902090331e+05, 1.7704592729488880e+05, 1.7705582561913438e+05, 1.7706572399363754e+05, 1.7707562241839571e+05, 1.7708552089340644e+05, 1.7709541941866715e+05, 1.7710531799417530e+05, 1.7711521661992843e+05, 1.7712511529592396e+05, 1.7713501402215939e+05, 1.7714491279863223e+05, 1.7715481162533988e+05, 1.7716471050227989e+05, 1.7717460942944966e+05, 1.7718450840684667e+05, 1.7719440743446851e+05, 1.7720430651231253e+05, 1.7721420564037628e+05, 1.7722410481865722e+05, 1.7723400404715282e+05, 1.7724390332586059e+05, 1.7725380265477797e+05, 1.7726370203390249e+05, 1.7727360146323158e+05, 1.7728350094276271e+05, 1.7729340047249341e+05, 1.7730330005242114e+05, 1.7731319968254335e+05, 1.7732309936285755e+05, 1.7733299909336117e+05, 1.7734289887405181e+05, 1.7735279870492683e+05, 1.7736269858598377e+05, 1.7737259851722012e+05, 1.7738249849863333e+05, 1.7739239853022087e+05, 1.7740229861198028e+05, 1.7741219874390896e+05, 1.7742209892600449e+05, 1.7743199915826431e+05, 1.7744189944068587e+05, 1.7745179977326674e+05, 1.7746170015600434e+05, 1.7747160058889617e+05, 1.7748150107193968e+05, 1.7749140160513239e+05, 1.7750130218847177e+05, 1.7751120282195535e+05, 1.7752110350558054e+05, 1.7753100423934488e+05, 1.7754090502324581e+05, 1.7755080585728088e+05, 1.7756070674144756e+05, 1.7757060767574332e+05, 1.7758050866016565e+05, 1.7759040969471203e+05, 1.7760031077937997e+05, 1.7761021191416695e+05, 1.7762011309907038e+05, 1.7763001433408784e+05, 1.7763991561921683e+05, 1.7764981695445481e+05, 1.7765971833979926e+05, 1.7766961977524767e+05, 1.7767952126079754e+05, 1.7768942279644634e+05, 1.7769932438219158e+05, 1.7770922601803078e+05, 1.7771912770396136e+05, 1.7772902943998086e+05, 1.7773893122608677e+05, 1.7774883306227654e+05, 1.7775873494854773e+05, 1.7776863688489777e+05, 1.7777853887132418e+05, 1.7778844090782446e+05, 1.7779834299439608e+05, 1.7780824513103656e+05, 1.7781814731774337e+05, 1.7782804955451400e+05, 1.7783795184134596e+05, 1.7784785417823674e+05, 1.7785775656518381e+05, 1.7786765900218472e+05, 1.7787756148923692e+05, 1.7788746402633793e+05, 1.7789736661348524e+05, 1.7790726925067633e+05, 1.7791717193790872e+05, 1.7792707467517987e+05, 1.7793697746248730e+05, 1.7794688029982845e+05, 1.7795678318720093e+05, 1.7796668612460213e+05, 1.7797658911202964e+05, 1.7798649214948091e+05, 1.7799639523695339e+05, 1.7800629837444468e+05, 1.7801620156195221e+05, 1.7802610479947351e+05, 1.7803600808700602e+05, 1.7804591142454729e+05, 1.7805581481209485e+05, 1.7806571824964613e+05, 1.7807562173719867e+05};

double stb_factorial(int n)
{
    if (n > STB_FACT_MAX) {
        fprintf(stderr, "Overflow!!! Use stb_log_factorial for large n.\n");
        return exp(log_fact_table[STB_FACT_MAX]);
    }

    return exp(log_fact_table[n]);

    // Old function

	// double factorial = 1;

	// while (n > 1) {
	// 	factorial = factorial * n;
	// 	n = n - 1;
	// }

	// return factorial;
}

double stb_log_factorial(int n)
{
    if (n > STB_LOG_FACT_MAX) {
        double factorial = 0;

        while (n > 1) {
            factorial = factorial + log(n);
            n = n - 1;
        }

        return factorial;
    } else {
        return log_fact_table[n];
    }
}


/* Returns and array with the ranks of the set of observations. This is needed by stb_spearman */ 
double *stb_rank(double *x, int n) 
{   
    double *rank_x = malloc(n * sizeof(double));

    for(int i = 0; i < n; i++) { 
        int r = 1, s = 1; 
          
        /* Count no of smaller elements in 0 to i-1 */ 
        for(int j = 0; j < i; j++) { 
            if (x[j] < x[i] ) {
            	r++; 
            }
            if (x[j] == x[i] ) {
            	s++;
            } 
        } 
      
        /* Count no of smaller elements in i+1 to n-1 */ 
        for (int j = i+1; j < n; j++) { 
            if (x[j] < x[i] ) {
            	r++;
            } 
            if (x[j] == x[i] ) {
            	s++;
            } 
        } 
  
        /* Use Fractional Rank formula: fractional_rank = r + (n-1)/2 */ 
        rank_x[i] = r + (s-1) * 0.5;         
    } 
      
    /* Return Rank Vector */ 
    return rank_x; 
} 
  
/* Returns Spearman's Rank correlation of two vectors x and y, each of size n */
double stb_spearman(double *x, double *y, int n) 
{  
    double sigma_x = 0, sigma_y = 0, sigma_xy = 0; 
    double sigma_xsq = 0, sigma_ysq = 0; 
	
	double *xrank = stb_rank(x, n);
	double *yrank = stb_rank(y, n);

    for (int i = 0; i < n; i++) { 
        /* sum of elements of array x and y */ 
        sigma_x += xrank[i]; 
        sigma_y += yrank[i]; 
  
        /* sum of x[i] * y[i] */ 
        sigma_xy += xrank[i] * yrank[i]; 
  
        /* sum of square of array elements */ 
        sigma_xsq += xrank[i] * xrank[i]; 
        sigma_ysq += yrank[i] * yrank[i]; 
    } 
  
    // Calculate spearman correlation coefficient.
    double num = n * sigma_xy -  sigma_x * sigma_y;
    double den = sqrt((n * sigma_xsq - sigma_x * sigma_x) *  (n * sigma_ysq - sigma_y * sigma_y)); 

	free(xrank);
	free(yrank);

    return num/den; 
} 

/* Returns Kendall's Rank correlation and the probability of two vectors x and y, each of size n */
double stb_kendall(double *x, double *y, int n, double *tau, double *z, double *prob)
{
	//double epsilon = 1e-6;
	int i, j;
	double yties = 0, xties = 0;
	double numerator=0;
	double dxdy, dx, dy;
	double var;

	for (i=0; i < n; i++){
		for (j=(i+1); j < n; j++) { 
			dx = x[i]- x[j];
			dy = y[i]- y[j];
			if (fabs(dx) < EPSILON) {
				dx = 0;
			}
			if (fabs(dy) < EPSILON) {
				dy = 0;
			}
			dxdy = dx*dy;

			if (dxdy) { // Neither array has a tie.
				yties++;
				xties++;
				dxdy > 0.0 ? numerator++ : numerator--;
			} else { //One or both arrays have ties.
				if(dy) yties++;	//an extra y event
				if(dx) xties++;	//an extra x event
			}
		}
	}

	*tau = (double) numerator / (sqrt(yties) * sqrt(xties));
	var = (double)(4.0 * n + 10.0) / (double)(9.0 * n * (n - 1.0));
	*z = (*tau) / sqrt(var);

	/* Note one sided so if testing for not equal multiply by 2! And a continuity correction is applied (0.5)*/
	//*prob = 1 - stb_phi(fabs(*z) + 0.5);

	/* two sided without continuity correction */
	*prob = 2*(1 - stb_phi(fabs(*z)));

	return *tau;
}

/* Calculate the liniear regression
 * x,y  = arrays of data
 * n = number of data points
 * a = output slope
 * b = output intercept
 * r = output correlation coefficient (can be NULL if you don't want it)
 */
int stb_linfit(const double *x, const double *y, int n, double *a, double *b, double *r)
{
	double   sumx = 0.0;                        /* sum of x                      */
	double   sumx2 = 0.0;                       /* sum of x**2                   */
	double   sumxy = 0.0;                       /* sum of x * y                  */
	double   sumy = 0.0;                        /* sum of y                      */
	double   sumy2 = 0.0;                       /* sum of y**2                   */
	int i;

	for (i = 0; i < n; i++) {
		sumx  += x[i];
		sumx2 += stb_sqr(x[i]);
		sumxy += x[i] * y[i];
		sumy  += y[i];
		sumy2 += stb_sqr(y[i]);
	}

	double denom = (n * sumx2 - stb_sqr(sumx));

	if (denom == 0) {
		/* singular matrix. can't solve the problem. */
		*a = 0;
		*b = 0;

		if (r) { *r = 0; }

		return 1;
	}

	*a = (n * sumxy  -  sumx * sumy) / denom;
	*b = (sumy * sumx2  -  sumx * sumxy) / denom;

	if (r != NULL) {
		/* compute correlation coefficient */
		*r = stb_sqr((sumxy - sumx * sumy / n)) / ((sumx2 - stb_sqr(sumx) / n) * (sumy2 - stb_sqr(sumy) / n));
	}

	return 0;
}

/* Calculate the exponential regression
 * x,y  = arrays of data
 * n = number of data points
 * a = output base
 * b = output exponent
 */
int stb_expfit(const double *x, const double *y, int n, double *a, double *b, double *r)
{
	double   sumx = 0.0;                        /* sum of x                      */
	double   sumx2 = 0.0;                       /* sum of x**2                   */
	double   sumxy = 0.0;                       /* sum of x * y                  */
	double   sumy = 0.0;                        /* sum of y                      */
	double   sumy2 = 0.0;                       /* sum of y**2                   */
	double   *Y = NULL;
	Y = calloc(n, sizeof(double));
	int i;

	for (i = 0; i < n; i++) {
		Y[i]  = log(y[i]);
	}

	for (i = 0; i < n; i++) {
		sumx  += x[i];
		sumx2 += stb_sqr(x[i]);
		sumxy += x[i] * Y[i];
		sumy  += Y[i];
		sumy2 += stb_sqr(Y[i]);
	}

	free(Y);

	double denom = (n * sumx2 - stb_sqr(sumx));

	if (denom == 0) {
		/* singular matrix. can't solve the problem. */
		*a = 0;
		*b = 0;

		if (r) { *r = 0; }

		return 1;
	}

	*a = (sumy * sumx2  -  sumx * sumxy) / denom;
	*a = exp(*a);
	*b = (n * sumxy  -  sumx * sumy) / denom;

	if (r != NULL) {
		/* compute correlation coefficient */
		*r = stb_sqr((sumxy - sumx * sumy / n)) / ((sumx2 - stb_sqr(sumx) / n) * (sumy2 - stb_sqr(sumy) / n));
	}

	return 0;
}

/* Calculate the power regression
 * y = a * x^b
 *
 * x,y  = arrays of data
 * n = number of data points
 * a = output base
 * b = output exponent
 */
int stb_powfit(const double *x, const double *y, int n, double *a, double *b, double *r)
{
	double   sumx = 0.0;                        /* sum of x                      */
	double   sumx2 = 0.0;                       /* sum of x**2                   */
	double   sumxy = 0.0;                       /* sum of x * y                  */
	double   sumy = 0.0;                        /* sum of y                      */
	double   sumy2 = 0.0;                       /* sum of y**2                   */
	double   *Y = NULL;
	double   *X = NULL;
	Y = calloc(n, sizeof(double));
	X = calloc(n, sizeof(double));
	int i;

	for (i = 0; i < n; i++) {
		Y[i] = log(y[i]);
		X[i] = log(x[i]);
	}

	for (i = 0; i < n; i++) {
		sumx  += X[i];
		sumx2 += stb_sqr(X[i]);
		sumxy += X[i] * Y[i];
		sumy  += Y[i];
		sumy2 += stb_sqr(Y[i]);
	}

	free(Y);
	free(X);

	*b = (sumxy - (sumx * sumy) / n) / (sumx2 - (stb_sqr(sumx) / n));
	*a = exp((sumy / n) - *b * (sumx / n));

	if (r != NULL) {
		/* compute correlation coefficient */
		*r = stb_sqr((sumxy - sumx * sumy / n)) / ((sumx2 - stb_sqr(sumx) / n) * (sumy2 - stb_sqr(sumy) / n));
	}

	return 0;
}

/* This Function performs Gauss-Elimination and "returns" the Upper triangular matrix and solution of 
 * equations. It passes the augmented matrix (a) as the parameter, and calculate and store the 
 * upperTriangular (Gauss-Eliminated Matrix) in it.
 */
void stb_gaussian_elimination(int m, int n, double a[m][n], double x[n-1]){
    int i,j,k;
    for(i = 0; i < (m - 1); i++){
        /* Partial Pivoting */
        for(k = i + 1; k < m; k++){
            /* If diagonal element(absolute vallue) is smaller than any of the terms below it */
            if(fabs(a[i][i]) < fabs(a[k][i])){
                /* Swap the rows */
                for(j = 0; j < n; j++){                
                    double temp;
                    temp = a[i][j];
                    a[i][j] = a[k][j];
                    a[k][j] = temp;
                }
            }
        }
        /* Begin Gauss Elimination */
        for(k = i + 1; k < m; k++){
            double  term = a[k][i] / a[i][i];
            for(j = 0; j < n; j++){
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }
         
    }
    /* Begin Back-substitution */
    for(i = m - 1; i >= 0; i--){
        x[i] = a[i][n-1];
        for(j = i + 1; j < (n - 1); j++){
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }
             
}

/* x = x-axis values
 * y = y-axis values
 * N = number of data-points
 * n = degree of polynomial
 *
 * Thanks Manas Sharma
 * Example:
 *
 * double x[5] = {0, 1, 2, 3, 4};
 * double y[5] = {1, 1.8, 1.3, 2.5, 6.3};
 * int degrees = 2;
 * double *fit;
 *
 * fit = stb_polyfit(x, y, 5, degrees);
 *
 * for (int i = 0; i < (degrees + 1); i++ {
 *     printf("%lfx^%i+", fit[i], i)	
 * }
 * printf("\n");
 */
double *stb_polyfit(const double *x, const double *y, int N, int n)
{
	int i, j;
    /* An array of size 2*n+1 for storing N, Sig xi, Sig xi^2, ...., etc. which 
     * are the independent components of the normal matrix
     */
    double X[2*n+1];
    for(i = 0; i <= 2 * n; i++){
        X[i] = 0;
        for(j = 0;j < N; j++){
            X[i] = X[i] + pow(x[j], i);
        }
    }

    /* The normal augmented matrix */
    double B[n + 1][n + 2];  
    /* Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi) */
    double Y[n + 1];      
    for(i = 0; i <= n; i++){
        Y[i] = 0;
        for(j = 0; j < N; j++){
        	/* Consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi) */
            Y[i] = Y[i] + pow(x[j], i) * y[j];
        }
    }
    
    for(i = 0; i <= n; i++){
        for(j = 0; j <= n; j++){
        	/* Build the Normal matrix by storing the corresponding coefficients at the right positions 
        	 * except the last column of the matrix
        	 */
            B[i][j]=X[i+j]; 
        }
    }
    
    for(i = 0; i <= n; i++){
    	/* Load the values of Y as the last column of B(Normal Matrix but augmented) */
        B[i][n+1]=Y[i];
    }

    double *A = calloc(n + 1, sizeof(double));
    stb_gaussian_elimination(n + 1, n + 2, B, A);

    return A;
}

/* Lagrange interpolation 
 * x = x values
 * y = y values
 * n = number of data points
 * xp = value of x to interpolate y -> f(xp) = sum (the interpolated value)
 */
double stb_lagrange(double *x, double *y, int n, double xp)
{
	double sum = 0.0;
	double product = 1.0;

	for(int i = 0; i < n; i++) {
    	product = 1.0;
    	for(int j = 0; j < n; j++) {
    		if(j != i) {
    			product = product * (xp - x[j]) / (x[i] - x[j]);
    		}
    	}
    	sum += y[i] * product;
	}

	return sum;
}

/* This function returns the area under a curve of a given function f
 * a = x[0]
 * b = x[n]
 * n = number of intervals
 * f = the curve function that is used to calculate y (the height of each slice)
 */
double stb_trap(double a, double b, int n, stb_function f)
{
	double area = 0.0;
	double sum = 0.0;
	/* Width of the intervals */
	double h = fabs(b - a) / (double) n;

	for (int x = 1; x < n; ++x) {
		sum += f(x * h + a);
	}
	area = (h / 2.0) * (2.0 * sum + f(a) + f(b));

	return area;
}

/* This function perform integration by trapezoidal rule until converged to the given accuracy
 * or till LIMIT iterations has been reached and returns the area under a curve of a given function f
 * a = x[0]
 * b = x[n]
 * n = current number of intervals
 * accuracy = the desired accuracy
 * f = the curve function that is used to calculate y (the height of each slice)
 *
 * Example:
 * double func1(double x)
 * {
 *     return (1.0 / (1.0 + stb_sqr(x)));
 * }
 * 
 * int main()
 * {
 *     int intervals = 10;
 *     double accuracy = 0.000001;
 *     double x0 = 1;
 *     double xn = 5;
 *     area = stb_trapezoidal(x0, xn, &intervals, accuracy, func1);
 *     printf("Area = %lf (required: %i intervals)\n", area, intervals);
 *     return 0;
 * }
 */
double stb_trapezoidal(double a, double b, int *n, double accuracy, stb_function f)
{
	double area;
	double new_area = stb_trap(a, b, *n, f);
	int i = *n;
	const int LIMIT = 10000;

	do { 
		area = new_area;
		i++;
		new_area = stb_trap(a, b, i, f);
		if (i == LIMIT) {
			printf("stb_trapezoidal: did not converge!\n");
			break;
		}
	} while (fabs(new_area - area) > accuracy);

	*n = i;
	return new_area;
}

/* for the initialization of xoshiro512** */
uint64_t stb_splitmix64(uint64_t seed)
{
	uint64_t z = (seed += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

/* Seed xoshiro512** */
uint64_t *stb_sxoshiro512(uint64_t seed)
{
	static uint64_t s[8];
	uint64_t t;

	s[0] = t = stb_splitmix64(seed);
	s[1] = t = stb_splitmix64(t);
	s[2] = t = stb_splitmix64(t);
	s[3] = t = stb_splitmix64(t);
	s[4] = t = stb_splitmix64(t);
	s[5] = t = stb_splitmix64(t);
	s[6] = t = stb_splitmix64(t);
	s[7] = t = stb_splitmix64(t);

	return s;
}

/* xoshiro512** Thanks David Blackman and Sebastiano Vigna */
uint64_t stb_xoshiro512(uint64_t *s)
{
	const uint64_t result_starstar = stb_rotate64(s[1] * 5, 7) * 9;
	unsigned int rotate = 21;
	const uint64_t t = s[1] << 11;

	s[2] ^= s[0];
	s[5] ^= s[1];
	s[1] ^= s[2];
	s[7] ^= s[3];
	s[3] ^= s[4];
	s[4] ^= s[5];
	s[0] ^= s[6];
	s[6] ^= s[7];

	s[6] ^= t;

	s[7] = stb_rotate64(s[7], rotate);

	return result_starstar;
}

/* stb_xoshiro512_bounded returns a uniformly distributed interger, r, in the range [0, n). */
uint64_t stb_xoshiro512_bounded(uint64_t *s, uint64_t n)
{
	if (n == UINT64_MAX - 1) {
		return stb_xoshiro512(s);
	} else {
		/* We avoid all values that would cause a skew */
		uint64_t end = UINT64_MAX / n;
		end *= n;

		/* We ignore results from rand() that fall above the limit. */
		uint64_t r;

		while ((r = stb_xoshiro512(s)) >= end);

		return r % n;
	}
}

/* PCG-XSH-RR with 64-bit state and 32-bit output */
uint32_t stb_pcg32(uint64_t *s)
{
	unsigned int rotate = (unsigned int)(*s >> 59);
	*s = *s * 6364136223846793005u + (1442695040888963407u | 1u); /* LCG see: https://en.wikipedia.org/wiki/Linear_congruential_generator (Knuth) */
	*s ^= *s >> 18; /* XSH */
	return stb_rotate32((uint32_t)(*s >> 27), rotate); /* RR */
}

/* Seed the pcg32 PRNG */
uint64_t stb_spcg32(uint64_t seed)
{
	return stb_splitmix64(seed);
}

/* stb_pcg32_bounded returns a uniformly distributed interger, r, in the range [0, n). */
uint32_t stb_pcg32_bounded(uint64_t *s, uint32_t n)
{
	if (n == UINT32_MAX - 1) {
		return stb_pcg32(s);
	} else {
		/* We avoid all values that would cause a skew */
		uint32_t end = UINT32_MAX / n;
		end *= n;

		/* We ignore results that fall above the limit. */
		uint32_t r;

		while ((r = stb_pcg32(s)) >= end);

		return r % n;
	}
}

double stb_pcg32_uniform(uint64_t *seed)
{
    return stb_uint32_to_double(stb_pcg32(seed));
}

// Gaussian (normal) random sample with mean 0 and standard deviation 1 from
// Knuth and Marsaglia and Bray, ``A Convenient Method for Generating Normal Variables''
double stb_pcg32_gauss(uint64_t *seed)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            double U1 = stb_pcg32_uniform(seed);;
            double U2 = stb_pcg32_uniform(seed);;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

// Gaussian (normal) random sample with specified mean and standard deviation
double stb_pcg32_gauss_msd(uint64_t *seed, double mean, double stdev)
{
    return mean + stdev * stb_pcg32_gauss(seed);
}

// Implementation based on "A Simple Method for Generating Gamma Variables"
// by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
// Vol 26, No 3, September 2000, pages 363-372.
// shape (alpha)  and scale (lambda)
double stb_pcg32_gamma(uint64_t *seed, double shape, double scale)
{
    double x, d, c, Z, V, U;
    int flag;
    if (shape >= 1.0) {
        d = shape - 1.0 / 3.0;
        c = 1.0 / sqrt(9.0 * d);
        flag = 1;
        while (flag) {
            Z = stb_pcg32_gauss(seed);
            //printf("Z:%lf\n", Z);
            if (Z > (-1.0 / c)) {
                V = pow(1.0 + c * Z, 3);
                U = stb_pcg32_uniform(seed);
                flag = ((log(U) > (0.5 * stb_sqr(Z) + d - d * V + d * log(V))) || (U > 1 - 0.0331 * pow(Z, 4.0)));
            }
        }
        x = d * V * scale;
    } else {
        x = stb_pcg32_gamma(seed, shape + 1, 1.0);
        x = scale * x * pow(stb_pcg32_uniform(seed), 1.0 / shape);
    }

    return x;
}

double stb_pcg32_exponential(uint64_t *seed)
{
    // exponential random sample with mean 1
    return -log(stb_pcg32_uniform(seed));
}

// exponential random sample with specified mean
double stb_pcg32_exponential_m(uint64_t *seed, double mean)
{
    // exponential random sample with mean 1
    return mean * stb_pcg32_exponential(seed);
}


// Knuth: mean (lambda)
double stb_pcg32_poisson(uint64_t *seed, const double mean)
{
    if (mean < 30) {
        double k = 0.0;
        double p = stb_pcg32_uniform(seed);
        const double target = exp(-mean);
        while (p > target) {
            p *= stb_pcg32_uniform(seed);
            k += 1.0;
        }
        return k;
    } else {
        // "Rejection method PA" from "The Computer Generation of Poisson Random Variables" by A. C. Atkinson
        // Journal of the Royal Statistical Society Series C (Applied Statistics) Vol. 28, No. 1. (1979)
        // The article is on pages 29-35. The algorithm given here is on page 32.
        // The algorithm is also described on https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
        double c = 0.767 - 3.36 / mean;
        double beta = M_PI / sqrt(3.0 * mean);
        double alpha = beta * mean;
        double k = log(c) - mean - log(beta);

        while(1)
        {
            double u = stb_pcg32_uniform(seed);
            double x = (alpha - log((1.0 - u) / u)) / beta;
            double n = floor(x + 0.5);
            if (n < 0) {
                continue;
            }
            double v = stb_pcg32_uniform(seed);
            double y = alpha - beta*x;
            double lhs = y + log(v / stb_sqr(1.0 + exp(y)));
            double rhs = k + n * log(mean) - stb_log_gamma(n + 1);
            if (lhs <= rhs) {
                return n;
            }
        }        
    }
}

double stb_pcg32_nbinom(uint64_t *seed, double size, double prob)
{
    if (prob == 1.0) {
        return 0.0;
    }

    return stb_pcg32_poisson(seed, stb_pcg32_gamma(seed, size, (1 - prob) / prob));
}

double stb_pcg32_chisquare(uint64_t *seed, double degrees_of_freedom)
{
    return stb_pcg32_gamma(seed, degrees_of_freedom / 2.0, 2.0);
}

double stb_pcg32_invgamma(uint64_t *seed, double shape, double scale)
{
    return 1.0 / stb_pcg32_gamma(seed, shape, 1.0 / scale);
}

double stb_pcg32_beta(uint64_t *seed, double a, double b)
{
    double u = stb_pcg32_gamma(seed, a, 1.0);
    double v = stb_pcg32_gamma(seed, b, 1.0);
    return u / (u + v);
}

double stb_pcg32_nbinom_mu(uint64_t *seed, double size, double mu)
{
    if (mu == 0.0) {
        return 0.0;
    }
    
    return stb_pcg32_poisson(seed, stb_pcg32_gamma(seed, size, mu / size));
}
/* -log(sqrt(2*pi)) */
#define MINUS_LOG_SQRT_2PI (-0.9189385332046727)

/* Robbins’s factorial approximation given n trials and a successes 
 * see: https://pvk.ca/Blog/2018/07/06/testing-slo-type-properties-with-the-confidence-sequence-method/
 */
double stb_robbins(uint64_t n, uint64_t a)
{
    /* Simple cases no/all succes*/
    if (a == 0 || a == n) {
        return 0.0;
    }

    /* almost all/no succes */
    if (a == 1 || a == n - 1) {
        return log(n);
    }

    /*Robbin;s factorial approximation */
    uint64_t b = n - a;
    double values[7] = {
        MINUS_LOG_SQRT_2PI,
        (n + 0.5) * log(n),
        //-n,
        1.0 / (12 * n),
        -(a + 0.5) * log(a),
        //a,
        -1.0 / (12 * a + 1),
        -(b + 0.5) * log(b),
        //b,
        -1.0 / (12 * b + 1)
    };
    // -n + a + b = 0

    return stb_sum(values, 7);
}

/**
 * Confidence Sequence Method.
 *
 * See: https://pvk.ca/Blog/2018/07/06/testing-slo-type-properties-with-the-confidence-sequence-method/
 * 
 * And
 *
 * See "A simple method for implementing Monte Carlo tests,"
 * Ding, Gandy, and Hahn, 2017 (https://arxiv.org/abs/1611.01675).
 *
 * Given n trials and a successes, can we conclude that the success rate
 * differs from alpha with exp(log_eps) false positive rate?
 *
 * Output the current log confidence level in OUT_log_level if non-NULL and
 * returns 1 if confidence is high enough to conclude that with n trails and a successes
 * the e-value is sufficientely robust
 */
int stb_csm(uint64_t n, double alpha, uint64_t a, double log_eps, double *OUT_log_level)
{
    uint64_t b = n - a;
    double values[4] = {
        log(n + 1),
        stb_robbins(n, a),
        a * log(alpha),
        b * log1p(-alpha)
    };
    double log_level = stb_sum(values, 4);

    if (OUT_log_level != NULL) {
        *OUT_log_level = log_level;
    }

    return log_level < log_eps;
}

/* Allocates a matrix of rowsize*colsize with a single call to malloc, so
 * only a single call to free is necessary
 */
void **stb_allocmat(int rowsize, int colsize, size_t item_size)
{
	void **matrix = NULL;
	int datasize = colsize * rowsize * item_size;
	int ptrsize = rowsize;
	int *rowptr;

	int totalsize = ptrsize * sizeof(void *) + datasize;

	if((matrix = (void **) calloc(totalsize, sizeof(int))) == NULL) {
		fprintf(stderr, "Can't allocate %d bytes\n.", totalsize);
		return NULL;
	}

	rowptr = (void *)(matrix + ptrsize);

	/* set pointers */
	for (int i = 0; i < rowsize; i++) {
		(matrix)[i] = rowptr + i * colsize * item_size;
	}

	return matrix;
}

/* Helper function for stb_combinations */
void stb_visit(int *c, int t, int ***comb, int which)
{
	int **mat = *comb;

	int i = 0;

	for (int j = t; j > 0; j--, i++) { 
		mat[which][i] = c[j]; 
	}

}

/* returns the number of combinations (nCr) of t out of n and returns a matrix (combinations) of them
 * Uses Knuth Algorithm T (section 2.7.1.3 from The Art of Computer Programming)
 */
int stb_combinations(int n, int t, int ***combinations)
{
	int **comb;
	int *c = calloc(t + 3, sizeof(int));
	int current = 0;
	/* The number of combinations = n! / (t! * (n-t)!) */
	int nCr = stb_factorial(n) / (stb_factorial(t) * stb_factorial(n - t));
	/* Allocate matrix */
	comb = (int **)stb_allocmat(nCr, t, sizeof(int));

	int j, x;
	/* Initialize */
	for (j = 1; j <= t; j++) {
		c[j] = j - 1;
	}

	c[t + 1] = n;
	c[t + 2] = 0;
	j = t;

	for (;;) {
		stb_visit(c, t, &comb, current++);

		while (j > 0) {
			x = j;
			/* Increase */
			c[j] = x;
			j = j - 1;
			stb_visit(c, t, &comb, current++);
		}

		while ((c[1] + 1) < c[2]) {
			c[1] = c[1] + 1;
			stb_visit(c, t, &comb, current++);

			while (j > 0) {
				x = j;
				/* Increase */
				c[j] = x;
				j = j - 1;
				stb_visit(c, t, &comb, current++);
			}
		}

		j = 2;

		/* Find J */
		c[j - 1] = j - 2;
		x = c[j] + 1;

		while (x == c[j + 1]) {
			j = j + 1;
			c[j - 1] = j - 2;
			x = c[j] + 1;
		}

		/* Done? */
		if (j > t) {
			*combinations = comb;
			free(c);
			return nCr;
		}

		c[j] = x;
		j = j - 1;
	}
}

/* Reads an entire file into an array of strings, needs only a single call to free */
char **stb_fgetlns(char *filename, size_t *number_of_lines)
{
	size_t count = 1;
	char **sfile = NULL;
	char *p;

	FILE *f = NULL;
	f = fopen(filename, "r");
	if (!f) {
		fprintf(stderr, "Unable to open: %s\n", filename);
		return NULL;
	}

	/* Determine the file size */
	fseek(f, 0, SEEK_END);
	size_t fsize = ftell(f);
	fseek(f, 0, SEEK_SET);

	/* Read the file into a temporary buffer */
	char *buffer = malloc(fsize + 1);
	fread(buffer, fsize, 1, f);
	if (!buffer) {
    	fprintf(stderr, "Failed to read %s into memory\n", filename);
    	return NULL;
    }
	buffer[fsize] = 0;
	
	/* Close the file */
	fclose(f);

	/* Count the number of new lines */
	p = buffer;
	size_t i = 0;
	while (p[i]) {
		if (p[i] == '\n') {
			if ( p[i+1] == '\r') {
				count++;
				i++;
			} else {
				count++;
			}
		} else if (*p == '\r') {
			count++;
		}
		i++;
	}

	if (number_of_lines) {
		*number_of_lines = count;
	}

	/* Allocate space to keep the entire file */
	sfile = (char **) malloc(sizeof(char *) * (count + 1) + fsize + 1);
	if (!sfile) {
		fprintf(stderr, "Could not copy the data\n");
		return NULL;
	}
	sfile[count] = NULL;
	/* Copy in the original data */
	memcpy(&sfile[count + 1], buffer, fsize + 1);
	
	free(buffer);
    buffer = (char *) &sfile[count + 1];

    /* Go over everything again and set the pointers */
	p = buffer;
	i = 0;
	count = 0;
	sfile[count] = &p[i];
	while (p[i]) {
		if (p[i] == '\n') {
			if ( p[i+1] == '\r') {
				p[i] = '\0';
				p[i+1] = '\0';
				count++;
				i++;
				if (p[i+1]) {
					sfile[count] = &p[i+1];
				}
			} else {
				p[i] = '\0';
				count++;
				if (p[i+1]) {
					sfile[count] = &p[i+1];
				}
			}
		} else if (*p == '\r') {
			p[i] = '\0';
			count++;
			if (p[i+1]) {
				sfile[count] = &p[i+1];
			}
		}
		i++;
	}

	return sfile;
}

/* Dynamic allocation version of fgets(), capable of reading "unlimited" line lengths. */
char *stb_fgetln(char **buf, int *n, FILE *fp)
{
	char *s;
	int   len;
	int   location;
	
	/* Initial cache size */
	size_t cache_sz = 1024;
	/* Thereafter, each time capacity needs to be increased,
	 * multiply the increment by this factor. 
	 */
	const size_t cache_sz_inc_factor = 2;

	if (*n == 0) {
		*buf = malloc(sizeof(char) * cache_sz);
		*n   = cache_sz;
	}

	/* We're sitting at EOF, or there's an error.*/
	if (fgets(*buf, *n, fp) == NULL)
		return NULL;

	/* We got a string AND it reached EOF. (last string without an '\n')*/
	if (feof(fp))
		return *buf;

	/* We got a complete string */
	len = strlen(*buf);
	if ((*buf)[len-1] == '\n')
		return *buf;

	/* We have an incomplete string and we have to extend the buffer.
	 * We make sure we overwrite the previous fgets \0 we got in the 
	 * first step (and subsequent steps (location - 1) and append 
	 * the new buffer until we find the \n or EOF.
	 */
	location = len;
	while (1) {
		*n  *= cache_sz_inc_factor;
		*buf = realloc(*buf, sizeof(char) * (*n));
		/* Append to previous buf */
		s = *buf + location;
		
		if (fgets(s, (*n - location), fp) == NULL)
			return *buf;
		
		if (feof(fp))
			return *buf;

		len = strlen(s);
		if (s[len - 1] == '\n')
			return *buf;

		location = *n - 1;
	}
}

void stb_free_csv_line(char **parsed)
{
	char **ptr;

	for (ptr = parsed; *ptr; ptr++) {
		free(*ptr);
	}

	free(parsed);
}

int stb_count_fields(const char *line, const char delim)
{
	const char *ptr;
	int cnt, fQuote;

	for (cnt = 1, fQuote = 0, ptr = line; *ptr; ptr++) {
		if (fQuote) {
			if (*ptr == CSV_QUOTE) {
				if (ptr[1] == CSV_QUOTE) {
					ptr++;
					continue;
				}
				fQuote = 0;
			}
			continue;
		}

		if (*ptr == CSV_QUOTE) {
			fQuote = 1;
			continue;
		} else if (*ptr == delim) {
			cnt++;
			continue;
		} else {
			continue;
		}
	}

	if (fQuote) {
		return -1;
	}

	return cnt;
}

/*
 *  Given a string containing no linebreaks, or containing line breaks
 *  which are escaped by "double quotes", extract a NULL-terminated
 *  array of strings, one for every cell in the row.
 */
char **stb_parse_csv(const char *line, const char delim, int *nrfields)
{
	char **buf, **bptr, *tmp, *tptr;
	const char *ptr;
	int fieldcnt = 0, fQuote, fEnd;

	// If we did not get a count for the number of fields determine it
	if (!*nrfields) {
		fieldcnt = stb_count_fields(line, delim);
		*nrfields = fieldcnt;
	}

	if (fieldcnt <= -1) {
		fprintf(stderr, "No tokens found!\n");
		return NULL;
	}

	buf = calloc((fieldcnt + 1), sizeof(char*));
	if (!buf) {
		return NULL;
	}

	tmp = calloc(strlen(line) + 1, sizeof(char));
	if (!tmp) {
		free(buf);
		return NULL;
	}

	bptr = buf;
	for (ptr = line, fQuote = 0, tptr = tmp, fEnd = 0; !fEnd; ptr++) {
		if (fQuote) {
			if (!*ptr) {
				break;
			}

			if (*ptr == CSV_QUOTE) {
				if (ptr[1] == CSV_QUOTE) {
					*tptr++ = CSV_QUOTE;
					ptr++;
					continue;
				}
				fQuote = 0;
			} else {
				*tptr++ = *ptr;
			}

			continue;
		}

		if (*ptr == CSV_QUOTE) {
			fQuote = 1;
			continue;
		} else if ((*ptr == '\0') || (*ptr == delim)) {
			if (*ptr == '\0') {
				fEnd = 1;
			}
			*tptr = '\0';
			*bptr = strdup(tmp);
			tptr = tmp;

			if (!*bptr) {
				for (bptr--; bptr >= buf; bptr--) {
					free(*bptr);
				}
				free(buf);
				free(tmp);

				return NULL;
			}

			bptr++;

			if (fEnd) {
				break;
			} else {
				continue;
			}
		} else {
			*tptr = *ptr;
			tptr++;
			continue;
		}
	}

	*bptr = NULL;
	free(tmp);
	return buf;
}

int stb_double_almost_equal(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

int stb_double_equal(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

int stb_double_greater(double a, double b, double epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

int stb_double_less(double a, double b, double epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

int dcmp(const void *a, const void *b)
{
#define _D(x) *(double*)x
	if (stb_double_greater(_D(a), _D(b), EPSILON)) {
		return 1;
	} else if (stb_double_less(_D(a), _D(b), EPSILON)) {
		return -1;
	} else {
		return 0;
	}
#undef _D
}

int icmp(const void *a, const void *b)
{
#define _I(x) *(int*)x
	if (_I(a) > _I(b)) {
		return 1;
	} else if (_I(a) < _I(b)) {
		return -1;
	} else {
		return 0;
	}
#undef _I
}

/* Moment estimation of the parameters of a scaled F-distribution (the prior) 
 * The first degrees of freedom is given 
 */
void stb_fit_f_dist(double *var, int len, int df1, double *pvar, double *pdf2)
{
    double median = stb_median(var, len);
    // More than half of residual variances are exactly zero
    if (stb_double_almost_equal(median, 0.0, 1e-5)) {
        printf("More than half of residual variances are exactly zero\n");
        median = 1.0;
    }

    // Zero sample variances will be offset
    double *z = malloc(len * sizeof(double));
    double *e = malloc(len * sizeof(double));
    double c1 = stb_polygamma(0, df1 / 2.0);
    double c2 = log(df1 / 2.0);
    for (int i = 0; i < len; i++) {
        // Better to work on with log(F)
        z[i] = log(fmax(var[i], 1e-5 * median));
        e[i] = z[i] - c1 + c2;
    }

    double tevar, evar, emean = 0.0;
    stb_meanvar(e, len, &emean, &tevar);
    double *etvar =  malloc(len * sizeof(double));
    c1 = (double) len / (double) (len - 1);
    c2 = stb_polygamma(1, df1 / 2.0);
    for (int i = 0; i < len; i++) {
        etvar[i] = c1 * stb_sqr(e[i] - emean) - c2;
    }
    stb_meanvar(etvar, len, &evar, &tevar);
    if (evar > 0.0) {
        *pdf2 = 2 * stb_trigamma_inverse(evar);
        *pvar = exp(emean + stb_polygamma(0, *pdf2 / 2.0) - log(*pdf2 / 2.0));
    } else {
        *pdf2 = stb_pos_inf;
        *pvar = exp(emean);
    }

    free(z);
    free(e);
    free(etvar);
}

/* Filter items in place and return list of unique numbers, their count and the tumber of unique numbers.
 * NOTE: if a separate list is desired, duplicate it before calling this function 
 */
int stb_dunique(double *values, double **counts, int len)
{
	int i, j;
	double count = 0;
	double *tcounts = calloc(len, sizeof(double));

	qsort(values, len, sizeof(double), dcmp);
	
	tcounts[0] = 0;
	for (i = j = 0; i < len; i++) {
		if (dcmp(&values[j], &values[i])) {
			//printf("%f != %f %f\n", values[j], values[i], count);
			values[++j] = values[i];
			tcounts[j - 1] = count;
			count = 1;
		} else {
			count++;
		}
	}
	tcounts[j] = count;

	*counts = tcounts;

	return j + 1;
}

/* Filter items in place and return list of unique numbers, their count and the tumber of unique numbers.
 * NOTE: if a separate list is desired, duplicate it before calling this function 
 */
int stb_iunique(int *values, int **counts, int len)
{
	int i, j;
	int count = 0;
	int *tcounts = calloc(len, sizeof(int));

	qsort(values, len, sizeof(int), icmp);
	
	tcounts[0] = 0;
	for (i = j = 0; i < len; i++) {
		if (icmp(&values[j], &values[i])) {
			//printf("%f != %f %f\n", values[j], values[i], count);
			values[++j] = values[i];
			tcounts[j - 1] = count;
			count = 1;
		} else {
			count++;
		}
	}
	tcounts[j] = count;

	*counts = tcounts;

	return j + 1;
}

void stb_matrix_print(STB_MAT *matrix)
{
	for (int i = 0; i < matrix->rows; i++) {
		for (int n = 0; n < matrix->columns; n++) {
			printf("%E\t", matrix->data[i][n]);
		}
		printf("\n");
	}
}

void stb_free_matrix(STB_MAT *matrix)
{
	free(matrix->data);
	free(matrix);
}

STB_MAT *stb_new_matrix(int rows, int columns)
{
	STB_MAT *matrix;
	matrix = (STB_MAT *) calloc(1, sizeof(STB_MAT));
	matrix->rows = rows;
	matrix->columns = columns;
	matrix->data = (double **) stb_allocmat(rows, columns, sizeof(double));

	return matrix;
}

STB_MAT *stb_identity_matrix(int rows, int columns)
{
	if (rows != columns) {
		printf("An identity matrix is square! So the rows should be equal to the columns\n");
		exit(1);
	}

	STB_MAT *matrix;
	matrix = (STB_MAT *) calloc(1, sizeof(STB_MAT));
	matrix->rows = rows;
	matrix->columns = columns;
	matrix->data = (double **) stb_allocmat(rows, columns, sizeof(double));

	for (int i = 0; i < matrix->columns; i++) { 
		matrix->data[i][i] = 1.0;
	}

	return matrix;
}

STB_MAT *stb_dup_matrix(STB_MAT *matrix)
{
	STB_MAT *dup;
	dup = (STB_MAT *) calloc(1, sizeof(STB_MAT));
	dup->rows = matrix->rows;
	dup->columns = matrix->columns;
	dup->data = (double **) stb_allocmat(matrix->rows, matrix->columns, sizeof(double));

	for (int i = 0; i < matrix->rows; i++) {
		for (int n = 0; n < matrix->columns; n++) {
			dup->data[i][n] = matrix->data[i][n];
		}
	}

	return dup;
}

void stb_fill_matrix(STB_MAT *matrix, double value)
{
	for (int i = 0; i < matrix->rows; i++) {
		for (int n = 0; n < matrix->columns; n++) {
			matrix->data[i][n] = value;
		}
	}
}

/* A |   B   =    C 
 *
 * 1 | 1 0 0   1 1 0 0
 * 1 | 0 1 0 = 1 0 1 0
 * 1 | 0 0 1   1 0 0 1
 */
void stb_join_matrix(STB_MAT *A, STB_MAT *B, STB_MAT **D)
{
	if (A->rows != B->rows) {
		printf("Unable to join to matrixes with unequal row sizes!\n");
		exit(1);
	}

	STB_MAT *C = stb_new_matrix(A->rows, A->columns + B->columns);

	/* Copy in A */
	for (int i = 0; i < A->rows; i++) { 
		for (int n = 0; n < A->columns; n++) { 
			C->data[i][n] = A->data[i][n];
		}
	}

	/* Copy in B */
	for (int i = 0; i < B->rows; i++) { 
		for (int n = 0; n < B->columns; n++) { 
			C->data[i][A->columns + n] = B->data[i][n];
		}
	}

	*D = C;
}

/* D = A * B */
void stb_matrix_multiply(STB_MAT *A, STB_MAT *B, STB_MAT **D)
{
	if (A->columns != B->rows) {
		printf("Cannot multiply matrix: rows of A need to be equal to the columns of B!\n");
		exit(1);
	}

	STB_MAT *C = stb_new_matrix(A->rows, B->columns);

	for (int i = 0; i < A->rows; i++) {
		for (int n = 0; n < A->columns; n++) {
			for (int j = 0; j < B->columns; j++) {
				C->data[i][j] += A->data[i][n] * B->data[n][j];
			}
		}
	}

	*D = C;
}

/* D = A - B */
void stb_matrix_subs(STB_MAT *A, STB_MAT *B, STB_MAT **D)
{
	if (A->columns != B->columns) {
		printf("Cannot substract matrix: columns of A need to be equal to the columns of B!\n");
		exit(1);
	} else 	if (A->rows != B->rows) {
		printf("Cannot substract matrix: rows of A need to be equal to the rows of B!\n");
		exit(1);
	}

	STB_MAT *C = stb_new_matrix(A->rows, A->columns);
	
	for (int i = 0; i < A->rows; i++) {
		for (int n = 0; n < A->columns; n++) {
			C->data[i][n] = A->data[i][n] - B->data[i][n];
		}
	}

	*D = C;
}

/* D = A + B */
void stb_matrix_add(STB_MAT *A, STB_MAT *B, STB_MAT **D)
{
	if (A->columns != B->columns) {
		printf("Cannot substract matrix: columns of A need to be equal to the columns of B!\n");
		exit(1);
	} else 	if (A->rows != B->rows) {
		printf("Cannot substract matrix: rows of A need to be equal to the rows of B!\n");
		exit(1);
	}

	STB_MAT *C = stb_new_matrix(A->rows, A->columns);
	
	for (int i = 0; i < A->rows; i++) {
		for (int n = 0; n < A->columns; n++) {
			C->data[i][n] = A->data[i][n] + B->data[i][n];
		}
	}

	*D = C;
}

/* A => A transposed (B)*/
void stb_transpose_matrix(STB_MAT *A, STB_MAT **Atransposed)
{ 
	STB_MAT *B = stb_new_matrix(A->columns, A->rows);

	for (int i = 0; i < A->columns; i++) { 
		for (int n = 0; n < A->rows; n++) { 
			B->data[i][n] = A->data[n][i];
		}
	}

	*Atransposed = B;
}

/* Invert matrix */
void stb_invert_matrix(STB_MAT *A, STB_MAT **Ainverted)
{
	STB_MAT *temp = NULL;
	STB_MAT *B = stb_identity_matrix(A->rows, A->columns);

	stb_join_matrix(A, B, &temp);

    for (int i = 0; i < A->rows; i++) {
        /* Search for maximum in this column */
        double maxEl = fabs(temp->data[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < A->rows; k++) {
            if (fabs(temp->data[k][i]) > maxEl) {
                maxEl = temp->data[k][i];
                maxRow = k;
            }
        }

        /* Swap maximum row with current row (column by column) */
        for (int k = i; k < 2 * A->rows;k++) {
            double tmp = temp->data[maxRow][k];
            temp->data[maxRow][k] = temp->data[i][k];
            temp->data[i][k] = tmp;
        }

        /* Make all rows below this one 0 in current column */
        for (int k = i + 1; k < A->rows; k++) {
            double c = -temp->data[k][i] / temp->data[i][i];
            for (int j=i; j < 2 * A->rows; j++) {
                if (i == j) {
                    temp->data[k][j] = 0;
                } else {
                    temp->data[k][j] += c * temp->data[i][j];
                }
            }
        }
    }

    /* Solve equation Ax=b for an upper triangular matrix A */
    for (int i = A->rows - 1; i >= 0; i--) {
        for (int k = A->rows; k < 2 * A->rows; k++) {
            temp->data[i][k] /= temp->data[i][i];
        }

        for (int rowModify = i - 1; rowModify >= 0; rowModify--) {
            for (int columModify = A->rows; columModify < 2 * A->rows; columModify++) {
                temp->data[rowModify][columModify] -= temp->data[i][columModify] * temp->data[rowModify][i];
            }
        }
    }

    /* Copy the final inverted matrix */
	for (int i = 0; i < A->rows; i++) {
		for (int k = 0; k < A->rows; k++) {
			B->data[i][k] = temp->data[i][k + A->rows];
		}
	}

	stb_free_matrix(temp);

	*Ainverted = B;
}

void stb_matrix_multiply_by_value(STB_MAT *A, double value)
{
	for (int i = 0; i < A->rows; i++) {
		for (int n = 0; n < A->columns; n++) {
				A->data[i][n] = A->data[i][n] * value;
		}
	}
}

void stb_matrix_divide_by_value(STB_MAT *A, double value)
{
	for (int i = 0; i < A->rows; i++) {
		for (int n = 0; n < A->columns; n++) {
				A->data[i][n] = A->data[i][n] / value;
		}
	}
}

double stb_matrix_sum(STB_MAT *A)
{
	double sum = 0.0;

	for (int i = 0; i < A->rows; i++) {
		for (int n = 0; n < A->columns; n++) {
				sum += A->data[i][n];
		}
	}

	return sum;
}

STB_MAT *stb_matrix_from_double(double **data, int rows, int columns)
{
	STB_MAT *dup;
	dup = (STB_MAT *) calloc(1, sizeof(STB_MAT));
	dup->rows = rows;
	dup->columns = columns;
	dup->data = (double **) stb_allocmat(rows, columns, sizeof(double));

	for (int i = 0; i < rows; i++) {
		for (int n = 0; n < columns; n++) {
			dup->data[i][n] = data[i][n];
		}
	}

	return dup;
}

STB_MAT *stb_matrix_from_file(char *filename)
{
	FILE *fin = NULL;
	if ((fin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Unable to open: %s\n", filename);
		exit(1);
	}

	/* first line indicates rows and columns */
	int rows, columns;
	if (fscanf(fin, "%i %i", &rows, &columns) != 2) {
		fprintf(stderr, "First line should indicate the rows and columns\n");
		exit(1);
	}

	STB_MAT *matrix = NULL;
	matrix = stb_new_matrix(rows, columns);
    if (!matrix) {
        fprintf(stderr, "Error creating matrix of %ix%i\n", rows, columns);
        return NULL;
    }

	/* the rest of the lines are the matrix */
	for (int i = 0; i < rows; i++) {
		for (int n = 0; n < columns; n++) {
			fscanf(fin, "%lf", &matrix->data[i][n]);
		}
	}

	/* We are done, close the file */
	fclose(fin);
	return matrix;
}


#define QNORM_SWAP(A,B,SIZE)                        \
    do {                                             \
        char *a_byte = A;                            \
        char *b_byte = B;                            \
                                                     \
        int idxA = ((double *)a_byte - values);        \
        int idxB = ((double *)b_byte - values);        \
                                                     \
        double dtemp;                                \
        dtemp = values[idxA];                          \
        values[idxA] = values[idxB];                     \
        values[idxB] = dtemp;                          \
                                                     \
        int itemp;                                   \
        itemp = tmpidx[idxA];                         \
        tmpidx[idxA] = tmpidx[idxB];                   \
        tmpidx[idxB] = itemp;                         \
                                                     \
    } while (0)

int stb_didx(double *values, int **idx, int **ranks, int **counts, int len)
{
        int current, prev, i;
        double count = 0;
        int *tmpcounts = calloc(len, sizeof(double));
        int *tmpranks = calloc(len, sizeof(double));
        int *tmpidx = calloc(len, sizeof(int));
        for (i = 0; i < len; i++) {
                tmpidx[i] = i;
                //printf("%lf, %i\n", values[i], tmpidx[i] );
        }

// use the swap function to co-sort a dependent array
#undef STB_QSORT_SWAP
#define STB_QSORT_SWAP(A,B,SIZE) QNORM_SWAP(A, B, SIZE)
    stb_qsort(values, len, sizeof(double), dcmp);
#undef STB_QSORT_SWAP

        tmpcounts[0] = 0;
        int curr_rank = 0;
        tmpranks[0] = curr_rank;
        tmpcounts[curr_rank]++;
        for (current = 1, prev = 0; current < len; current++, prev++) {
                if (dcmp(&values[prev], &values[current])) {
                        // printf("curr_rank counts: %i\n", tmpcounts[curr_rank]);
                        // printf("%lf != %lf\n", values[prev], values[current]);
                        curr_rank += tmpcounts[curr_rank];
                        tmpranks[current] = curr_rank;
                        tmpcounts[curr_rank]++;
                        // printf("curr_rank: %i\n", curr_rank);
                } else {
                        tmpcounts[curr_rank]++;
                        tmpranks[current] = curr_rank;
                }
        }

        *counts = tmpcounts;
        *ranks = tmpranks;
        *idx = tmpidx;

        return tmpcounts[curr_rank] + 1;
}

/* Perform quantile normalization between columns without a reference
In statistics, quantile normalization is a technique for making two distributions identical in statistical properties.
see: https://en.wikipedia.org/wiki/Quantile_normalization
Arrays 1 to 3, genes A to D
A    5    4    3
B    2    1    4
C    3    4    6
D    4    2    8
Converted to STB_MAT file test.mat:
4 3
5 4 3
2 1 4
3 4 6
4 2 8

will become:
5.666667E+00 5.166667E+00 2.000000E+00 
2.000000E+00 2.000000E+00 3.000000E+00 
3.000000E+00 5.166667E+00 4.666667E+00 
4.666667E+00 3.000000E+00 5.666667E+00

statistics:
Min.   :2.000   Min.   :2.000   Min.   :2.000  
1st Qu.:2.750   1st Qu.:2.750   1st Qu.:2.750  
Median :3.833   Median :4.083   Median :3.833  
Mean   :3.833   Mean   :3.833   Mean   :3.833  
3rd Qu.:4.917   3rd Qu.:5.167   3rd Qu.:4.917  
Max.   :5.667   Max.   :5.167   Max.   :5.667  
 */
void stb_qnorm_matrix(STB_MAT *original)
{
    STB_MAT *transposed = NULL;
    stb_transpose_matrix(original, &transposed);

    int **ranks = malloc(sizeof(int *) * original->columns);
    int **counts = malloc(sizeof(int *) * original->columns);
    int **idx = malloc(sizeof(int *) * original->columns);

    /* For each row determine a rank from lowest to highest taking into consideration possible ties
     * at the same time, we rearrange that first set of row values so each row is in order going lowest to highest value
     */
    for (int i = 0; i < transposed->rows; i++) {
        stb_didx(transposed->data[i], &idx[i], &ranks[i], &counts[i], transposed->columns);
    }

    /* Now that everything is sorted, we find the mean for each row to determine the value of the ranks (this is for quantile normalization without reference,
     * otherwise it is the reference value)
     */
    double *means = malloc(sizeof(double) * original->rows);
    double var;
    STB_MAT *sorted = NULL;
    stb_transpose_matrix(transposed, &sorted);
    stb_free_matrix(transposed);
    for (int i = 0; i < sorted->rows; i++) {
        stb_meanvar(sorted->data[i], sorted->columns, &means[i], &var);
    }

    stb_free_matrix(sorted);

    // Now take the ranking order and substitute in new values at the original locations. 
    for (int i = 0; i < original->columns; i++) {
        for (int n = 0; n < original->rows; n++) {
            // find the index in the original dataset
            original->data[idx[i][n]][i] = 0.0;
            //Note that when, values are tied in rank, they should instead be assigned the mean of the values.
            for (int j = 0; j < counts[i][ranks[i][n]]; j++) {
                original->data[idx[i][n]][i] += means[ranks[i][n] + j]; 
            }
            original->data[idx[i][n]][i] /=  (double) counts[i][ranks[i][n]];
        }
    }

    // Cleanup
    for (int i = 0; i < original->columns; i++) {
        free(idx[i]);
        free(ranks[i]);
        free(counts[i]);
    }

    free(idx);
    free(ranks);
    free(counts);

    free(means);
}

void stb_qnorm(double **data, int rows, int columns)
{
    // Transform data into STB_)MAT
    STB_MAT *original = malloc(sizeof(STB_MAT));
    original->data = data;
    original->rows = rows;
    original->columns = columns;

    // Call the matrix qnorm function
    stb_qnorm_matrix(original);

    // Free the temproary STB_MAT struct
    free(original);
}

/* Perform quantile normalization between columns without a reference
In statistics, quantile normalization is a technique for making two distributions identical in statistical properties.

NOTE:
Given a reference distribution, the original distribution is normalized by replacing each of its values by the value of the variable with the same rank 
in the reference distribution. If the reference distribution contains multiple samples, the original and reference distributions will only be identical if 
the reference distribution is first quantile normalized across all samples!

see: https://en.wikipedia.org/wiki/Quantile_normalization
Arrays 1 to 3, genes A to D
A    5    4    3
B    2    1    4
C    3    4    6
D    4    2    8
Converted to STB_MAT file test.mat:
4 3
5 4 3
2 1 4
3 4 6
4 2 8

Using the QUNATILE NORMALIZED reference:
4 3
5.666667E+00 5.166667E+00 2.000000E+00 
2.000000E+00 2.000000E+00 3.000000E+00 
3.000000E+00 5.166667E+00 4.666667E+00 
4.666667E+00 3.000000E+00 5.666667E+00

will become:
5.666667E+00 5.166667E+00 2.000000E+00 
2.000000E+00 2.000000E+00 3.000000E+00 
3.000000E+00 5.166667E+00 4.666667E+00 
4.666667E+00 3.000000E+00 5.666667E+00

statistics:
Min.   :2.000   Min.   :2.000   Min.   :2.000  
1st Qu.:2.750   1st Qu.:2.750   1st Qu.:2.750  
Median :3.833   Median :4.083   Median :3.833  
Mean   :3.833   Mean   :3.833   Mean   :3.833  
3rd Qu.:4.917   3rd Qu.:5.167   3rd Qu.:4.917  
Max.   :5.667   Max.   :5.167   Max.   :5.667  
 */
void stb_qnorm_matrix_with_reference(STB_MAT *original, STB_MAT *reference)
{
    STB_MAT *transposed = NULL;
    stb_transpose_matrix(original, &transposed);

    int **ranks = malloc(sizeof(int *) * original->columns);
    int **counts = malloc(sizeof(int *) * original->columns);
    int **idx = malloc(sizeof(int *) * original->columns);
    /* For each row determine a rank from lowest to highest taking into consideration possible ties
     * at the same time, we rearrange that first set of row values so each row is in order going lowest to highest value
     */
    for (int i = 0; i < transposed->rows; i++) {
        stb_didx(transposed->data[i], &idx[i], &ranks[i], &counts[i], transposed->columns);
    }

    stb_free_matrix(transposed);

    STB_MAT *reference_transposed = NULL;
    stb_transpose_matrix(original, &reference_transposed);

    // Sort the reference:
    // Note: Given a reference distribution, a target distribution is normalized by replacing each of its values by the value of the variable with the same rank 
    // in the reference distribution. If the reference distribution contains multiple samples, the target and reference distributions will only be identical if 
    // the reference distribution is first quantile normalized across all samples.
    for (int i = 0; i < reference_transposed->rows; i++) {
        // use the default swap function
    #undef STB_QSORT_SWAP
    #define STB_QSORT_SWAP(A,B,SIZE) STB_DEFAULT_SWAP(A, B, SIZE)
        stb_qsort(reference_transposed->data[i], reference_transposed->columns, sizeof(double), dcmp);
    }

    /* Now that everything is sorted, we find the mean for each row to determine the value of the ranks (this is for quantile normalization without reference,
     * otherwise it is the reference value)
     */
    double *means = malloc(sizeof(double) * reference->rows);
    double var;
    STB_MAT *reference_sorted = NULL;
    stb_transpose_matrix(reference_transposed, &reference_sorted);
    stb_free_matrix(reference_transposed);
    for (int i = 0; i < reference_sorted->rows; i++) {
        stb_meanvar(reference_sorted->data[i], reference_sorted->columns, &means[i], &var);
    }

    stb_free_matrix(reference_sorted);

    // Now take the ranking order and substitute in new values at the original locations. 
    for (int i = 0; i < original->columns; i++) {
        for (int n = 0; n < original->rows; n++) {
            // find the index in the original dataset
            original->data[idx[i][n]][i] = 0.0;
            //Note that when, values are tied in rank, they should instead be assigned the mean of the values.
            for (int j = 0; j < counts[i][ranks[i][n]]; j++) {
                original->data[idx[i][n]][i] += means[ranks[i][n] + j]; 
            }
            original->data[idx[i][n]][i] /=  (double) counts[i][ranks[i][n]];
        }
    }

    // Cleanup
    for (int i = 0; i < original->columns; i++) {
        free(idx[i]);
        free(ranks[i]);
        free(counts[i]);
    }
    free(idx);
    free(ranks);
    free(counts);

    free(means);
}

void stb_qnorm_with_reference(double **data, double **reference, int rows, int columns)
{
    // Transform data into STB_)MAT
    STB_MAT *original = malloc(sizeof(STB_MAT));
    original->data = data;
    original->rows = rows;
    original->columns = columns;

    // Transform data into STB_)MAT
    STB_MAT *target = malloc(sizeof(STB_MAT));
    target->data = reference;
    target->rows = rows;
    target->columns = columns;


    // Call the matrix qnorm function
    stb_qnorm_matrix_with_reference(original, target);

    // Free the temproary STB_MAT struct
    free(original);
    free(target);
}

/* Perform a simple linear regression and return a vector containing the Beta values, the T-test values 
 * and the corresponding P-values. The formula determined using the least squared method is:
 * Y = Beta[0] + Beta[1] * X[0] + Beta[2] * X[1] + Beta[n] * X[n-1]
 * 
 * Note: This can also be calculated using a design matrix (1 on first column and X values for the rest)
 */
void stb_multi_linear_regression(STB_MAT *A, STB_MAT *Y, double **beta, double **tvalue, double **pvalue)
{
	/* Make a matrix containing 1s that we can add to A to make the design matrix X */
	STB_MAT *temp = NULL;
	temp = stb_new_matrix(A->rows, 1);
	stb_fill_matrix(temp, 1.0);

	/* Add the temp matrix (for the intercepts (column of 1s)) to the matrix with the X values (A),
	 * this makes the design matrix X.
	 */
	STB_MAT *X = NULL;
	stb_join_matrix(temp, A, &X);

	/* Transpose X for: Beta = (Xtranspose*X)inverse * Xtranspose * Y */
	STB_MAT *Xt = NULL;
	stb_transpose_matrix(X, &Xt);

	/* Xm = (Xtranspose*X) */
	STB_MAT *Xm = NULL;
	stb_matrix_multiply(Xt, X, &Xm);

	/* Xinv = (Xtranspose*X)inverse */
	STB_MAT *Xinv = NULL;
	stb_invert_matrix(Xm, &Xinv);

	/* P = Xtrans * Y */
	STB_MAT *P = NULL;
	stb_matrix_multiply(Xt, Y, &P);

	/* Beta = (Xtranspose*X)inverse * Xtranspose * Y == Xinv * P */
	STB_MAT *Beta = NULL;
	stb_matrix_multiply(Xinv, P, &Beta);

	/* Allocate space to return beta */
	double *betaret = NULL;
	betaret = calloc(X->columns, sizeof(double));
	/* copy in the beta values */
	for (int i = 0; i < X->columns; i++) {
		betaret[i] = Beta->data[i][0];
	}
	*beta = betaret;

	/* Calculate the residuals == real Y minus the predicted Y = (Y - Y') == Y - X * Beta */
	STB_MAT *Ypredicted;
	stb_matrix_multiply(X, Beta, &Ypredicted);

	STB_MAT *Residuals;
	stb_matrix_subs(Y, Ypredicted, &Residuals);
	
	/* Determine the sum of squares */
	double sumofsquares = 0.0;
	for (int i = 0; i < Residuals->rows; i++) {
		for (int n = 0; n < Residuals->columns; n++) {
			sumofsquares += stb_sqr(Residuals->data[i][n]);
		}
	}

	/* The Means Square Error = SSE / (number of data values) - (number of predictor values + 1) */
	double MSE = sumofsquares / (A->rows - (A->columns + 1));

	/* Multiply the (Xtranspose*X)inverse with the MSE to get the variance-covariance matrix 
	 * of the estimated regression coefficients. The diagonal will have the std.err.
	 */
	stb_matrix_multiply_by_value(Xinv, MSE);

	/* Determine the estimated std.err to normalize the regression coefficients */
	STB_MAT *stdev;
	stdev = stb_new_matrix(Beta->rows, Beta->columns);

	for (int i = 0; i < Xinv->rows; i++) {
		for (int n = 0; n < Xinv->columns; n++) {
			if (i == n) {
				/* std.dev = sqrt(std.err.) <- the diagonal of the variance-covariance matrix */
				stdev->data[i][0] = sqrt(Xinv->data[i][n]);
			}
		}
	}

	/* Allocate space to return the T-value and the P-value */
	double *tvalret = NULL;
	tvalret = calloc(X->columns, sizeof(double));
	double *pvalret = NULL;
	pvalret = calloc(X->columns, sizeof(double));

	/* Normalize values (e.a. divide by the standard deviation) and lookup the T-value to get the p-value */
	/* if the p-values are all extremely small, it indicates that the predictor are significantly related to the response */
	for (int i = 0; i < Beta->rows; i++) {
		tvalret[i] = Beta->data[i][0] / stdev->data[i][0];
		if (tvalret[i] > 0) {
			pvalret[i] = 2 * (1 - stb_cdf_student_t(Beta->data[i][0] / stdev->data[i][0], (A->rows - (A->columns + 1))));	
		} else {
			pvalret[i] = 2 * (1 - stb_cdf_student_t(-Beta->data[i][0] / stdev->data[i][0], (A->rows - (A->columns + 1))));		
		}
	}

	*tvalue = tvalret;
	*pvalue = pvalret;

	/* Not needed anymore */
	stb_free_matrix(temp);
	stb_free_matrix(X);
	stb_free_matrix(Xt);
	stb_free_matrix(Xm);
	stb_free_matrix(Xinv);
	stb_free_matrix(P);
	stb_free_matrix(Beta);
	stb_free_matrix(Ypredicted);
	stb_free_matrix(Residuals);
	stb_free_matrix(stdev);
}

/* Perform a simple logistic regression and return a vector containing the Beta values, the Z-test values 
 * , the corresponding P-values and return the log-likelihood. The formula determined using the newton method is:
 * ln(1 / 1 - Y) = Beta[0] + Beta[1] * X[0] + Beta[2] * X[1] + Beta[n] * X[n-1]
 *
 * Note: This can also be calculated using a design matrix (1 on first column and X values for the rest)
 */
double stb_multi_logistic_regression(STB_MAT *A, STB_MAT *Y, double **beta, double **zvalue, double **pvalue)
{
	STB_MAT *Pi = NULL;
	STB_MAT *V = NULL;
	STB_MAT *Xt = NULL;
	STB_MAT *Xm = NULL;
	STB_MAT *Xmm = NULL;
	STB_MAT *Xinv = NULL;
	STB_MAT *Ym;
	STB_MAT *P = NULL;
	STB_MAT *Beta1 = NULL;
	STB_MAT *Beta = NULL;
	double loglikelihood = 0.0;
	double oldloglikelihood = 0.0;
	double delta = 0.0;
	/* Make a matrix containing 1s that we can add to A to make the design matrix X */
	STB_MAT *temp = NULL;
	temp = stb_new_matrix(A->rows, 1);
	stb_fill_matrix(temp, 1.0);

	/* Make the design matrix including intercepts (column of 1s; temp) */
	STB_MAT *X = NULL;
	stb_join_matrix(temp, A, &X);

	/* Not needed anymore */
	stb_free_matrix(temp);

	/* Initialize Beta0 to values close to 0 (0.01) except for the intercept (1.0), this hopefully leads to a
	 * lower amount of iterations needed to converge.
	 */
	STB_MAT *Beta0 = NULL;
	Beta0 = stb_new_matrix(X->columns, 1);
	stb_fill_matrix(Beta0, 0.01);
	Beta0->data[0][0] = 1.0;

	/* Transpose X for: Beta = (Xtranspose*V*X)inverse * Xtranspose * (Y - Pm)
	 * We only have to do this once, so it is outside the loop.
	 */
	stb_transpose_matrix(X, &Xt);

	/* Perform the Newton iterations to approximate the values of Beta 
	 * Beta = (Xtranspose*V*X)inverse * Xtranspose * (Y - Pm)
	 */
	for (int i = 1; i <= 1024; i++) {
		/* Calculate the predicted probabilities */
		stb_matrix_multiply(X, Beta0, &Pi);
		for (int i = 0; i < Pi->rows; i++) {
			for (int n = 0; n < Pi->columns; n++) {
				Pi->data[i][n] = 1.0 / (1.0 + exp(-1.0 * Pi->data[i][n]));
			}
		}

		/* make the V matrix with Pi * (1- Pi) on the diagonal */
		V = stb_new_matrix(Pi->rows, Pi->rows);
		for (int i = 0; i < V->rows; i++) {
			V->data[i][i] = Pi->data[i][0] * (1 - Pi->data[i][0]);
		}

		/* (Xtranspose*V) */
		stb_matrix_multiply(Xt, V, &Xm);

		/* (Xtranspose*V*X) */
		stb_matrix_multiply(Xm, X, &Xmm);

		/* Calculate: (Xtranspose*V*X)inverse */
		stb_invert_matrix(Xmm, &Xinv);

		/* Calculate: Xtrans * (Y - Pm) -> Ym = Y - Pi = real probability - predicted probability */
		stb_matrix_subs(Y, Pi, &Ym);

		/* Calculate: Xtrans * (Y - Pm)*/
		stb_matrix_multiply(Xt, Ym, &P);

		/* Beta = (Xtranspose*X)inverse * Xtranspose * Y <= Xinv * P */
		stb_matrix_multiply(Xinv, P, &Beta1);

		stb_matrix_add(Beta0, Beta1, &Beta);

		/* Determine the new probabilities */
		stb_matrix_multiply(X, Beta, &Pi);
		for (int i = 0; i < Pi->rows; i++) {
			for (int n = 0; n < Pi->columns; n++) {
				Pi->data[i][n] = 1.0 / (1.0 + exp(-1.0 * Pi->data[i][n]));
			}
		}

		/* Calculate the new log-likelihood 
		 * Pi = X * B
		 * log-likelihood = sum(-1ln(1+exp(Pi))+YiPi) 
		 */
		loglikelihood = 0.0;
		for (int i = 0; i < Y->rows; i++) {
			for (int n = 0; n < Y->columns; n++) {
				loglikelihood += -log(1 + exp(Pi->data[i][n])) + Y->data[i][n] * Pi->data[i][n];
			}
		}

		/* The first iteration, we do not have a old log-likelihood yet */
		if (i == 1)	{
			oldloglikelihood = 2 * loglikelihood;
		}	

		delta = loglikelihood - oldloglikelihood;

		printf("delta: %lf \n", delta);

		/* Needs to be freed before the break statement */
		stb_free_matrix(Pi);
		stb_free_matrix(V);
		stb_free_matrix(Xm);
		stb_free_matrix(Xmm);
		stb_free_matrix(Ym);
		stb_free_matrix(Beta1);

		if (delta < EPSILON) {
			/* Found it */
			break;
		}

		oldloglikelihood = loglikelihood;

		/* Update the Beta and free some matrix that has to be freed after the break statement */
		stb_free_matrix(Beta0);
		Beta0 = stb_dup_matrix(Beta);
		stb_free_matrix(Xinv);
	}	

	if (delta > EPSILON) {
		printf("We did not converge!\n");
	}

	/* Allocate space to return beta */
	double *betaret = NULL;
	betaret = calloc(X->columns, sizeof(double));
	
	/* copy in the beta values */
	for (int i = 0; i < X->columns; i++) {
		betaret[i] = Beta->data[i][0];
	}
	*beta = betaret;

	/* Since we converged (I hope), we will use the "old" values to determine the z-value */

	// /* new Pi*/
	// stb_matrix_multiply(X, Beta, &Pi);

	// /* The new predicted values */
	// for (int i = 0; i < Pi->rows; i++) {
	// 	for (int n = 0; n < Pi->columns; n++) {
	// 		Pi->data[i][n] = 1.0 / (1.0 + exp(-1.0 * Pi->data[i][n]));
	// 	}
	// }	

	// /* Construct the new V matrix with Pi->data[i][0] * (1 - Pi->data[i][0]) on the diagonal */
	// V = stb_new_matrix(Pi->rows, Pi->rows);
	// for (int i = 0; i < V->rows; i++) {
	// 	V->data[i][i] = Pi->data[i][0] * (1 - Pi->data[i][0]);
	// }

	// /* (Xtranspose*V) */
	// stb_matrix_multiply(Xt, V, &Xm);

	// /* (Xtranspose*V*X) */
	// stb_matrix_multiply(Xm, X, &Xmm);

	// /* Calculate the variance-covariance matrix of the estimated regression coefficients: 
	//  * (Xtranspose*V*X)inverse 
	//  */
	// stb_invert_matrix(Xmm, &Xinv);

	STB_MAT *stdev;
	stdev = stb_new_matrix(Beta->rows, Beta->columns);
	for (int i = 0; i < Xinv->rows; i++) {
		for (int n = 0; n < Xinv->columns; n++) {
			if (i == n) {
				/* Store the standard deviations for the estimated values */
				stdev->data[i][0] = sqrt(Xinv->data[i][n]);
			}
		}
	}

	/* Allocate space to return the T-value and the P-value */
	double *zvalret = NULL;
	zvalret = calloc(X->columns, sizeof(double));
	double *pvalret = NULL;
	pvalret = calloc(X->columns, sizeof(double));

	/* Normalize values (e.a. divide by the standard deviation) and lookup the T-value to get the p-value */
	/* if the p-values are all extremely small, it indicates that the predictor are significantly related to the response */
	for (int i = 0; i < Beta->rows; i++) {
		zvalret[i] = Beta->data[i][0] / stdev->data[i][0];
		if (zvalret[i] > 0) {
			pvalret[i] = 2 * (1 - stb_phi(Beta->data[i][0] / stdev->data[i][0]));	
		} else {
			pvalret[i] = 2 * (1 - stb_phi(-Beta->data[i][0] / stdev->data[i][0]));		
		}
	}

	*zvalue = zvalret;
	*pvalue = pvalret;

	/* Not needed anymore */
	stb_free_matrix(Xt);
	stb_free_matrix(Beta0);
	stb_free_matrix(Beta);
	stb_free_matrix(X);
	stb_free_matrix(stdev);
	/* Free here since we still need it to calculate the std.dev */
	stb_free_matrix(Xinv);

	return loglikelihood;
}

double stb_euclidean_distance(const double *a, const double *b, const int size)
{
    double dist = 0.0;

    for (int i = 0; i < size; i++) {
        dist += stb_sqr(a[i] - b[i]);
    }

    return sqrt(dist);
}

double stb_euclidean_distance_sqr(const double *a, const double *b, const int size)
{
    double dist = 0.0;

    for (int i = 0; i < size; i++) {
        dist += stb_sqr(a[i] - b[i]);
    }

    return dist;
}

double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}

// Not that fast ?!
double fast_sigmoid(double value)
{
    //double e = 2.718281828;
    //return 1.0 / (1.0 + pow(e, -x));
    //return atan(M_PI*x/2)*2/M_PI;
    //return atan(x);
    //return 1/(1+exp(-x));
    //return x/sqrt(1+x*x);
    //return erf(sqrt(M_PI)*x/2);
    //return tanh(x);
    //return x / (1.0 + fabs(x));

    double x = fabs(value);
    double x2 = x*x;
    double e = 1.0f + x + x2*0.555f + x2*x2*0.143f;
    return 1.0f / (1.0f + (value > 0 ? 1.0f / e : e));
}

// simple logistic regression
// NOTE: Unlike stb_multi_logistic_regression, to get the intercept (bias), add a column of 1.0s to matrix A
// NOTE: Unlike stb_multi_logistic_regression, Y = 1 or -1 instead of 1 or 0!
void stb_logistic_regression(STB_MAT *A, STB_MAT *Y, double **beta, double **zvalue, double **pvalue)
{
    // the convergence rate
    double delta = 1e-8;
    // the learning rate
    double epsilon = 0.001;
    int max_iters = 1e6;
    int iter = 0;

    // init
    double *beta_old = calloc(A->columns, sizeof(double));
    /* Initialize Beta0 to values close to 0 (0.01) except for the intercept (1.0), this hopefully leads to a
     * lower amount of iterations needed to converge.
     */
    beta_old[0] = 1.0;
    for (size_t k=1; k<A->columns; ++k) {
        beta_old[k] = 0.01;
    }
    double *beta_new = calloc(A->columns, sizeof(double));     

    while (iter <= max_iters) {
        // update each weight
        for (size_t k=0; k<A->columns; ++k) {
            double gradient = 0;
            for (size_t i=0; i<A->rows; ++i) {
                double z_i = 0;
                for (size_t j=0; j<A->columns; ++j) {
                    z_i += beta_old[j] * A->data[i][j];
                }
                gradient += Y->data[i][0] * A->data[i][k] * sigmoid((-Y->data[i][0]) * z_i);
            }            
            beta_new[k] = beta_old[k] + epsilon * gradient;
        }

        double dist = stb_euclidean_distance(beta_new, beta_old, A->columns);
        if (dist < delta && iter > 2) {
            break;
        } else {
        	for (int i = 0; i < A->columns; i++) {
            	beta_old[i] = beta_new[i];
            }
        }

        iter++;
        //printf("%i delta: %f\n", iter, dist);
    }

    if (iter >= max_iters) {
        printf("Did not converge in %i iterations!\n", iter);
    }

    //variance-covariance matrix
    double *var = calloc(A->columns, sizeof(double));
    for (size_t k=0; k<A->columns; ++k) {
        for (size_t i=0; i<A->rows; ++i) {
            double z_i = 0;
            for (size_t j=0; j<A->columns; ++j) {
                z_i += beta_new[j] * A->data[i][j];
            }
            var[k] += A->data[i][k] * A->data[i][k] * sigmoid(z_i) * (1 - sigmoid(z_i));
        }
        var[k] = sqrt(1.0/var[k]);
    }

    double *zvalret = NULL;
    zvalret = calloc(A->columns, sizeof(double));
    double *pvalret = NULL;
    pvalret = calloc(A->columns, sizeof(double));

    /* Normalize values (e.a. divide by the standard deviation) and lookup the T-value to get the p-value */
    /* if the p-values are all extremely small, it indicates that the predictor are significantly related to the response */
    for (int i = 0; i < A->columns; i++) {
        zvalret[i] = beta_new[i] / var[i];
        if (zvalret[i] > 0) {
            pvalret[i] = 2 * (1 - stb_phi(beta_new[i] / var[i]));   
        } else {
            pvalret[i] = 2 * (1 - stb_phi(-beta_new[i] / var[i]));      
        }
    }

    printf("zvalue\n");
    for (size_t k=0; k<A->columns; ++k) {
        printf("%f\n", zvalret[k]);
    }
    printf("pvalue\n");
    for (size_t k=0; k<A->columns; ++k) {
        printf("%f\n", pvalret[k]);
    }

    printf("the best weight:\n");
    for (int i = 0; i < A->columns; i++) {
        printf("%f\t", beta_new[i]);
    }
    printf("\n");
    
    free(beta_old);
    free(var);
    
    *beta = beta_new;
    *zvalue = zvalret;
    *pvalue = pvalret;
}

// L2-regularized logistic regression
// NOTE: Unlike stb_multi_logistic_regression, to get the intercept (bias), add a column of 1.0s to matrix A
// NOTE: Unlike stb_multi_logistic_regression, Y = 1 or -1 instead of 1 or 0!
void stb_logistic_regression_L2(STB_MAT *A, STB_MAT *Y, double **beta, double **zvalue, double **pvalue)
{
    // the convergence rate
    double delta = 1e-8;
    // the learning rate
    double epsilon = 0.001;
    double cost = 1.0;
    int max_iters = 1e6;
    int iter = 0;

    // init
    double *beta_old = calloc(A->columns, sizeof(double));
    /* Initialize Beta0 to values close to 0 (0.01) except for the intercept (1.0), this hopefully leads to a
     * lower amount of iterations needed to converge.
     */
    beta_old[0] = 1.0;
    for (size_t k=1; k<A->columns; ++k) {
        beta_old[k] = 0.01;
    }
    double *beta_new = calloc(A->columns, sizeof(double));     
    
    while (iter <= max_iters) {
        // update each weight
        for (size_t k=0; k<A->columns; ++k) {
            double gradient = 0;
            for (size_t i=0; i<A->rows; ++i) {
                double z_i = 0;
                for (size_t j=0; j<A->columns; ++j) {
                    z_i += beta_old[j] * A->data[i][j];
                }
                gradient += Y->data[i][0] * A->data[i][k] * sigmoid((-Y->data[i][0]) * z_i);
            }            
            beta_new[k] = beta_old[k] + epsilon * gradient - epsilon * cost * beta_old[k];
        }

        double dist = stb_euclidean_distance(beta_new, beta_old, A->columns);
        if (dist < delta && iter > 2) {
            break;
        } else {
            for (int i = 0; i < A->columns; i++) {
                beta_old[i] = beta_new[i];
            }
        }

        iter++;
        //printf("%i delta: %f\n", iter, dist);
    }

    if (iter >= max_iters) {
        printf("Did not converge in %i iterations!\n", iter);
    }

    //variance-covariance matrix
    double *var = calloc(A->columns, sizeof(double));
    for (size_t k=0; k<A->columns; ++k) {
        for (size_t i=0; i<A->rows; ++i) {
            double z_i = 0;
            for (size_t j=0; j<A->columns; ++j) {
                z_i += beta_new[j] * A->data[i][j];
            }
            var[k] += A->data[i][k] * A->data[i][k] * sigmoid(z_i) * (1 - sigmoid(z_i));
        }
        var[k] = sqrt(1.0/var[k]);
    }

    double *zvalret = NULL;
    zvalret = calloc(A->columns, sizeof(double));
    double *pvalret = NULL;
    pvalret = calloc(A->columns, sizeof(double));

    /* Normalize values (e.a. divide by the standard deviation) and lookup the T-value to get the p-value */
    /* if the p-values are all extremely small, it indicates that the predictor are significantly related to the response */
    for (int i = 0; i < A->columns; i++) {
        zvalret[i] = beta_new[i] / var[i];
        if (zvalret[i] > 0) {
            pvalret[i] = 2 * (1 - stb_phi(beta_new[i] / var[i]));   
        } else {
            pvalret[i] = 2 * (1 - stb_phi(-beta_new[i] / var[i]));      
        }
    }

    printf("zvalue\n");
    for (size_t k=0; k<A->columns; ++k) {
        printf("%f\n", zvalret[k]);
    }
    printf("pvalue\n");
    for (size_t k=0; k<A->columns; ++k) {
        printf("%f\n", pvalret[k]);
    }

    printf("the best weight:\n");
    for (int i = 0; i < A->columns; i++) {
        printf("%f\t", beta_new[i]);
    }
    printf("\n");
    
    free(beta_old);
    free(var);

    *beta = beta_new;
    *zvalue = zvalret;
    *pvalue = pvalret;
}

struct Jenks {
        double *cumulValues_value;
        double *cumulValues_weight;
        size_t cumulValues_idx;
        size_t numValues; // number of data points
        size_t numBreaks; // k
        size_t bufferSize; // = number of data points - (numBreaks - 1)
        double *previousSSM; // bufferSize long
        double *currentSSM; // bufferSize long
        int *classBreaks; // bufferSize * (numBreaks - 1) long
        int classBreaksIndex; // 0
        int completedRows; // 0
};

struct Jenks *Init_JenksFisher(const double *values, const double *counts, size_t m, size_t k)
{
	struct Jenks *JenksFisher = NULL;
	JenksFisher = calloc(1, sizeof(struct Jenks));

	JenksFisher->numValues = m;
	JenksFisher->numBreaks = k;
	JenksFisher->bufferSize = JenksFisher->numValues - (JenksFisher->numBreaks - 1);
	JenksFisher->previousSSM = calloc(JenksFisher->bufferSize, sizeof(double));
	JenksFisher->currentSSM = calloc(JenksFisher->bufferSize, sizeof(double));
	JenksFisher->classBreaks = calloc(JenksFisher->bufferSize * (JenksFisher->numBreaks - 1), sizeof(int));
	JenksFisher->classBreaksIndex = 0;
	JenksFisher->completedRows = 0;
	JenksFisher->cumulValues_value = calloc(JenksFisher->numValues, sizeof(double));
	JenksFisher->cumulValues_weight = calloc(JenksFisher->numValues, sizeof(double));
	JenksFisher->cumulValues_idx = 0;

	double w = 0.0, cw = 0.0, cwv = 0.0;
	for (int i = 0; i < JenksFisher->numValues; ++i) {
        w = counts[i];
        cw += w;

        if (cw < w) {
        	fprintf(stderr, "Overflow detected!\n");
        	return NULL;
        }

        cwv += w * values[i];

        JenksFisher->cumulValues_weight[JenksFisher->cumulValues_idx] = cw;
        JenksFisher->cumulValues_value[JenksFisher->cumulValues_idx] = cwv;
        JenksFisher->cumulValues_idx++;

        if (i < JenksFisher->bufferSize) {
            JenksFisher->previousSSM[i] = cwv * cwv / cw; // prepare sum of squared means for first class. Last (k-1) values are omitted
        }
		
	}

	return JenksFisher;
}

/**
 * Gets sum of weighs for elements with index b..e.
 *
 * @param b index of begin element
 * @param e index of end element
 * @return sum of weights.
 */
double getSumOfWeights(struct Jenks *JenksFisher, size_t b, size_t e) {
    //assert (b != 0);    // First element always belongs to class 0, thus queries should never include it.
    assert (b <= e);
    assert (e < JenksFisher->numValues);

    double res = JenksFisher->cumulValues_weight[e];
    res -= JenksFisher->cumulValues_weight[b - 1];

    return res;
}

/**
 * Gets sum of weighed values for elements with index b..e
 *
 * @param b index of begin element
 * @param e index of end element
 * @return the cumul. sum of the values*weight
 */
double getSumOfWeightedValues(struct Jenks *JenksFisher, size_t b, size_t e) {
    //assert (b != 0);    // First element always belongs to class 0, thus queries should never include it.
    assert (b <= e);
    assert (e < JenksFisher->numValues);

    double res = JenksFisher->cumulValues_value[e];
    res -= JenksFisher->cumulValues_value[b - 1];

    return res;
}

/**
 * Gets the Squared Mean for elements within index b..e, multiplied by weight. Note that
 * n*mean^2 = sum^2/n when mean := sum/n
 *
 * @param b index of begin element
 * @param e index of end element
 * @return the sum of squared mean
 */
double getSSM(struct Jenks *JenksFisher, int b, int e) {
    double res = getSumOfWeightedValues(JenksFisher, b, e);
    return res * res / getSumOfWeights(JenksFisher, b, e);
}

/**
 * Finds CB[i+completedRows] given that the result is at least
 * bp+(completedRows-1) and less than ep+(completedRows-1)
 * Complexity: O(ep-bp) <= O(m) @
 *
 * @param i startIndex
 * @param bp endindex
 * @param ep
 *
 * @return the index
 */
int findMaxBreakIndex(struct Jenks *JenksFisher, int i, int bp, int ep) {
    assert (bp < ep);
    assert (bp <= i);
    assert (ep <= i + 1);
    assert (i < JenksFisher->bufferSize);
    assert (ep <= JenksFisher->bufferSize);

    double minSSM = JenksFisher->previousSSM[bp] + getSSM(JenksFisher, bp + JenksFisher->completedRows, i + JenksFisher->completedRows);
    int foundP = bp;
    while (++bp < ep) {
        double currSSM = JenksFisher->previousSSM[bp] + getSSM(JenksFisher, bp + JenksFisher->completedRows, i + JenksFisher->completedRows);
        if (currSSM > minSSM) {
            minSSM = currSSM;

            foundP = bp;
        }
    }
    JenksFisher->currentSSM[i] = minSSM;
    return foundP;
}

/**
 * Find CB[i+completedRows] for all i>=bi and i<ei given that the
 * results are at least bp+(completedRows-1) and less than
 * ep+(completedRows-1)
 * Complexity: O(log(ei-bi)*Max((ei-bi),(ep-bp)))
 * <= O(m*log(m))
 *
 *
 * @param bi
 * @param ei
 * @param bp
 * @param ep
 * @return
 */
void calcRange(struct Jenks *JenksFisher, int bi, int ei, int bp, int ep) {
    assert (bi <= ei);

    assert (ep <= ei);
    assert (bp <= bi);

    if (bi == ei) {
        return;
    }
    assert (bp < ep);

    int mi = (int) floor((double)(bi + ei) / 2.0);
    int mp = findMaxBreakIndex(JenksFisher, mi, bp, fmin((double)ep, (double) mi + 1.0));

    assert (bp <= mp);
    assert (mp < ep);
    assert (mp <= mi);

    // solve first half of the sub-problems with lower 'half' of possible outcomes
    calcRange(JenksFisher, bi, mi, bp, fmin((double)mi, (double)mp + 1.0));

    JenksFisher->classBreaks[JenksFisher->classBreaksIndex + mi] = mp; // store result for the middle element.

    // solve second half of the sub-problems with upper 'half' of possible outcomes
    calcRange(JenksFisher, mi + 1, ei, mp, ep);
}

/**
 * Swaps the content of the two lists with each other.
 */
void swapArrays(struct Jenks *JenksFisher) {
    double temp;
    int i;

    for (i = 0; i < JenksFisher->bufferSize; i++) {
    	temp = JenksFisher->previousSSM[i];
    	JenksFisher->previousSSM[i] = JenksFisher->currentSSM[i];
    	JenksFisher->currentSSM[i] = temp;
    }
}

/**
 * Starting point of calculation of breaks.
 *
 * complexity: O(m*log(m)*k)
 */
void calcAll(struct Jenks *JenksFisher) {
    if (JenksFisher->numBreaks > 1) {
        JenksFisher->classBreaksIndex = 0;
        for (JenksFisher->completedRows = 1; JenksFisher->completedRows < JenksFisher->numBreaks - 1; ++JenksFisher->completedRows) {
            calcRange(JenksFisher, 0, JenksFisher->bufferSize, 0, JenksFisher->bufferSize); // complexity: O(m*log(m))
            
            swapArrays(JenksFisher);
            JenksFisher->classBreaksIndex += JenksFisher->bufferSize;
        }
    }
}

double getSSM2(const double *values, const double *counts, int start, int end)
//double weighted_meanvar(const double *values, const double *counts, int start, int end)
{
	double w_sum = 0, w_sum2 = 0, S = 0;
	double mean = 0, mean_old = 0;

	printf("start: %i end: %i\n", start, end);
	if (start == end) {
		mean = values[start];
		S = 0.0;
		return 0.0;
	}

	int i;
	for (i = start; i <= end; i++) {
		w_sum = w_sum + counts[i];
		w_sum2 = w_sum2 + counts[i] * counts[i];

		mean_old = mean;
		mean = mean_old + (counts[i] / w_sum) * (values[i] - mean_old);
		S = S + counts[i] * (values[i] - mean_old) * (values[i] - mean);
		//printf("value: %f (%f)\n", values[i], counts[i]);

	}
	double population_variance = S / w_sum;
	//double sample_variance = S / (w_sum  - 1);
	//double sample_reliability_variance = S / (w_sum - w_sum2 / w_sum);

	return population_variance;
}

// values is sorted stb_uniqueue values
// counts is weight of the values
// k number of breaks
// m number of stb_uniqueue values (= length array)
// return breaksArray 
double *ClassifyJenksFisherFromValueCountPairs(const double *values, const double *counts, size_t m, size_t k)
{
	double *breaksArray = NULL;
	struct Jenks *JenksFisher; // Helper struct
	size_t nbreaks = k;
	breaksArray = calloc(k, sizeof(double)); // The array holding the final results

	assert(k <= m); // PRECONDITION

	if (!k) {
		fprintf(stderr, "Need to have at least 1 class\n");
		return NULL;
	}

	// initialize struct
	JenksFisher = Init_JenksFisher(values, counts, m, k);

	if (k > 1) {
		calcAll(JenksFisher);

		size_t lastClassBreakIndex = findMaxBreakIndex(JenksFisher, JenksFisher->bufferSize - 1, 0, JenksFisher->bufferSize);

		while (--k) {
			breaksArray[k] = values[lastClassBreakIndex+k];
			assert(lastClassBreakIndex < JenksFisher->bufferSize);
			if (k > 1) {
				JenksFisher->classBreaksIndex -= JenksFisher->bufferSize;
				lastClassBreakIndex = JenksFisher->classBreaks[JenksFisher->classBreaksIndex + lastClassBreakIndex];
			}
		}
		assert(JenksFisher->classBreaksIndex == JenksFisher->classBreaks[0]);
	}

	breaksArray[0] = values[0]; // break for the first class is the minimum of the dataset.

	// the goodness of variance fit (GVF) is calculated. GVF is defined as (SDAM - SDCM) / SDAM. GVF ranges from 0 (worst fit) to 1 (perfect fit).
	int start = 0;
	int end = 0;
	double ssm = 0.0;
	int n , i;
	for (n = 1; n < nbreaks; n++) {
		for (i = start; i < m; ++i) {
			if (values[i] < breaksArray[n]) {
				continue;
			} else {
				end = i - 1;
				break;
			}
		}
		ssm += getSSM2(values, counts, start, end);
		++end;
		start = end;
	}
	end = m - 1;
	ssm += getSSM2(values, counts, start, end);

	double ssa = getSSM2(values, counts, 0, end);
	double gvf = (ssa - ssm)/ssa;
	fprintf(stderr, "Fit: %f\n", gvf);
	return breaksArray;
}

/**
 * Main entry point for creation of Jenks-Fisher natural breaks.
 * Port of Jenks/Fisher breaks originally created in C++ by Maarten Hilferink.
 * @param values array of the values, do not need to be sorted.
 * @param k number of breaks to create
 * @param len length of values array
 * @return Array with breaks
 */
double *stb_jenks(double *ivalues, int len, int k)
{
	// Make a copy of the input values
	double *values = calloc(len, sizeof(double));
	for (int i = 0; i < len; i++) {
		values[i] = ivalues[i];
	}

	double *counts;
	int unique_vals = stb_dunique(values, &counts, len);
    double *breaks_array = NULL;
    
    breaks_array = ClassifyJenksFisherFromValueCountPairs(values, counts, unique_vals, k);
    
    free(values);
    free(counts);

    return breaks_array;
}
// int main(int argc, char *argv[])
// {
// 	STB_MAT *X = stb_matrix_from_file("x.train");
// 	STB_MAT *Y = stb_matrix_from_file("y.train");
// 	double *beta;
// 	double *zvalue;
// 	double *pvalue;
//     // stb_multi_logistic_regression works for -1 - 1 values in Y
//     printf("==========================\n");
// 	stb_logistic_regression(X, Y, &beta, &zvalue, &pvalue);
//     free(beta);
//     free(zvalue);
//     free(pvalue);
//     printf("\n");

//     printf("==========================\n");
//     stb_logistic_regression_L2(X, Y, &beta, &zvalue, &pvalue);
//     free(beta);
//     free(zvalue);
//     free(pvalue);
//     printf("\n");

//     STB_MAT *BIAS = stb_new_matrix(X->rows, 1);
//     for (int i = 0; i < X->rows; i++) {
//         BIAS->data[i][0] = 1.0;
//     }

//     STB_MAT *XBIAS = NULL;
//     stb_join_matrix(BIAS, X, &XBIAS);
//     printf("==========================\n");
//     stb_logistic_regression(XBIAS, Y, &beta, &zvalue, &pvalue);
//     free(beta);
//     free(zvalue);
//     free(pvalue);
//     printf("\n");

//     printf("==========================\n");
//     stb_logistic_regression_L2(XBIAS, Y, &beta, &zvalue, &pvalue);
//     free(beta);
//     free(zvalue);
//     free(pvalue);
//     printf("\n");

//     // stb_multi_logistic_regression works for 0 - 1 values in Y
//     for (int i = 0; i < X->rows; i++) {
//         if (Y->data[i][0] < 0.5) {
//             Y->data[i][0] = 0;
//         }
//     }

//     printf("==========================\n");
//     stb_multi_logistic_regression(X, Y, &beta, &zvalue, &pvalue);     
//     printf("zvalue\n");
//     for (size_t k=0; k<XBIAS->columns; ++k) {
//         printf("%f\n", zvalue[k]);
//     }
//     printf("pvalue\n");
//     for (size_t k=0; k<XBIAS->columns; ++k) {
//         printf("%f\n", pvalue[k]);
//     }

//     printf("the best weight:\n");
//     for (int i = 0; i < XBIAS->columns; i++) {
//         printf("%f\t", beta[i]);
//     }
//     printf("\n");
//     free(beta);
//     free(zvalue);
//     free(pvalue);  

//     stb_free_matrix(X);
//     stb_free_matrix(Y);
//     stb_free_matrix(BIAS);
//     stb_free_matrix(XBIAS);

// 	return 0;
// }

/* K-Means++ data clustering
 * x[n][d] = the data points
 * n       = Number of points
 * d       = Dimension of the data (e.g. color, weight, etc.) 
 * k       = # clusters
*  c[k][d] = Center points of clusters
*  z[n]    = What cluster a point is in
*  wss[k]  = The within-cluster sum of square of each cluster (optional)
*
* Note: This algorithm does not scale very well to a large number of points, consider using k-Means||
*/
void stb_kmeans(double **x, int n, int d, int k, double ***cret, int **zret, double **wssret)
{
    /* Local variables */
    double *work = malloc(n * sizeof(double));
    int *z = calloc(n, sizeof(int));
    *zret = z;
    double **c = (double **) stb_allocmat(k, d, sizeof(double));
    *cret = c;
    int h, i, j, l;
    int current_center, chosen_center, closest_center;
    double u, dist;
    double total_dist, best_dist;
    bool change;

    /* Init the random number generator */
    uint64_t seed;
    seed = stb_spcg32(time(NULL));
    
    /* constants */
    const int ITER = 1000000;
    const double BIG = 1e33f;

    /* Begin. */
    if (k < 1 || k > n) {
        fprintf(stderr, "The value of %i for k is out of bounds and should be between 1 and %i!\n", k, n);
        exit(1);
    }

    /* Clear initial centers */
    for (i = 0; i < n; ++i) {
        work[i] = 1e33f;
    }
    
    /* The inital center is choosen randomly */
    chosen_center = stb_pcg32_bounded(&seed, n);
    for (i = 0; i < d; ++i) {
        c[0][i] = x[chosen_center][i];     /* copy */
    }
    
    /* initialize other centers */
    for (current_center = 1; current_center < k; ++current_center) {        
        total_dist = 0.f;

        /* measure from each point */
        for (i = 0; i < n; ++i) {            
            best_dist = work[i];

            /* Squared Euclidean distance to prior center */
            dist = stb_euclidean_distance_sqr(x[i], c[current_center - 1], d);

            /* shortest squared distance? */
            if (dist < best_dist) {
                best_dist = dist;
            }

            work[i] = best_dist;

            /* cumulative squared distance */
            total_dist += best_dist;
        } /* next data point */

        /* Choose center with probability proportional to its squared distance from existing centers. */
        u = stb_uint32_to_double(stb_pcg32(&seed));
        /* uniform at random over cumulative distance */
        u *= total_dist;                     
        total_dist = 0.f;
        for (i = 0; i < n; ++i) {
            chosen_center = i;
            total_dist += work[i];
            if (total_dist > u) {
                break;
            }
        } /* next i */

        for (j = 0; j < d; ++j) {
            c[current_center][j] = x[chosen_center][j];        /* assign center */
        }
    } /* next center to initialize */

    /* Main loop */
    for (h = 0; h < ITER; ++h) {
        change = false;

        /* Find nearest center for each point */
        for (i = 0; i < n; ++i) {
            current_center = z[i];
            closest_center = 0;
            best_dist = BIG;
            for (l = 0; l < k; ++l) {
                dist = stb_euclidean_distance_sqr(x[i], c[l], d);

                if (dist < best_dist) {
                    best_dist = dist;
                    closest_center = l;
                }
            }

            if (current_center != closest_center) {
                /* reassign point */
                z[i] = closest_center;         
                change = true;
            }
        }

        /* There are no more changes: Success!!! */
        if (!change) {
            break;
        }

        /* Find cluster centers */
        
        /* zero population */
        for (l = 0; l < k; ++l) {
            work[l] = 0.f;   
        }

         /* zero centers */
        for (j = 0; j < k; ++j) {       
            for (l = 0; l < d; ++l) {
                c[j][l] = 0.f;
            }
        }

        for (i = 0; i < n; ++i) {
            current_center = z[i];
            work[current_center] += 1.f;             /* count */
            for (j = 0; j < d; ++j) {
                c[current_center][j] += x[i][j];    /* add */
            }
        }
        
        for (current_center = 0; current_center < k; ++current_center) {
            /* empty cluster?*/
            if (work[current_center] < .5f) {
                fprintf(stderr, "Something strange happened. We have an empty cluster\n");
                exit(1);
            }

            u = 1.f / work[current_center];
            for (j = 0; j < d; ++j) {
                c[current_center][j] *= u; /* multiplication is faster than division*/
            }
        }
    } /* next h */

    /* Optionally compute the within-cluster sum of squares for each cluster. */
    if (wssret) {
        double *wss = calloc(k, sizeof(double));
        *wssret = wss;

        for (i = 0; i < n; ++i) {
            current_center = z[i];
            for (j = 0; j < d; ++j) {
                wss[current_center] += stb_sqr(x[i][j] - c[current_center][j]);;
            }
        }
    }

    if (h == ITER) {
        fprintf(stderr, "Not enough iterations! Did not converge\n");
        exit(1);
    }
}

/*----------------------------------------------------------------------
  TRED2
  DATE WRITTEN   760101   (YYMMDD)
  REVISION DATE  830518   (YYMMDD)
  CATEGORY NO.  D4C1B1
  AUTHOR  Smith, B. T., et al.
  PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal
           matrix using and accumulating orthogonal transformation
  DESCRIPTION

    This subroutine is a translation of the ALGOL procedure TRED2,
    NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
    Handbook for Auto. Comp., Vol.II-Linear Algebra, 212-226(1971).

    This subroutine reduces a REAL SYMMETRIC matrix to a
    symmetric tridiagonal matrix using and accumulating
    orthogonal similarity transformations.

    On Input

       N is the order of the matrix.

       A contains the real symmetric input matrix.  Only the
         lower triangle of the matrix need be supplied.

    On Output

       D contains the diagonal elements of the tridiagonal matrix.

       E contains the subdiagonal elements of the tridiagonal
         matrix in its last N-1 positions.  E(0) is set to zero.

       Z contains the orthogonal transformation matrix
         produced in the reduction.

       A and Z may coincide.  If distinct, A is unaltered.

    Questions and comments should be directed to B. S. Garbow,
    Applied Mathematics Division, Argonne National Laboratory
    ------------------------------------------------------------------
  REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
                Y. Ikebe, V. C. Klema, C. B. Moler, *Matrix Eigen-
                system Routines - EISPACK Guide*, Springer-Verlag,
                1976.
  ROUTINES CALLED  (NONE)
----------------------------------------------------------------------- */
void tred2 (double **a, int n, double *d, double *e, double **z)
{
/* Local variables */
    double f, g, h;
    int i, j, k, l;
    double hh;
    int ii, jp1;
    double scale;
    double tmp;
/* Function Body */
    for (i = 0; i < n; ++i) {
        for (j = 0; j <= i; ++j) {
            z[i][j] = a[i][j];
        }
    }
    if (1 == n) goto L320;
/* .......... FOR I=N STEP -1 UNTIL 2 DO -- .......... */
    for (ii = 2; ii <= n; ++ii) {
        i = n + 2 - ii;
        l = i - 1;
        h = 0.f;
        scale = 0.f;
        if (l < 2) goto L130;
/* .......... Scale row (ALGOL TOL then not needed) .......... */
        for (k = 1; k <= l; ++k) scale += fabsf(z[i-1][k-1]);
        if (scale != 0.f)  goto L140;
L130:
        e[i-1] = z[i-1][l-1];
        goto L290;
L140:
        for (k = 1; k <= l; ++k) {
            z[i-1][k-1] /= scale;
            tmp = z[i-1][k-1];
            h += tmp * tmp;
        }
        f = z[i-1][l-1];
        g = -copysignf(sqrtf(h), f);
        e[i-1] = scale * g;
        h -= f * g;
        z[i-1][l-1] = f - g;
        f = 0.f;
        for (j = 1; j <= l; ++j) {
            z[j-1][i-1] = z[i-1][j-1] / h;
            g = 0.f;
/* .......... Form element of A*U .......... */
            for (k = 1; k <= j; ++k) g += z[j-1][k-1] * z[i-1][k-1];
            jp1 = j + 1;
            if (l < jp1) goto L220;
            for (k = jp1; k <= l; ++k) g += z[k-1][j-1] * z[i-1][k-1];
/* .......... Form element of P .......... */
L220:
            e[j-1] = g / h;
            f += e[j-1] * z[i-1][j-1];
        }
        hh = f / (h + h);
/* .......... Form reduced A .......... */
        for (j = 1; j <= l; ++j) {
            f = z[i-1][j-1];
            g = e[j-1] - hh * f;
            e[j-1] = g;
            for (k = 1; k <= j; ++k)
                z[j-1][k-1] = z[j-1][k-1] - f * e[k-1] - g * z[i-1][k-1];
        }
L290:
        d[i-1] = h;
    }
L320:
    d[0] = 0.f;
    e[0] = 0.f;
/* .......... Accumulation of transformation matrices .......... */
    for (i = 1; i <= n; ++i) {
        l = i - 1;
        if (0.f == d[i-1]) goto L380;
        for (j = 1; j <= l; ++j) {
            g = 0.f;
            for (k = 1; k <= l; ++k) g += z[i-1][k-1] * z[k-1][j-1];
            for (k = 1; k <= l; ++k) z[k-1][j-1] -= g * z[k-1][i-1];
        }
L380:
        d[i-1] = z[i-1][i-1];
        z[i-1][i-1] = 1.f;
        if (l < 1) continue;
        for (j = 1; j <= l; ++j) {
            z[i-1][j-1] = 0.f;
            z[j-1][i-1] = 0.f;
        }
    }
    return;
} /* end of tred2 */


/*----------------------------------------------------------------------
  TQL2
  DATE WRITTEN   760101   (YYMMDD)
  REVISION DATE  830518   (YYMMDD)
  CATEGORY NO.  D4A5,D4C2A
  AUTHOR  Smith, B. T., et al.
  PURPOSE  Compute eigenvalues and eigenvectors of symmetric
           tridiagonal matrix.
  DESCRIPTION

    This subroutine is a translation of the ALGOL procedure TQL2,
    Num. Math. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
    Wilkinson.
    Handbook for Auto. Comp., Vol.II-Linear Algebra, 227-240(1971).

    This subroutine finds the eigenvalues and eigenvectors
    of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
    The eigenvectors of a FULL SYMMETRIC matrix can also
    be found if  TRED2  has been used to reduce this
    full matrix to tridiagonal form.

    On Input

       N is the order of the matrix.

       D contains the diagonal elements of the input matrix.

       E contains the subdiagonal elements of the input matrix
         in its last N-1 positions.  E(0) is arbitrary.

       Z contains the transformation matrix produced in the
         reduction by  TRED2, if performed.  If the eigenvectors
         of the tridiagonal matrix are desired, Z must contain
         the identity matrix.

     On Output

       D contains the eigenvalues in ascending order.  If an
         error exit is made, the eigenvalues are correct but
         unordered for indices 1,2,...,IERR-1.

       E has been destroyed.

       Z contains orthonormal eigenvectors of the symmetric
         tridiagonal (or full) matrix.  If an error exit is made,
         Z contains the eigenvectors associated with the stored
         eigenvalues.

       Return value is set to
         Zero       for normal return,
         J          if the J-th eigenvalue has not been
                    determined after 30 iterations.

    Questions and comments should be directed to B. S. Garbow,
    Applied Mathematics Division, Argonne National Laboratory
    ------------------------------------------------------------------
  REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
                Y. Ikebe, V. C. Klema, C. B. Moler, *Matrix Eigen-
                system Routines - EISPACK Guide*, Springer-Verlag,
                1976.
----------------------------------------------------------------------*/
int tql2 (int n, double *d, double *e, double **z)
{
/* Local variables */
    double b, c, f, g, h;
    int i, j, k, l, m;
    double p, r, s, c2, c3;
    int l1, l2;
    double s2;
    int ii;
    double dl1, el1;
    int mml;
    double one;
/* Function Body */
    one = 1.f;
    if (1 == n) return 0;
    for (i = 2; i <= n; ++i) e[i-2] = e[i-1];
    f = 0.f;
    b = 0.f;
    e[n-1] = 0.f;
    s2 = 0.f;
    c3 = 0.f;
    for (l = 1; l <= n; ++l) {
        j = 0;
        h = fabsf(d[l-1]) + fabsf(e[l-1]);
        if (b < h) b = h;
/* .......... Look for small sub-diagonal element .......... */
        for (m = l; m <= n; ++m) {
            if (b + fabsf(e[m-1]) == b) break;
/* .......... E(N) is always zero, so there is no exit */
/*            through the bottom of the loop .......... */
/* added 2011/06/13 AJA */
            if (m > n) return -1;
        }
        if (m == l)  goto L220;
L130:
/* .......... Set error -- No convergence to an */
/*            eigenvalue after 30 iterations .......... */
        if (30 == j) return l;
        ++j;
/* .......... Form shift .......... */
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l-1];
        p = (d[l1-1] - g) / (e[l-1] * 2.f);
        r = hypotf(p, one);
        d[l-1] = e[l-1] / (p + copysignf(r, p));
        d[l1-1] = e[l-1] * (p + copysignf(r, p));
        dl1 = d[l1-1];
        h = g - d[l-1];
        if (l2 > n) goto L145;
        for (i = l2; i <= n; ++i) d[i-1] -= h;
L145:
        f += h;
/* .......... QL transformation .......... */
        p = d[m-1];
        c = 1.f;
        c2 = c;
        el1 = e[l1-1];
        s = 0.f;
        mml = m - l;
/* .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... */
        for (ii = 1; ii <= mml; ++ii) {
            c3 = c2;
            c2 = c;
            s2 = s;
            i = m - ii;
            g = c * e[i-1];
            h = c * p;
            if (fabsf(p) < fabsf(e[i-1])) goto L150;
            c = e[i-1] / p;
            r = sqrtf(c * c + 1.f);
            e[i] = s * p * r;
            s = c / r;
            c = 1.f / r;
            goto L160;
L150:
            c = p / e[i-1];
            r = sqrtf(c * c + 1.f);
            e[i] = s * e[i-1] * r;
            s = 1.f / r;
            c *= s;
L160:
            p = c * d[i-1] - s * g;
            d[i] = h + s * (c * g + s * d[i-1]);
/* .......... Form vector .......... */
            for (k = 1; k <= n; ++k) {
                h = z[k-1][i];
                z[k-1][i] = s * z[k-1][i-1] + c * h;
                z[k-1][i-1] = c * z[k-1][i-1] - s * h;
            }
        }
        p = -s * s2 * c3 * el1 * e[l-1] / dl1;
        e[l-1] = s * p;
        d[l-1] = c * p;
        if (b + fabsf(e[l-1]) > b) goto L130;
L220:
        d[l-1] += f;
    }
/* .......... Order eigenvalues and eigenvectors .......... */
    for (ii = 2; ii <= n; ++ii) {
        i = ii - 1;
        k = i;
        p = d[i-1];
        for (j = ii; j <= n; ++j) {
            if (d[j-1] >= p) continue;
            k = j;
            p = d[j-1];
        }
        if (k == i) continue;
        d[k-1] = d[i-1];
        d[i-1] = p;
        for (j = 1; j <= n; ++j) {
            p = z[j-1][i-1];
            z[j-1][i-1] = z[j-1][k-1];
            z[j-1][k-1] = p;
        }
    }
    return 0;
} /* end of tql2 */

/* Compute eigenvalues and eigenvectors of a symmetric matrix
 * a[n][n] = The matrix
 * n       = Order of a
 * w[n]    = Eigenvalues
 * z[n][n] = Eigenvectors
 */
int stb_eigenv (double **a, int n, double **wret, double ***zret)
{
/* Local variables */
    int ierr;
/* Function Body */
    
    /* contains the subdiagonal elements of the tridiagonal
     * matrix in its last N-1 positions 
     */
    double *t = calloc(n, sizeof(double));
    double *w = calloc(n, sizeof(double));
    *wret = w;
    double **z = (double **) stb_allocmat(n, n, sizeof(double));
    *zret = z;
    /* Reduce real symmetric matrix to symmetric tridiagonal
     * matrix using and accumulating orthogonal transformation
     */    
    tred2 (a, n, w, t, z);

    ierr = tql2(n, w, t, z);

    free(t);
    return ierr;
} /* end of seigv */

/* This function returns the median */
double stb_median(double *data, int n)
{
    double median;

    if (n < 2) {
        fprintf(stderr, "Need at least 2 elements\n");
        exit(1);
    }
    
    if (n == 2) {
        return ((data[0] + data[1]) / 2.0);
    }

    /* Copy the data */
    double *sorted_data = (double *) calloc(n, sizeof(double));
    for(int i = 0; i < n; i++) {
        sorted_data[i] = data[i];
    }

    /* Sort it */
// use the default swap function
#undef STB_QSORT_SWAP
#define STB_QSORT_SWAP(A,B,SIZE) STB_DEFAULT_SWAP(A, B, SIZE)
    stb_qsort(sorted_data, n, sizeof(double), dcmp);

    if (IS_ODD(n)) {
        int split = (int) floor((double) n / 2.0);
        median = sorted_data[split];
        free(sorted_data);
        return median;
    } 


    /* If there are an even number of data points in the original ordered data set, split this data set exactly in half.
     * The lower quartile value is the median of the lower half of the data. The upper quartile value is the median of the upper half of the data.
     */
    int split = n / 2;
    median = (sorted_data[split - 1] + sorted_data[split]) / 2.0;
    free(sorted_data);
    return median;
}

/* Check whether division will overflow */
int stb_safdiv (double numer, double denom)
{
    int retval;

    retval = fabs(denom) > fmax(0.f, fabs(numer / DBL_MAX));

    return retval;
} /* end of safdiv */

/* stb_pca Principal Component Analysis
 * x[n][p]    = Data matrix
 * nx[n][p]   = 1 if data exists this point, 0 otherwise (optional (no missing data), can be NULL)
 * n          = Number of objects
 * p          = Variables each object
 * weights[p] = Weight of each variable (optional, can be NULL)
 * eret[p]    = Eigenvalues of covariance matrix
 * vret[p]    = Eigenvectors of covariance matrix
 * rret[n][m] = Projected data
 * m          = # of dimensions to project
 * level      = Level of robustness:
 *      -1 => flimsy statistics, Chebyshev codeviation
 *       0 => regular statistics, covariance matrix
 *       1 => semi-robust statistics, Manahattan codeviation
 *       2 => robust statistics, comedian matrix
 */
void stb_pca(double **x, int n, int p, int **nx, double *weights, int m, int level, double **eret, double ***vret, double ***rret)
{
/* Local variables */
    int h, i, j, k, ifault;
    double xj, xk;
    double xdif, vall, rtot, vtot, xsum;

    /* Averages each variable */
    double *avg = calloc(p, sizeof(double));
    /* arrays to hold intermediate values */
    double *t = calloc(p, sizeof(double));
    double *u = calloc(n, sizeof(double));

    /* Eigenvalues and Eigenvectors of covariance matrix */
    double *eval;
    double **evec;
    
    /* Eigenvectors of covariance matrix */
    double **var = (double **) stb_allocmat(p, p, sizeof(double));

    /* Projected data */
    double **r = (double **) stb_allocmat(n, m, sizeof(double));
    *rret = r;

/* Function Body */
    if (m > p) {
        fprintf(stderr, "Cannot reduce %i dimensions to %i\n", p, m);
        exit(1);
    }

    if (level < -1 || level > 2) {
        fprintf(stderr, "Invalid level %i! Should be between [-1,2] \n", level);
        exit(1);
    }

    /* Calculate empirical mean along each column */
    for (j = 0; j < p; ++j) {
        h = 0;
        for (i = 0; i < n; ++i) {
            if (nx) {
                if (nx[i][j] > 0) {
                    u[h++] = x[i][j];
                }
            } else {
                u[h++] = x[i][j];
            }
        }
        if (h > 0) {
            if (0 == level) {
                rtot = stb_sum(u, h);
                avg[j] = rtot / (double) h;
            } else if (-1 == level) {
                xj = u[0];
                xk = u[0];
                for (i = 1; i < h; ++i) {
                    xj = fmin(xj, u[i]);
                    xk = fmax(xk, u[i]);
                }
                avg[j] = (xj + xk) * .5f;
            } else {
                avg[j] = stb_median(u, h);
            }
        /* No values */
        } else {
            avg[j] = 0.f;
        }
    }

    /*----------------------------------------------------------------------
    Overall codeviation / covariance matrix.
    A criterion for the orthogonality of two vectors u and v is that
    the distance from u to v should equal the distance from u to -v.

    Applying the Linfinity metric, this means that
      max|u - v| = max|u - (-v)|
      (1/2) (max|u + v| - max|u - v|) = 0

    Applying the L2 metric, this means that
      sum(u - v)^2 = sum(u - (-v))^2
      sum(uv) = 0

    Applying the L1 metric, this means that
      sum|u - v| = sum|u - (-v)|
      (1/2) sum(|u + v| - |u - v|) = 0

    which is identical to the expression
      max(min(u, v), min(-u, -v))

    Thus, the L2 orthogonality criterion yields the formula for covariance,
    and the L1 orthogonality criterion yields the formula for codeviation.
    For more information on applications of L1 orthogonality, request the
    technical report "A direct method for L1-norm codeviation" [2013] from
      andy@13olive.net

    The comedian matrix is a highly robust estimate of covariance
    that is easy to compute.  Refer to:
         Michael Falk
         "On MAD and Comedians"
         Ann. Inst. Statist. Math., v.49 n.4 pp.615-644, 1997

    The Linfinity codeviation matrix generalizes the range.
    The covariance matrix generalizes the variance.
    The L1 codeviation matrix generalizes the mean absolute deviation.
    The comedian matrix generalizes the median absolute deviation.

    It should be acknowledged that there are many ways to handle the missing
    data problem, many robust variants of the covariance matrix have been
    proposed, and there are many ways to approximately factor a matrix.
    ----------------------------------------------------------------------*/
    /* Calculate the deviations from the mean */
    for (k = 0; k < p; ++k) {
        for (j = k; j < p; ++j) {
            xsum = 0.f;
            xdif = 0.f;
            h = 0;
            for (i = 0; i < n; ++i) {
                if (nx) {
                    if (nx[i][j] > 0 && nx[i][k] > 0) {
                        if (weights) {
                            xj = (x[i][j] - avg[j]) * weights[j];
                            xk = (x[i][k] - avg[k]) * weights[k];
                        } else {
                            xj = (x[i][j] - avg[j]);
                            xk = (x[i][k] - avg[k]);                            
                        }
                        if (-1 == level) {
                            xsum = fmax(xsum, fabsf(xj + xk));
                            xdif = fmax(xdif, fabsf(xj - xk));
                        } else if (0 == level || 2 == level) {
                            u[h] = xj * xk;
                        } else {
                            u[h] = fmax(fmin(xj, xk), fmin(-xj, -xk));
                        }
                        ++h;
                    }
                } else {
                    if (weights) {
                        xj = (x[i][j] - avg[j]) * weights[j];
                        xk = (x[i][k] - avg[k]) * weights[k];
                    } else {
                        xj = (x[i][j] - avg[j]);
                        xk = (x[i][k] - avg[k]);                        
                    }
                    if (-1 == level) {
                        xsum = fmax(xsum, fabsf(xj + xk));
                        xdif = fmax(xdif, fabsf(xj - xk));
                    } else if (0 == level || 2 == level) {
                        u[h] = xj * xk;
                    } else {
                        u[h] = fmax(fmin(xj, xk), fmin(-xj, -xk));
                    }
                    ++h;                    
                }
            }
            if (-1 == level) {
                var[j][k] = (xsum - xdif) * .5f;
            } else if (0 == level || 1 == level) {
                if (h > 1) {
                    rtot = stb_sum(u, h);
                    var[j][k] = rtot / (double) (h - 1);
                } else {
                    var[j][k] = 0.f;
                }
            } else {
                if (h > 0) {
                    var[j][k] = stb_median(u, h);
                } else {
                    var[j][k] = 0.f;
                }
            }
        }
    }

/*----------------------------------------------------------------------
         Factor covariance matrix
----------------------------------------------------------------------*/
    ifault = stb_eigenv(var, p, &eval, &evec);
    if (ifault != 0) {
        fprintf(stderr, "Error caluclating Eigenvalues and Eigenvectors\n");
        exit(1);
    }
    /* Eigenvalues and Eigenvectors of covariance matrix */
    *eret = eval;
    *vret = evec;

/*----------------------------------------------------------------------
The singular value decomposition is X = U S V^T
The symmetric eigenproblem is X^T X = V L V^T
Step through columns of V in reverse order so that the first column
of R has the greatest variance.
----------------------------------------------------------------------*/
/* R[i,k] = X[j,i] V[j,p-k+1] */
    for (k = 0; k < m; ++k) {
        for (j = 0; j < p; ++j) {
            if (weights) {
                t[j] = evec[j][p-k-1] * weights[j];
            } else {
                t[j] = evec[j][p-k-1];
            }
        }
        vall = stb_sum(t, p);
        for (i = 0; i < n; ++i) {
            rtot = 0.f;
            vtot = 0.f;
            for (j = 0; j < p; ++j) {
                if (nx) {
                    if (nx[i][j] > 0) {
                        rtot += (x[i][j] - avg[j]) * t[j];
                        vtot += t[j];
                    }
                } else {
                    rtot += (x[i][j] - avg[j]) * t[j];
                    vtot += t[j];                    
                }
            }
            r[i][k] = stb_safdiv(vall, vtot) ? rtot * (vall / vtot) : 0.f;
        }
    }

    /* Clean */
    free(var);
    free(avg);
    free(t);
    free(u);
} /* end of pca */

#define CUSTOM_SWAP(A,B,SIZE)                        \
    do {                                             \
        char *a_byte = A;                            \
        char *b_byte = B;                            \
                                                     \
        int idxA = ((double *)a_byte - dist);        \
        int idxB = ((double *)b_byte - dist);        \
                                                     \
        double dtemp;                                \
        dtemp = dist[idxA];                          \
        dist[idxA] = dist[idxB];                     \
        dist[idxB] = dtemp;                          \
                                                     \
        int itemp;                                   \
        itemp = index[idxA];                         \
        index[idxA] = index[idxB];                   \
        index[idxB] = itemp;                         \
                                                     \
    } while (0)

/* Arrange the N elements of data in random order. */
void stb_shuffle(void *data, size_t n, size_t size, uint64_t *seed) {
    char tmp[size];
    char *arr = data;
    size_t stride = size * sizeof(char);

    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; ++i) {
            size_t j = i + stb_pcg32_bounded(seed, n - i + 1);

            memcpy(tmp, arr + j * stride, size);
            memcpy(arr + j * stride, arr + i * stride, size);
            memcpy(arr + i * stride, tmp, size);
        }
    }
}

/* returns a with n non-repeating random integers in [0,high). This 
 * is useful when n is small, but the range of numbers is big. Otherwise
 * it might be best to use stb_shuffle
 */
int *stb_unique_random(int n, int high, uint64_t *seed) {
    int i, j, duplicate;

    int *a = calloc(n, sizeof(int));

    for (i = 0; i < n; i++) {
        do {
            duplicate = 0;
            a[i] = stb_pcg32_bounded(seed, high);
            for (j = i - 1; j >= 0; j--) {
                if (a[j] == a[i]) {
                    duplicate = 1;
                    break;
                }
            }
        } while (duplicate);
    }

    return a;
}

/* stb_neugas, a data neural gas clustering algorithm
 * See: www.demogng.de/JavaPaper/node16.html
 *
 * x[n][p] = Data Matrix
 * n       = Number of objects
 * p       = Measurements per object
 * k       = Number of clusters
 * c[k][p] = Cluster centers
 * z[n]    = What cluster a point is in (optional)
 * wss[k]  = The within-cluster sum of square of each cluster (optional only possible in combination with z!)
 *
 * Note: Neural gas was developed with a focus on learning a representation of the data space, rather than partitioning a data set
 */
void stb_neugas(double **x, int n, int p, int k, double ***cret, int **zret, double **wssret)
{
/* constants */
    double lami = 10.f; /* range of adjustment, initial */
    double lamf = 0.01f; /* range of adjustment, final */
    double epsi = 0.3f; /* adjustment, initial */
    double epsf = 0.005f; /* adjustment, final */
    int tau = 40000; /* iterations */

    if (k > n) {
        fprintf(stderr, "Cannot make %i clusters with only %i data points\n", k, n);
        exit(1);
    }

    double **c = (double **) stb_allocmat(k, p, sizeof(double));
    *cret = c;

    /* Distances to each center */
    double *dist = calloc(k, sizeof(double));
    /* Index array */
    int *index = calloc(k, sizeof(int));

    /* Local variables */
    int h, i, j, l;
    double r, y;
    double eps, efact, lfact, lambda;
    
    /* Function Body */
    lambda = lami;
    eps = epsi;
    lfact = pow( lamf / lami, 1.f / (double) (tau - 1) );
    efact = pow( epsf / epsi, 1.f / (double) (tau - 1) );

    /* Init the random number generator */
    uint64_t seed;
    seed = stb_spcg32(time(NULL));

    /* The inital centers are choosen randomly but unique*/
    int *chosen_centers = stb_unique_random(k, n, &seed);
    for (j = 0; j < k; ++j) {
        for (i = 0; i < p; ++i) {
            c[j][i] = x[chosen_centers[j]][i];     /* copy */
        }
    }
    free(chosen_centers);

    /* main loop */
    for (h = 0; h < tau; ++h) {
        /* choose a datum at random */
        r = stb_uint32_to_double(stb_pcg32(&seed));;
        i = (int) (r * n);
        if (i > (n - 1)) {
            i = n - 1;
        }
        /* find distance to each center */
        for (l = 0; l < k; ++l) {
            dist[l] = stb_euclidean_distance_sqr(x[i], c[l], p);
            index[l] = l;
        }
        /* sort distances */

// use the swap function to co-sort a dependent array
#undef STB_QSORT_SWAP
#define STB_QSORT_SWAP(A,B,SIZE) CUSTOM_SWAP(A, B, SIZE)
        stb_qsort(dist, k, sizeof(double), dcmp);

        for (l = 0; l < k; ++l) { 
            dist[(index[l] - 1)] = (double) l; /* rank */
        }
        /* adapt the cluster centers */
        for (l = 0; l < k; ++l) {
            y = eps * exp(-dist[l] / lambda);
            for (j = 0; j < p; ++j) {
                c[l][j] += y * (x[i][j] - c[l][j]);
            }
        }
        /* adjust parameters */
        lambda *= lfact;
        eps *= efact;
    }


    /* Optionally find the nearest center for each point */
    if (zret) {
        int current_center, closest_center;
        double curr_dist, best_dist;
        int *z = calloc(n, sizeof(int));
        *zret = z;
        const double BIG = 1e33f;

        /* Find nearest center for each point */
        for (i = 0; i < n; ++i) {
            current_center = z[i];
            closest_center = 0;
            best_dist = BIG;
            for (l = 0; l < k; ++l) {
                curr_dist = stb_euclidean_distance_sqr(x[i], c[l], p);

                if (curr_dist < best_dist) {
                    best_dist = curr_dist;
                    closest_center = l;
                }
            }

            if (current_center != closest_center) {
                /* reassign point */
                z[i] = closest_center;         
            }
        }

        /* Optionally compute the within-cluster sum of squares for each cluster. */
        if (wssret) {
            double *wss = calloc(k, sizeof(double));
            *wssret = wss;

            for (i = 0; i < n; ++i) {
                current_center = z[i];
                for (j = 0; j < p; ++j) {
                    wss[current_center] += stb_sqr(x[i][j] - c[current_center][j]);;
                }
            }
        }
    }
    return;
} /* end of neugas */

// Murmer One At A Time 32 bit
uint32_t inline stb_murmer32(const char *key)
{
    uint32_t h = 3323198485ul;
    for (;*key;++key) {
        h ^= *key;
        h *= 0x5bd1e995;
        h ^= h >> 15;
  }

  return h;
}

// Murmer One At A Time 64 bit
uint64_t inline stb_murmer64(const char *key)
{
    uint64_t h = 525201411107845655ull;
    for (;*key;++key) {
        h ^= *key;
        h *= 0x5bd1e9955bd1e995;
        h ^= h >> 47;
    }

    return h;
}

/* Calculate the Shannon Index (a.k.a. Shannon's diversity index, the Shannon–Wiener index, 
 * the Shannon–Weaver index, the Shannon entropy). Pilou evenness compares the actual diversity
 * value (such as the Shannon Index, H′) to the maximum possible diversity value (when all species
 * are equally common, Hmax=ln S where S is the total number of species). 
 * The higher the Shannon index, the more diverse the species are in the habitat (alpha diversity).
 * Evenness gives you a value between 0 and 1. The lower the evenness, the higher the diversity.
 */
void stb_shannon(double *data, size_t n, double *index, double *evenness)
{
    double sum = stb_sum(data, n);
    double H = 0.0, value;
    
    int i;
    for (i = 0; i < n; i++) {
        value = data[i]/sum;
        H += value * log(value);
    }
    H *= -1;

    *index = H;

    *evenness = H/log((double) n);
}

/* Simpson's Diversity Index is a measure of diversity which takes into account the number of species
 * present, as well as the relative abundance of each species. As species richness and evenness 
 * increase, so diversity increases. The value ranges between 0 and 1. One represents infinite diversity
 * and 0, no diversity.
 */
void stb_simpson(double *data, size_t n, double *index)
{
    double value = 0.0;
    double sum = 0.0;

    int i;
    for (i = 0; i < n; i++) {
        value += data[i] * (data[i] - 1);
        sum += data[i];
    }

    *index = 1 - value / (sum * (sum - 1));
}

/* The Jaccard similarity index (sometimes called the Jaccard similarity coefficient) compares members 
 * for two sets to see which members are shared and which are distinct. It’s a measure of similarity for 
 * the two sets of data, with a range from 0% to 100%. The higher the percentage, the more similar the 
 * two populations. Although it’s easy to interpret, it is extremely sensitive to small samples sizes and 
 * may give erroneous results, especially with very small samples or data sets with missing observations.
 */

int stb_strcmp(char * l, char *r)
{
    if (!l || !r)
        return 1;

    return (strcmp(l, r));
}

stb_create_htable(stb_JHT_, char *, double, stb_murmer64(it), stb_strcmp(l, r) == 0)

double stb_jaccard(char **setA, size_t n_setA, char **setB, size_t n_setB)
{
    stb_JHT_t JHT = stb_JHT_new();
    int i;
    for (i = 0; i < n_setA; i++) {
        stb_JHT_put(JHT, setA[i], 1.0);
    }

    double ret;
    double shared = 0.0;
    for (i = 0; i < n_setB; i++) {
        ret = 0.0;
        stb_JHT_get(JHT, setB[i], &ret);
        if (ret) {
            shared++;
        }
    }

    stb_JHT_free(JHT);

    return shared / (double) (n_setA + n_setB - shared);
}

/* The Bray–Curtis dissimilarity is a measure used to quantify the compositional dissimilarity between two 
 * different sites, based on counts (c_setA and c_setB) at each site. The Bray–Curtis dissimilarity is bounded between 0 and 1, 
 * where 0 means the two sites have the same composition and 1 means the two sites do not share any species.
 */
double stb_bray_curtis(char **setA, double *c_setA, size_t n_setA, char **setB, double *c_setB, size_t n_setB)
{
    stb_JHT_t JHT = stb_JHT_new();
    int i;
    double sumA = 0.0, sumB = 0.0;

    for (i = 0; i < n_setA; i++) {
        // Add species and counts to hash table
        stb_JHT_put(JHT, setA[i], c_setA[i]);
        sumA += c_setA[i];
    }

    double ret;
    double lesser_sum = 0.0;
    for (i = 0; i < n_setB; i++) {
        ret = 0.0;
        stb_JHT_get(JHT, setB[i], &ret);
        sumB += c_setB[i];
        if (ret) {
            // Use lesser sum
            if (ret < c_setB[i]) {
                lesser_sum += ret;
            } else {
                lesser_sum += c_setB[i];
            }
        }
    }

    stb_JHT_free(JHT);

    return 1 - (2 * lesser_sum) / (sumA + sumB);
}

#endif //STB_STATS_DEFINE
#endif //STB__STATS__H