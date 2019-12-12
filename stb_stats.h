/* stb_stats.h - v1.09 - Statistics Tool Box -- public domain
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
 Manas Sharma (Polynomial fitting, Public domain)

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

/* Compute (numerical stable) sum of the data using the Neumaier summation algorithm */
STB_EXTERN double stb_sum(double *data, int n);

/* Returns Spearman's Rank correlation of two vectors x and y, each of size n */
STB_EXTERN double stb_spearman(double *x, double *y, int n);

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

/**
 * Main entry point for creation of Jenks-Fisher natural breaks.
 * Port of Jenks/Fisher breaks originally created in C++ by Maarten Hilferink.
 * @param values array of the values, do not need to be sorted.
 * @param k number of breaks to create
 * @param len length of values array
 * @return Array with breaks
 */
STB_EXTERN double *stb_jenks(double *values, int len, int k);

#ifdef STB_STATS_DEFINE

/* The incomplete gamma function (Thanks Jacob Wells) */
static double stb_incgamma(double S, double Z)
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

/* Given  an  array  of data (length n), this routine  returns  its  mean and standard  deviation standard_dev calculated using the Welford’s method. */
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
	*sample_variance = sq_sum / (n - 1);
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
	if (n < 4) {
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
		*median = (sorted_data[split - 1] + data[split]) / 2;

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

double stb_factorial(int n)
{
	double factorial = 1;

	while (n > 1) {
		factorial = factorial * n;
		n = n - 1;
	}

	return factorial;
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

	if (fieldcnt == -1) {
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
			printf("%E ", matrix->data[i][n]);
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

double euclidian_dist(const double *v1, const double *v2, const int size) {
    double sum = 0.0;

    for (size_t i=0; i<size; ++i) {
        double minus = v1[i] - v2[i];
        sum += (minus * minus);
    }

    return sqrt(sum);
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

        double dist = euclidian_dist(beta_new, beta_old, A->columns);
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

        double dist = euclidian_dist(beta_new, beta_old, A->columns);
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
	double sample_variance = S / (w_sum  - 1);
	double sample_reliability_variance = S / (w_sum - w_sum2 / w_sum);

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

#endif //STB_STATS_DEFINE
#endif //STB__STATS__H