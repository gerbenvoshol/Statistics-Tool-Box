# Statistics-Tool-Box

This is a single header file inspired by stb.h by Sean Barrett with a bunch of useful statistical functions

 ============================================================================
 
	 You MUST

			#define STB_STATS_DEFINE

	 in EXACTLY _one_ C or C++ file that includes this header, BEFORE the
	 include, like this:

			#define STB_STATS_DEFINE
			#include "stb_stats.h"

	 All other files should just #include "stb_stats.h" without the #define.
 ============================================================================

Functions included are:
* stb_moderated_ttest, stb_cosine_similarity, RSE Normalization (stb_calc_geometric_scaling_factors and stb_meanvar_counts_to_common_scale)
* stb_shannon (Shannon's diversity index, Pilou evenness, stb_simpson (Simpson's Diversity Index), stb_jaccard (Jaccard similarity index), stb_bray_curtis (Bray–Curtis dissimilarity) and stb_create_htable a simple basic hash table
* stb_pdf_hypgeo hypergeometric distribution probability density function, speedup stb_log_factorial using lookup table
* stb_fisher2x2 simple fisher exact test for 2x2 contigency tables
* stb_pdf_binom and stb_pdf_pois, the binomial and poison probability density functions
* stb_polygamma, stb_trigamma_inverse gamme functions and stb_fit_f_dist for moment estimation of the scaled F-distribution
* stb_qnorm and stb_qnorm_with_reference (also matrix variants) quantile normalization between columns with and without a reference
* stb_neugas Neural gas clustering algorithm
* stb_pca Principal Component Analysis
* stb_csm (confident sequence method) for monte-carlo simulations
* stb_kmeans k-means++ classical data clustering
* stb_qsort (Quicksort), could be used to replace current sorting method
* stb_cdf_gumbel, stb_pdf_gumbel, stb_icdf_gumbel and stb_est_gumbel, the (inverse) cumulative/probability 
 		      density functions for the gumbel distribution and the ML estimator of the gumbel parameters
* stb_kendall (Kendall's Rank correlation)
* stb_jenks Initial port of O(k×n×log(n)) Jenks-Fisher algorithm originally created by Maarten Hilferink
* stb_logistic_regression_L2 simple L2-regularized logistic regression
* stb_spearman (Spearman's Rank correlation)
* stb_invert_matrix, stb_transpose_matrix, stb_matrix_multiply, ..., stb_multi_linear_regression and stb_multi_logistic_regression 
* stb_ksample_anderson_darling, stb_2sample_anderson_darling, (one sample) stb_anderson_darling
* stb_expfit (Exponential fitting), stb_polyfit (Polynomial fitting), stb_powfit (Power curve fitting), stb_linfit (Liniear fitting)
* stb_trap, stb_trapezoidal (returns the integral (area under the cruve) of a given function and interval)
* stb_lagrange (polynomial interpolation), stb_sum (Neumaier summation algorithm)
* stb_mann_whitney, stb_kruskal_wallis (Unfinished, needs a better way to handle Dunn's post-hoc test)
* stb_combinations
* stb_allocmat (simple allocation of 2d array, but might not work on all systems?!)
* stb_fgetln, stb_fgetlns
* stb_pcg32 (PCG-XSH-RR) and stb_xoshiro512 (xoshiro512**) Pseudo Random Number Generators
* stb_anova (One-Way Anova with Tukey HSD test and Scheffe T-statistics method (post-hoc) (Unfinished))
* stb_quartiles
* stb_histogram (very simple histogram), stb_print_histogram, ...
* stb_factorial
* stb_meanvar
* stb_ttest, stb_uttest
* stb_ftest, 
* stb_benjamini_hochberg
* stb_chisqr, stb_chisqr_matrix, stb_gtest, stb_gtest_matrix, 

CITATION

If you use this Tool-Box in a publication, please reference:

Voshol, G.P. (2022). STB: A simple Statistics Tool Box (Version 1.23) [Software]. 
Available from https://github.com/gerbenvoshol/Statistics-Tool-Box
