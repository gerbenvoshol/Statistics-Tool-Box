# Statistics-Tool-Box

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

Functions included are:
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
