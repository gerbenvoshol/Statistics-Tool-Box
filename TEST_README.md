# STB_STATS Test Suite

## Overview

This test suite provides comprehensive testing for the `stb_stats.h` statistical library. It tests 97+ individual test cases covering the major statistical functions.

## Test Coverage

### ✅ Basic Statistics (6 tests)
- `stb_meanvar` - Mean and sample variance calculation
- `stb_sum` - Neumaier summation algorithm
- `stb_quartiles` - Quartile and median calculation

### ✅ Statistical Tests (12 tests)
- `stb_ttest` - Student's t-test for equal variances
- `stb_uttest` - Welch's t-test for unequal variances  
- `stb_ftest` - F-test for variance comparison

### ✅ P-value Corrections (13 tests)
- `stb_benjamini_hochberg` - Benjamini-Hochberg FDR threshold
- `stb_adjust_pvalues_bh` - Batch BH correction
- `stb_bonferroni` - Bonferroni correction threshold
- `stb_sidak` - Sidak correction
- `stb_log2_fold_change` - Log2 fold change with pseudocount

### ✅ Distribution Functions (12 tests)
- Gumbel distribution (PDF, CDF, ICDF, parameter estimation)
- `stb_pdf_binom` - Binomial PMF
- `stb_pdf_pois` - Poisson PMF

### ✅ Factorial & Combinatorics (5 tests)
- `stb_factorial` - Factorial calculation
- `stb_log_factorial` - Log factorial

### ✅ Contingency Tables (5 tests)
- `stb_fisher2x2` - Fisher exact test for 2×2 tables
- `stb_chisqr` - Chi-square test

### ✅ Correlation (6 tests)
- `stb_spearman` - Spearman's rank correlation
- `stb_kendall` - Kendall's tau correlation

### ✅ Regression (4 tests)
- `stb_linfit` - Linear regression (y = ax + b)

### ✅ Matrix Operations (8 tests)
- Matrix creation, addition, transpose
- Matrix sum operations
- Memory management

### ✅ Diversity Metrics (7 tests)
- `stb_shannon` - Shannon diversity & Pielou evenness
- `stb_simpson` - Simpson diversity index
- `stb_jaccard` - Jaccard similarity
- `stb_bray_curtis` - Bray-Curtis dissimilarity

### ✅ Distance Functions (4 tests)
- `stb_euclidean_distance` - Euclidean distance
- `stb_euclidean_distance_sqr` - Squared Euclidean distance
- `stb_cosine_similarity` - Cosine similarity

### ✅ Random Number Generation (5 tests)
- `stb_pcg32` - PCG32 RNG
- Uniform, Gaussian, exponential distributions
- Bounded integer generation

### ⚠️ Known Limitations

Some functions are not currently tested due to technical issues:
- `stb_xoshiro512` family - Uses static state causing segfaults in test suite
- `stb_dunique` / `stb_iunique` - In-place modification causes memory issues
- Some CDF functions - Require additional investigation

## Building and Running

### Build Test Suite
```bash
make test_stb_stats
```

### Run Tests
```bash
make test
```

Or directly:
```bash
./test_stb_stats
```

### Clean Build Artifacts
```bash
make clean
```

## Continuous Integration

The test suite is integrated with GitHub Actions CI/CD:
- **Ubuntu** and **macOS** testing
- **Code coverage** reporting
- Automatic execution on push/PR to main branches

See `.github/workflows/test.yml` for configuration.

## Test Output

The test suite provides colored output:
- ✓ **Green** - Test passed
- ✗ **Red** - Test failed
- **Blue** - Section headers

Example output:
```
╔════════════════════════════════════════════════════════════╗
║       COMPREHENSIVE STB_STATS TEST SUITE                  ║
╚════════════════════════════════════════════════════════════╝

==== stb_meanvar - Mean and Variance ====
  ✓ Mean of [1,2,3,4,5] is 3.0 (3.0000000000 ≈ 3.0000000000)
  ✓ Sample variance of [1,2,3,4,5] is 2.5 (2.5000000000 ≈ 2.5000000000)
  ...

╔════════════════════════════════════════════════════════════╗
║                      TEST SUMMARY                          ║
╚════════════════════════════════════════════════════════════╝
  Total tests run:    97
  Tests passed:       97
  Tests failed:       0

OVERALL: ALL TESTS PASSED ✓
```

## Adding New Tests

To add tests for additional functions:

1. Create a test function following the naming pattern `test_<function_name>()`
2. Use the provided assertion macros:
   - `TEST_ASSERT(condition, message)` - Boolean assertion
   - `TEST_ASSERT_DOUBLE_EQ(actual, expected, message)` - Exact float comparison
   - `TEST_ASSERT_DOUBLE_NEAR(actual, expected, tolerance, message)` - Tolerance-based comparison
3. Add test function call to `main()`
4. Rebuild and run

Example:
```c
void test_my_function(void) {
    TEST_SECTION("stb_my_function - Description");
    
    double result = stb_my_function(input);
    TEST_ASSERT_DOUBLE_NEAR(result, expected, EPSILON_RELAXED, "Test description");
}
```

## Troubleshooting

### Compilation Warnings

Some warnings from `stb_stats.h` are expected and don't affect functionality:
- Sign comparison warnings in loops
- `fabsf` vs `fabs` type warnings
- Unused parameter warnings

### Segmentation Faults

If you encounter segfaults:
1. Check if testing functions with static state (e.g., `stb_sxoshiro512`)
2. Verify input arrays aren't modified in-place unexpectedly
3. Ensure proper memory cleanup with `free()` calls

### Test Failures

Tests may fail due to:
- Floating-point precision differences across platforms
- Incorrect expected values
- Bugs in the implementation

Adjust `EPSILON` or `EPSILON_RELAXED` tolerances as needed.

## Future Enhancements

Planned improvements:
- [ ] Test remaining CDF functions
- [ ] Resolve xoshiro512 static state issues
- [ ] Add tests for advanced analytics (PCA, k-means, t-SNE, UMAP)
- [ ] Test normalization functions
- [ ] Test logistic regression variants
- [ ] Add performance benchmarks
- [ ] Increase code coverage to 90%+

## License

This test suite is provided under the same license as `stb_stats.h`.

## Contributing

Contributions are welcome! Please:
1. Follow existing code style
2. Add tests for new functions
3. Ensure all tests pass before submitting PR
4. Update this README with new test coverage
