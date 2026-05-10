# EMgLASSO News

## Version 0.1.1 (Development)

### Public Readiness
- Updated package metadata and repository links to the live GitHub repository.
- Aligned license declaration with current LICENSE text.
- Corrected citation metadata (year, URL, maintainer contact).
- Added `CONTRIBUTING.md` and `CODE_OF_CONDUCT.md`.
- Added CI workflow for automated `R CMD check` on push and pull request.
- Added docs index file at `docs/README.md`.

### Reliability Improvements
- Added input validation to `emglasso()`, `rmmnorm()`, and `rspdmatrix()`.
- Normalized mixture proportions (`iTau` / `Tau`) consistently before use.
- Replaced fragile sample allocation in `rmmnorm()` with multinomial allocation,
  ensuring exactly `N` output rows.
- Improved warning behavior for non-convergence and degenerate E-step outputs.

### Testing
- Expanded unit tests to cover argument validation and normalization behavior.

## Version 0.1.0 (Initial Release)

### Major Changes
- First release of EMgLASSO package
- Complete R package structure with DESCRIPTION, NAMESPACE, and roxygen documentation
- Mathematical corrections to EM algorithm derivations:
  - Fixed complete likelihood notation
  - Corrected log-likelihood indicator weighting
  - Corrected posterior responsibility formula
  - Fixed tau update indices in M-step
- Added comprehensive roxygen2 documentation for all exported functions
- Added unit tests in `tests/testthat/`
- Added examples for precision network visualization and EM algorithm demo

### Features
- `emglasso()`: Main function for fitting Gaussian mixture models with sparse covariance estimation
- `rmmnorm()`: Generate samples from Gaussian mixture models
- `rspdmatrix()`: Generate random sparse positive definite matrices
- `estep()`: E-step of EM algorithm (internal function)
- `mstep()`: M-step of EM algorithm (internal function)
- `lstep()`: Graphical LASSO with BIC-based penalty selection (internal function)

### Documentation
- Complete user-focused README.md with quick start guide
- Technical manuscript source in `main.tex` with supporting documentation in `docs/`
- Roxygen-generated help pages for all functions
- Example scripts demonstrating usage

### Infrastructure
- DESCRIPTION file with package metadata and dependencies
- NAMESPACE file specifying exports and imports
- .Rbuildignore for clean package builds
- Unit tests for core functions
- .gitignore for version control

### Dependencies
- mvtnorm (multivariate normal distributions)
- glasso (graphical LASSO estimation)
- MASS (multivariate sampling)
- stats (base R functions)

### Bug Fixes & Corrections
- Removed encoding corruption in original manuscript
- Corrected mathematical notation throughout documentation
- Fixed missing figure references
- Cleaned up obsolete code artifacts

### Known Limitations
- Cluster assignment probabilities require separate computation after model fitting
- No automatic selection of number of components (must be specified via initial `iTau`)
- Performance may degrade for very high dimensions or small sample sizes
