# Contributing to EMgLASSO

Thanks for your interest in contributing.

## How to Contribute

1. Fork the repository and create a feature branch.
2. Make focused changes with clear commit messages.
3. Add or update tests in `tests/testthat/` when behavior changes.
4. Update documentation (`README.md`, roxygen comments, and `NEWS.md`) as needed.
5. Open a pull request describing motivation, changes, and validation steps.

## Development Checklist

- Keep changes scoped and backward-compatible when possible.
- Ensure examples are still runnable.
- Run package checks locally before submitting:

```r
devtools::document()
devtools::test()
devtools::check()
```

## Reporting Issues

Use GitHub Issues with:

- Reproducible code example
- Session information (`sessionInfo()`)
- Expected vs. observed behavior

## Code Style

- Follow existing style in `R/` files.
- Prefer clear variable names and explicit input checks.
- Keep public function documentation in roxygen blocks.