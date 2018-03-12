# Partial correlation

This package provides a function for computing the partial correlation directly and
another for testing the hypothesis of zero partial correlation using the HypothesisTests
framework.

Partial correlation is a generalization of correlation to conditional distributions;
it's akin to Pearson's correlation, substituting conditional variances and covariances in
the computation.
It can be thought of as the correlation between two random variables for the subpopulation
defined by the conditional variable(s).

## API

```@docs
MultivariateTests.partialcor
MultivariateTests.PartialCorTest
```
