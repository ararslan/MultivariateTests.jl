# Measures of dispersion

Two summary statistics are provided, both of which are generalizations of the sample
variance to multivariate data.

The first is the *generalized variance* of Wilks (1960), which provides a scalar
measure of multidimensional scatter.
For a random vector ``X``, the generalized variance is defined as the determinant of the
sample variance-covariance matrix of ``X``, i.e. ``|\text{Cov}(X)|``.
Note that this can be 0 if there is linear dependence among the columns of ``X``.

The second is the *total variance*.
For a random vector ``X``, the total variance is the matrix trace of the sample
variance-covariance matrix of ``X``, i.e. the sum of the sample variances of the columns
of ``X``.
This can only be 0 if there is only a single observation.

When ``X`` has only one column, both of these measures are equivalent to the sample
variance of the column.

## API

```@docs
MultivariateTests.genvar
MultivariateTests.totalvar
```

## References

Wilks, S.S. (1960). "Multidimensional Statistical Scatter." In *Contributions to
Probability and Statistics*, I. Olkin et al., ed. Stanford University Press, Stanford,
CA, pp. 486-503.
