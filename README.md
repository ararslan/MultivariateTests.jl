# MultivariateTests.jl

This package provides utilities for multivariate statistical hypothesis testing in Julia.
It extends the framework provided by
[HypothesisTests](https://github.com/JuliaStats/HypothesisTests.jl)
and aims to provide APIs consistent with those familiar from univariate testing.

Note that the functions here are defined in terms of matrices rather than tabular structures
such as [`DataFrame`s](https://github.com/JuliaData/DataFrames.jl).
Depending on your data format, this means that some conversion may be necessary before the
functions can be used.
