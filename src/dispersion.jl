# Summary statistics

"""
    genvar(X)

Compute the generalized sample variance of `X`. If `X` is a vector, one-column matrix,
or other one-dimensional iterable, this is equivalent to the sample variance.
Otherwise if `X` is a matrix, this is equivalent to the determinant of the covariance
matrix of `X`.

!!! note
    The generalized sample variance will be 0 if the columns of the matrix of deviations
    are linearly dependent.
"""
genvar(X::AbstractMatrix) = size(X, 2) == 1 ? var(vec(X)) : det(cov(X))
genvar(itr) = var(itr)

"""
    totalvar(X)

Compute the total sample variance of `X`. If `X` is a vector, one-column matrix,
or other one-dimensional iterable, this is equivalent to the sample variance.
Otherwise if `X` is a matrix, this is equivalent to the sum of the diagonal elements
of the covariance matrix of `X`.
"""
totalvar(X::AbstractMatrix) = sum(var(X, dims=1))
totalvar(itr) = var(itr)
