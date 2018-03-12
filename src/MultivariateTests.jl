__precompile__()

module MultivariateTests

using Compat
using Compat.LinearAlgebra
using StatsBase
using Distributions
using HypothesisTests

using HypothesisTests: HypothesisTest

export
    # General measures of dispersion
    genvar,
    totalvar,
    # Partial correlation
    PartialCorTest,
    partialcor,
    # Equality of mean vectors
    OneSampleHotellingT2,
    EqualCovHotellingT2,
    UnequalCovHotellingT2,
    # Equality of covariances
    BartlettsTest

## Common utility functions

function checkdims(X::AbstractMatrix, Y::AbstractMatrix)
    nx, p = size(X)
    ny, q = size(Y)
    p == q || throw(DimensionMismatch("Inconsistent number of variables: $p, $q"))
    (nx > 0 && ny > 0) || throw(ArgumentError("Inputs must be non-empty"))
    return (p, nx, ny)
end

# Pool the given covariance matrices, trashing them both
function _poolcov!(Sx::AbstractMatrix, nxm1::Int, Sy::AbstractMatrix, nym1::Int)
    Sx .*= nxm1
    Sy .*= nym1
    Sx .+= Sy
    Sx ./= nxm1 + nym1
    return Sx
end

# Benchmarks suggest that this is faster than using pivoted Cholesky when the input
# is not positive definite
function trychol(X::AbstractMatrix)
    S = Symmetric(X)
    return isposdef(S) ? cholfact(S) : factorize(S)
end

# Helper function for computing A'B⁻¹A when B is a covariance matrix, taking advantage
# of the fact that covariance matrices are symmetric, positive (semi-)definite. We'll
# check for positive definiteness and use Cholesky, otherwise let Julia pick Bunch-
# Kaufman or LU in the case of failure.
At_Binv_A(A::AbstractArray, B::AbstractArray) = A'*(trychol(B) \ A)

include("dispersion.jl")
include("partialcor.jl")
include("hotelling.jl")

end # module
