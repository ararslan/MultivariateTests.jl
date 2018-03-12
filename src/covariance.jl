# Tests for equality of covariance matrices

abstract type CovarianceEqualityTest <: HypothesisTest end

HypothesisTests.population_param_of_interest(::CovarianceEqualityTest) =
    ("Equality of covariance matrices", NaN, NaN)

# Finite population correction factor
@inline _correction(p::Int, nxm1::Int, nym1::Int) =
    1 - ((2p^2 + 3p - 1) * (inv(nxm1) + inv(nym1) - inv(nxm1 + nym1)) / 6(p + 1))

## Bartlett's test

struct BartlettsTest <: CovarianceEqualityTest
    L′::Real
    p::Int
    nx::Int
    ny::Int
end

function BartlettsTest(X::AbstractMatrix, Y::AbstractMatrix)
    nx, p = size(X)
    ny, q = size(Y)
    p == q || throw(DimensionMismatch("Inconsistent number of variables"))
    a = nx - 1
    b = ny - 1
    Sx = cov(X)
    Sy = cov(Y)
    L′ = -a * logdet(Sx) - b * logdet(Sy)
    L′ += (a + b) * logdet(_poolcov!(Sx, a, Sy, b))
    L′ *= _correction(p, a, b)
    return BartlettsTest(L′, p, nx, ny)
end

StatsBase.nobs(B::BartlettsTest) = (B.nx, B.ny)
StatsBase.dof(B::BartlettsTest) = div(B.p * (B.p + 1), 2)

HypothesisTests.testname(::BartlettsTest) =
    "Bartlett's Test for Equality of Covariance Matrices"
HypothesisTests.pvalue(B::BartlettsTest) = pvalue(Chisq(dof(B)), B.L′)

function HypothesisTests.show_params(io::IO, B::BartlettsTest, indent="")
    println(io, indent, "number of observations: ", nobs(B))
    println(io, indent, "number of variables:    ", B.p)
    println(io, indent, "χ² statistic:           ", B.L′)
    println(io, indent, "degrees of freedom:     ", dof(B))
end
