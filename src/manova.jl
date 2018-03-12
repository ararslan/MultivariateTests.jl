# MANOVA

abstract type MANOVATest <: HypothesisTest end

for (Tn, f) in [(:WilksLambdaTest, :_wilks),
                (:PillaiTraceTest, :_pillai),
                (:LawleyTraceTest, :_lawley),
                (:RoyMaxRootTest,  :_roy)]
    @eval struct $Tn <: MANOVATest
        stat::Real
        F::Real
        df::NTuple{2,Int}
        SS::SSCP
    end
    @eval ($Tn)(X::AbstractMatrix, group::CategoricalArray) = ($f)(SSCP(X, group))
end

StatsBase.nobs(M::MANOVATest) = M.SS.dfT + 1
StatsBase.dof(M::MANOVATest) = M.df

HypothesisTests.pvalue(M::MANOVATest) = pvalue(FDist(dof(M)...), M.F)
HypothesisTests.population_param_of_interest(M::MANOVATest) =
    ("Equality of group means", NaN, NaN)

function HypothesisTests.show_params(io::IO, M::MANOVATest, indent="")
    println(io, indent, "number of observations: ", nobs(M))
    println(io, indent, "number of variables:    ", size(M.SS.H, 2))
    println(io, indent, "number of groups:       ", M.SS.dfH + 1)
    println(io, indent, "test statistic:         ", M.stat)
    println(io, indent, "transformed statistic:  ", M.F)
    println(io, indent, "degrees of freedom:     ", dof(M))
end

## Sums of squares and cross products

struct SSCP{T}
    H::Symmetric{T,<:AbstractMatrix{T}}
    E::Factorization{T}
    T::Factorization{T}
    dfH::Int
    dfE::Int
    dfT::Int
end

"""
    SSCP(X::AbstractMatrix, group::CategoricalArray)

Compute the sums of squares and cross products (SSCP) matrices based on the data matrix
`X` and the grouping `group`.
"""
function SSCP(X::AbstractMatrix, group::CategoricalArray)
    nobs, p = size(X)
    nobs == length(group) || throw(DimensionMismatch("Inconsistent number of observations"))
    groups = levels(group)
    ngroups = length(groups)
    isempty(groups) && throw(ArgumentError("At least one non-missing group is required"))
    g1 = first(groups)
    E = cov(view(X, group .== g1, :))
    @inbounds for i = 2:ngroups
        Xᵢ = view(X, group .== groups[i], :)
        E .+= (size(Xᵢ, 1) - 1) .* cov(Xᵢ)
    end
    Z = X .- mean(X, 1)
    T = Z'Z
    H = T - E
    return SSCP{eltype(H)}(Symmetric(H), trychol(E), trychol(T),
                           ngroups - 1, nobs - ngroups, nobs - 1)
end

## Wilk's lambda

HypothesisTests.testname(::WilksLambdaTest) = "Wilks Lambda Test"

function _wilks(SS::SSCP)
    Λ = det(SS.E) / det(SS.T)
    p = size(SS.H, 2)
    q = SS.dfH
    r = SS.dfE - (p - q + 1) / 2
    u = (p * q - 2) / 4
    tmp = p^2 + q^2 - 5
    t = tmp > 0 ? sqrt((p^2 * q^2 - 4) / tmp) : one(sqrt(middle(tmp)))
    iΛ = Λ^inv(t)
    df2 = trunc(Int, p * q)
    df1 = trunc(Int, r * t - 2u)
    F = ((1 - iΛ) * df2) / (iΛ * df1)
    return WilksLambdaTest(Λ, F, (df1, df2), SS)
end

## Pilai's trace

HypothesisTests.testname(::PillaiTraceTest) = "Pillai Trace Test"

function _pillai(SS::SSCP)
    V = trace(SS.H * inv(SS.T))
    p = size(SS.H, 2)
    q = SS.dfH
    s = min(p, q)
    m = (abs(p - q) - 1) / 2
    n = (SS.dfE - p - 1) / 2
    num = 2n + s + 1
    den = 2m + s + 1
    F = num * V / (den * (s - V))
    return PillaiTraceTest(V, F, (s * den, s * num), SS)
end

## Hotelling-Lawley trace

HypothesisTests.testname(::LawleyTraceTest) = "Hotelling-Lawley Trace Test"

function _lawley(SS::SSCP)
    U = trace(SS.H * inv(SS.E))
    p = size(SS.H, 2)
    q = SS.dfH
    s = min(p, q)
    m = (abs(p - q) - 1) / 2
    n = (SS.dfE - p - 1) / 2
    if n > 0
        b = (p + 2n) * (q + 2n) / (2 * (2n + 1) * (n - 1))
        c = (2 + (p * q + 2) / (b - 1)) / 2n
        F = (U / c) * ((4 + (p * q + 2) / (b - 1)) / (p * q))
        df1 = p * q
        df2 = div(4 + (p * q + 2), b - 1)
    else
        F = 2 * (s * n + 1) * U / (s^2 * (2m + s + 1))
        df1 = s * (2m + s + 1)
        df2 = 2 * (s * n + 1)
    end
    return LawleyTraceTest(U, F, (df1, df2), SS)
end

## Roy's maximum root

HypothesisTests.testname(::RoyMaxRootTest) = "Roy's Maximum Root Test"

function _roy(SS::SSCP)
    Θ = eigmax(SS.H * inv(SS.E))
    p = size(SS.H, 2)
    q = SS.dfH
    r = min(p, q)
    den = SS.dfE - r + q
    F = Θ * den / r
    return RoyMaxRootTest(Θ, F, (r, den), SS)
end
