# Partial correlation

"""
    partialcor(X, Y, Z)

Compute the partial correlation of the vectors `X` and `Y` given `Z`, which can be
a vector or matrix.
"""
function partialcor(x::AbstractVector, y::AbstractVector, Z::AbstractVecOrMat)
    length(x) == length(y) == size(Z, 1) ||
        throw(DimensionMismatch("Inputs must have the same number of observations"))
    length(x) > 0 || throw(ArgumentError("Inputs must be non-empty"))
    return Base.clampcor(_partialcor(x, mean(x), y, mean(y), Z))
end

function _partialcor(x::AbstractVector, μx, y::AbstractVector, μy, Z::AbstractMatrix)
    p = size(Z, 2)
    p == 1 && return _partialcor(x, μx, y, μy, vec(Z))
    z₀   = view(Z, :, 1)
    Zmz₀ = view(Z, :, 2:p)
    μz₀ = mean(z₀)
    rxz = _partialcor(x,  μx,  z₀, μz₀, Zmz₀)
    rzy = _partialcor(z₀, μz₀, y,  μy,  Zmz₀)
    rxy = _partialcor(x,  μx,  y,  μy,  Zmz₀)::typeof(rxz)
    return (rxy - rxz * rzy) / (sqrt(1 - rxz^2) * sqrt(1 - rzy^2))
end

function _partialcor(x::AbstractVector, μx, y::AbstractVector, μy, z::AbstractVector)
    μz = mean(z)

    # Initialize all of the accumulators to 0 of the appropriate types
    Σxx = abs2(zero(eltype(x)) - zero(μx))
    Σyy = abs2(zero(eltype(y)) - zero(μy))
    Σzz = abs2(zero(eltype(z)) - zero(μz))
    Σxy = zero(Σxx * Σyy)
    Σxz = zero(Σxx * Σzz)
    Σzy = zero(Σzz * Σyy)

    # We only want to make one pass over all of the arrays
    @inbounds begin
        @simd for i in eachindex(x, y, z)
            xi = x[i] - μx
            yi = y[i] - μy
            zi = z[i] - μz

            Σxx += abs2(xi)
            Σyy += abs2(yi)
            Σzz += abs2(zi)

            Σxy += xi * yi
            Σxz += xi * zi
            Σzy += zi * yi
        end
    end

    # Individual pairwise correlations
    rxy = Σxy / sqrt(Σxx * Σyy)
    rxz = Σxz / sqrt(Σxx * Σzz)
    rzy = Σzy / sqrt(Σzz * Σyy)

    return (rxy - rxz * rzy) / (sqrt(1 - rxz^2) * sqrt(1 - rzy^2))
end

"""
    PartialCorTest(X, Y, Z)

Perform a t-test for the hypothesis that the partial correlation of `X` and `Y` given
`Z` is zero against the alternative that it's nonzero.

Implements `pvalue` and `confint`. See also [`partialcor`](@ref).
"""
struct PartialCorTest <: HypothesisTest
    r::Real
    n::Int
    k::Int
    t::Real

    # Error checking is done in `partialcor`
    function PartialCorTest(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
        r = partialcor(x, y, Z)
        n, k = size(Z)
        t = r * sqrt((n - 2 - k) / (1 - r^2))
        return new(r, n, k, t)
    end

    function PartialCorTest(x::AbstractVector, y::AbstractVector, z::AbstractVector)
        r = partialcor(x, y, z)
        n = length(z)
        t = r * sqrt((n - 3) / (1 - r^2))
        return new(r, n, 1, t)
    end
end

HypothesisTests.testname(::PartialCorTest) = "Test for partial correlation"
HypothesisTests.population_param_of_interest(p::PartialCorTest) =
    ("Partial correlation", zero(p.r), p.r)

StatsBase.nobs(p::PartialCorTest) = p.n
StatsBase.dof(p::PartialCorTest) = p.n - 2 - p.k

function StatsBase.confint(test::PartialCorTest, alpha::Float64=0.05)
    dof(test) > 1 || return (-one(test.r), one(test.r)) # Otherwise we can get NaNs
    q = quantile(Normal(), 1 - alpha / 2)
    fisher = (log1p(test.r) - log1p(-test.r)) / 2
    bound = q / sqrt(test.n - 3 - test.k)
    lo = 2 * (fisher - bound)
    hi = 2 * (fisher + bound)
    elo = Base.clampcor(expm1(lo) / (exp(lo) + 1))
    ehi = Base.clampcor(expm1(hi) / (exp(hi) + 1))
    return (elo, ehi)
end

HypothesisTests.default_tail(::PartialCorTest) = :both
HypothesisTests.pvalue(test::PartialCorTest; tail=:both) =
    pvalue(TDist(dof(test)), test.t, tail=tail)

function HypothesisTests.show_params(io::IO, test::PartialCorTest, indent="")
    println(io, indent, "number of observations:          ", nobs(test))
    println(io, indent, "number of conditional variables: ", test.k)
    println(io, indent, "t-statistic:                     ", test.t)
    println(io, indent, "degrees of freedom:              ", dof(test))
end
