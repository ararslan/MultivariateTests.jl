using Documenter
using MultivariateTests

makedocs(modules=[MultivariateTests],
         clean=false,
         format=:html,
         sitename="MultivariateTests.jl",
         authors="Alex Arslan and other contributors",
         pages=Any["Home" => "index.md",
                   "Measures of Dispersion" => "dispersion.md",
                   "Partial Correlation" => "partialcor.md",
                   "Equality of Mean Vectors" => "hotelling.md",
                   "Equality of Covariance Matrices" => "covariance.md"])

deploydocs(repo="github.com/ararslan/MultivariateTests.jl.git",
           target="build",
           deps=nothing,
           make=nothing,
           julia="0.6")
