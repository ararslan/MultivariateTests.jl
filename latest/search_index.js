var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MultivariateTests.jl-1",
    "page": "Home",
    "title": "MultivariateTests.jl",
    "category": "section",
    "text": "MultivariateTests is a Julia package that extends the framework provided by the HypothesisTests to include functionality for multivariate statistics."
},

{
    "location": "index.html#Functionality-1",
    "page": "Home",
    "title": "Functionality",
    "category": "section",
    "text": "Pages = [\n    \"dispersion.md\",\n    \"partialcor.md\",\n    \"hotelling.md\",\n    \"covariance.md\",\n]\nDepth = 1"
},

{
    "location": "index.html#Acknowledgements-1",
    "page": "Home",
    "title": "Acknowledgements",
    "category": "section",
    "text": "Initial work on this package by Alex Arslan was part of a research project for the master\'s program in applied statistics at Penn State University. Alex thanks Andreas Noack, Jiahao Chen, Moritz Schauer, and Simon Byrne for initial feedback and discussion in the early stages of development."
},

{
    "location": "dispersion.html#",
    "page": "Measures of Dispersion",
    "title": "Measures of Dispersion",
    "category": "page",
    "text": ""
},

{
    "location": "dispersion.html#Measures-of-dispersion-1",
    "page": "Measures of Dispersion",
    "title": "Measures of dispersion",
    "category": "section",
    "text": "Two summary statistics are provided, both of which are generalizations of the sample variance to multivariate data.The first is the generalized variance of Wilks (1960), which provides a scalar measure of multidimensional scatter. For a random vector X, the generalized variance is defined as the determinant of the sample variance-covariance matrix of X, i.e. textCov(X). Note that this can be 0 if there is linear dependence among the columns of X.The second is the total variance. For a random vector X, the total variance is the matrix trace of the sample variance-covariance matrix of X, i.e. the sum of the sample variances of the columns of X. This can only be 0 if there is only a single observation.When X has only one column, both of these measures are equivalent to the sample variance of the column."
},

{
    "location": "dispersion.html#MultivariateTests.genvar",
    "page": "Measures of Dispersion",
    "title": "MultivariateTests.genvar",
    "category": "function",
    "text": "genvar(X)\n\nCompute the generalized sample variance of X. If X is a vector, one-column matrix, or other one-dimensional iterable, this is equivalent to the sample variance. Otherwise if X is a matrix, this is equivalent to the determinant of the covariance matrix of X.\n\nnote: Note\nThe generalized sample variance will be 0 if the columns of the matrix of deviations are linearly dependent.\n\n\n\n"
},

{
    "location": "dispersion.html#MultivariateTests.totalvar",
    "page": "Measures of Dispersion",
    "title": "MultivariateTests.totalvar",
    "category": "function",
    "text": "totalvar(X)\n\nCompute the total sample variance of X. If X is a vector, one-column matrix, or other one-dimensional iterable, this is equivalent to the sample variance. Otherwise if X is a matrix, this is equivalent to the sum of the diagonal elements of the covariance matrix of X.\n\n\n\n"
},

{
    "location": "dispersion.html#API-1",
    "page": "Measures of Dispersion",
    "title": "API",
    "category": "section",
    "text": "MultivariateTests.genvar\nMultivariateTests.totalvar"
},

{
    "location": "dispersion.html#References-1",
    "page": "Measures of Dispersion",
    "title": "References",
    "category": "section",
    "text": "Wilks, S.S. (1960). \"Multidimensional Statistical Scatter.\" In Contributions to Probability and Statistics, I. Olkin et al., ed. Stanford University Press, Stanford, CA, pp. 486-503."
},

{
    "location": "partialcor.html#",
    "page": "Partial Correlation",
    "title": "Partial Correlation",
    "category": "page",
    "text": ""
},

{
    "location": "partialcor.html#Partial-correlation-1",
    "page": "Partial Correlation",
    "title": "Partial correlation",
    "category": "section",
    "text": "This package provides a function for computing the partial correlation directly and another for testing the hypothesis of zero partial correlation using the HypothesisTests framework.Partial correlation is a generalization of correlation to conditional distributions; it\'s akin to Pearson\'s correlation, substituting conditional variances and covariances in the computation. It can be thought of as the correlation between two random variables for the subpopulation defined by the conditional variable(s)."
},

{
    "location": "partialcor.html#MultivariateTests.partialcor",
    "page": "Partial Correlation",
    "title": "MultivariateTests.partialcor",
    "category": "function",
    "text": "partialcor(x, y, Z)\n\nCompute the partial correlation of the vectors x and y given Z, which can be a vector or matrix.\n\n\n\n"
},

{
    "location": "partialcor.html#MultivariateTests.PartialCorTest",
    "page": "Partial Correlation",
    "title": "MultivariateTests.PartialCorTest",
    "category": "type",
    "text": "PartialCorTest(x, y, Z)\n\nPerform a t-test for the hypothesis that textCor(xyZ=z) = 0, i.e. the partial correlation of vectors x and y given the matrix Z is zero.\n\nImplements pvalue and confint. See also partialcor.\n\n\n\n"
},

{
    "location": "partialcor.html#API-1",
    "page": "Partial Correlation",
    "title": "API",
    "category": "section",
    "text": "MultivariateTests.partialcor\nMultivariateTests.PartialCorTest"
},

{
    "location": "hotelling.html#",
    "page": "Equality of Mean Vectors",
    "title": "Equality of Mean Vectors",
    "category": "page",
    "text": ""
},

{
    "location": "hotelling.html#Hotelling\'s-T2-test-1",
    "page": "Equality of Mean Vectors",
    "title": "Hotelling\'s T^2 test",
    "category": "section",
    "text": "The HypothesisTests framework is extended here to include Hotelling\'s T^2 test for one and two samples."
},

{
    "location": "hotelling.html#MultivariateTests.OneSampleHotellingT2",
    "page": "Equality of Mean Vectors",
    "title": "MultivariateTests.OneSampleHotellingT2",
    "category": "type",
    "text": "OneSampleHotellingT2(X::AbstractMatrix, μ₀=<zero vector>)\n\nPerform a one sample Hotelling\'s T^2 test of the hypothesis that the vector of column means of X is equal to μ₀.\n\n\n\nOneSampleHotellingT2(X::AbstractMatrix, Y::AbstractMatrix, μ₀=<zero vector>)\n\nPerform a paired Hotelling\'s T^2 test of the hypothesis that the vector of mean column differences between X and Y is equal to μ₀.\n\n\n\n"
},

{
    "location": "hotelling.html#MultivariateTests.EqualCovHotellingT2",
    "page": "Equality of Mean Vectors",
    "title": "MultivariateTests.EqualCovHotellingT2",
    "category": "type",
    "text": "EqualCovHotellingT2(X::AbstractMatrix, Y::AbstractMatrix)\n\nPerform a two sample Hotelling\'s T^2 test of the hypothesis that the difference in the mean vectors of X and Y is zero, assuming that X and Y have equal covariance matrices.\n\n\n\n"
},

{
    "location": "hotelling.html#MultivariateTests.UnequalCovHotellingT2",
    "page": "Equality of Mean Vectors",
    "title": "MultivariateTests.UnequalCovHotellingT2",
    "category": "type",
    "text": "UnequalCovHotellingT2(X::AbstractMatrix, Y::AbstractMatrix)\n\nPerform a two sample Hotelling\'s T^2 test of the hypothesis that the difference in the mean vectors of X and Y is zero, without assuming that X and Y have equal covariance matrices.\n\n\n\n"
},

{
    "location": "hotelling.html#API-1",
    "page": "Equality of Mean Vectors",
    "title": "API",
    "category": "section",
    "text": "MultivariateTests.OneSampleHotellingT2\nMultivariateTests.EqualCovHotellingT2\nMultivariateTests.UnequalCovHotellingT2"
},

{
    "location": "covariance.html#",
    "page": "Equality of Covariance Matrices",
    "title": "Equality of Covariance Matrices",
    "category": "page",
    "text": ""
},

{
    "location": "covariance.html#Bartlett\'s-test-1",
    "page": "Equality of Covariance Matrices",
    "title": "Bartlett\'s test",
    "category": "section",
    "text": "Bartlett\'s test for equality of two variance-covariance matrices is provided. This is equivalent to Box\'s M-test for two groups."
},

{
    "location": "covariance.html#MultivariateTests.BartlettsTest",
    "page": "Equality of Covariance Matrices",
    "title": "MultivariateTests.BartlettsTest",
    "category": "type",
    "text": "BartlettsTest(X::AbstractMatrix, Y::AbstractMatrix)\n\nPerform Bartlett\'s test of the hypothesis that the covariance matrices of X and Y are equal.\n\nnote: Note\nBartlett\'s test is sensitive to departures from multivariate normality.\n\n\n\n"
},

{
    "location": "covariance.html#API-1",
    "page": "Equality of Covariance Matrices",
    "title": "API",
    "category": "section",
    "text": "MultivariateTests.BartlettsTest"
},

]}
