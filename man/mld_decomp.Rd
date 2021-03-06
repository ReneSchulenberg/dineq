% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mld_decomp.R
\name{mld_decomp}
\alias{mld_decomp}
\title{Decomposition of the mean log deviation}
\source{
Plat, D. (2012). IC2: Inequality and Concentration Indices and Curves. R package
version 1.0-1. https://CRAN.R-project.org/package=IC2
}
\usage{
mld_decomp(x, z, weights = NULL)
}
\arguments{
\item{x}{a numeric vector containing at least non-negative elements.}

\item{z}{a factor containing the population subgroups.}

\item{weights}{an optional vector of weights of x to be used in the computation of
the decomposition. Should be NULL or a numeric vector.}
}
\value{
a list with the results of the decomposition and the parts used for the
decomposition, containing the following components:
   \item{mld_decomp}{a list containing the decomposition: mld_total (value of the mean log
   deviation index of x) mld_within (value of within-group inequality) and mld_between (value
   of between-group inequality)}
   \item{mld_group}{a list containing mld_group (the mean log deviations of the different
   subgroups) and mld_group_contribution(the contribution of the subgroups to the total
   within-group inequality: adds up to mld_within)}
   \item{mld_decomp}{a list containing the means of x: mean_total (value of the mean of x of
   all subgroups combined) and mean_group (value of the mean of x of the individual subgroups)
   inequality) and mld_between (value of between-group inequality)}
   \item{share_groups}{the distribution of the subgroups z }
   \item{share_income_groups}{the distribution of vector x by subgroups z }
   \item{number_cases}{a list containing the number of cases in total, by subgroup (weighted
   and unweighted): n_unweighted (total number of unweighted x), n_weighted (total number of
   weighted x), n_group_unweighted (number of unweighted x by subgroup z), n_group_unweighted
   (number of weighted x by subgroup z)}
   \item{note}{number of zero or negative observations. The mean log deviation uses a logarithmic
   transformation of x. Therefore these observations are deleted from the analysis}
}
\description{
Decomposes the mean log deviation into non overlapping population subgroups.
Distinction is made by between and within group inequality.
}
\details{
The decomposition of the mean log deviation by between and within group
inequality. Within group inequality is calculated by using the mean log deviation for each sub
group. Between group inequality by the mean log deviation of the average of both sub groups.

It uses a logarithmic transformation of the values of the distribution. Therefore it
cannot handle negative or zero values. Those are excluded from the computation in this function.

Based on calcGEI function in IC2 package. Handles missings.
}
\examples{
#Decomposition of mean log deviation by level of education using Mexican Income data set
data(mex_inc_2008)
education_decomp <- mld_decomp(x=mex_inc_2008$income,z=mex_inc_2008$education,
weights=mex_inc_2008$factor)

#complete output
education_decomp

#Selected output: decomposition into between- and within-group inequality
education_decomp["mld_decomp"]

}
\references{
Mookherjee, D. and A. Shorrocks (1982) A decomposition analysis of the trend in UK income
inequality, \emph{Economic Journal}, 92 (368), p. 886-902.

Brewer M., and L. Wren-Lewis (2016) Accounting for Changes in Income Inequality:
Decomposition Analyses for the UK, 1978–2008. \emph{Oxford Bulletin of economics and
statistics}, 78 (3), p. 289-322,

Haughton, J. and S. Khandker. (2009) \emph{Handbook on poverty and inequality},
Washington, DC: World Bank.
}
\seealso{
\code{\link{mld_change}}
         \code{\link{gini_decomp}}
}
