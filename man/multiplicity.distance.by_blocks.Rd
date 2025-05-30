% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiplicity.distance.R
\name{multiplicity.distance.by_blocks}
\alias{multiplicity.distance.by_blocks}
\title{Distance based Multiplicity by Blocks}
\usage{
multiplicity.distance.by_blocks(ids, abs, diss_frame, clust, sigma = 1)
}
\arguments{
\item{ids}{1 x n vector of the ids of the elements.}

\item{diss_frame}{a three columns data frame, with the distances between
elements. Columns are assumed to be in the order: ID1, ID2 and distance.}

\item{clust}{List or vector  element of size n with the unit's corresponding
cluster}

\item{ab}{1 x n vector of abundances before clustering, in the same order
as the ids.}

\item{sig}{float value determining that two units are different}
}
\value{
\eqn{FM_{\sigma}}.
}
\description{
Distance Multiplicity defined as the ratio between the functional
diversities of the unclustered set and the cluster set. Computed over a
block distance matrix. Here we assume that the distance between elements of
different clusters is maximum (i.e. sigma)
}
