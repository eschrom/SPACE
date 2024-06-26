% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measure_mutual_info.R
\name{random_census}
\alias{random_census}
\title{Generate a Random Census}
\usage{
random_census(pls, osp)
}
\arguments{
\item{pls}{List: for each object map or independent scalar, patches of contiguous pixels, at just one radius}

\item{osp}{List: for each object map with attached scalars, a matrix giving the combinations of interest}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{ Data frame: randomized census simulated by scrambling contiguous patches
}\if{html}{\out{</div>}}
}
\description{
Generate a random census to simulate expectations in the absence of patterning
}
\details{
Generate a randomized census from a census and a patch list at a single radius.
}
