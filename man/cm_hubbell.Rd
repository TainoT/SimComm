% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CM_Models.R
\docType{class}
\name{cm_hubbell}
\alias{cm_hubbell}
\title{Hubbell's Neutral Theory Integration Class}
\description{
Hubbell's Neutral Theory Integration Class

Hubbell's Neutral Theory Integration Class
}
\details{
A \link{community_matrixmodel} where each cell contains an individual of one species.
At each generation, an individual may survive or not and empty cells be filled with a new individual.

The survival and generation rules are fixed by the mortality and reproduction rate of each individual.
Edge effects are eliminated by a toroidal correction. (??? check that out)
}
\section{Super classes}{
\code{\link[SimComm:community_model]{SimComm::community_model}} -> \code{\link[SimComm:community_matrixmodel]{SimComm::community_matrixmodel}} -> \code{cm_hubbell}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{death_rate}}{The mortality rate of an individual.
Default is \code{0.1}}

\item{\code{disturbance_rate}}{The disturbance rate of a community.
Default is \code{1}}

\item{\code{migration_rate}}{The migration rate of an individual.
Default is \code{0.005}}

\item{\code{speciation_rate}}{The speciation rate of an individual.
Default is \code{0.0001}}

\item{\code{local_cm}}{The matrix local community}

\item{\code{meta_cm}}{The matrix meta community}

\item{\code{model}}{The model to use}

\item{\code{species_mapping}}{A mapping of meta species to local species}

\item{\code{next_local_id}}{The next avilaible local id}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-cm_hubbell-new}{\code{cm_hubbell$new()}}
\item \href{#method-cm_hubbell-show_values}{\code{cm_hubbell$show_values()}}
\item \href{#method-cm_hubbell-clone}{\code{cm_hubbell$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_model" data-id="along_time"><a href='../../SimComm/html/community_model.html#method-community_model-along_time'><code>SimComm::community_model$along_time()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_model" data-id="run"><a href='../../SimComm/html/community_model.html#method-community_model-run'><code>SimComm::community_model$run()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_model" data-id="saved_pattern"><a href='../../SimComm/html/community_model.html#method-community_model-saved_pattern'><code>SimComm::community_model$saved_pattern()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_matrixmodel" data-id="autoplot"><a href='../../SimComm/html/community_matrixmodel.html#method-community_matrixmodel-autoplot'><code>SimComm::community_matrixmodel$autoplot()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_matrixmodel" data-id="graph"><a href='../../SimComm/html/community_matrixmodel.html#method-community_matrixmodel-graph'><code>SimComm::community_matrixmodel$graph()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_matrixmodel" data-id="neighbors"><a href='../../SimComm/html/community_matrixmodel.html#method-community_matrixmodel-neighbors'><code>SimComm::community_matrixmodel$neighbors()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_matrixmodel" data-id="plot"><a href='../../SimComm/html/community_matrixmodel.html#method-community_matrixmodel-plot'><code>SimComm::community_matrixmodel$plot()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_matrixmodel" data-id="prepare_buffer"><a href='../../SimComm/html/community_matrixmodel.html#method-community_matrixmodel-prepare_buffer'><code>SimComm::community_matrixmodel$prepare_buffer()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_matrixmodel" data-id="synthesis"><a href='../../SimComm/html/community_matrixmodel.html#method-community_matrixmodel-synthesis'><code>SimComm::community_matrixmodel$synthesis()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cm_hubbell-new"></a>}}
\if{latex}{\out{\hypertarget{method-cm_hubbell-new}{}}}
\subsection{Method \code{new()}}{
community evolution rules for the neutral theory model


disturb the community at given time and replace the dead cells with the neighbors


Create a new instance of this \link[R6:R6Class]{R6} class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cm_hubbell$new(
  local_pattern = NULL,
  timeline = 0,
  type = "Species",
  neighborhood = "Moore 1",
  global_pattern = NULL,
  model = self$model,
  death_rate = self$death_rate,
  disturbance_rate = self$disturbance_rate,
  event = NULL,
  migration_rate = self$migration_rate,
  speciation_rate = self$speciation_rate
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{timeline}}{A numeric vector that contains the points of time of interest.}

\item{\code{type}}{The type of individuals. Information only.}

\item{\code{neighborhood}}{A character string defining what is the neighborhood of a cell:
"von Neumann 1" or "4" for the closest four neighbors (North, West, South, East);
"Moore 1" or "8" for all adjacent cells (the first four and North-West, etc.);
"Moore 2" or "24" for two rings of neighbors.
"Global" or "0" for the whole matrix as neighborhood.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cm_hubbell-show_values"></a>}}
\if{latex}{\out{\hypertarget{method-cm_hubbell-show_values}{}}}
\subsection{Method \code{show_values()}}{
Show the values of the model
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cm_hubbell$show_values(values = "basic")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cm_hubbell-clone"></a>}}
\if{latex}{\out{\hypertarget{method-cm_hubbell-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cm_hubbell$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
