% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CP_Models.R
\name{cp_hubbell}
\alias{cp_hubbell}
\title{Hubbell's Neutral Theory Integration Class}
\description{
Hubbell's Neutral Theory Integration Class

Hubbell's Neutral Theory Integration Class
}
\details{
A \link{community_scpmodel} where each point contains an individual of one species.
}
\section{Super classes}{
\code{\link[SimComm:community_model]{SimComm::community_model}} -> \code{\link[SimComm:community_spcmodel]{SimComm::community_spcmodel}} -> \code{cp_hubbell}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{death_rate}}{The mortality rate of an individual.
Default is \code{0.1}}

\item{\code{kill_rate}}{The disturbance rate of the community.
Default is \code{0}}

\item{\code{migration_rate}}{The migration rate of an individual.
Default is '0.005'}

\item{\code{speciation_rate}}{The speciation rate of an individual.
Default is '0.001'}

\item{\code{model}}{The wmppp local community.}

\item{\code{model}}{The wmppp meta community.}

\item{\code{model}}{The model to use}

\item{\code{n_neighbors}}{The number of nearest neighbors to take into account.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-cp_hubbell-new}{\code{cp_hubbell$new()}}
\item \href{#method-cp_hubbell-clone}{\code{cp_hubbell$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_model" data-id="along_time"><a href='../../SimComm/html/community_model.html#method-community_model-along_time'><code>SimComm::community_model$along_time()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_model" data-id="autoplot"><a href='../../SimComm/html/community_model.html#method-community_model-autoplot'><code>SimComm::community_model$autoplot()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_model" data-id="run"><a href='../../SimComm/html/community_model.html#method-community_model-run'><code>SimComm::community_model$run()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_model" data-id="saved_pattern"><a href='../../SimComm/html/community_model.html#method-community_model-saved_pattern'><code>SimComm::community_model$saved_pattern()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_spcmodel" data-id="neighbor_types_r"><a href='../../SimComm/html/community_spcmodel.html#method-community_spcmodel-neighbor_types_r'><code>SimComm::community_spcmodel$neighbor_types_r()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_spcmodel" data-id="neighbors_n"><a href='../../SimComm/html/community_spcmodel.html#method-community_spcmodel-neighbors_n'><code>SimComm::community_spcmodel$neighbors_n()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_spcmodel" data-id="neighbors_r"><a href='../../SimComm/html/community_spcmodel.html#method-community_spcmodel-neighbors_r'><code>SimComm::community_spcmodel$neighbors_r()</code></a></span></li>
<li><span class="pkg-link" data-pkg="SimComm" data-topic="community_spcmodel" data-id="plot"><a href='../../SimComm/html/community_spcmodel.html#method-community_spcmodel-plot'><code>SimComm::community_spcmodel$plot()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cp_hubbell-new"></a>}}
\if{latex}{\out{\hypertarget{method-cp_hubbell-new}{}}}
\subsection{Method \code{new()}}{
disturb the community at given time and replace the dead cells with the neighbours
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cp_hubbell$new(
  local_pattern = NULL,
  timeline = 0,
  type = "Species",
  n_neighbors = 6,
  global_pattern = NULL,
  model = self$model,
  death_rate = self$death_rate,
  disturbance_rate = self$kill_rate,
  event = NULL,
  migration_rate = self$migration_rate,
  speciation_rate = self$speciation_rate
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cp_hubbell-clone"></a>}}
\if{latex}{\out{\hypertarget{method-cp_hubbell-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{cp_hubbell$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
