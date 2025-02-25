% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulations.R
\docType{class}
\name{model_simulation}
\alias{model_simulation}
\title{Simulation generation attribute class}
\description{
This class generate mass simulation of communities
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-model_simulation-new}{\code{model_simulation$new()}}
\item \href{#method-model_simulation-simulate}{\code{model_simulation$simulate()}}
\item \href{#method-model_simulation-new_data}{\code{model_simulation$new_data()}}
\item \href{#method-model_simulation-graph}{\code{model_simulation$graph()}}
\item \href{#method-model_simulation-clone}{\code{model_simulation$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-model_simulation-new"></a>}}
\if{latex}{\out{\hypertarget{method-model_simulation-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{model_simulation$new(
  count = 15,
  timeline = 0:15,
  neighborhood = "Moore 1",
  model = "local",
  local_pattern = NULL,
  global_pattern = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-model_simulation-simulate"></a>}}
\if{latex}{\out{\hypertarget{method-model_simulation-simulate}{}}}
\subsection{Method \code{simulate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{model_simulation$simulate(
  count = self$count,
  timeline = self$timeline,
  neighborhood = self$neighborhood,
  model = self$model,
  fashion = "matrix"
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-model_simulation-new_data"></a>}}
\if{latex}{\out{\hypertarget{method-model_simulation-new_data}{}}}
\subsection{Method \code{new_data()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{model_simulation$new_data(summary = "highest")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-model_simulation-graph"></a>}}
\if{latex}{\out{\hypertarget{method-model_simulation-graph}{}}}
\subsection{Method \code{graph()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{model_simulation$graph(graph = "highest")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-model_simulation-clone"></a>}}
\if{latex}{\out{\hypertarget{method-model_simulation-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{model_simulation$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
