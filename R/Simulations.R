#' Simulation generation attribute class
#'
#' @description
#' This class generate mass simulation of communities
#'
#' @docType class
#'
#' @param time integer number of time steps
#' @param n integer number of communities
#'
#' @export

model_simulation <- R6::R6Class("model_simulation",
  private = list(
    sim_loop = function(run, timeline, neighborhood, model,
                       fashion = "matrix", save = T) {
      if (fashion == "matrix")
        running <- SimComm::cm_hubbell$new(timeline = timeline,
                              neighborhood = neighborhood, model = model)
      else if (fashion == "wmppp")
        running <- SimComm::cp_hubbell$new(timeline = timeline,
                              neighborhood = neighborhood, model = model)
      running$run(save = save)
      df <- running$graph(save = save)
      df$simulation <- run
      print("1")
      return(df)
    }
  ),
  public = list(
    count = NULL,
    # run = NULL,
    timeline = NULL,
    neighborhood = NULL,
    model = NULL,
    data_sim = NULL,
    data_synth = NULL,
    data_graph = list(),

    initialize = function(
        count = 15,
        timeline = 0:15,
        neighborhood = "Moore 1",
        model = "local") {
      self$count <- count
      self$timeline <- timeline
      self$neighborhood <- neighborhood
      self$model <- model
      self$data_sim <- vector("list", length = count)
    },

    simulate = function(count = self$count, timeline = self$timeline,
                        neighborhood = self$neighborhood, model = self$model,
                        fashion = "matrix") {
      for (i in 1:self$count)
        self$data_sim[[i]] <- private$sim_loop(run = i, timeline = timeline,
                                               neighborhood = neighborhood,
                                               model = model, save = T,
                                               fashion = fashion)

      self$data_synth <- bind_rows(self$data_sim)
    },

    new_data = function(summary = "highest") {
      if (summary == "highest") {
        highest_ranked <- self$data_synth %>%
          group_by(simulation, species) %>%
          top_n(1, desc(rank)) %>%
          ungroup()
        mean_abundance <- highest_ranked %>%
          group_by(time) %>%
          summarise(mean_abundance = mean(count))
        # store highest_ranked and mean_abundance in a named list called self$data_list
        self$data_graph <- list(highest_ranked, mean_abundance)
        names(self$data_graph) <- c("highest_ranked", "mean_abundance_highest")
      }
    },

    graph = function(graph = "highest") {
      if (graph == "highest") {
        gg <- ggplot() +
          geom_line(data = self$data_graph$highest_ranked,
                    aes(x = time, y = count, group = interaction(simulation, species)),
                    alpha = 0.5, color = "grey") +
          geom_line(data = self$data_graph$mean_abundance_highest,
                    aes(x = time, y = mean_abundance, group = 1),
                    color = "red") +
          labs(title = "Highest ranked species",
               x = "Time",
               y = "Abundance") +
          # facet_wrap(~simulation) +
          theme_minimal()
        print(gg)
      }
    }
  )
)
