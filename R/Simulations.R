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
        running <- cm_hubbell$new(timeline = timeline,
                              neighborhood = neighborhood, model = model,
                              local_pattern = self$local_pattern,
                              global_pattern = self$global_pattern)
      else if (fashion == "wmppp")
        running <- cp_hubbell$new(timeline = timeline,
                              neighborhood = neighborhood, model = model,
                              local_pattern = self$local_pattern,
                              global_pattern = self$global_pattern)
      running$run(save = save)
      df <- running$synthesis(save = save, calc_neighbors = F)
      df$simulation <- run
      print(paste("Simulation", run, "completed"))
      return(df)
    }
  ),
  public = list(
    count = NULL,
    timeline = NULL,
    neighborhood = NULL,
    model = NULL,
    local_pattern = NULL,
    global_pattern = NULL,
    data_sim = NULL,
    data_synth = NULL,
    data_graph = list(),

    initialize = function(
        count = 15,
        timeline = 0:15,
        neighborhood = "Moore 1",
        model = "local",
        local_pattern = NULL,
        global_pattern = NULL) {
      self$count <- count
      self$timeline <- timeline
      self$neighborhood <- neighborhood
      self$model <- model
      self$data_sim <- vector("list", length = count)

      if (is.null(local_pattern))
        self$local_pattern <- local_pc$new()
      else
        self$local_pattern <- local_pattern
      if (is.null(global_pattern))
        self$global_pattern <- meta_pc$new()
      else
        self$global_pattern <- global_pattern
    },

    simulate = function(count = self$count, timeline = self$timeline,
                        neighborhood = self$neighborhood, model = self$model,
                        fashion = "matrix") {
      self$count <- count
      self$timeline <- timeline
      self$neighborhood <- neighborhood
      self$model <- model

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
          filter(rank == 1) %>%
          ungroup()
        mean_abundance <- highest_ranked %>%
          group_by(time) %>%
          summarise(mean_abundance = mean(count))
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
          geom_smooth(data = self$data_graph$mean_abundance_highest,
                      aes(x = time, y = mean_abundance),
                      method = "lm", se = FALSE, color = "blue", linewidth = 0.5) +
          labs(title = "Highest ranked species",
               x = "Time",
               y = "Abundance") +
          ylim(0, self$local_pattern$ny * self$local_pattern$nx) +
          theme_minimal()
        print(gg)
      }
    }
  )
)
