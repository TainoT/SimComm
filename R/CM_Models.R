#' Community Drift
#'
#' A [community_matrixmodel] where each cell contains an individual.
#' Marks are species.
#' At each generation, each individual is replaced by one of its neighbors.
#'
#' Edge effects are eliminated by a toroidal correction.
#'
#' @docType class
#' @export
cm_drift <- R6::R6Class("cm_drift",
  inherit = community_matrixmodel,
  private = list(
    evolve =  function(time, save) {
      # Prepare the buffer
      self$prepare_buffer()

      # Draw a neighbor
      for(row in seq(nrow(self$pattern))) {
        for(col in seq(ncol(self$pattern))) {
          self$pattern[row, col] <- sample(self$neighbors(row, col), size = 1)
        }
      }

      if(save) {
        # Save the new pattern
        self$run_patterns[, , which(self$timeline == time)] <- self$pattern
      }
    }
  )
)




#' Conway's game of life
#'
#' A [community_matrixmodel] where each cell contains or not an individual.
#' At each generation, an individual may survive or not and empty cells be filled by a new individual.
#'
#' The survival and generation rules are fixed by the number of neighbors of each cell.
#' Edge effects are eliminated by a toroidal correction.
#'
#' @docType class
#' @param neighborhood A character string defining what is the neighborhood of a cell:
#' "von Neumann 1" or "4" for the closest four neighbors (North, West, South, East);
#' "Moore 1" or "8" for all adjacent cells (the first four and North-West, etc.);
#' "Moore 2" or "24" for two rings of neighbors.
#' @param pattern The pattern which describes the location of agents.
#' @param time The point of the timeline considered.
#' Its value should be in `timeline`.
#' @param timeline A numeric vector that contains the points of time of interest.
#' @param type The type of individuals. Informational only.
#' @export
cm_Conway <- R6::R6Class("cm_Conway",
                         inherit = community_matrixmodel,
                         private = list(
                           evolve = function(time, save) {
      # Prepare the buffer
      self$prepare_buffer()

      # Change cells
      for(row in seq(nrow(self$pattern))) {
        for(col in seq(ncol(self$pattern))) {
          # Count the neighbors
          n_neighbors <- sum(self$neighbors(row, col))
          # Apply the rule
          self$pattern[row, col] <- (self$pattern[row, col] & (n_neighbors %in% self$to_survive)) | (!self$pattern[row, col] & (n_neighbors %in% self$to_generate))
        }
      }

      if(save) {
        # Save the new pattern
        self$run_patterns[, , which(self$timeline == time)] <- self$pattern
      }
    }
  ),
  public = list(
    #' @field to_survive The number of neighbors necessary for an individual to survive.
    #' Default is `2:3`.
    to_survive = c(2, 3),
    #' @field to_generate The number of neighbors necessary for an empty cell to generate an individual.
    #' Default is 3.
    to_generate = 3,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function(
        pattern = pattern_matrix_individuals(),
        timeline = 0,
        type = "Alive",
        neighborhood = "Moore 1") {
      super$initialize(
        pattern = pattern,
        timeline = timeline,
        type = type,
        neighborhood = neighborhood
      )
    }
  )
)


#' Hubbell's Neutral Theory Integration Class
#'
#' A [community_matrixmodel] where each cell contains an individual of one species.
#' At each generation, an individual may survive or not and empty cells be filled with a new individual.
#'
#' The survival and generation rules are fixed by the mortality and reproduction rate of each individual.
#' Edge effects are eliminated by a toroidal correction. (??? check that out)
#'
#' @docType class
#' @param pattern The pattern which describe the location of agents.
#' @param time The point of the timeline considered.
#' Its value should be in `timeline`
#' @param timeline A numeric vector that contains the points of time of interest.
#' @param type The type of individuals. Information only.
#' @export
cm_hubbell <- R6::R6Class("cm_hubbell",
  inherit = community_matrixmodel,
  private = list(
    evolve = function(time, save) {
      # Prepare the buffer
      self$prepare_buffer()

      # Change cells
      for(row in seq(nrow(self$pattern))) {
        for(col in seq(ncol(self$pattern))) {
          # Roll a dice between 0 and 1, if the roll is lower than drate, then run it
          if (isFALSE(runif(1, 0, 1) < self$death_rate))
            self$pattern[row, col] <- sample(self$neighbors(row, col), size = 1)
          else
            self$pattern[row, col] <- self$pattern[row, col]

          # Maybe too simple, does not count in the death state
          # death state = empty cell (perturbation, dead tree, etc)
        }
      }

      if(save){
        #Save the new pattern
        self$run_patterns[, , which(self$timeline == time)] <- self$pattern
      }
    }
  ),
  public = list(
    #' @field death_rate The mortality rate of an individual.
    #' Default is `0.1`
    death_rate = 0.1,
    #' @field birth_rate The reproduction rate of an individual.
    #' Default is `0.2`
    birth_rate = 0.2,
    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function(
        pattern = NULL,
        timeline = 0,
        type = "Species",
        neighborhood = "Moore 1") {
      super$initialize(
        pattern = pattern,
        timeline = timeline,
        type = type,
        neighborhood = neighborhood
      )#,
      # I can add new functions on initialization here, including the meta com
      # its resetting the thing, meaning im not even taking in the object that exist
      # how about i could do a new() on the model, that take in a previously created community class? 
      # basically, I would have to do a class = NULL in the list, this way, i would be able to access the same community elsewhere, without that silly defaulting 
      commu <- local_pc$new(death_rate = self$death_rate, birth_rate = self$birth_rate)
      if(is.null(self$pattern))
      {
        self$pattern <- commu$local_matrix_class()
        print("sum")
      }
      else
        warning("Something broke, no pattern returned")
    },

    #' #' @description
    #' #' Change the demographic rates
    #' #' `death_rate`, `birth_rate` and `migration_rate`(unused for now)
    #' set_rate = function(
    #' death_rate = 0.1,
    #' birth_rate = 0.2) {
    #'   death_rate = death_rate
    #'   birth_rate = birth_rate
    #' },

    #' @description
    #' Draw the abundance of each species over time
    plot_line = function() {
      plot(self$along_time(entropart::Richness, Correction = "None", CheckArguments = F), type = "l")
      #the_plotline <- data.frame(self$timeline, self$saved_pattern(self$timeline))

      # ggplot(the_plotline, aes(x = self$timeline, y = self$saved_pattern(timeline))) +
      #   geom_line()
      }

    #' @description
    #' call for community ?
 )
)
