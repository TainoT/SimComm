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
    #' TODO add an option to change the rules on the fly ?
    #' Realm of Cellular Automata here
    evolve = function(time, save) {
      # Prepare the buffer
      self$prepare_buffer()

      # Change cells

          ## RULE FOR DEMOGRAPHIC DRIFT
      for(row in seq(nrow(self$pattern))) {
        for(col in seq(ncol(self$pattern))) {
          # Roll a dice between 0 and 1, if the roll is lower than drate, then run it
          if (isFALSE(runif(1, 0, 1) < self$death_rate))
            self$pattern[row, col] <- sample(self$neighbors(row, col), size = 1)
          else
            self$pattern[row, col] <- self$pattern[row, col]

          ## RULE TO ADD FOR SPECIATION
          #'
          #' check if alive, we cant have a death and speciation at the same time
          #' when its dead, let the neighbour colonize the space
          #' so for 0.1 --> check on daily 10


          ## RULE TO ADD FOR MIGRATION

          # Maybe too simple, does not count in the death state, NOT A PRIORITY
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
    death_rate = 0.9,
    #' @field birth_rate The reproduction rate of an individual.
    #' Default is `0.2`
    birth_rate = 0.2,
    #' @field community The local community
    local_cm = NULL,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    #' What it does : default the values
    #'                check if pattern is empty (no argument given)
    #'                else
    #'                create new community (local_cm), by default, draw it
    #'                send the resulting matrix to pattern
    #' TODO
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
      )
      # TODO N/A : add whatever I want here when I create a neutral theory model
      #            meaning : the model that will plot, the local community, the meta community

      # We're creating the local community model
      self$local_cm <- local_pc$new(death_rate = self$death_rate, birth_rate = self$birth_rate, draw = TRUE)
      print("its ok here")

      self$pattern <- self$local_cm$the_matrix
      print("its ok here too")
      # TODO Check if there's no given pattern first - Generally from a community created prior the model
      # if(is.null(self$pattern) == FALSE)
      # {
      #   ## check typing, we need something that fit $pattern and send an error otherwise
      #   ## if.matrix == TRUE maybe
      #   self$pattern <- pattern
      #   print("pattern set from outside model")
      # }
    },

    #' @description
    #' Redraw the model with new values of communities
    #' TODO : it should call set_values and draw_matrix
    #'        call make_local(draw = T) if the option is available
    #'
    #' TODO : add option to modify value of model ?
    # redraw = function(
    #   nx = 20,
    #   ny = nx,
    #   S = 2,
    #   Distribution = "lnorm",
    #   sd = 1,
    #   prob = 0.1,
    #   alpha = 40) {
    #   self$local_cm$set_values(
    #     nx = self$local_cm$nx,
    #     ny = self$local_cm$ny,
    #     S = self$local_cm$S,
    #     Distribution = self$local_cm$Distribution,
    #     sd = self$local_cm$sd,
    #     prob = self$local_cm$prob,
    #     alpha = self$local_cm$alpha)
    #   self$pattern <- self$local_cm$the_matrix
    #   },

    ## TODO : delete ?
    #' #' @description
    #' #' Change the demographic rates
    #' #' `death_rate`, `birth_rate` and `migration_rate`(unused for now)
    #' set_rate = function(
    #' death_rate = 0.1,
    #' birth_rate = 0.2) {
    #'   death_rate = death_rate
    #'   birth_rate = birth_rate
    #' },


    ## TODO NA : once all is clean, try it
    ## Not the priority, first finish classes, then meta community
    ## Only then, we can explore how to display results
    #' @description
    #' Draw the abundance of each species over time
    plot_line = function() {
      plot(self$along_time(entropart::Richness, Correction = "None", CheckArguments = F), type = "l")
      #the_plotline <- data.frame(self$timeline, self$saved_pattern(self$timeline))

      # ggplot(the_plotline, aes(x = self$timeline, y = self$saved_pattern(timeline))) +
      #   geom_line()
      }
 )
)
