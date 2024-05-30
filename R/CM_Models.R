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
          self$pattern[row, col] <- (self$pattern[row, col] &
                                       (n_neighbors %in% self$to_survive)) |
            (!self$pattern[row, col] & (n_neighbors %in% self$to_generate))
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
#' @param neighborhood A character string defining what is the neighborhood of a cell:
#' "von Neumann 1" or "4" for the closest four neighbors (North, West, South, East);
#' "Moore 1 or 8" for all adjacent cells (the first four and North-West, etc.);
#' "Moore 2 or 24" for two rings of neighbors.
#' "Global or 0" for the whole matrix as neighborhood.
#' @param pattern The pattern which describe the location of agents.
#' @param time The point of the timeline considered.
#' Its value should be in `timeline`
#' @param timeline A numeric vector that contains the points of time of interest.
#' @param type The type of individuals. Information only.
#'
#'
#' @export
cm_hubbell <- R6::R6Class("cm_hubbell",
  inherit = community_matrixmodel,
  private = list(
    #' @description
    #' community evolution rules for the neutral theory model
    evolve = function(time, save) {
      # Prepare the buffer
      self$prepare_buffer()

      # Prep the model
      if (isTRUE(self$model == "local"))
        self$migration_rate <- 0
      # Change cells
      for(row in seq(nrow(self$pattern))) {
        for(col in seq(ncol(self$pattern))) {
          r <- runif(1, 0, 1)
          # Roll a dice between 0 and 1, if the roll is lower than rate, then run it
          #TODO verify if meta community is not empty
          if ((self$model == "default")
              && r < self$migration_rate) {
            # Sample from the meta community matrix
            self$pattern[row, col] <- sample(self$meta_cm$the_matrix, size = 1)
          }
          # death rate is removed from the spc and migration rate as we already assume death here
          else if ((self$model == "default" || self$model == "local")
                   && r < self$death_rate - self$migration_rate) {
            self$pattern[row, col] <- sample(self$neighbors(row, col), size = 1)
          }
          else
            self$pattern[row, col] <- self$pattern[row, col]

          if ((self$model == "default" || self$model == "meta")
              && r < self$speciation_rate) {
            #TODO
            # add a new species to the meta community
            print("Speciation")
          }
        }
      }

      if(save){
        #Save the new pattern
        self$run_patterns[, , which(self$timeline == time)] <- self$pattern
      }
    },

    #' @description
    #' distrub the community at given time and replace the dead cells with the neighbors
    disturb = function(time, save) {
      # dead_c <- 0
      # Prepare the buffer
      self$prepare_buffer()

      # Store matrix in pattern
      self$pattern <- self$local_cm$the_matrix

      # Make n cell in pattern as NA
      # TODO change up 36 to something the user can modify
      # why is it 36 again ? is thhat the disturbance rat
      self$pattern[sample(1:length(self$pattern), 36)] <- 0
      # dead_c <- sample(self$pattern, size = 36)

      # Change cells to fill the 0 with the neighbors
      for (row in seq(nrow(self$pattern))) {
        for (col in seq(ncol(self$pattern))) {
          if (isTRUE(self$pattern[row, col] == 0)) {
            self$pattern[row, col] <- sample(self$neighbors(row, col), size = 1)
          }
        }
      }

      if(save) {
        # Save the new pattern
        self$run_patterns[, , which(self$timeline == time)] <- self$pattern
      }
    }
  ),
  public = list(
    #' @field death_rate The mortality rate of an individual.
    #' Default is `0.1`
    death_rate = 0.1,
    #' @field migration_rate The reproduction rate of an individual.
    #' Default is `0.005`
    migration_rate = 0.005,
    #' @field speciation_rate The reproduction rate of an individual.
    #' Default is `0.0001`
    speciation_rate = 0.0001,
    #' @field local_cm The matrix local community
    local_cm = NULL,
    #' @field meta_cm The matrix meta community
    meta_cm = NULL,
    #' @field model The model to use
    model = "default",

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function(
        local_pattern = NULL,
        timeline = 0,
        type = "Species",
        neighborhood = "Moore 1",
        global_pattern = NULL,
        model = model,
        death_rate = self$death_rate,
        migration_rate = self$migration_rate,
        speciation_rate = self$speciation_rate
        ) {
      super$initialize(
        pattern = local_pattern,
        timeline = timeline,
        type = type,
        neighborhood = neighborhood
      )
      self$model <- model
      self$death_rate <- death_rate
      self$migration_rate <- migration_rate
      self$speciation_rate <- speciation_rate
      # TODO : add whatever I want here when I create a neutral theory model
      #        - the local community
      #        - the meta community
      #        - the species count
      #        - the speciation rate
      #        - the migration rate
      #        - the death rate
      #        - the birth rate
      #        - the disturbance rate

      #TODO BY DEFAULT, no pattern is given, we create communities
      # We're creating the local community model
      if (isTRUE(self$model == "local")) {
        self$local_cm <- local_pc$new(death_rate = self$death_rate)
      } else if (isTRUE(self$model == "meta")) {
        self$meta_cm <- meta_pc$new(migration_rate = self$migration_rate,
                                    speciation_rate = self$speciation_rate)
      } else {
        if (isTRUE(self$model == "default")) {
          self$local_cm <- local_pc$new(death_rate = self$death_rate)
          self$meta_cm <- meta_pc$new(migration_rate = self$migration_rate,
                                      speciation_rate = self$speciation_rate)
        } else {
          stop("Model not recognized. Please use 'local', 'meta' or 'default'")
        }
      }

      # We store the result of that community in pattern for further control
      self$pattern <- self$local_cm$the_matrix

      # We're creating the meta community model
      #TODO if patterns are given, we use them

    },

    #' @description
    #' Show the values of the model
    show_values = function(values = "basic") {
      if (values == "basic") {
        print(paste("Model: ", self$model))
        print(paste("Death rate: ", self$death_rate))
        print(paste("Migration rate: ", self$migration_rate))
        print(paste("Speciation rate: ", self$speciation_rate))
      }
      else if (values == "all"){
        if (is.null(self$local_cm) == FALSE)
          self$local_cm$show_values()
        print(paste("\n====================\n"))
        if (is.null(self$meta_cm) == FALSE)
          self$meta_cm$show_values()
      }
      else
        warning("Values not recognized. Please use 'basic' or 'all'")
    },

    #' @description
    #' DRAFT
    #TODO
    count_species = function() {
      self$species_count <- table(unlist(self$run_patterns))
    }
)
)

