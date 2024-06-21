#' Community Drift
#'
#' A [community_spcmodel].
#' At each generation, each individual is replaced by one of its neighbors.
#'
#'
#' @param n_neighbors Number of neighbors.
#' @param pattern The pattern which describes the location of agents.
#' @param time The point of the timeline considered.
#' Its value should be in `timeline`.
#' @param type The type of individuals. Informational only.
#' @param timeline A numeric vector that contains the points of time of interest.
#' @export
cp_drift <- R6::R6Class("cp_drift",
  inherit = community_spcmodel ,
  private = list(
    neighbors = NULL,
    evolve =  function(time, save) {
      # Choose one neighbor for each point
      chosen_neighbors <- apply(private$neighbors, 1, sample, size=1)
      # Change the types
      self$pattern$marks$PointType <- self$pattern$marks$PointType[chosen_neighbors]

      if(save) {
        # Save the new pattern
        self$run_patterns[[which(self$timeline == time)]] <- self$pattern
      }
    }
  ),

  public = list (
    #' @field n_neighbors The number of nearest neighbors to take into account.
    n_neighbors = NULL,

    #' @description
    #' Create a new object.
    initialize = function(
        pattern = SpatDiv::rSpCommunity(n = 1, size = 2, CheckArguments = FALSE),
        timeline = 0,
        type = "Species",
        n_neighbors = 6) {
      super$initialize(pattern = pattern, timeline = timeline, type = type)
      self$n_neighbors <- n_neighbors
      # Store the neighbors
      private$neighbors <- self$neighbors_n(self$n_neighbors)
    }
  )
)

#TODO redo that documentation properly
#' Hubbell's Neutral Theory Integration Class
#'
#' A [community_scpmodel] where each point contains an individual of one species.
#'
#' @export
cp_hubbell <- R6::R6Class("cp_hubbell",
  inherit = community_spcmodel,
  private = list(
    neighbors = NULL,
    evolve = function(time, save) {

      print("in")
      chosen_neighbors <- apply(private$neighbors, 1, sample, size=1)

      if (isTRUE(self$model == "local"))
        self$migration_rate <- 0
      r <- runif(1, 0, 1)
      # Migration happens if both communities are present
      if ((self$model == "default")
          && r < self$migration_rate) {
        self$pattern$marks$PointType <- self$meta_cp$the_points$marks$PointType[chosen_neighbors]
      }
      # Death happens only in the local community, substraced from the migration rate
      else if ((self$model == "default" || self$model == "local")
                 && r < self$death_rate - self$migration_rate) {
        self$pattern$marks$PointType <- self$pattern$marks$PointType[chosen_neighbors]
      }
      else
        self$pattern$marks$PointType <- self$pattern$marks$PointType
      # Speciation happens only in the meta community
      if ((self$model == "default" || self$model == "meta")
          && r < self$speciation_rate) {
        #TODO
        # add a new species to the meta community
        print("Speciation")
      }

      if(save){
        #Save the new pattern
        self$run_patterns[[which(self$timeline == time)]] <- self$pattern
      }    },

    #' @description
    #' disturb the community at given time and replace the dead cells with the neighbours
    disturb = function(time, save) {

      # Store community in pattern
      self$pattern <- self$saved_pattern(time)
      # self$pattern$marks$PointType <- NULL # or 0 ?

      self$pattern$marks$PointType[sample(1:length(self$pattern), self$kill_rate)] <- 0

      # change cells to fill the 0 with the neighbors
      chosen_neighbors <- apply(private$neighbors, 1, sample, size=1)
      if (isTRUE(self$pattern$marks$PointType == 0))
        self$pattern$marks$PointType <- self$pattern$marks$PointType[chosen_neighbors]

      if (save){
        # Save the new pattern
        self$run_patterns[[which(self$timeline == time)]] <- self$pattern
      }
    }

  ),
  public = list(
    #' @field death_rate The mortality rate of an individual.
    #' Default is `0.1`
    death_rate = 0.1,
    #' @field kill_rate The disturbance rate of the community.
    #' Default is `0`
    kill_rate = 0,
    #' @field migration_rate The migration rate of an individual.
    #' Default is '0.005'
    migration_rate = 0.005,
    #' @field speciation_rate The speciation rate of an individual.
    #' Default is '0.001'
    speciation_rate = 0.001,
    #' @field model The wmppp local community.
    local_cp = NULL,
    #' @field model The wmppp meta community.
    meta_cp = NULL,
    #' @field model The model to use
    model = "default",
    new_sp = 1,

    #' @field n_neighbors The number of nearest neighbors to take into account.
    n_neighbors = NULL,

    initialize = function(
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
        speciation_rate = self$speciation_rate) {
      if (is.null(local_pattern)) {
        super$initialize(
          pattern = local_pattern,
          timeline = timeline,
          type = type)
      } else {
        super$initialize(
          pattern = local_pattern$the_wmppp,
          timeline = timeline,
          type = type)
      }
      self$model <- model
      self$death_rate <- death_rate
      self$kill_rate <- disturbance_rate
      self$event <- event
      self$migration_rate <- migration_rate
      self$speciation_rate <- speciation_rate

      if (isTRUE(self$model == "local")) {
        self_local_cp <- local_pc$new(death_rate = self$death_rate,
                                      fashion = "wmppp")
      } else if (isTRUE(self$model == "meta")) {
        self$meta_cp <- meta_pc$new(migration_rate = self$migration_rate,
                                    speciation_rate = self$speciation_rate)
      } else {
        if (isTRUE(self$model == "default")) {
          self$local_cp <- local_pc$new(death_rate = self$death_rate,
                                        fashion = "wmppp")
          self$meta_cp <- meta_pc$new(migration_rate = self$migration_rate,
                                      speciation_rate = self$speciation_rate)
        } else {
          stop("Model not recognized. Please use 'local', 'meta' or 'default'.")
        }
      }

      # We store the result of that community in pattern for further control
      self$pattern <- self$local_cp$the_wmppp

      #' @description
      #' Show the values of the model
      show_values = function(values = "basic") {
        if (values == "basic") {
          print(paste("Model: ", self$model))
          print(paste("Death rate: ", self$death_rate))
          print(paste("Migration rate: ", self$migration_rate))
          print(paste("Speciation rate: ", self$speciation_rate))
        }
        else if (values == "all") {
          if (is.null(self$local_cp) == FALSE)
            self$local_cp$show_values()
          if (is.null(self$meta_cp) == FALSE)
            self$meta_cp$show_values()
        }
        else
          warning("Values not recognized. Please use 'basic' or 'all'.")
      }
    }
  )
)


# }
#       # print(self$pattern)
#       self$n_neighbors <- n_neighbors
#       # Store the neighbors
#       private$neighbors <- self$neighbors_n(self$n_neighbors)
#
#       # We're creating a local point community
#
#       self$local_cp <- local_pc$new(death_rate = self$death_rate, fashion = "wmppp")
#
#       # self$pattern <- self$local_cp$the_wmppp
#
#       print("local_pc works with wmppp")
#
#       #Meta to do later on
#       # self$meta_cm <- meta_pc$new(migration_rate = self$migration_rate)
#       # print("meta_cm works")
#     }
#  )
# )
