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
        pattern = SpatDiv::rSpCommunity(n = 1, size = 100, CheckArguments = FALSE),
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

      # # cp drift is picking one nhood at random, but thats where we're
      # # rolling with our algo so it should look like
      # if (runif(1, 0, 1) < self$death_rate)
      #   chosen_neighbors <- apply(private$neigbors, 1, sample, size = 1)
      # # that look like a function tho, sorta, that does not make sense
      #
      # #how about :
      # if (runif(1, 0, 1) < self$death_rate)
      #   self$pattern$marks$PointType <- self$pattern$marks$PointType[chosen_neighbors]
      # # so what im doing here is
      # #if the death check is okay--->
      # # ----> pick from one of the neighbors thanks to the call of chosen_neighbors

      # so removing the double 'for', we have
      if (runif(1, 0, 1) > 1 - self$speciation_rate) {
        # self$pattern[row, col] <- self$pattern[row, col] + self$new_sp
        #     should look like ....
        # actually
        self$pattern$marks$PointType <- self$pattern$marks$PointType[]
        self$new_sp <- self$new_sp + 1
        # nothing change
      }
      else if (runif(1, 0, 1) < self$migration_rate) {
        # Sample from the meta community matrix
        #       hould that look like ...
        self$pattern$marks$PointType <- self$meta_cp$the_points$marks$PointType[chosen_neighbors]
        # or somehting like that, since pattern is local, we jsut need to look for the class creating our pattern ?
      }
      # death rate is removed from the spc and migration rate as we already assume death here
      else if (runif(1, 0, 1) < self$death_rate - (self$speciation_rate + self$migration_rate))
        self$pattern$marks$PointType <- self$pattern$marks$PointType[chosen_neighbors]
        # that one is fine, its just above
      else
        # self$pattern[row, col] <- self$pattern[row, col]


      if(save){
        #Save the new pattern
        self$run_patterns[[which(self$timeline == time)]] <- self$pattern
      }
    }
  ),
  #removing the commentaries here while I do the port
  public = list(
    death_rate = 0.1,
    migration_rate = 0.005,
    speciation_rate = 0.001,
    local_cm = NULL,
    meta_cm = NULL,
    new_sp = 1,

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
      self$local_cp <- local_pc$new(death_rate = self$death_rate, draw = TRUE)
      self$pattern <- self$local_cm$the_matrix
      print("local_cm works")

      #Meta to do later on
      # self$meta_cm <- meta_pc$new(migration_rate = self$migration_rate)
      # print("meta_cm works")
    }
 )
)
