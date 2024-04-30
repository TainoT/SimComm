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

      # Change cells
      if (runif(1, 0, 1) > 1 - self$speciation_rate) {
        self$pattern$marks$PointType <- self$pattern$marks$PointType['self$new_sp']
        self$new_sp <- self$new_sp + 1
        # unclear in bracket, study it up, we have more species, so spnames is bigger
        # its not spXXXX anymore, we can handle this right
        }
      else if (runif(1, 0, 1) < self$migration_rate) {
        # Sample from the meta community matrix
        self$pattern$marks$PointType <- self$meta_cp$the_points$marks$PointType[chosen_neighbors]
        # or somehting like that, since pattern is local, we jsut need to look for the class creating our pattern ?
      }
      else if (runif(1, 0, 1) < self$death_rate - (self$speciation_rate + self$migration_rate))
        self$pattern$marks$PointType <- self$pattern$marks$PointType[chosen_neighbors]
        # that one is fine, its just above
      else
        # self$pattern[row, col] <- self$pattern[row, col]


      if(save){
        #Save the new pattern
        self$run_patterns[[which(self$timeline == time)]] <- self$pattern
      }
    },
    disturb = function(time, save){
      # keep a counter, we keep deleting cells til its over, we ought to track that
      # each time we're onto a cell we roll a dice, if its less than 1/J, we can kill the cell
      # so even if we dont kill the sufficient amount of cells on the first loop of the entire area
      # we still can keep going and start all over again

      self$pattern$marks$PointType <- NULL # or 0 ?

      if (save){
        # Save the new pattern
        self$run_patterns[[which(self$timeline == time)]] <- self$pattern
      }
    }

    ),

    # As per the UNTB, there's a component of distrubance to be applied upon the community
    # the Disturbance rate is removing D individuals in the community
    # on the next cycle, the dead individuals are replaced by new ones
    # the new individuals can be drawn from nearby cells or from the meta community


  # ),

  #removing the commentaries here while I do the porting
  public = list(
    death_rate = 0.1,
    migration_rate = 0.005,
    speciation_rate = 0.001,
    local_cp = NULL,
    meta_cp = NULL,
    new_sp = 1,

    disturbance_rate = 10,

    #' @field n_neighbors The number of nearest neighbors to take into account.
    n_neighbors = NULL,

    initialize = function(
        pattern = NULL,
        timeline = 0,
        type = "Species",
        n_neighbors = 6) {
      super$initialize(
        pattern = pattern,
        timeline = timeline,
        type = type)
      self$n_neighbors <- n_neighbors
      # Store the neighbors
      private$neighbors <- self$neighbors_n(self$n_neighbors)

      # Change local_cp to include wmppp
      self$local_cp <- local_pc$new(death_rate = self$death_rate,
                                    draw = TRUE, fashion = "wmppp")
      self$pattern <- self$local_cm$the_matrix
      print("local_cm works with wmppp")

      #Meta to do later on
      # self$meta_cm <- meta_pc$new(migration_rate = self$migration_rate)
      # print("meta_cm works")
    }
 )
)
