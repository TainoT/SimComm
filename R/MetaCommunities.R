#' Meta Community Class
#'
#' @description
#' The class generate a community of a large number of species.
#'
#' @docType class
#' @param J The abundance of individuals
#' @param S The number of species.
#' @param Distribution The distribution of species frequencies. May be "lnorm" (log-normal), "lseries" (log-series), "geom" (geometric) or "bstick" (broken stick).
#' @param sd The simulated distribution standard deviation. For the log-normal distribution, this is the standard deviation on the log scale.
#' @param prob The proportion of ressources taken by successive species in the geometric model.
#' @param alpha Fisher's alpha in the log-series model.
#'
#' @param abundance The amount of individual in the community (community size)
#' @param death_rate The mortality rate of individuals
#' @param birth_rate The reproduction rate of individuals
#' @param migration_rate The migration rate of individuals from the meta community to the local community
#' @param speciation_rate The speciation rate of individual in the local community
#' @param community_type The type of community (local or meta)
#' @param theta The compound value of community size and speciation rate
#'
#' @export
community_param <- R6::R6Class("community_param",
  private = list(
  ),
  public = list(
    S = 300,
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    alpha = 40

  )
)

#' Local
#'
#'
local_pc <- R6::R6Class("local_pc",
  inherit = community_param,
  private = list(

  ),
  public = list(
    #' @field death_rate The mortality rate of an individual.
    #' Default is `0.1`
    death_rate = 0.1,
    #' @field birth_rate The reproduction rate of an individual.
    #' Default is `0.2`
    birth_rate = 0.2
    #' @field speciation_rate The speciation rate of an individual
    #' Default is `0.00001`
  )
)

local_pc <- R6::R6Class("local_pc",
  inherit = community_param,
  private = list(

  ),
  public = list(
    #' @field migration_rate The migration rate of a species of a meta community
    #' Default is `0.005`
    migration_rate = 0.005

  )
)
