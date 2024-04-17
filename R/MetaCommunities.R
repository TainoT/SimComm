#' Meta Community Class
#'
#' @description
#' The class generate a community of a large number of species.
#'
#' ### this part isn't complete yet
#' @docType class
#' @param J The abundance of individuals
#' @param S The number of species.
#'
#' @param nx The size of horizontal axis
#' @param ny The size of vertical axis
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
    nx = NULL,
    ny = NULL,
    S = NULL,
    Distribution = NULL,
    sd = NULL,
    prob = NULL,
    alpha = NULL,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$set_values()
      print("instanciation has been called : community_param")
    },

    #' @description
    #' Set values of the community generation parameters (set to default if no arguments)
##### it should probably be on the cm_hubbell -> actually no, user can call default community
##### but for finer details, they're free to create a whole community object and set its values
    set_values = function(
    nx = 20,
    ny = nx,
    S = 2,
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    alpha = 40) {
      self$nx = nx
      self$ny = ny
      self$S = S
      self$Distribution = Distribution
      self$sd = sd
      self$prob = prob
      self$alpha = alpha
      if (nx == 20)
        print("values are defaulted : community_param")
      else
        print("new values are in : community_param")
    },

    #' @description
    #' Draw a community with no adjectives
    pattern_matrix_class = function(){
      the_community <- entropart::rCommunity(
      1,
      size = 100 * self$nx * self$ny,
      S = self$S,
      Distribution = self$Distribution,
      sd = self$sd,
      prob = self$prob,
      alpha = self$alpha,
      CheckArguments = FALSE
    )
      # Names are numbers, find a way to make it better, red and white shouldnt become orange when one dies
      spNames <- seq(length(the_community))
      # Make a matrix, thats the matrix we sending and that should be studied
      the_matrix <- matrix(
        sample(spNames, size = self$nx * self$ny, replace = TRUE,
               prob = the_community / sum(the_community)),
        nrow = self$ny,
        ncol = self$nx
      )
      # are we sure we need to add that name ?
      # i can send the_matrix right away and let the next function deal with it
      class(the_matrix) <- c("pattern_matrix_class", class(the_matrix))
      print("community drawn : community_param")
      return(the_matrix)
    }

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
    death_rate = 0.11,
    #' @field birth_rate The reproduction rate of an individual.
    #' Default is `0.2`
    birth_rate = 0.22,
    #' @field speciation_rate The speciation rate of an individual
    #' Default is `0.00001`

    #' keep the value as default, its not great if its NULL
    #' but allow modification
    #' maybe create a brand new set_values (maybe not overwriting) for ecological values
    initialize = function(death_rate = self$death_rate, birth_rate = self$birth_rate){
      self$death_rate <- death_rate
      self$birth_rate <- birth_rate
      super$set_values()
      print("instanciation has been called : local_pc")
    },

    #' @description
    #' Local community default community drawing
    local_matrix_class = function(
    nx = 10,
    ny = nx,
    S = 2,
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    theta = 40,
    death_rate = 0.1,
    birth_rate = 0.2) {
      # not ideal, im calling the superfunction, but i gotta actually rewrite it because it's
      # using brand new values
      # remember : the superclass is only creating a community
      # in the local, we're specifying default values closer to a local community
      # and adding death/birth rate
      the_matrix <- super$pattern_matrix_class()
      # remove the first str from the inherit
      class(the_matrix) <- c("local_matrix_class", class(the_matrix), death_rate, birth_rate)
      print("community drawn : local_pc")
      return(the_matrix)
    }
  )
)


# TODO for 18/04, keep working on making sure this code is rock solid, then open up to meta community
# with options and all that


#' meta_pc <- R6::R6Class("local_pc",
#'   inherit = community_param,
#'   private = list(
#'
#'   ),
#'   public = list(
#'     nx = 10,
#'     ny = nx,
#'     S = 2,
#'     Distribution = "lnorm",
#'     sd = 1,
#'     prob = 0.1,
#'     alpha = 40,
#'
#'     #' @field migration_rate The migration rate of a species of a meta community
#'     #' Default is `0.005`
#'     migration_rate = 0.005
#'     #' @field theta The diversity number of the meta community
#'
#'
#'   )
#' )
