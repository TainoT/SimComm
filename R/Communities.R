#' Community Attributes Class
#'
#' @description
#' The class generate a community of a large number of species.
#'
#' ### this part isn't complete yet
#' @docType class
#' @param J The abundance of individuals
#' @param S The number of species.
#'
#' For the community drawing with entropart
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
    J = NULL, # TODO use it to replace nx/ny
    nx = NULL,
    ny = NULL,
    S = NULL,
    Distribution = NULL,
    sd = NULL,
    prob = NULL,
    alpha = NULL,

    the_matrix = NULL,
    the_wmppp = NULL,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$set_values()
      print("instanciation has been called : community_param")
    },

    #' @description
    #' Set values of the community generation parameters (set to default if no arguments)
    #TODO set_values should be defaulting t
    #'
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
    },

    #' @description
    #' Show the values of the community generation parameters
    show_values = function() {
    values <- paste("Type :", class(self)[1],
                "\n> Landscape ::\nnx = ", self$nx,
                "| ny = ", self$ny,
                "\n> Community ::\nS = ", self$S,
                "| size = ", self$nx * self$ny,
                "| alpha = ", self$alpha,
                "\n> Simulation ::\nDistribution : ", self$Distribution,
                "| sd : ", self$sd,
                "| prob : ", self$prob)
    cat(values)
    },

    #' @description
    #' Draw a community with no adjectives
    draw_matrix = function(){
      the_community <- entropart::rCommunity(
      1,
      #bit64::as.integer64 is an option here if we want to handle bigger numbers
      size = 100 * self$nx * self$ny,
      S = self$S,
      Distribution = self$Distribution,
      sd = self$sd,
      prob = self$prob,
      alpha = self$alpha,
      CheckArguments = FALSE
    )
      # Names are numbers, find a way to make it better, red and white shouldnt become orange when one dies
      # while it affects graphics, verify affects on plots, thats more important
      spNames <- seq(length(the_community))
      # Make a matrix, thats the matrix we sending and that should be studied
      # TODO : scrap that out for the meta com, i think, check performances
      self$the_matrix <- NULL
      self$the_matrix <- matrix(
        sample(spNames, size = self$nx * self$ny, replace = TRUE,
               prob = the_community / sum(the_community)),
        nrow = self$ny,
        ncol = self$nx
      )
      # are we sure we need to add that name ? this is confusing rn
      # i can send the_matrix right away and let the next function deal with it
#-----class(self$the_matrix) <- c("draw_matrix", class(self$the_matrix))
      return(self$the_matrix)
    },

    #' @description
    #' Draw a wmppp community with no adjective
    draw_wmppp = function(){
      spNames <- seq(self$S)
      the_community <- SpatDiv::rSpCommunity(
        1,
        size = self$nx * self$ny,#100 * self$nx * self$ny,
        S = self$S,
        Distribution = self$Distribution,
        sd = self$sd,
        prob = self$prob,
        alpha = self$alpha,
        CheckArguments = FALSE
      )

      # TODO I need to SAMPLE the community from 100 * x * y to a size of 100
      self$the_wmppp <- the_community
      # self$the_wmpp <- matrix()
      return(self$the_wmppp)
      # self$the_wmppp <- the_community
    }


  )
)

#' Local community
#'
#' @description The class generate a local community of a large number of species.
#'
#' @docType class
#' @param death_rate The mortality rate of an individual.
#' @param fashion The fashion of the community (matrix or wmppp)
#'
#' @export
local_pc <- R6::R6Class("local_pc",
  inherit = community_param,
  private = list(

  ),
  public = list(
    #' @field death_rate The mortality rate of an individual.
    #' Default is `0.1`
    death_rate = 0.11,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function(fashion = "matrix",
        death_rate = self$death_rate){
      self$death_rate <- death_rate
      if (fashion == "matrix"){
        self$the_matrix <- self$make_local(fashion = fashion)
      }
      else if (fashion == "wmppp") {
        self$the_wmppp <- self$make_local(fashion = fashion)
      }
      else
        stop("No defined fashion")
      print("instanciation has been called : local_pc")
    },

    #' @description
    #' Local community default community drawing - add death/birth_rate in the matrix
    #TODO : Verify if defaulting values in argument is right
    #              I could do an option and make a set_value ?
    #             strengthen the code, verify if the rates are on same level IN the matrix
    make_local = function(
    nx = 10,
    ny = nx,
    S = 10,
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    theta = 40,
    death_rate = self$death_rate,
    fashion = "matrix") {
      # We're NOT drawing the matrix
      self$set_values(nx = nx, ny = ny, S = S,
                      Distribution = Distribution, sd = sd,
                      prob = prob, alpha = theta)

      if (isTRUE(fashion == "matrix")){
        self$the_matrix <- NULL
        self$the_matrix <- self$draw_matrix()
        return(self$the_matrix)
      }
      else if (isTRUE(fashion == "wmppp")){
        self$the_wmppp <- NULL
        self$the_wmppp <- self$draw_wmppp()
        return(self$the_wmppp)
      }
      else
        stop("No defined fashion")
     }
  )
)

# TODO for 19/04, keep working on making sure this code is rock solid, then open up to meta community
# with options and all that

#' Meta community
#'
#' @description The class generate a meta community of a large number of species.
#'
#' @docType class
#' @param migration_rate The migration rate of an individual.
#' @param speciation_rate The speciation rate of an individual.
#'
#' @export
meta_pc <- R6::R6Class("meta_pc",
  inherit = community_param,
  private = list(

  ),
  public = list(
    #' @field migration_rate The migration of an individual.
    #' Default is `0.005`
    migration_rate = 0.0055,
    #' @field speciation_rate The speciation rate of an individual
    #' Default is `0.00001`
    speciation_rate = 0.000011,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function(
        migration_rate = self$migration_rate,
        speciation_rate = self$speciation_rate
        ){
      self$migration_rate <- migration_rate
      self$speciation_rate <- speciation_rate
      self$the_matrix <- self$make_meta()
      print("instanciation has been called : meta_pc")
    },

    #' @description
    #' Meta community default community drawing - add migr/spec_rate in the matrix
    make_meta = function(
    nx = 100,
    ny = nx,
    S = 100, #TODO : use Fisher's to get a proper N number here
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    theta = 40,
    migration_rate = self$migration_rate,
    speciation_rate = self$speciation_rate) {
      if (nx > 4000 || ny > 4000) {
        nx = 4000
        ny = nx
        warning("Too large number, defaulted back to the max allowed")
      }
      self$set_values(nx = nx, ny = ny, S = S,
                      Distribution = Distribution, sd = sd,
                      prob = prob, alpha = theta)
      self$the_matrix <- NULL
      self$the_matrix <- self$draw_matrix()

      #TODO remove the first str from the inherit but is that necessary to have that ?
#-----class(self$the_matrix) <- c("make_meta", class(self$the_matrix), migration_rate)

      return(self$the_matrix)
    }
    )
)
