#' Community Attributes Class
#'
#' @description
#' The class generate a community of a large number of species.
#'
#' @docType class
#' @param S The number of species.
#' @param nx The size of horizontal axis
#' @param ny The size of vertical axis
#'
#' @param Distribution The distribution of species frequencies. May be "lnorm" (log-normal), "lseries" (log-series), "geom" (geometric) or "bstick" (broken stick).
#' @param sd The simulated distribution standard deviation. For the log-normal distribution, this is the standard deviation on the log scale.
#' @param prob The proportion of ressources taken by successive species in the geometric model.
#' @param alpha Fisher's alpha in the log-series model.
#'
#' @param distribute The type of distribution. May be "entropart" (entropart package), "random" or "uniform".
#' @param style The style of the community. May be "random", "checkerboard" or "uniform".
#'
#' @param the_matrix The matrix of the community.
#' @param the_wmppp The wmppp of the community.
#'
#' @param death_rate The mortality rate of individuals
#' @param migration_rate The migration rate of individuals from the meta community to the local community
#' @param speciation_rate The speciation rate of individual in the local community
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

    distribute = NULL,
    style = NULL,

    the_matrix = NULL,
    the_wmppp = NULL,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      self$set_values()
    },

    #' @description
    #' Set values of the community generation parameters (set to default if no arguments)
    set_values = function(
    nx = 20,
    ny = nx,
    S = 20,
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    alpha = 40,
    distribute = "entropart",
    style = "random") {
      self$nx = nx
      self$ny = ny
      self$S = S
      self$Distribution = Distribution
      self$sd = sd
      self$prob = prob
      self$alpha = alpha
      self$distribute = distribute
      self$style = style

      if (!is.null(self$the_matrix)){
        temp <- self$draw_matrix()
        self$the_matrix <- temp
      }
      else if (!is.null(self$the_wmppp)){
        temp <- self$draw_wmppp()
        self$the_wmppp <- temp
      }
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
    #' Draw a community
    draw_matrix = function(){
      if (self$distribute != "entropart" && self$distribute != "random" && self$distribute != "uniform") {
        self$distribute <- "entropart"
        warning("No defined distribution, defaulting to entropart")
      }

      if (self$distribute == "entropart") {
        the_community <- entropart::rCommunity(
          1,
          #bit64::as.integer64 is an option here if we want to handle bigger numbers
          size = 100 * self$nx * self$ny,
          S = self$S,
          Distribution = self$Distribution,
          sd = self$sd,
          prob = self$prob,
          alpha = self$alpha,
          CheckArguments = FALSE)
      }
      else if (self$distribute == "random") {
        the_community <- rnorm(self$S, mean = 100, sd = 10)
      }
      else if (self$distribute == "uniform") {
        indiv_species <- (self$nx * self$ny) / self$S
        if (indiv_species %% 1 != 0) {
          #numbers::divisors is an option to get the divisors
          stop("The number of species is not a multiple of the number of individuals")}
        species <- rep(1:self$S, each = indiv_species)
        the_community <- species
      }
      else
        stop("No defined distribution")

      spNames <- seq(length(the_community))
      self$the_matrix <- NULL
      if (self$style != "checkerboard" && self$style != "random" && self$style != "uniform") {
        self$style <- "random"
        warning("No defined style, defaulting to random")
      }
      if (self$style == "random") {
        spNames <- seq(length(the_community))
        self$the_matrix <- matrix(
          sample(spNames, size = self$nx * self$ny, replace = TRUE,
                 prob = the_community / sum(the_community)),
          nrow = self$ny,
          ncol = self$nx
        )
      }
      else if (self$style == "uniform")
        self$the_matrix <- matrix(
          the_community,
          # rep(sample(spNames, size = self$nx * self$ny, replace = TRUE),
          # each = self$nx * self$ny / self$S),
          nrow = self$ny,
          ncol = self$nx
        )
      else if (self$style == "checkerboard") {
        self$the_matrix <- matrix(nrow = self$ny, ncol = self$nx)
        for (i in 1:self$nx) {
          for (j in 1:self$ny) {
            individual = (i + j) %% self$S
            self$the_matrix[j, i] <- individual
          }
        }
      }
      return(self$the_matrix)
    },

    #' @description
    #' Draw a wmppp community with no adjective
    draw_wmppp = function(){
      spNames <- seq(self$S)
      the_community <- SpatDiv::rSpCommunity(
        1,
        size = self$nx * self$ny,
        S = self$S,
        Distribution = self$Distribution,
        sd = self$sd,
        prob = self$prob,
        alpha = self$alpha,
        CheckArguments = FALSE
      )
      self$the_wmppp <- the_community
      return(self$the_wmppp)
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
    death_rate = 0.1,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function(fashion = "matrix",
        death_rate = self$death_rate,
        distribute = "entropart",
        style = "random"){
      self$death_rate <- death_rate
      self$distribute <- distribute
      self$style <- style
      if (fashion == "matrix"){
        self$the_matrix <- self$make_local(fashion = fashion, distribute = distribute, style = style)
      }
      else if (fashion == "wmppp") {
        self$the_wmppp <- self$make_local(fashion = fashion)
      }
      else
        stop("No defined fashion")
    },

    #' @description
    #' Local community default community drawing - add death_rate
    make_local = function(
    nx = 10,
    ny = nx,
    S = 10,
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    alpha = 40,
    death_rate = self$death_rate,
    fashion = "matrix",
    distribute = "entropart",
    style = "random") {
      self$set_values(nx = nx, ny = ny, S = S,
                      Distribution = Distribution, sd = sd,
                      prob = prob, alpha = alpha, distribute = distribute, style = style)

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
    migration_rate = 0.005,
    #' @field speciation_rate The speciation rate of an individual
    #' Default is `0.00001`
    speciation_rate = 0.00001,

    #' @description
    #' Create a new instance of this [R6][R6::R6Class] class.
    initialize = function(
        migration_rate = self$migration_rate,
        speciation_rate = self$speciation_rate,
        distribute = "entropart",
        style = "random") {
      self$migration_rate <- migration_rate
      self$speciation_rate <- speciation_rate
      self$distribute <- distribute
      self$style <- style
      self$the_matrix <- self$make_meta()
    },

    #' @description
    #' Meta community default community drawing - add migr/spec_rate in the matrix
    make_meta = function(
    nx = 100,
    ny = nx,
    S = 100,
    Distribution = "lnorm",
    sd = 1,
    prob = 0.1,
    alpha = 40,
    migration_rate = self$migration_rate,
    speciation_rate = self$speciation_rate) {
      if (nx > 4000 || ny > 4000) {
        nx = 4000
        ny = nx
        warning("Too large number, defaulted back to the max allowed")
      }
      self$set_values(nx = nx, ny = ny, S = S,
                      Distribution = Distribution, sd = sd,
                      prob = prob, alpha = alpha)
      self$the_matrix <- NULL
      self$the_matrix <- self$draw_matrix()

      return(self$the_matrix)
    }
    )
)
