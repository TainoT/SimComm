#' Graphs
#'
#' @description
#' Graphs describe the structure of the community.
#'
#' @details
#' Several graphs are available by default.
#'
#' Available graphs are:
#'
#'   - `timed_abundance`
#'
#'   - `timed_relative_abundance`
#'
#'   - `overall_abundance_distribution`
#'
#'   - `overall_abundance_rank`
#'
#'
#' @name graphs
NULL


#' @rdname timed_abundance
#'
#' @export
timed_abundance <- function(
    time = NULL,
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")
  if (is.null(time) || isTRUE(time < 1) || isTRUE(time > nrow(df))) {
    warning("Time is unfit, using the first step")
    time <- 1
  }

  df <- data[data$time == time, ]
  df <- df[order(-df$count), ]

  ggplot(df, aes(x = factor(species, levels = species), y = count)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = paste("Abundance of species at time", time),
         x = "Species", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' @rdname timed_relative_abundance
#'
#' @export
timed_relative_abundance <- function(
    time = NULL,
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")
  if (is.null(time) || isTRUE(time < 1) || isTRUE(time > nrow(df))) {
    warning("Time is unfit, using the first step")
    time <- 1
  }

  df <- data[data$time == time, ]
  total_abu <- sum(df$count)
  df$rela_abu <- df$count / total_abu
  df <- df[order(-df$rela_abu), ]

  ggplot(df, aes(x = reorder(factor(species), -rela_abu), y = rela_abu, fill = factor(species))) +
    geom_bar(stat = "identity") +
    labs(title = paste("Relative Abundance of Species at Time", time),
         x = "Species", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    scale_y_continuous(labels = scales::percent_format())

  # ggplot(df, aes(x = factor(species, levels = species), y = rela_abu)) +
  #   geom_bar(stat = "identity", fill = "skyblue") +
  #   labs(title = paste("Relative abundance of species at time", time),
  #        x = "Species", y = "Relative Abundance") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

#' @rdname overall_abundance_distribution
#'
#' @export
overall_abundance_distribution <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  abu_distr <- NULL
  abu_distr <- aggregate(count ~ species, data, sum)
  abu_distr <- abu_distr[order(-abu_distr$count), ]
  ggplot(data, aes(x = factor(time), y = count, fill = factor(species))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = rainbow(length(abu_distr$species))) +
    labs(title = "Distribution of Abundance of Species Over Time",
         x = "Time", y = "Count",
         fill = "Species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")
}


#' @rdname overall_abundance_rank
#'
#' @export
overall_abundance_rank <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  abu_rank <- aggregate(count ~ rank, data, function(x) c(mean = mean(x), sd = sd(x)))
  abu_rank <- do.call(data.frame, abu_rank)
  colnames(abu_rank) <- c("rank", "mean", "sd")
  abu_rank$rank <- reorder(abu_rank$rank, -abu_rank$mean)

  ggplot(abu_rank, aes(x = rank, y = mean, fill = as.factor(rank))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(0.9)) +
    labs(title = "Mean and SD of Abundance per Rank of Species Abundance Over Time",
         x = "Rank of Species Abundance", y = "Count",
         fill = "Rank") +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position = "none"
    )
}

