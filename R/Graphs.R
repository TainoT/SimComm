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
#'   - `relative_abundance_rank`
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

  g <- ggplot(df, aes(x = factor(species, levels = species), y = count)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = paste("Abundance of species at time", time),
         x = "Species", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
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

  g <- ggplot(df, aes(x = reorder(factor(species), -rela_abu),
                 y = rela_abu, fill = factor(species))) +
    geom_bar(stat = "identity") +
    labs(title = paste("Relative Abundance of Species at Time", time),
         x = "Species", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    scale_y_continuous(labels = scales::percent_format())
  print(g)
}

#' @rdname overall_abundance_distribution
#'
#' @export
overall_abundance_distribution <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  # abu_distr <- NULL
  abu_distr <- aggregate(count ~ species, data, sum)
  abu_distr <- abu_distr[order(-abu_distr$count), ]

  g <- ggplot(data, aes(x = factor(time), y = count, fill = factor(species))) +
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
  print(g)
}


#' @rdname overall_abundance_rank
#'
#' @export
overall_abundance_rank <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  abu_rank <- aggregate(count ~ rank, data,
                        function(x) c(mean = mean(x), sd = sd(x)))
  abu_rank <- do.call(data.frame, abu_rank)
  colnames(abu_rank) <- c("rank", "mean", "sd")
  abu_rank$rank <- reorder(abu_rank$rank, -abu_rank$mean)

  g <- ggplot(abu_rank, aes(x = rank, y = mean, fill = as.factor(rank))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  width = 0.2, position = position_dodge(0.9)) +
    labs(title = "Mean/SD of Abundance per Rank of Species Abundance Over Time",
         x = "Rank of Species Abundance", y = "Count",
         fill = "Rank") +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position = "none"
    )
  print(g)
}

#' @rdname relative_abundance_rank
#'
#' @export
relative_abundance_rank <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  sp_count <- aggregate(count ~ species, data, sum)
  sp_count <- sp_count[order(-sp_count$count), ]
  sp_count$rank <- 1:nrow(sp_count)

  sp_count$rela_abu <- sp_count$count / sum(sp_count$count)

  g <- ggplot(sp_count, aes(x = rank, y = rela_abu)) +
    geom_line() +
    labs(title = "Relative Abundance of Species by Rank",
         x = "Species Rank",
         y = "Relative Abundance") +
    theme_minimal() +
    scale_y_log10()
  print(g)
}

#' @rdname neighors_rank
#'
#' @export
# make a function that plot the neighbors for all species by rank over time
neighbors_rank <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  g <- ggplot(data, aes(x = rank, y = neighbors, color = factor(species))) +
    geom_point(size = 2) +
    geom_line(aes(group = species), linetype = "dashed") +
    labs(title = "Rank and Neighbors over time",
         x = "Rank", y = "Neighbors",
         color = "Species") +
    theme_minimal()+
    theme(legend.position = "none")
  print(g)
}

#' @rdname neighbors_time
#'
#' @export
neighbors_time <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  g <- ggplot(data, aes(x = time, y = neighbors, color = factor(species))) +
    geom_point(size = 2) +
    geom_line(aes(group = species)) +
    labs(title = "Time and Neighbors per species",
         x = "Time", y = "Neighbors",
         color = "Species") +
    theme_minimal()+
    theme(legend.position = "none")
  print(g)
}

#' @rdname neighbors_abundance
#'
#' @export
neighbors_abundance <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  g <- ggplot(data, aes(x = count, y = neighbors, color = factor(species))) +
    geom_point(size = 2) +
    geom_line(aes(group = species), linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Abundance and Neighbors",
         x = "Abundance", y = "Neighbors",
         color = "Species") +
    theme_minimal()+
    theme(legend.position = "none")
  print(g)
}

#' @rdname neighbors_abundance_time
#'
#' @export
neighbors_abundance_time <- function(
    data = NULL) {
  if (is.null(data))
    stop("Data frame is missing")

  g <- ggplot(data, aes(x = time)) +
    geom_line(aes(y = count, color = factor(species)), size = 1) +
    geom_point(aes(y = count, color = factor(species)), size = 1) +
    scale_y_continuous(name = "Abundance") +
    labs(title = "Abundance over time",
         x = "Time",
         color = "Species") +
    theme_minimal()+
    theme(legend.position = "none")

  g <- g +
    geom_line(aes(y = neighbors * max(count) / max(neighbors),
                 color = factor(species)), linetype = "dashed") +
    geom_point(aes(y = neighbors * max(count) / max(neighbors),
                  color = factor(species)), shape = 1, size = 1) +
    scale_y_continuous(
      sec.axis = sec_axis(~ . * max(data$neighbors) / max(data$count),
                          name = "Neighbors")
    )
  print(g)
}


plot_neighbors_over_time <- function(data) {
  if (is.null(data)) {
    stop("Data frame is missing")
  }

  # Plot: Relationship between Neighbors and Time, grouped by Species
  ggplot(data, aes(x = time, y = neighbors, color = factor(species))) +
    geom_point(size = 3) +
    geom_line(aes(group = species)) +
    labs(
      title = "Relationship between Neighbors and Time, grouped by Species",
      x = "Time",
      y = "Average Number of Same-species Neighbors",
      color = "Species"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

plot_abundance_and_neighbors <- function(data) {
  if (is.null(data)) {
    stop("Data frame is missing")
  }

  max_count <- max(data$count, na.rm = TRUE)
  max_neighbors <- max(data$neighbors, na.rm = TRUE)

  # Plot with dual y-axes
  ggplot(data, aes(x = time)) +
    geom_line(aes(y = count, color = factor(species)), size = 1) +
    geom_point(aes(y = count, color = factor(species)), size = 3) +
    scale_y_continuous(
      name = "Abundance (Count)",
      sec.axis = sec_axis(~ . * max_neighbors / max_count, name = "Average Number of Same-species Neighbors")
    ) +
    geom_line(aes(y = neighbors * max_count / max_neighbors, color = factor(species)), linetype = "dashed") +
    geom_point(aes(y = neighbors * max_count / max_neighbors, color = factor(species)), shape = 1, size = 3) +
    labs(
      title = "Abundance and Neighbors over Time",
      x = "Time",
      color = "Species"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}


cumul_new_extinct_species <- function(data) {
  if (is.null(data)) {
    stop("Data frame is missing")
  }

new_species <- numeric()
extinct_species <- numeric()
seen_species <- c()

for (t in unique(data$time)) {
  current_species <- data %>% filter(time == t) %>% pull(species)
  if (t == min(unique(data$time))) {
    new_species_count <- length(current_species)
    extinct_species_count <- 0
  } else {
    new_species_count <- sum(!current_species %in% seen_species)
    previous_species <- data %>% filter(time == (t - 1)) %>% pull(species)
    extinct_species_count <- sum(!previous_species %in% current_species)
  }

  new_species <- c(new_species, new_species_count)
  extinct_species <- c(extinct_species, extinct_species_count)
  seen_species <- union(seen_species, current_species)
}

summary <- data.frame(  time = unique(data$time),  new_species = new_species,  extinct_species = extinct_species)
summary <- summary %>% mutate(cumulative_new = cumsum(new_species),
                              cumulative_extinct = cumsum(extinct_species))

ggplot(summary, aes(x = time)) +
  geom_line(aes(y = cumulative_new, color = 'New Species')) +
  geom_line(aes(y = cumulative_extinct, color = 'Extinct Species')) +
  labs(title = 'Cumulative New and Extinct Species Over Time',
       x = 'Time',
       y = 'Cumulative Count') +
  scale_color_manual(name = 'Legend',
                     values = c('New Species' = 'blue', 'Extinct Species' = 'red')) +
  theme_minimal()
}

#' @rdname map_of_graphs
#'
#' @export
map_of_graphs <- list(
  "timed_abundance" = timed_abundance,
  "timed_relative_abundance" = timed_relative_abundance,
  "overall_abundance_distribution" = overall_abundance_distribution,
  "overall_abundance_rank" = overall_abundance_rank,
  "relative_abundance_rank" = relative_abundance_rank,
  "neighbors_rank" = neighbors_rank,
  "neighbors_time" = neighbors_time,
  "neighbors_abundance" = neighbors_abundance,
  "neighbors_abundance_time" = neighbors_abundance_time,
  "plot_neighbors_over_time" = plot_neighbors_over_time,
  "plot_abundance_and_neighbors" = plot_abundance_and_neighbors,
  "cumul_new_extinct_species" = cumul_new_extinct_species
)
