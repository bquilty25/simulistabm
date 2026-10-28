#' Create a geography for a spatial outbreak simulation
#'
#' @description A geography describes the places in which an outbreak is
#' simulated and how people move between them. Supplying a geography to
#' [sim_outbreak()], [sim_linelist()] or [sim_contacts()] gives every
#' individual a location of residence, drawn along the transmission chain:
#' each contact lives in a location sampled from the mobility probabilities of
#' the location of the individual who named them.
#'
#' @details
#' The probability that a contact of an individual in location \eqn{i} lives
#' in location \eqn{j} is given by the mobility matrix \eqn{M_{ij}}, whose
#' rows sum to 1. When `mobility` is a function (see [mobility_gravity()] and
#' [mobility_radiation()]) the off-diagonal entries are taken from it and the
#' diagonal is set to `prob_stay`. When `mobility` is a matrix it is used as
#' supplied and `prob_stay` is ignored.
#'
#' `contact_multiplier` scales the number of contacts made by individuals in
#' each location. Multipliers are divided by their maximum, so the
#' `contact_distribution` supplied to the simulation function describes the
#' location with the most contacts, and the number of contacts elsewhere is
#' thinned binomially. Binomial thinning of a negative binomial offspring
#' distribution keeps its dispersion, so superspreading is preserved.
#'
#' @param locations A `<data.frame>` with one row per location and columns:
#' * `name`: unique location names (`character`).
#' * `population`: population size (positive `numeric`).
#' * `x` and `y`: planar coordinates in kilometres (`numeric`). Required when
#'   `mobility` is a function.
#' * `contact_multiplier`: optional relative contact rate (positive
#'   `numeric`), default 1 in every location.
#' @param mobility A `function` that takes `locations` and returns a square
#' matrix of non-negative weights for moving between locations (see
#' [mobility_gravity()] and [mobility_radiation()]), or a square row-stochastic
#' `matrix` of mobility probabilities. Default is [mobility_gravity()].
#' @param prob_stay A single `numeric` between 0 and 1 for the probability that
#' a contact lives in the same location as the individual who named them. Only
#' used when `mobility` is a function. Default is 0.8.
#' @param seed_location A single `character` naming the location of the index
#' case, or `NULL` (default) to draw it in proportion to population.
#'
#' @return A `<simulist_geography>` object: a list with the checked
#' `locations` (with `contact_multiplier` scaled to a maximum of 1), the
#' `mobility` matrix and `seed_location`.
#' @export
#'
#' @examples
#' locations <- data.frame(
#'   name = c("North", "Central", "South"),
#'   population = c(40000, 120000, 60000),
#'   x = c(0, 3, 7),
#'   y = c(5, 0, -2)
#' )
#' geography <- create_geography(locations)
#' geography
#'
#' # a geography with denser contact in the centre and a short-range kernel
#' locations$contact_multiplier <- c(0.8, 1.2, 1)
#' create_geography(
#'   locations,
#'   mobility = mobility_gravity(kernel_km = 2),
#'   seed_location = "Central"
#' )
create_geography <- function(locations,
                             mobility = mobility_gravity(),
                             prob_stay = 0.8,
                             seed_location = NULL) {
  checkmate::assert_data_frame(locations, min.rows = 1)
  checkmate::assert_names(colnames(locations), must.include = c("name", "population"))
  checkmate::assert_character(locations$name, any.missing = FALSE, unique = TRUE)
  checkmate::assert_numeric(locations$population, lower = 0, any.missing = FALSE)
  stopifnot("Every location must have a positive population" = all(locations$population > 0))
  if (is.null(locations$contact_multiplier)) {
    locations$contact_multiplier <- 1
  }
  checkmate::assert_numeric(locations$contact_multiplier, lower = 0, any.missing = FALSE)
  stopifnot("`contact_multiplier` must be positive in at least one location" = max(locations$contact_multiplier) > 0)
  locations$contact_multiplier <- locations$contact_multiplier / max(locations$contact_multiplier)

  n <- nrow(locations)
  if (is.function(mobility)) {
    checkmate::assert_number(prob_stay, lower = 0, upper = 1)
    if (n > 1) {
      checkmate::assert_names(colnames(locations), must.include = c("x", "y"))
      checkmate::assert_numeric(locations$x, any.missing = FALSE)
      checkmate::assert_numeric(locations$y, any.missing = FALSE)
    }
    weights <- if (n > 1) mobility(locations) else matrix(0, 1, 1)
    checkmate::assert_matrix(weights, mode = "numeric", nrows = n, ncols = n, any.missing = FALSE)
    stopifnot("Mobility weights must be non-negative" = all(weights >= 0))
    diag(weights) <- 0
    row_totals <- rowSums(weights)
    # a location with no outward weight keeps all of its contacts
    mobility_matrix <- weights / ifelse(row_totals > 0, row_totals, 1) * (1 - prob_stay)
    diag(mobility_matrix) <- ifelse(row_totals > 0, prob_stay, 1)
  } else {
    checkmate::assert_matrix(mobility, mode = "numeric", nrows = n, ncols = n, any.missing = FALSE)
    stopifnot(
      "Mobility probabilities must be between 0 and 1" = all(mobility >= 0 & mobility <= 1),
      "Each row of the mobility matrix must sum to 1" = isTRUE(all.equal(unname(rowSums(mobility)), rep(1, n)))
    )
    mobility_matrix <- mobility
  }
  dimnames(mobility_matrix) <- list(locations$name, locations$name)

  if (!is.null(seed_location)) {
    checkmate::assert_choice(seed_location, choices = locations$name)
  }

  structure(
    list(
      locations = locations,
      mobility = mobility_matrix,
      seed_location = seed_location
    ),
    class = "simulist_geography"
  )
}

#' @export
print.simulist_geography <- function(x, ...) {
  writeLines(c(
    sprintf("<simulist_geography> with %d locations", nrow(x$locations)),
    sprintf("Total population: %s", format(sum(x$locations$population), big.mark = ",")),
    sprintf("Mean probability of staying in a location: %.2f", mean(diag(x$mobility))),
    sprintf("Index case location: %s", if (is.null(x$seed_location)) "drawn by population" else x$seed_location)
  ))
  invisible(x)
}

#' Gravity model for mobility between locations
#'
#' @description Returns a function for the `mobility` argument of
#' [create_geography()]. The weight for moving from location \eqn{i} to
#' location \eqn{j} is the population of \eqn{j} multiplied by a kernel of the
#' distance between them, either \eqn{\exp(-d_{ij} / \kappa)} (exponential) or
#' \eqn{(1 + d_{ij} / \kappa)^{-\gamma}} (power law).
#'
#' @param kernel_km A single positive `numeric` for the distance scale
#' \eqn{\kappa} of the kernel, in kilometres. Default is 5.
#' @param kernel A `character` string, either `"exponential"` (default) or
#' `"power"`.
#' @param power A single positive `numeric` for the exponent \eqn{\gamma} of
#' the power-law kernel. Default is 2.
#'
#' @return A `function` that takes a `locations` `<data.frame>` and returns a
#' matrix of mobility weights.
#' @export
#'
#' @examples
#' locations <- data.frame(
#'   name = c("A", "B", "C"),
#'   population = c(1000, 5000, 2000),
#'   x = c(0, 1, 10),
#'   y = c(0, 0, 0)
#' )
#' mobility_gravity(kernel_km = 3)(locations)
mobility_gravity <- function(kernel_km = 5,
                             kernel = c("exponential", "power"),
                             power = 2) {
  checkmate::assert_number(kernel_km, lower = 0)
  stopifnot("`kernel_km` must be positive" = kernel_km > 0)
  kernel <- match.arg(kernel)
  checkmate::assert_number(power, lower = 0)
  function(locations) {
    distance <- .location_distances(locations)
    decay <- if (kernel == "exponential") {
      exp(-distance / kernel_km)
    } else {
      (1 + distance / kernel_km)^-power
    }
    matrix(locations$population, nrow(locations), nrow(locations), byrow = TRUE) * decay
  }
}

#' Radiation model for mobility between locations
#'
#' @description Returns a function for the `mobility` argument of
#' [create_geography()]. The radiation model (Simini et al. 2012) has no free
#' parameters: the weight for moving from location \eqn{i} to location \eqn{j}
#' is \eqn{m_i n_j / ((m_i + s_{ij})(m_i + n_j + s_{ij}))}, where \eqn{m_i}
#' and \eqn{n_j} are the populations of \eqn{i} and \eqn{j}, and
#' \eqn{s_{ij}} is the population living closer to \eqn{i} than \eqn{j} is,
#' excluding \eqn{i} and \eqn{j}.
#'
#' @return A `function` that takes a `locations` `<data.frame>` and returns a
#' matrix of mobility weights.
#' @export
#'
#' @references Simini, F., González, M. C., Maritan, A. and Barabási, A.-L.
#' (2012) A universal model for mobility and migration patterns. Nature 484,
#' 96–100. \doi{10.1038/nature10856}
#'
#' @examples
#' locations <- data.frame(
#'   name = c("A", "B", "C"),
#'   population = c(1000, 5000, 2000),
#'   x = c(0, 1, 10),
#'   y = c(0, 0, 0)
#' )
#' mobility_radiation()(locations)
mobility_radiation <- function() {
  function(locations) {
    distance <- .location_distances(locations)
    pop <- locations$population
    n <- nrow(locations)
    weights <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)[-i]) {
        closer <- distance[i, ] < distance[i, j]
        closer[c(i, j)] <- FALSE
        s <- sum(pop[closer])
        weights[i, j] <- pop[i] * pop[j] / ((pop[i] + s) * (pop[i] + pop[j] + s))
      }
    }
    weights
  }
}

#' Euclidean distances between locations, in the units of `x` and `y`
#' @keywords internal
.location_distances <- function(locations) {
  as.matrix(stats::dist(cbind(locations$x, locations$y)))
}

#' Probability of infection that changes with time and location
#'
#' @description Builds a function for the `prob_infection` argument of the
#' simulation functions, for interventions that reduce transmission in some
#' locations over some period. Reductions from overlapping interventions
#' multiply.
#'
#' @param base A single `numeric` between 0 and 1, or a named `numeric`
#' vector with one probability per location, for the probability of infection
#' without interventions.
#' @param interventions A `<data.frame>` with columns `location` (a location
#' name, or `"all"`), `start` and `end` (days since the start of the
#' outbreak; `end` can be `Inf`), and `reduction` (proportional reduction in
#' the probability of infection, between 0 and 1). Can have zero rows.
#'
#' @return A `function` with arguments `time` and `location`, returning the
#' probability of infection for contacts at those times and locations.
#' @export
#'
#' @examples
#' prob_infection <- intervention_schedule(
#'   base = 0.4,
#'   interventions = data.frame(
#'     location = c("Central", "all"),
#'     start = c(30, 60),
#'     end = c(Inf, 90),
#'     reduction = c(0.5, 0.3)
#'   )
#' )
#' prob_infection(time = c(10, 45, 70), location = c("Central", "Central", "North"))
intervention_schedule <- function(base,
                                  interventions = data.frame(
                                    location = character(),
                                    start = numeric(),
                                    end = numeric(),
                                    reduction = numeric()
                                  )) {
  checkmate::assert_numeric(base, lower = 0, upper = 1, any.missing = FALSE, min.len = 1)
  if (length(base) > 1) {
    checkmate::assert_names(names(base), type = "unique")
  }
  checkmate::assert_data_frame(interventions)
  checkmate::assert_names(colnames(interventions), must.include = c("location", "start", "end", "reduction"))
  checkmate::assert_numeric(interventions$reduction, lower = 0, upper = 1, any.missing = FALSE)
  stopifnot("Intervention `end` must not be before `start`" = all(interventions$end >= interventions$start))
  function(time, location) {
    p <- if (length(base) == 1) rep(base, length(time)) else unname(base[location])
    if (anyNA(p)) {
      stop("`base` has no probability of infection for some locations.", call. = FALSE)
    }
    for (k in seq_len(nrow(interventions))) {
      applies <- time >= interventions$start[k] & time < interventions$end[k] &
        (interventions$location[k] == "all" | location == interventions$location[k])
      p[applies] <- p[applies] * (1 - interventions$reduction[k])
    }
    p
  }
}
