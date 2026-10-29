#' Turn a transmission chain from any model into a line list and contacts table
#'
#' @description [sim_outbreak()] simulates a transmission chain with a
#' branching process and turns it into a line list and contacts table.
#' `chain_to_outbreak()` does the second step for a chain simulated by another
#' model, for example an agent-based model with households and finite
#' populations, so that its output has the same structure as [sim_outbreak()]
#' and works with [messy_linelist()], [truncate_linelist()] and
#' [censor_linelist()].
#'
#' @details
#' The chain has one row per individual who was infected or was a contact of
#' an infected individual. Infected individuals appear once, with the
#' individual who infected them. A contact who was not infected appears once
#' for each infected individual who named them. Times are days since
#' `outbreak_start_date`: for infected individuals, the time of symptom onset
#' (which [sim_outbreak()] takes to be the time of infection); for contacts who
#' were not infected, the time of contact.
#'
#' Reporting delays, contact dates, hospitalisation, outcomes, names, case
#' types, Ct values and contact tracing status are added as in
#' [sim_outbreak()]. Sex and age are sampled as in [sim_outbreak()] unless the
#' chain has `sex` or `age` columns, in which case they are kept, so that
#' age-specific risks apply to the ages the model simulated.
#'
#' @param chain A `<data.frame>` with columns:
#' * `id`: unique positive integer for each row.
#' * `infector`: the `id` of the infected individual who infected or named
#'   this individual, or `NA` for index cases. Must refer to an infected row.
#' * `infected`: `TRUE` (or `"infected"`) if the individual was infected,
#'   `FALSE` (or `"contact"`) if they were a contact who was not infected.
#' * `time`: non-negative `numeric` time in days, as described in Details.
#' * Optional `location` (`character`), `age` (whole numbers) and `sex`
#'   (`"m"` or `"f"`).
#' @inheritParams sim_outbreak
#'
#' @return A list with a line list `<data.frame>` and a contacts
#' `<data.frame>`, with the same columns as [sim_outbreak()], plus `location`
#' in the line list and `from_location` and `to_location` in the contacts
#' table when the chain has a `location` column.
#' @export
#'
#' @examples
#' # two index cases; case 1 infects case 3 and names a contact who is not
#' # infected; case 3 infects case 5
#' chain <- data.frame(
#'   id = 1:5,
#'   infector = c(NA, NA, 1, 1, 3),
#'   infected = c(TRUE, TRUE, TRUE, FALSE, TRUE),
#'   time = c(0, 1.5, 4.2, 3.1, 9.8),
#'   location = c("North", "South", "North", "North", "South")
#' )
#' outbreak <- chain_to_outbreak(chain)
#' outbreak$linelist
#' outbreak$contacts
chain_to_outbreak <- function(chain,
                              onset_to_hosp = function(n) stats::rlnorm(n = n, meanlog = 1.5, sdlog = 0.5), # nolint start: line_length_linter.
                              onset_to_death = function(n) stats::rlnorm(n = n, meanlog = 2.5, sdlog = 0.5), # nolint end
                              onset_to_recovery = NULL,
                              reporting_delay = NULL,
                              hosp_risk = 0.2,
                              hosp_death_risk = 0.5,
                              non_hosp_death_risk = 0.05,
                              outbreak_start_date = as.Date("2023-01-01"),
                              anonymise = FALSE,
                              population_age = c(0, 90),
                              case_type_probs = c(
                                suspected = 0.2,
                                probable = 0.3,
                                confirmed = 0.5
                              ),
                              contact_tracing_status_probs = c(
                                under_followup = 0.7,
                                lost_to_followup = 0.2,
                                unknown = 0.1
                              ),
                              config = create_config()) {
  chain <- .check_chain(chain)

  # the transmission arguments of as_function() and .check_sim_input() are not
  # used here, so placeholders satisfy the shared conversions and checks
  funcs <- as_function(
    list(
      contact_distribution = function(x) x,
      infectious_period = function(n) n,
      onset_to_hosp = onset_to_hosp,
      onset_to_death = onset_to_death,
      onset_to_recovery = onset_to_recovery,
      reporting_delay = reporting_delay
    )
  )

  .check_sim_input(
    sim_type = "outbreak",
    contact_distribution = function(x) x,
    infectious_period = function(n) n,
    prob_infection = 0.5,
    outbreak_start_date = outbreak_start_date,
    outbreak_size = c(1, 1),
    onset_to_hosp = funcs$onset_to_hosp,
    onset_to_death = funcs$onset_to_death,
    onset_to_recovery = funcs$onset_to_recovery,
    reporting_delay = funcs$reporting_delay,
    anonymise = anonymise,
    case_type_probs = case_type_probs,
    contact_tracing_status_probs = contact_tracing_status_probs,
    hosp_risk = hosp_risk,
    hosp_death_risk = hosp_death_risk,
    non_hosp_death_risk = non_hosp_death_risk,
    population_age = population_age
  )
  .cross_check_sim_input(
    onset_to_hosp = funcs$onset_to_hosp,
    onset_to_death = funcs$onset_to_death,
    hosp_risk = hosp_risk,
    hosp_death_risk = hosp_death_risk,
    non_hosp_death_risk = non_hosp_death_risk
  )

  if (is.data.frame(population_age)) {
    population_age <- .check_df(population_age, df_type = "age")
    age_range <- c(
      lower = min(population_age[, "min_age"]),
      upper = max(population_age[, "max_age"])
    )
  } else {
    population_age <- sort(population_age)
    names(population_age) <- c("lower", "upper")
    age_range <- population_age
  }
  # risks are checked against the population age range, as in sim_outbreak(),
  # not the ages that happen to be in the chain
  for (risk in c("hosp_risk", "hosp_death_risk", "non_hosp_death_risk")) {
    value <- get(risk)
    if (is.data.frame(value)) {
      assign(risk, .check_df(value, df_type = "risk", age_range = age_range))
    }
  }

  .chain_to_outputs(
    .data = chain,
    sim_type = "outbreak",
    onset_to_hosp = funcs$onset_to_hosp,
    onset_to_death = funcs$onset_to_death,
    onset_to_recovery = funcs$onset_to_recovery,
    reporting_delay = funcs$reporting_delay,
    hosp_risk = hosp_risk,
    hosp_death_risk = hosp_death_risk,
    non_hosp_death_risk = non_hosp_death_risk,
    outbreak_start_date = outbreak_start_date,
    anonymise = anonymise,
    population_age = population_age,
    case_type_probs = case_type_probs,
    contact_tracing_status_probs = contact_tracing_status_probs,
    config = config
  )
}

#' Check a transmission chain and convert it to the internal structure
#'
#' @inheritParams chain_to_outbreak
#'
#' @return The chain as a `<data.frame>` with `infected` coded as
#' `"infected"` or `"contact"`, ordered by `id`.
#' @keywords internal
.check_chain <- function(chain) {
  checkmate::assert_data_frame(chain, min.rows = 1)
  checkmate::assert_names(colnames(chain), must.include = c("id", "infector", "infected", "time"))
  checkmate::assert_integerish(chain$id, lower = 1, any.missing = FALSE, unique = TRUE)
  checkmate::assert_numeric(chain$time, lower = 0, any.missing = FALSE)

  infected <- chain$infected
  if (is.character(infected)) {
    checkmate::assert_subset(infected, c("infected", "contact"))
    infected <- infected == "infected"
  }
  checkmate::assert_logical(infected, any.missing = FALSE)

  roots <- is.na(chain$infector)
  stopifnot(
    "The chain must have at least one index case (`infector` is NA)" = any(roots),
    "Index cases (`infector` is NA) must be infected" = all(infected[roots]),
    "Every `infector` must be the `id` of an infected row" =
      all(chain$infector[!roots] %in% chain$id[infected])
  )
  # every row must trace back to an index case, which rules out individuals
  # infecting themselves and cycles; each pass reaches one more generation
  reached <- roots
  repeat {
    newly <- !reached & chain$infector %in% chain$id[reached]
    if (!any(newly)) break
    reached <- reached | newly
  }
  stopifnot(
    "Every row must trace back to an index case through `infector` (no cycles)" = all(reached)
  )
  if ("location" %in% colnames(chain)) {
    checkmate::assert_character(chain$location, any.missing = FALSE)
  }
  if ("age" %in% colnames(chain)) {
    checkmate::assert_integerish(chain$age, lower = 0, any.missing = FALSE)
  }
  if ("sex" %in% colnames(chain)) {
    checkmate::assert_subset(chain$sex, c("m", "f"))
  }

  chain$id <- as.integer(chain$id)
  chain$infector <- as.integer(chain$infector)
  chain$infected <- ifelse(infected, "infected", "contact")
  chain <- chain[order(chain$id), , drop = FALSE]
  row.names(chain) <- NULL
  chain
}
