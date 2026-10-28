#' Simulate a line list and contacts table with an agent-based model
#'
#' @description Simulates an outbreak in a finite, structured population with
#' the \pkg{individual} agent-based modelling framework, and returns a line
#' list and contacts table with the same structure as [sim_outbreak()] (with
#' a geography), so [messy_linelist()], [truncate_linelist()] and
#' [censor_linelist()] apply unchanged.
#'
#' @details
#' Every person in the `geography` is an agent. Agents live in households
#' within their location, and are susceptible, latent, infectious or
#' recovered (SEIR). Each day, a susceptible agent is infected with
#' probability \eqn{1 - \exp(-(h_H + h_C))}, where the household hazard
#' \eqn{h_H} is `household_hazard` times the total infectiousness of
#' infectious members of their household, and the community hazard \eqn{h_C}
#' is `community_hazard` times the infectiousness per person, averaged over
#' locations weighted by the agent's row of the mobility matrix. Both are
#' scaled by the contact multiplier of the agent's location and by any active
#' `interventions`. Because the population is finite, the outbreak slows as
#' susceptible agents are used up.
#'
#' The infector of each new case is recorded: the infection is attributed to
#' the household or the community in proportion to the two hazards, and the
#' infector is drawn from the infectious agents there in proportion to their
#' infectiousness. Each case's contacts are its uninfected household members
#' and a Poisson number (mean `contacts_named`) of uninfected agents from its
#' location, contacted at a uniform time in its infectious period. Symptom
#' onset is taken as the start of infectiousness, as [sim_outbreak()] takes
#' onset as the time of infection.
#'
#' The simulation runs on a daily time step until no agent is latent or
#' infectious, `max_days` have passed, or cumulative infections exceed the
#' maximum of `outbreak_size` (with a warning, as in [sim_outbreak()]). Cases
#' with onset after the last simulated day are left out, so the line list is
#' complete up to the date in the `"complete_until"` attribute. Simulations
#' with fewer cases than the minimum of `outbreak_size` are rerun.
#'
#' @param geography A `<simulist_geography>` object from [create_geography()].
#' Its populations are the numbers of agents.
#' @param household_hazard A single non-negative `numeric`: the daily
#' infection hazard from one infectious household member of average
#' infectiousness.
#' @param community_hazard A single non-negative `numeric`: the daily community
#' infection hazard if everyone an agent mixes with were infectious with
#' average infectiousness. In a fully susceptible population, an infectious
#' agent of average infectiousness causes about `community_hazard` times the
#' mean infectious period community infections.
#' @param latent_period,infectious_period A `function` or an `<epiparameter>`
#' object for the latent and infectious period distributions, in days. A
#' function must have a single argument and return that many random numbers.
#' @param infectiousness A `function` with a single argument returning that
#' many relative infectiousness values, one per agent, with mean 1. The
#' default Gamma distribution with shape 1 gives offspring dispersion close to
#' 1; smaller shapes give more superspreading.
#' @param household_size A `function` with a single argument returning that
#' many household sizes (whole numbers of at least 1).
#' @param interventions A `<data.frame>` with columns `location` (a location
#' name, or `"all"`), `start` and `end` (days since the start of the outbreak;
#' `end` can be `Inf`) and `reduction` (proportional reduction in
#' susceptibility, between 0 and 1), or `NULL` (default) for none.
#' @param n_seeds A single whole number of agents infected at the start, in
#' the geography's `seed_location` (or a location drawn by population).
#' @param contacts_named A single non-negative `numeric`: the mean number of
#' uninfected community contacts each case names to contact tracers, in
#' addition to their uninfected household members.
#' @param max_days A single whole number: the longest the simulation runs.
#' @param outbreak_size A `numeric` vector of length 2 with the minimum and
#' maximum number of infections. Simulations with fewer cases than the
#' minimum are rerun; simulations stop at the end of the day on which
#' infections exceed the maximum.
#' @inheritParams sim_outbreak
#'
#' @return A list with a line list `<data.frame>` and a contacts
#' `<data.frame>`, as from [sim_outbreak()] with a geography, with
#' attributes:
#' * `"daily"`: a `<data.frame>` of new infections and infectious agents by
#'   date and location.
#' * `"complete_until"`: the last date with complete onsets.
#' * `"attempts"`: the number of simulations run to reach the minimum size.
#' @export
#'
#' @examples
#' geography <- create_geography(
#'   data.frame(
#'     name = c("North", "Central", "South"),
#'     population = c(4000, 12000, 6000),
#'     x = c(0, 3, 7),
#'     y = c(5, 0, -2)
#'   ),
#'   seed_location = "Central"
#' )
#' set.seed(1)
#' outbreak <- sim_outbreak_abm(
#'   geography,
#'   household_hazard = 0.08,
#'   community_hazard = 0.3,
#'   outbreak_size = c(50, 2000)
#' )
#' head(outbreak$linelist)
#' head(outbreak$contacts)
#' head(attr(outbreak, "daily"))
sim_outbreak_abm <- function(geography,
                             household_hazard = 0.05,
                             community_hazard = 0.25,
                             latent_period = function(n) stats::rlnorm(n = n, meanlog = 1, sdlog = 0.3),
                             infectious_period = function(n) stats::rlnorm(n = n, meanlog = 1.5, sdlog = 0.4),
                             infectiousness = function(n) stats::rgamma(n = n, shape = 1, rate = 1),
                             household_size = function(n) 1L + stats::rpois(n = n, lambda = 3.5),
                             interventions = NULL,
                             n_seeds = 1,
                             contacts_named = 2,
                             max_days = 730,
                             onset_to_hosp = function(n) stats::rlnorm(n = n, meanlog = 1.5, sdlog = 0.5), # nolint start: line_length_linter.
                             onset_to_death = function(n) stats::rlnorm(n = n, meanlog = 2.5, sdlog = 0.5), # nolint end
                             onset_to_recovery = NULL,
                             reporting_delay = NULL,
                             hosp_risk = 0.2,
                             hosp_death_risk = 0.5,
                             non_hosp_death_risk = 0.05,
                             outbreak_start_date = as.Date("2023-01-01"),
                             anonymise = FALSE,
                             outbreak_size = c(10, 1e4),
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
  checkmate::assert_class(geography, "simulist_geography")
  checkmate::assert_number(household_hazard, lower = 0)
  checkmate::assert_number(community_hazard, lower = 0)
  latent_period <- as.function(latent_period, func_type = "generate")
  infectious_period <- as.function(infectious_period, func_type = "generate")
  .check_func_req_args(latent_period, func_name = "latent_period")
  .check_func_req_args(infectious_period, func_name = "infectious_period")
  .check_func_req_args(infectiousness, func_name = "infectiousness")
  .check_func_req_args(household_size, func_name = "household_size")
  checkmate::assert_count(n_seeds, positive = TRUE)
  checkmate::assert_number(contacts_named, lower = 0)
  checkmate::assert_count(max_days, positive = TRUE)
  checkmate::assert_integerish(outbreak_size, lower = 1, len = 2)
  if (!is.null(interventions)) {
    checkmate::assert_data_frame(interventions)
    checkmate::assert_names(colnames(interventions), must.include = c("location", "start", "end", "reduction"))
    checkmate::assert_subset(interventions$location, c("all", geography$locations$name))
    checkmate::assert_numeric(interventions$reduction, lower = 0, upper = 1, any.missing = FALSE)
  }
  # chain_to_outbreak() checks the argument as supplied; the agents' ages
  # come from the checked version
  population_age_input <- population_age
  if (is.data.frame(population_age)) {
    population_age <- .check_df(population_age, df_type = "age")
  } else {
    checkmate::assert_numeric(population_age, len = 2)
    population_age <- sort(population_age)
  }
  stopifnot(
    "`prob_male` in `config` must be >= 0 and <=1" =
      config$prob_male >= 0 && config$prob_male <= 1
  )

  # the population is built once and reused by every attempt, since rebuilding
  # millions of agents costs more than a failed introduction
  agents <- .abm_population(geography, household_size, population_age, config$prob_male, infectiousness)
  for (attempt in seq_len(1e3)) {
    run <- .sim_abm(
      geography = geography,
      agents = agents,
      household_hazard = household_hazard,
      community_hazard = community_hazard,
      latent_period = latent_period,
      infectious_period = infectious_period,
      interventions = interventions,
      n_seeds = n_seeds,
      contacts_named = contacts_named,
      max_days = max_days,
      max_infections = max(outbreak_size)
    )
    if (sum(run$chain$infected) >= min(outbreak_size)) break
    if (attempt == 1e3) {
      stop(
        "Exceeded maximum number of iterations for simulating outbreak. \n",
        "Change input parameters or outbreak_size.",
        call. = FALSE
      )
    }
  }
  if (run$hit_max) {
    warning(
      "Number of cases exceeds maximum outbreak size. \n",
      "Returning data early with ", run$n_infected, " infections.",
      call. = FALSE
    )
  }

  outbreak <- chain_to_outbreak(
    run$chain[c("id", "infector", "infected", "time", "location", "age", "sex")],
    onset_to_hosp = onset_to_hosp,
    onset_to_death = onset_to_death,
    onset_to_recovery = onset_to_recovery,
    reporting_delay = reporting_delay,
    hosp_risk = hosp_risk,
    hosp_death_risk = hosp_death_risk,
    non_hosp_death_risk = non_hosp_death_risk,
    outbreak_start_date = outbreak_start_date,
    anonymise = anonymise,
    population_age = population_age_input,
    case_type_probs = case_type_probs,
    contact_tracing_status_probs = contact_tracing_status_probs,
    config = config
  )
  run$daily$date <- outbreak_start_date + run$daily$day
  attr(outbreak, "daily") <- run$daily[c("date", "location", "new_infections", "infectious")]
  attr(outbreak, "complete_until") <- outbreak_start_date + run$last_day
  attr(outbreak, "attempts") <- attempt
  outbreak
}

#' Build the agents of an agent-based simulation
#'
#' @inheritParams sim_outbreak_abm
#' @param prob_male Probability an agent is male.
#'
#' @return A list with each agent's location, household, age, sex and
#' relative infectiousness, and the population of each location.
#' @keywords internal
.abm_population <- function(geography, household_size, population_age, prob_male, infectiousness) {
  population <- round(geography$locations$population)
  n <- sum(population)
  # agents are stored in location order, so each location is a contiguous block
  location <- rep(seq_along(population), population)
  relative_infectiousness <- infectiousness(n)
  stopifnot(
    "`infectiousness` must return non-negative numbers" =
      checkmate::test_numeric(relative_infectiousness, lower = 0, any.missing = FALSE, len = n)
  )
  list(
    population = population,
    location = location,
    household = .make_households(location, household_size),
    age = .sample_ages(n, population_age),
    sex = ifelse(stats::runif(n) < prob_male, "m", "f"),
    relative_infectiousness = relative_infectiousness
  )
}

#' Run the agent-based model once and build the transmission chain
#'
#' @inheritParams sim_outbreak_abm
#' @param agents The agents, from [.abm_population()].
#' @param max_infections Stop at the end of the first day on which cumulative
#' infections exceed this.
#'
#' @return A list with the `chain` (for [chain_to_outbreak()]), the `daily`
#' series, the `last_day` simulated, the number of infections and whether the
#' maximum was reached.
#' @keywords internal
.sim_abm <- function(geography,
                     agents,
                     household_hazard,
                     community_hazard,
                     latent_period,
                     infectious_period,
                     interventions,
                     n_seeds,
                     contacts_named,
                     max_days,
                     max_infections) {
  locations <- geography$locations
  mobility <- geography$mobility
  n_locations <- nrow(locations)
  population <- agents$population
  n <- sum(population)
  location <- agents$location
  household <- agents$household
  n_households <- max(household)
  age <- agents$age
  sex <- agents$sex
  relative_infectiousness <- agents$relative_infectiousness

  exposure_day <- rep(NA_real_, n)
  onset_day <- rep(NA_real_, n)
  recovery_day <- rep(NA_real_, n)
  infector <- rep(NA_integer_, n)
  daily_new <- matrix(0L, max_days, n_locations)
  daily_infectious <- matrix(0L, max_days, n_locations)

  seed_location <- if (is.null(geography$seed_location)) {
    sample.int(n_locations, 1, prob = population)
  } else {
    match(geography$seed_location, locations$name)
  }
  stopifnot("`n_seeds` is larger than the population of the seed location" = n_seeds <= population[seed_location])
  seeds <- which(location == seed_location)[sample.int(population[seed_location], n_seeds)]
  initial <- rep("S", n)
  initial[seeds] <- "E"
  health <- individual::CategoricalVariable$new(c("S", "E", "I", "R"), initial)

  # agents are stored in location order, so each location's agents are a
  # contiguous block of indices; household members are indexed for lookup
  is_susceptible <- rep(TRUE, n)
  is_susceptible[seeds] <- FALSE
  is_candidate <- rep(FALSE, n)
  susceptible_count <- tabulate(location[is_susceptible], nbins = n_locations)
  location_first <- cumsum(population) - population + 1L
  members_by_household <- order(household)
  household_count <- tabulate(household, nbins = n_households)
  household_first <- cumsum(household_count) - household_count + 1L

  to_recovered <- individual::TargetedEvent$new(n)
  to_recovered$add_listener(function(t, target) health$queue_update("R", target))
  to_infectious <- individual::TargetedEvent$new(n)
  to_infectious$add_listener(function(t, target) health$queue_update("I", target))

  # someone exposed on day d is infectious from day ceiling(d + latent) and
  # recovered from day ceiling(d + latent + infectious period); events fire
  # after the day's transmission and updates apply at the end of the day, so
  # the delays are one day shorter than those day counts
  schedule_progression <- function(uids, day) {
    # a Bitset holds its members in ascending order and TargetedEvent matches
    # a vector of delays to that order, so the agents must be sorted first
    uids <- sort(uids)
    latent_days <- latent_period(length(uids))
    infectious_days <- infectious_period(length(uids))
    delay_infectious <- pmax(ceiling(latent_days) - 1, 0)
    delay_recovered <- pmax(ceiling(latent_days + infectious_days) - 1, delay_infectious + 1)
    exposure_day[uids] <<- day
    onset_day[uids] <<- day + delay_infectious + 1
    recovery_day[uids] <<- day + delay_recovered + 1
    to_infectious$schedule(individual::Bitset$new(n)$insert(uids), delay = delay_infectious)
    to_recovered$schedule(individual::Bitset$new(n)$insert(uids), delay = delay_recovered)
  }
  schedule_progression(seeds, 0)

  susceptibility <- function(day) {
    rel <- locations$contact_multiplier
    for (k in seq_len(if (is.null(interventions)) 0 else nrow(interventions))) {
      if (day >= interventions$start[k] && day < interventions$end[k]) {
        hit <- if (interventions$location[k] == "all") seq_len(n_locations) else match(interventions$location[k], locations$name)
        rel[hit] <- rel[hit] * (1 - interventions$reduction[k])
      }
    }
    rel
  }

  stopped <- FALSE
  last_day <- max_days - 1
  infection <- function(t) {
    day <- t - 1
    if (stopped) return(invisible(NULL))
    infectious_now <- health$get_index_of("I")$to_vector()
    daily_infectious[t, ] <<- tabulate(location[infectious_now], nbins = n_locations)
    if (length(infectious_now) == 0) {
      if (health$get_size_of("E") == 0) {
        stopped <<- TRUE
        last_day <<- day
      }
      return(invisible(NULL))
    }
    weight <- relative_infectiousness[infectious_now]
    household_load <- numeric(n_households)
    by_household <- rowsum(weight, household[infectious_now])
    household_load[as.integer(rownames(by_household))] <- by_household[, 1]
    location_load <- numeric(n_locations)
    by_location <- rowsum(weight, location[infectious_now])
    location_load[as.integer(rownames(by_location))] <- by_location[, 1]
    location_load <- location_load / population
    # community hazard for a resident of each location, from where they mix
    community_by_location <- community_hazard * as.numeric(mobility %*% location_load)
    rel_location <- susceptibility(day)

    # susceptible agents with an infectious household member face both
    # hazards, so they are drawn one by one
    loaded <- which(household_load > 0)
    candidates <- members_by_household[sequence(household_count[loaded], from = household_first[loaded])]
    candidates <- candidates[is_susceptible[candidates]]
    rel <- rel_location[location[candidates]]
    h_household <- rel * household_hazard * household_load[household[candidates]]
    h_community <- rel * community_by_location[location[candidates]]
    infected_candidates <- stats::runif(length(candidates)) < 1 - exp(-(h_household + h_community))

    # every other susceptible agent faces only their location's community
    # hazard, so the number infected in each location is binomial and they
    # are a uniform sample; this gives each agent the same daily probability
    # of infection as drawing them one by one, at a cost that grows with the
    # epidemic rather than the population
    is_candidate[candidates] <<- TRUE
    others <- susceptible_count - tabulate(location[candidates], nbins = n_locations)
    n_community_only <- stats::rbinom(n_locations, others, 1 - exp(-rel_location * community_by_location))
    community_only <- .sample_susceptible(
      n_community_only, location_first, population, is_susceptible, is_candidate
    )
    is_candidate[candidates] <<- FALSE

    newly <- c(candidates[infected_candidates], community_only)
    if (length(newly) > 0) {
      is_susceptible[newly] <<- FALSE
      susceptible_count <<- susceptible_count - tabulate(location[newly], nbins = n_locations)
      # attribute each infection to the household or the community in
      # proportion to the two hazards, then draw the infector
      share_household <- c(
        h_household[infected_candidates] / (h_household[infected_candidates] + h_community[infected_candidates]),
        rep(0, length(community_only))
      )
      from_household <- stats::runif(length(newly)) < share_household
      who <- integer(length(newly))
      if (any(from_household)) {
        who[from_household] <- .pick_infectors(household[newly[from_household]], infectious_now, household[infectious_now], weight)
      }
      if (any(!from_household)) {
        targets <- which(!from_household)
        source <- integer(length(targets))
        for (a in unique(location[newly[targets]])) {
          here <- location[newly[targets]] == a
          source[here] <- sample.int(n_locations, sum(here), replace = TRUE, prob = mobility[a, ] * location_load)
        }
        who[targets] <- .pick_infectors(source, infectious_now, location[infectious_now], weight)
      }
      infector[newly] <<- who
      daily_new[t, ] <<- tabulate(location[newly], nbins = n_locations)
      health$queue_update("E", newly)
      schedule_progression(newly, day)
    }
    if (sum(!is.na(exposure_day)) > max_infections) {
      stopped <<- TRUE
      last_day <<- day
    }
  }

  individual::simulation_loop(
    variables = list(health),
    events = list(to_infectious, to_recovered),
    processes = list(infection),
    timesteps = max_days
  )

  days_run <- seq_len(last_day + 1)
  list(
    chain = .abm_chain(location, household, age, sex, infector, exposure_day, onset_day, recovery_day,
                       locations$name, contacts_named, last_day),
    daily = data.frame(
      day = rep(days_run - 1, n_locations),
      location = rep(locations$name, each = length(days_run)),
      new_infections = as.vector(daily_new[days_run, , drop = FALSE]),
      infectious = as.vector(daily_infectious[days_run, , drop = FALSE])
    ),
    last_day = last_day,
    n_infected = sum(!is.na(exposure_day)),
    hit_max = sum(!is.na(exposure_day)) > max_infections,
    infections = data.frame(
      agent = which(!is.na(exposure_day)),
      infector = infector[!is.na(exposure_day)],
      exposure_day = exposure_day[!is.na(exposure_day)],
      onset_day = onset_day[!is.na(exposure_day)],
      recovery_day = recovery_day[!is.na(exposure_day)]
    )
  )
}

#' Sample agents uniformly from those still eligible in each location
#'
#' Locations are contiguous blocks of agent indices. Agents are drawn by
#' rejection sampling within the block, which is fast while a fair share of
#' the block is eligible; when few are left, the eligible agents are listed
#' directly. Rejection sampling repeats until enough are accepted, so it is a
#' loop.
#'
#' @param n_per_location Number to draw in each location.
#' @param location_first First agent index of each location.
#' @param population Number of agents in each location.
#' @param is_susceptible,is_candidate Logical vectors over all agents; agents
#' are eligible if susceptible and not a candidate. They are checked only for
#' proposed agents, so the cost does not grow with the population.
#' @return Agent indices.
#' @keywords internal
.sample_susceptible <- function(n_per_location, location_first, population, is_susceptible, is_candidate) {
  drawn <- vector("list", length(n_per_location))
  for (a in which(n_per_location > 0)) {
    k <- n_per_location[a]
    offset <- location_first[a] - 1L
    accepted <- integer(0)
    tries <- 0
    while (length(accepted) < k && tries < 20) {
      proposal <- offset + sample.int(population[a], ceiling(2 * (k - length(accepted))) + 10L, replace = TRUE)
      accepted <- unique(c(accepted, proposal[is_susceptible[proposal] & !is_candidate[proposal]]))
      tries <- tries + 1
    }
    if (length(accepted) < k) {
      block <- offset + seq_len(population[a])
      pool <- block[is_susceptible[block] & !is_candidate[block]]
      accepted <- pool[sample.int(length(pool), min(k, length(pool)))]
    }
    drawn[[a]] <- accepted[seq_len(min(k, length(accepted)))]
  }
  unlist(drawn, use.names = FALSE)
}

#' Household of each agent, with sizes from `household_size`, drawn within
#' each location
#' @keywords internal
.make_households <- function(location, household_size) {
  household <- integer(length(location))
  offset <- 0L
  for (a in sort(unique(location))) {
    members <- which(location == a)
    sizes <- household_size(length(members))
    stopifnot("`household_size` must return whole numbers of at least 1" =
                checkmate::test_integerish(sizes, lower = 1, any.missing = FALSE, len = length(members)))
    sizes <- sizes[seq_len(which(cumsum(sizes) >= length(members))[1])]
    ids <- offset + rep(seq_along(sizes), sizes)[seq_along(members)]
    household[members[sample.int(length(members))]] <- ids
    offset <- offset + length(sizes)
  }
  household
}

#' Ages of agents from a population age range or structure, as in
#' [sim_outbreak()]
#' @keywords internal
.sample_ages <- function(n, population_age) {
  if (is.data.frame(population_age)) {
    age_groups <- apply(population_age, MARGIN = 1, function(x) x[1]:x[2])
    sample_weight <- rep(population_age$proportion, times = lengths(age_groups)) /
      rep(lengths(age_groups), times = lengths(age_groups))
    sample(unlist(age_groups), size = n, replace = TRUE, prob = sample_weight)
  } else {
    sample(population_age[1]:population_age[2], size = n, replace = TRUE)
  }
}

#' Draw one infector for each target from groups of infectious agents,
#' weighted by infectiousness
#'
#' @param group_of_target Group (household or location) each target is
#' infected from.
#' @param infectious Agent indices of infectious agents.
#' @param group_of_infectious Group of each infectious agent.
#' @param weight Infectiousness of each infectious agent.
#' @keywords internal
.pick_infectors <- function(group_of_target, infectious, group_of_infectious, weight) {
  ord <- order(group_of_infectious)
  infectious <- infectious[ord]
  group_sorted <- group_of_infectious[ord]
  cum_weight <- cumsum(weight[ord])
  first <- match(group_of_target, group_sorted)
  last <- length(group_sorted) - match(group_of_target, rev(group_sorted)) + 1
  base <- c(0, cum_weight)[first]
  target <- base + stats::runif(length(group_of_target)) * (cum_weight[last] - base)
  infectious[pmin(findInterval(target, cum_weight, left.open = TRUE) + 1, last)]
}

#' Transmission chain for [chain_to_outbreak()] from an agent-based run
#'
#' Cases are agents with onset by `last_day`, whose whole line of infection
#' back to a seed are cases; onset is the start of infectiousness, with a
#' uniform time within the day. Each case names its uninfected household
#' members and a Poisson number of uninfected agents from its location,
#' contacted at a uniform time in its infectious period.
#' @keywords internal
.abm_chain <- function(location, household, age, sex, infector, exposure_day, onset_day, recovery_day,
                       location_names, contacts_named, last_day) {
  is_case <- !is.na(onset_day) & onset_day <= last_day
  # each pass drops cases whose infector was dropped
  repeat {
    orphan <- is_case & !is.na(infector) & !is_case[pmax(infector, 1L)]
    if (!any(orphan)) break
    is_case[orphan] <- FALSE
  }
  cases <- which(is_case)
  case_row <- match(seq_along(location), cases)

  # uninfected agents sorted by household and by location, so each case's
  # household members and community pool are contiguous blocks
  uninfected <- which(is.na(exposure_day))
  by_household <- uninfected[order(household[uninfected])]
  household_count <- tabulate(household[by_household], nbins = max(household))
  household_start <- cumsum(household_count) - household_count + 1L
  by_location <- uninfected[order(location[uninfected])]
  location_count <- tabulate(location[by_location], nbins = length(location_names))
  location_start <- cumsum(location_count) - location_count

  n_household <- household_count[household[cases]]
  household_case <- rep(seq_along(cases), n_household)
  household_agent <- by_household[sequence(n_household, from = household_start[household[cases]])]
  # community contacts are drawn with replacement; repeats are rare in a
  # large pool
  n_community <- stats::rpois(length(cases), contacts_named) * (location_count[location[cases]] > 0)
  community_case <- rep(seq_along(cases), n_community)
  pool_size <- location_count[location[cases]][community_case]
  community_agent <- by_location[location_start[location[cases]][community_case] + 1L +
                                   floor(stats::runif(length(community_case)) * pool_size)]

  contact_of <- c(household_case, community_case)
  contact_agent <- c(household_agent, community_agent)
  start <- onset_day[cases][contact_of]
  end <- recovery_day[cases][contact_of]
  contact_time <- start + stats::runif(length(contact_of)) * pmax(end - start, 1)

  rbind(
    data.frame(
      id = seq_along(cases),
      infector = case_row[infector[cases]],
      infected = rep(TRUE, length(cases)),
      time = onset_day[cases] + stats::runif(length(cases)),
      location = location_names[location[cases]],
      age = age[cases],
      sex = sex[cases],
      agent = cases
    ),
    data.frame(
      id = length(cases) + seq_along(contact_agent),
      infector = contact_of,
      infected = rep(FALSE, length(contact_agent)),
      time = contact_time,
      location = location_names[location[contact_agent]],
      age = age[contact_agent],
      sex = sex[contact_agent],
      agent = contact_agent
    )
  )
}
