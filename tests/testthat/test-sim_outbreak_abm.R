abm_geography <- create_geography(
  data.frame(
    name = c("North", "Central", "South"),
    population = c(1500, 4000, 2500),
    x = c(0, 3, 7),
    y = c(5, 0, -2)
  ),
  seed_location = "Central"
)

test_that("sim_outbreak_abm has the same structure as sim_outbreak", {
  set.seed(1)
  outbreak <- suppressWarnings(sim_outbreak_abm(
    abm_geography, household_hazard = 0.1, community_hazard = 0.35, outbreak_size = c(20, 500)
  ))
  set.seed(1)
  reference <- suppressWarnings(sim_outbreak(outbreak_size = c(10, 50), geography = abm_geography))
  expect_identical(colnames(outbreak$linelist), colnames(reference$linelist))
  expect_identical(colnames(outbreak$contacts), colnames(reference$contacts))
  expect_gte(nrow(outbreak$linelist), 20)
  expect_true(all(outbreak$linelist$location %in% abm_geography$locations$name))
  expect_s3_class(attr(outbreak, "daily"), "data.frame")
  expect_true(all(outbreak$linelist$date_onset <= attr(outbreak, "complete_until") + 1))
  expect_s3_class(messy_linelist(outbreak$linelist), "data.frame")
  expect_s3_class(truncate_linelist(outbreak$linelist), "data.frame")
  expect_s3_class(censor_linelist(outbreak$linelist, interval = "daily"), "data.frame")
})

test_that("sim_outbreak_abm contacts follow the transmission chain", {
  set.seed(2)
  outbreak <- suppressWarnings(sim_outbreak_abm(
    abm_geography, household_hazard = 0.1, community_hazard = 0.35, outbreak_size = c(20, 500)
  ))
  linelist <- outbreak$linelist
  contacts <- outbreak$contacts
  # every case except the seed is named exactly once, by a case
  case_contacts <- contacts[contacts$was_case, ]
  expect_identical(nrow(case_contacts), nrow(linelist) - 1L)
  expect_true(all(case_contacts$from %in% linelist$case_name))
  expect_setequal(c(case_contacts$to, setdiff(linelist$case_name, case_contacts$to)), linelist$case_name)
  # the infector named in the contacts table had onset before the case
  onset <- stats::setNames(linelist$date_onset, linelist$case_name)
  expect_true(all(onset[case_contacts$from] <= onset[case_contacts$to]))
})

test_that("sim_outbreak_abm infections happen while the infector is infectious", {
  agents <- .abm_population(
    abm_geography, function(n) 1L + stats::rpois(n, 3.5), c(0, 90), 0.5, function(n) stats::rgamma(n, 1, 1)
  )
  set.seed(8)
  run <- .sim_abm(
    abm_geography, agents, household_hazard = 0.1, community_hazard = 0.35,
    latent_period = function(n) stats::rlnorm(n, 1, 0.3), infectious_period = function(n) stats::rlnorm(n, 1.5, 0.4),
    interventions = NULL, n_seeds = 3, contacts_named = 2, max_days = 365, max_infections = 2000
  )
  x <- run$infections
  infector <- match(x$infector, x$agent)
  has_infector <- !is.na(infector)
  expect_gt(sum(has_infector), 100)
  expect_true(all(x$exposure_day[has_infector] >= x$onset_day[infector[has_infector]]))
  expect_true(all(x$exposure_day[has_infector] < x$recovery_day[infector[has_infector]]))
})

test_that("sim_outbreak_abm outbreaks deplete and end", {
  set.seed(3)
  outbreak <- sim_outbreak_abm(
    abm_geography, household_hazard = 0.1, community_hazard = 0.35, outbreak_size = c(20, 1e6)
  )
  attack <- table(outbreak$linelist$location) / abm_geography$locations$population[match(names(table(outbreak$linelist$location)), abm_geography$locations$name)]
  expect_true(all(attack < 1))
  expect_lt(attr(outbreak, "complete_until"), as.Date("2023-01-01") + 729)
})

test_that("sim_outbreak_abm applies interventions by location", {
  # nobody moves, and North is fully protected from the start, so there are
  # no cases in North even though it is next to the seed location
  geography <- create_geography(
    abm_geography$locations,
    mobility = diag(3),
    seed_location = "Central"
  )
  set.seed(4)
  outbreak <- sim_outbreak_abm(
    geography, household_hazard = 0.1, community_hazard = 0.35, outbreak_size = c(20, 1e6),
    interventions = data.frame(location = "North", start = 0, end = Inf, reduction = 1)
  )
  expect_false(any(outbreak$linelist$location == "North"))
})

test_that("sim_outbreak_abm gives agents ages from a population age structure", {
  # as in simulist, the last age limit is the oldest age
  population_age <- data.frame(age_limit = c(1, 20, 60, 90), proportion = c(0.5, 0.4, 0.1, 0))
  set.seed(7)
  outbreak <- suppressWarnings(sim_outbreak_abm(
    abm_geography, household_hazard = 0.1, community_hazard = 0.35, outbreak_size = c(20, 300),
    population_age = population_age,
    hosp_risk = data.frame(age_limit = c(1, 60), risk = c(0.1, 0.5))
  ))
  expect_true(all(outbreak$linelist$age >= 1 & outbreak$linelist$age <= 90))
  expect_true(all(outbreak$contacts$age >= 1))
})

test_that("sim_outbreak_abm with no transmission gives only the seeds", {
  set.seed(5)
  outbreak <- sim_outbreak_abm(
    abm_geography, household_hazard = 0, community_hazard = 0, n_seeds = 3, outbreak_size = c(1, 100)
  )
  expect_identical(nrow(outbreak$linelist), 3L)
  expect_true(all(outbreak$linelist$location == "Central"))
  expect_false(any(outbreak$contacts$was_case))
})

test_that("sim_outbreak_abm warns at the maximum outbreak size", {
  set.seed(6)
  expect_warning(
    sim_outbreak_abm(abm_geography, household_hazard = 0.1, community_hazard = 0.5, outbreak_size = c(5, 50)),
    regexp = "maximum outbreak size"
  )
})

test_that("sim_outbreak_abm errors with invalid input", {
  expect_error(sim_outbreak_abm(list()), regexp = "simulist_geography")
  expect_error(sim_outbreak_abm(abm_geography, household_hazard = -1))
  expect_error(
    sim_outbreak_abm(abm_geography, interventions = data.frame(location = "East", start = 0, end = 1, reduction = 0.5))
  )
  expect_error(sim_outbreak_abm(abm_geography, n_seeds = 1e5), regexp = "seed location")
})
