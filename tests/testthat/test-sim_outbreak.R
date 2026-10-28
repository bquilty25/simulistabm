test_that("sim_outbreak works as expected with defaults", {
  set.seed(1)
  expect_snapshot(sim_outbreak())
})


contact_distribution <- function(x) stats::dpois(x = x, lambda = 2)
infectious_period <- function(n) stats::rgamma(n = n, shape = 1, scale = 1)
onset_to_hosp <- function(n) {
  stats::rlnorm(n = n, meanlog = 0.947, sdlog = 1.628)
}
onset_to_death <- function(n) {
  stats::rlnorm(n = n, meanlog = 2.863, sdlog = 0.534)
}

test_that("sim_outbreak works as expected", {
  set.seed(1)
  expect_snapshot(
    sim_outbreak(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      onset_to_hosp = onset_to_hosp,
      onset_to_death = onset_to_death
    )
  )
})

test_that("sim_outbreak works as expected with anonymised", {
  set.seed(1)
  expect_snapshot(
    sim_outbreak(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      onset_to_hosp = onset_to_hosp,
      onset_to_death = onset_to_death,
      anonymise = TRUE
    )
  )
})

test_that("sim_outbreak works as expected with age-strat risks", {
  age_dep_hosp_risk <- data.frame(
    age_limit = c(0, 5, 80),
    risk = c(0.1, 0.05, 0.2)
  )
  age_dep_hosp_death_risk <- data.frame(
    age_limit = c(0, 5, 80),
    risk = c(0.1, 0.05, 0.2)
  )
  age_dep_non_hosp_death_risk <- data.frame(
    age_limit = c(0, 5, 80),
    risk = c(0.05, 0.025, 0.1)
  )
  set.seed(1)
  expect_snapshot(
    sim_outbreak(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      onset_to_hosp = onset_to_hosp,
      onset_to_death = onset_to_death,
      hosp_risk = age_dep_hosp_risk,
      hosp_death_risk = age_dep_hosp_death_risk,
      non_hosp_death_risk = age_dep_non_hosp_death_risk
    )
  )
})

test_that("sim_outbreak works as expected with age structure", {
  age_struct <- data.frame(
    age_limit = c(1, 5, 80, 90),
    proportion = c(0.1, 0.7, 0.2, 0)
  )
  set.seed(1)
  expect_snapshot(
    sim_outbreak(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      onset_to_hosp = onset_to_hosp,
      onset_to_death = onset_to_death,
      population_age = age_struct
    )
  )
})

geography <- create_geography(
  data.frame(
    name = c("North", "Central", "South"),
    population = c(40000, 120000, 60000),
    x = c(0, 3, 7),
    y = c(5, 0, -2)
  ),
  seed_location = "Central"
)

test_that("sim_outbreak adds location columns with a geography", {
  set.seed(1)
  outbreak <- suppressWarnings(sim_outbreak(
    contact_distribution = function(x) stats::dpois(x = x, lambda = 4),
    infectious_period = infectious_period,
    outbreak_size = c(10, 500),
    geography = geography
  ))
  linelist <- outbreak$linelist
  contacts <- outbreak$contacts
  expect_identical(colnames(linelist)[which(colnames(linelist) == "age") + 1], "location")
  expect_true(all(linelist$location %in% geography$locations$name))
  expect_true(all(c("from_location", "to_location") %in% colnames(contacts)))
  # the location of each case matches its location in the contacts table
  case_contacts <- contacts[contacts$was_case, ]
  expect_identical(
    linelist$location[match(case_contacts$to, linelist$case_name)],
    case_contacts$to_location
  )
  # each contact's from_location is where the case who named them lives
  expect_identical(
    linelist$location[match(contacts$from, linelist$case_name)],
    contacts$from_location
  )
})

test_that("sim_outbreak has no location columns without a geography", {
  set.seed(1)
  outbreak <- sim_outbreak(outbreak_size = c(10, 200))
  expect_false("location" %in% colnames(outbreak$linelist))
  expect_false(any(c("from_location", "to_location") %in% colnames(outbreak$contacts)))
})

test_that("sim_linelist and sim_contacts accept a geography", {
  set.seed(1)
  linelist <- suppressWarnings(sim_linelist(
    contact_distribution = function(x) stats::dpois(x = x, lambda = 4),
    outbreak_size = c(10, 200),
    geography = geography
  ))
  expect_true("location" %in% colnames(linelist))
  set.seed(1)
  contacts <- suppressWarnings(sim_contacts(
    contact_distribution = function(x) stats::dpois(x = x, lambda = 4),
    outbreak_size = c(10, 200),
    geography = geography
  ))
  expect_true(all(c("from_location", "to_location") %in% colnames(contacts)))
})

test_that("sim_outbreak errors with an invalid geography or prob_infection", {
  expect_error(sim_outbreak(geography = list()), regexp = "simulist_geography")
  expect_error(
    sim_outbreak(prob_infection = function(t) 0.5),
    regexp = "prob_infection supplied must have 2 arguments"
  )
})
