contact_distribution <- function(x) stats::dpois(x = x, lambda = 2)
infectious_period <- function(n) stats::rgamma(n = n, shape = 1, scale = 1)

test_that(".sim_network_bp works as expected", {
  set.seed(1)
  expect_snapshot(
    .sim_network_bp(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      max_outbreak_size = 1e4,
      config = create_config()
    )
  )
})

test_that(".sim_network_bp works as expected with no contacts", {
  contact_distribution <- function(x) stats::dpois(x = x, lambda = 1)
  set.seed(1)
  expect_snapshot(
    .sim_network_bp(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      max_outbreak_size = 1e4,
      config = create_config()
    )
  )
})

test_that(".sim_network_bp works as expected with unadjusted network", {
  set.seed(1)
  expect_snapshot(
    .sim_network_bp(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      max_outbreak_size = 1e4,
      config = create_config(network = "unadjusted")
    )
  )
})

test_that(".sim_network_bp warns as expected", {
  set.seed(2)
  expect_warning(
    .sim_network_bp(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      max_outbreak_size = 10,
      config = create_config()
    ),
    regexp = "(Number of cases exceeds maximum)*(Returning data early)"
  )
})

test_that(".sim_network_bp errors with negative infectious period", {
  infectious_period <- function(n) stats::rnorm(n = n, mean = 10, sd = 5)
  set.seed(3)
  expect_error(
    .sim_network_bp(
      contact_distribution = contact_distribution,
      infectious_period = infectious_period,
      prob_infection = 0.5,
      max_outbreak_size = 1e4,
      config = create_config()
    ),
    regexp = "(Negative infectious period sampled)"
  )
})

spatial_locations <- data.frame(
  name = c("North", "Central", "South"),
  population = c(40000, 120000, 60000),
  x = c(0, 3, 7),
  y = c(5, 0, -2)
)
fast_contacts <- function(x) stats::dpois(x = x, lambda = 4)
# fast-growing test outbreaks stop at the maximum size, which warns by design
run_bp <- function(...) suppressWarnings(.sim_network_bp(...))

test_that(".sim_network_bp adds locations with a geography", {
  geography <- create_geography(spatial_locations, seed_location = "Central")
  set.seed(1)
  chain <- run_bp(
    contact_distribution = fast_contacts,
    infectious_period = infectious_period,
    prob_infection = 0.5,
    max_outbreak_size = 500,
    config = create_config(),
    geography = geography
  )
  expect_true("location" %in% colnames(chain))
  expect_identical(chain$location[1], "Central")
  expect_true(all(chain$location %in% spatial_locations$name))
})

test_that(".sim_network_bp keeps everyone in place when nobody moves", {
  geography <- create_geography(spatial_locations, prob_stay = 1, seed_location = "South")
  set.seed(1)
  chain <- run_bp(
    contact_distribution = fast_contacts,
    infectious_period = infectious_period,
    prob_infection = 0.5,
    max_outbreak_size = 500,
    config = create_config(),
    geography = geography
  )
  expect_true(all(chain$location == "South"))
})

test_that(".sim_network_bp gives no contacts where the contact multiplier is 0", {
  geography <- create_geography(
    transform(spatial_locations, contact_multiplier = c(0, 1, 1)),
    prob_stay = 1,
    seed_location = "North"
  )
  set.seed(1)
  chain <- run_bp(
    contact_distribution = fast_contacts,
    infectious_period = infectious_period,
    prob_infection = 0.5,
    max_outbreak_size = 500,
    config = create_config(),
    geography = geography
  )
  expect_identical(nrow(chain), 1L)
})

test_that(".sim_network_bp uses a time-varying probability of infection", {
  # transmission stops completely after day 5
  prob_infection <- function(time, location) ifelse(time < 5, 0.5, 0)
  set.seed(1)
  chain <- run_bp(
    contact_distribution = fast_contacts,
    infectious_period = infectious_period,
    prob_infection = prob_infection,
    max_outbreak_size = 1e4,
    config = create_config()
  )
  expect_true(all(chain$time[chain$infected == "infected"] < 5))
  expect_false("location" %in% colnames(chain))
})

test_that(".sim_network_bp uses a location-dependent probability of infection", {
  geography <- create_geography(spatial_locations, seed_location = "Central")
  prob_infection <- intervention_schedule(base = c(North = 0.5, Central = 0.5, South = 0))
  set.seed(1)
  chain <- run_bp(
    contact_distribution = fast_contacts,
    infectious_period = infectious_period,
    prob_infection = prob_infection,
    max_outbreak_size = 500,
    config = create_config(),
    geography = geography
  )
  expect_false(any(chain$location[chain$infected == "infected"] == "South"))
  # contacts in South are still made, just never infected
  expect_true(any(chain$location == "South"))
})

test_that(".sim_network_bp errors when prob_infection returns invalid values", {
  set.seed(1)
  expect_error(
    .sim_network_bp(
      contact_distribution = fast_contacts,
      infectious_period = infectious_period,
      prob_infection = function(time, location) rep(2, length(time)),
      max_outbreak_size = 500,
      config = create_config()
    ),
    regexp = "between 0 and 1"
  )
  expect_error(
    .sim_network_bp(
      contact_distribution = fast_contacts,
      infectious_period = infectious_period,
      prob_infection = function(time, location) 0.5,
      max_outbreak_size = 500,
      config = create_config()
    ),
    regexp = "one probability"
  )
})

test_that(".sim_network_bp offspring by location follows the contact multiplier", {
  skip_on_cran()
  # with no movement, mean offspring in each location is R times its multiplier
  geography <- create_geography(
    transform(spatial_locations, contact_multiplier = c(0.5, 1, 1)),
    mobility = diag(3)
  )
  offspring_mean <- function(seed_location) {
    geography$seed_location <- seed_location
    set.seed(10)
    chains <- lapply(1:300, function(i) {
      .sim_network_bp(
        contact_distribution = function(x) stats::dpois(x = x, lambda = 1.6),
        infectious_period = infectious_period,
        prob_infection = 0.5,
        max_outbreak_size = 1e4,
        config = create_config(network = "unadjusted"),
        geography = geography
      )
    })
    # offspring of the index case only, so every count is fully observed
    mean(vapply(chains, function(x) sum(x$ancestor %in% 1 & x$infected == "infected"), numeric(1)))
  }
  expect_equal(offspring_mean("Central"), 0.8, tolerance = 0.15)
  expect_equal(offspring_mean("North"), 0.4, tolerance = 0.2)
})
