locations <- data.frame(
  name = c("North", "Central", "South"),
  population = c(40000, 120000, 60000),
  x = c(0, 3, 7),
  y = c(5, 0, -2)
)

test_that("create_geography works as expected", {
  geography <- create_geography(locations)
  expect_s3_class(geography, "simulist_geography")
  expect_identical(geography$locations$name, locations$name)
  expect_equal(unname(rowSums(geography$mobility)), rep(1, 3))
  expect_equal(unname(diag(geography$mobility)), rep(0.8, 3))
  expect_identical(dimnames(geography$mobility), list(locations$name, locations$name))
  expect_identical(geography$locations$contact_multiplier, rep(1, 3))
  expect_null(geography$seed_location)
})

test_that("create_geography scales contact multipliers to a maximum of 1", {
  locations$contact_multiplier <- c(1, 4, 2)
  geography <- create_geography(locations)
  expect_identical(geography$locations$contact_multiplier, c(0.25, 1, 0.5))
})

test_that("create_geography works with one location", {
  geography <- create_geography(locations[1, c("name", "population")])
  expect_identical(unname(geography$mobility), matrix(1, 1, 1))
})

test_that("create_geography accepts a mobility matrix", {
  mobility <- matrix(c(0.9, 0.1, 0, 0, 1, 0, 0.2, 0.3, 0.5), 3, byrow = TRUE)
  geography <- create_geography(locations, mobility = mobility)
  expect_equal(unname(geography$mobility), mobility)
})

test_that("create_geography errors with invalid input", {
  expect_error(create_geography(locations[, c("name", "x", "y")]), regexp = "population")
  expect_error(
    create_geography(transform(locations, name = c("A", "A", "B"))),
    regexp = "duplicated"
  )
  expect_error(
    create_geography(transform(locations, population = c(0, 1, 1))),
    regexp = "positive population"
  )
  expect_error(create_geography(locations[, c("name", "population")]), regexp = "x")
  expect_error(
    create_geography(locations, mobility = matrix(0.5, 3, 3)),
    regexp = "sum to 1"
  )
  expect_error(create_geography(locations, prob_stay = 1.5))
  expect_error(create_geography(locations, seed_location = "East"), regexp = "East")
})

test_that("print.simulist_geography works as expected", {
  expect_output(print(create_geography(locations)), regexp = "3 locations")
})

test_that("mobility_gravity weights fall with distance and rise with population", {
  equal_pop <- transform(locations, population = 1000)
  weights <- mobility_gravity(kernel_km = 2)(equal_pop)
  # North is nearer Central than South
  expect_gt(weights[1, 2], weights[1, 3])
  power <- mobility_gravity(kernel_km = 2, kernel = "power")(equal_pop)
  expect_gt(power[1, 2], power[1, 3])
  weights <- mobility_gravity(kernel_km = 2)(locations)
  expect_true(all(weights >= 0))
  expect_error(mobility_gravity(kernel_km = 0), regexp = "positive")
})

test_that("mobility_radiation gives non-negative weights", {
  weights <- mobility_radiation()(locations)
  expect_true(all(weights >= 0))
  expect_true(all(diag(weights) == 0))
  geography <- create_geography(locations, mobility = mobility_radiation())
  expect_equal(unname(rowSums(geography$mobility)), rep(1, 3))
})

test_that("intervention_schedule works as expected", {
  prob_infection <- intervention_schedule(
    base = 0.4,
    interventions = data.frame(
      location = c("Central", "all"),
      start = c(30, 60),
      end = c(Inf, 90),
      reduction = c(0.5, 0.25)
    )
  )
  expect_equal(
    prob_infection(
      time = c(10, 45, 70, 70, 100),
      location = c("Central", "Central", "North", "Central", "North")
    ),
    c(0.4, 0.2, 0.3, 0.15, 0.4)
  )
})

test_that("intervention_schedule works with a probability per location", {
  prob_infection <- intervention_schedule(base = c(North = 0.2, South = 0.6))
  expect_equal(prob_infection(time = c(1, 2), location = c("South", "North")), c(0.6, 0.2))
  expect_error(prob_infection(time = 1, location = "Central"), regexp = "no probability")
})

test_that("intervention_schedule errors with invalid input", {
  expect_error(intervention_schedule(base = 1.2))
  expect_error(
    intervention_schedule(
      base = 0.5,
      interventions = data.frame(location = "all", start = 10, end = 5, reduction = 0.5)
    ),
    regexp = "before"
  )
})
