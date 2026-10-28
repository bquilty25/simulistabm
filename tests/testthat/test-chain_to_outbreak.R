chain <- data.frame(
  id = 1:5,
  infector = c(NA, NA, 1, 1, 3),
  infected = c(TRUE, TRUE, TRUE, FALSE, TRUE),
  time = c(0, 1.5, 4.2, 3.1, 9.8),
  location = c("North", "South", "North", "North", "South")
)

test_that("chain_to_outbreak works as expected", {
  set.seed(1)
  outbreak <- chain_to_outbreak(chain)
  linelist <- outbreak$linelist
  contacts <- outbreak$contacts
  expect_identical(nrow(linelist), 4L)
  expect_identical(linelist$id, c(1L, 2L, 3L, 5L))
  expect_identical(linelist$location, c("North", "South", "North", "South"))
  expect_identical(linelist$date_onset, as.Date("2023-01-01") + c(0, 1.5, 4.2, 9.8))
  # every row except the two index cases is a contact
  expect_identical(nrow(contacts), 3L)
  expect_identical(contacts$was_case, c(TRUE, FALSE, TRUE))
  expect_identical(contacts$from_location, c("North", "North", "North"))
  expect_identical(contacts$to_location, c("North", "North", "South"))
  expect_identical(contacts$from[c(1, 2)], rep(linelist$case_name[1], 2))
  expect_identical(contacts$from[3], linelist$case_name[3])
})

test_that("chain_to_outbreak has the same columns as sim_outbreak", {
  set.seed(1)
  simulated <- suppressWarnings(sim_outbreak(outbreak_size = c(10, 100)))
  set.seed(1)
  converted <- chain_to_outbreak(chain[, c("id", "infector", "infected", "time")])
  expect_identical(colnames(converted$linelist), colnames(simulated$linelist))
  expect_identical(colnames(converted$contacts), colnames(simulated$contacts))
  expect_s3_class(messy_linelist(converted$linelist), "data.frame")
})

test_that("chain_to_outbreak keeps the chain's ages and sexes", {
  chain$age <- c(5L, 40L, 70L, 22L, 15L)
  chain$sex <- c("f", "m", "f", "m", "f")
  set.seed(1)
  outbreak <- chain_to_outbreak(chain)
  expect_identical(outbreak$linelist$age, c(5L, 40L, 70L, 15L))
  expect_identical(outbreak$linelist$sex, c("f", "m", "f", "f"))
  expect_identical(outbreak$contacts$age, c(70L, 22L, 15L))
})

test_that("chain_to_outbreak accepts character infection status", {
  chain$infected <- ifelse(chain$infected, "infected", "contact")
  set.seed(1)
  expect_identical(nrow(chain_to_outbreak(chain)$linelist), 4L)
})

test_that("chain_to_outbreak errors with an invalid chain", {
  expect_error(chain_to_outbreak(chain[, c("id", "time")]), regexp = "infector")
  expect_error(
    chain_to_outbreak(transform(chain, infector = c(1, NA, 1, 1, 3))),
    regexp = "trace back"
  )
  expect_error(
    chain_to_outbreak(transform(chain, infector = c(NA, NA, 5, 1, 3))),
    regexp = "trace back"
  )
  expect_error(
    chain_to_outbreak(transform(chain, infector = c(NA, NA, 1, 1, 4))),
    regexp = "infected row"
  )
  expect_error(
    chain_to_outbreak(transform(chain, infected = c(FALSE, TRUE, TRUE, FALSE, TRUE))),
    regexp = "must be infected"
  )
  expect_error(chain_to_outbreak(transform(chain, time = -1)))
  expect_error(chain_to_outbreak(transform(chain, sex = "x")))
})
