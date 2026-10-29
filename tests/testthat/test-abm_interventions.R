geography <- create_geography(
  data.frame(
    name = c("Ashford", "Brayton", "Cole"),
    population = c(8000, 4000, 2000),
    x = c(0, 12, 30),
    y = 0
  ),
  mobility = mobility_gravity(kernel_km = 10),
  prob_stay = 0.85,
  seed_location = "Ashford"
)

run_abm <- function(..., seed = 3) {
  set.seed(seed)
  suppressWarnings(sim_outbreak_abm(
    geography,
    household_hazard = 0.05,
    community_hazard = 0.3,
    latent_period = function(n) rep(2, n),
    reporting_delay = function(n) stats::rgamma(n, shape = 2, rate = 0.5),
    outbreak_size = c(50, 1e5),
    ...
  ))
}

test_that("vaccination() and contact_tracing() check their inputs", {
  expect_s3_class(vaccination("ring", efficacy = 0.9), "simulist_vaccination")
  expect_s3_class(contact_tracing(), "simulist_contact_tracing")
  expect_error(vaccination("mass", efficacy = 1.2, doses_per_day = 10))
  expect_error(vaccination("mass", efficacy = 0.5))
  expect_error(contact_tracing(prob_traced = -0.1))
  expect_error(sim_outbreak_abm(geography, vaccination = list(type = "ring")))
})

test_that("a fully effective vaccine prevents infection once protective", {
  outbreak <- run_abm(
    vaccination = vaccination("mass", efficacy = 1, start_after_reports = 30, doses_per_day = 300, delay_to_protection = 7)
  )
  cases <- attr(outbreak, "case_interventions")
  onset <- outbreak$linelist$date_onset[match(cases$id, outbreak$linelist$id)]
  # exposure is two days before onset; dates are floored, so allow a day
  expect_gt(sum(!is.na(cases$date_vaccination)), 0)
  expect_identical(sum(cases$date_vaccination + 7 < onset - 3, na.rm = TRUE), 0L)
})

test_that("mass vaccination respects its start, capacity and affected areas", {
  outbreak <- run_abm(
    vaccination = vaccination("mass", efficacy = 0.8, start_after_reports = 30, doses_per_day = 300)
  )
  doses <- attr(outbreak, "vaccinations")
  start <- attr(outbreak, "interventions")$start_date
  expect_identical(sum(doses$doses[doses$date < start]), 0L)
  expect_true(all(tapply(doses$doses, doses$date, sum) <= 300))
  # the start is the day after the 30th report
  report_day <- floor(as.numeric(outbreak$linelist$date_reporting))
  expect_identical(as.numeric(start), sort(report_day)[30] + 1)
  first_report <- tapply(report_day, outbreak$linelist$location, min)
  vaccinated <- doses[doses$doses > 0, ]
  expect_true(all(as.numeric(vaccinated$date) >= first_report[vaccinated$location]))
})

test_that("ring vaccination follows reported cases and reduces the outbreak", {
  base <- run_abm()
  ring <- run_abm(
    vaccination = vaccination("ring", efficacy = 0.95, start_after_reports = 20, coverage = 0.9, delay = 1)
  )
  doses <- attr(ring, "vaccinations")
  expect_gt(sum(doses$doses), 0)
  expect_identical(sum(doses$doses[doses$date < attr(ring, "interventions")$start_date]), 0L)
  expect_lt(nrow(ring$linelist), nrow(base$linelist))
})

test_that("contact tracing quarantines contacts and records tracing", {
  base <- run_abm()
  traced <- run_abm(
    contact_tracing = contact_tracing(prob_traced = 0.8, delay = 1, start_after_reports = 20)
  )
  expect_lt(nrow(traced$linelist), nrow(base$linelist))
  expect_true(all(c("date_traced", "date_vaccination") %in% colnames(traced$contacts)))
  expect_setequal(
    unique(traced$contacts$status),
    c("case", "under_followup", "lost_to_followup", "unknown")
  )
  # contacts with a tracing date are under follow-up
  uninfected <- traced$contacts[!traced$contacts$was_case, ]
  expect_true(all(uninfected$status[!is.na(uninfected$date_traced)] == "under_followup"))
  # cases in quarantine at onset are reported that day
  cases <- attr(traced, "case_interventions")
  expect_gt(sum(cases$quarantined_at_onset), 0)
  linelist <- traced$linelist[match(cases$id, traced$linelist$id), ]
  quarantined <- cases$quarantined_at_onset
  expect_true(all(linelist$date_reporting[quarantined] == linelist$date_onset[quarantined]))
})

test_that("the counterfactual matches until the interventions start", {
  outbreak <- run_abm(
    vaccination = vaccination("mass", efficacy = 0.8, start_after_reports = 30, doses_per_day = 300),
    contact_tracing = contact_tracing(prob_traced = 0.5, delay = 2, start_after_reports = 30),
    counterfactual = TRUE
  )
  plain <- run_abm()
  cf <- attr(outbreak, "counterfactual")
  start <- min(attr(outbreak, "interventions")$start_date)
  before <- cf$date < start
  expect_identical(cf$infections_with[before], cf$infections_without[before])
  expect_lt(sum(cf$infections_with, na.rm = TRUE), sum(cf$infections_without, na.rm = TRUE))
  # the rerun without interventions is the plain simulation
  expect_identical(sum(cf$onsets_without, na.rm = TRUE), nrow(plain$linelist))
  expect_null(attr(plain, "counterfactual"))
})

test_that("paired reruns branch at the start and are capped in length", {
  outbreak <- run_abm(
    contact_tracing = contact_tracing(prob_traced = 0.5, delay = 2, start_after_reports = 30),
    counterfactual = TRUE,
    counterfactual_pairs = 2,
    counterfactual_days = 60
  )
  pairs <- attr(outbreak, "counterfactual_pairs")
  kept <- attr(outbreak, "counterfactual")
  start <- min(attr(outbreak, "interventions")$start_date)
  expect_setequal(unique(pairs$pair), 1:2)
  for (k in 1:2) {
    pair <- pairs[pairs$pair == k & pairs$date < start, ]
    expect_identical(pair$infections_with, kept$infections_with[kept$date < start])
    expect_identical(pair$infections_without, pair$infections_with)
  }
  expect_true(all(is.na(pairs$onsets_without[pairs$date > start + 61])))
})

test_that("interventions keep the simulist structure", {
  outbreak <- run_abm(
    vaccination = vaccination("ring", efficacy = 0.9, start_after_reports = 20),
    contact_tracing = contact_tracing(start_after_reports = 20)
  )
  expect_identical(
    colnames(outbreak$linelist),
    c("id", "case_name", "case_type", "sex", "age", "location", "date_onset", "date_reporting",
      "date_admission", "outcome", "date_outcome", "date_first_contact", "date_last_contact", "ct_value")
  )
  expect_true(all(outbreak$linelist$date_reporting >= outbreak$linelist$date_onset))
  expect_s3_class(truncate_linelist(outbreak$linelist, truncation_day = 14), "data.frame")
  expect_setequal(attr(outbreak, "interventions")$intervention, c("vaccination", "contact_tracing"))
})
