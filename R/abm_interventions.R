#' Vaccination for the agent-based model
#'
#' @description Describes a vaccination programme for [sim_outbreak_abm()].
#' The vaccine is leaky: from `delay_to_protection` days after vaccination, a
#' vaccinated agent's daily infection hazard is multiplied by
#' `1 - efficacy`, so `efficacy` is the true vaccine effectiveness against
#' infection.
#'
#' A `"mass"` campaign vaccinates `doses_per_day` agents each day, shared
#' between locations with at least one reported case in proportion to their
#' populations, and drawn uniformly from agents there who are not yet
#' vaccinated and have not been reported as cases.
#'
#' `"ring"` vaccination follows each case reported once the programme has
#' started: the case's household members, the agents it infected before it was
#' reported, and their household members are each vaccinated with probability
#' `coverage`, `delay` days after the report.
#'
#' @param type `"mass"` or `"ring"`.
#' @param efficacy A single `numeric` between 0 and 1: the proportional
#' reduction in a protected agent's infection hazard.
#' @param start A single non-negative `numeric`: the day (since the start of
#' the outbreak) the programme starts. Ignored if `start_after_reports` is
#' given.
#' @param start_after_reports A single whole number, or `NULL` (default): if
#' given, the programme starts the day after this many cases have been
#' reported.
#' @param delay_to_protection A single non-negative `numeric`: days from
#' vaccination until protection.
#' @param doses_per_day A single non-negative whole number of doses a day, for
#' `"mass"` campaigns.
#' @param coverage A single `numeric` between 0 and 1: the probability each
#' agent in a ring is vaccinated, for `"ring"` vaccination.
#' @param delay A `function` with a single argument returning that many delays
#' in days from a report to the vaccination of its ring, or a single
#' non-negative `numeric` for a fixed delay.
#'
#' @return A `<simulist_vaccination>` object.
#' @export
#'
#' @examples
#' vaccination("ring", efficacy = 0.9, start_after_reports = 20, coverage = 0.8)
#' vaccination("mass", efficacy = 0.85, start = 30, doses_per_day = 500)
vaccination <- function(type = c("mass", "ring"),
                        efficacy,
                        start = 0,
                        start_after_reports = NULL,
                        delay_to_protection = 7,
                        doses_per_day = NULL,
                        coverage = 0.8,
                        delay = 2) {
  type <- match.arg(type)
  checkmate::assert_number(efficacy, lower = 0, upper = 1)
  checkmate::assert_number(start, lower = 0)
  checkmate::assert_count(start_after_reports, positive = TRUE, null.ok = TRUE)
  checkmate::assert_number(delay_to_protection, lower = 0)
  if (type == "mass") checkmate::assert_count(doses_per_day)
  checkmate::assert_number(coverage, lower = 0, upper = 1)
  structure(
    list(
      type = type,
      efficacy = efficacy,
      start = start,
      start_after_reports = start_after_reports,
      delay_to_protection = delay_to_protection,
      doses_per_day = if (is.null(doses_per_day)) 0 else doses_per_day,
      coverage = coverage,
      delay = .as_delay(delay, "delay")
    ),
    class = "simulist_vaccination"
  )
}

#' Contact tracing for the agent-based model
#'
#' @description Describes case isolation and contact tracing for
#' [sim_outbreak_abm()]. Once tracing has started, each reported case is
#' isolated from the day of its report, which multiplies its infectiousness by
#' `1 - isolation_reduction`. Its household members and the agents it infected
#' before it was reported are each traced with probability `prob_traced`,
#' `delay` days after the report, and quarantined for `quarantine_days`. An
#' infected agent in quarantine has its infectiousness multiplied by
#' `1 - quarantine_reduction`, and is reported on the day of onset rather than
#' after the usual reporting delay. Named community contacts in the contacts
#' table are traced with the same probability and delay.
#'
#' @param prob_traced A single `numeric` between 0 and 1: the probability each
#' contact is traced.
#' @param delay A `function` with a single argument returning that many delays
#' in days from a report to tracing its contacts, or a single non-negative
#' `numeric` for a fixed delay.
#' @param quarantine_days A single non-negative `numeric`: how long traced
#' contacts are quarantined.
#' @param quarantine_reduction,isolation_reduction Single `numeric`s between 0
#' and 1: the proportional reduction in the infectiousness of quarantined
#' contacts and isolated cases.
#' @inheritParams vaccination
#'
#' @return A `<simulist_contact_tracing>` object.
#' @export
#'
#' @examples
#' contact_tracing(prob_traced = 0.7, delay = 2, start_after_reports = 20)
contact_tracing <- function(prob_traced = 0.7,
                            delay = function(n) stats::rgamma(n, shape = 2, rate = 1),
                            quarantine_days = 14,
                            quarantine_reduction = 0.8,
                            isolation_reduction = 0.8,
                            start = 0,
                            start_after_reports = NULL) {
  checkmate::assert_number(prob_traced, lower = 0, upper = 1)
  checkmate::assert_number(quarantine_days, lower = 0)
  checkmate::assert_number(quarantine_reduction, lower = 0, upper = 1)
  checkmate::assert_number(isolation_reduction, lower = 0, upper = 1)
  checkmate::assert_number(start, lower = 0)
  checkmate::assert_count(start_after_reports, positive = TRUE, null.ok = TRUE)
  structure(
    list(
      prob_traced = prob_traced,
      delay = .as_delay(delay, "delay"),
      quarantine_days = quarantine_days,
      quarantine_reduction = quarantine_reduction,
      isolation_reduction = isolation_reduction,
      start = start,
      start_after_reports = start_after_reports
    ),
    class = "simulist_contact_tracing"
  )
}

#' A delay as a function returning that many non-negative delays
#' @keywords internal
.as_delay <- function(delay, name) {
  if (is.numeric(delay)) {
    checkmate::assert_number(delay, lower = 0, .var.name = name)
    value <- delay
    return(function(n) rep(value, n))
  }
  .check_func_req_args(delay, func_name = name)
  delay
}
