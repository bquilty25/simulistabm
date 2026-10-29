
<!-- README.md is generated from README.Rmd. Please edit that file. -->

<!-- The code to render this README is stored in .github/workflows/render-readme.yaml -->

<!-- Variables marked with double curly braces will be transformed beforehand: -->

<!-- `packagename` is extracted from the DESCRIPTION file -->

<!-- `gh_repo` is extracted via a special environment variable in GitHub Actions -->

# *simulistabm*: Agent-based line list simulation

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`{simulistabm}` is a fork of
[`{simulist}`](https://github.com/epiverse-trace/simulist) that
simulates outbreaks with an agent-based model built on
[`{individual}`](https://github.com/mrc-ide/individual) and returns the
same line list and contacts table as `{simulist}`. Every `{simulist}`
function is included, so truncating, messing up and censoring the output
works unchanged.

It was written to generate synthetic outbreaks for training and for
testing analysis tools, where the branching process in `{simulist}`
falls short: it has no depletion of susceptibles, no households and no
geography, so outbreaks never peak and control cannot take effect.

## Key features

`sim_outbreak_abm()` adds:

:house: Households within locations, with separate household and
community transmission <br> :busts_in_silhouette: Finite populations, so
outbreaks peak and end <br> :world_map: Locations with populations and a
mobility matrix (gravity or radiation model) for community mixing <br>
:fire: Individual variation in infectiousness, for superspreading <br>
:construction: Interventions that reduce transmission by location and
time <br> :link: A contacts table whose links follow the simulated
transmission chain <br>

`chain_to_outbreak()` turns a transmission chain from any model into
`{simulist}` output, adding delays, outcomes, names, case types and
contact tracing status.

From `{simulist}`, all functions are included, which simulate realistic
line list and contact tracing data with:

:hourglass_flowing_sand: Parameterised epidemiological delay
distributions <br> :hospital: Population-wide or age-stratified
hospitalisation and death risks <br> :bar_chart: Uniform or
age-structured populations <br> :chart_with_upwards_trend: Constant or
time-varying case fatality risk <br> :clipboard: Customisable
probability of case types and contact tracing follow-up <br>

Post-process simulated line list data for:

:date: Real-time outbreak snapshots with right-truncation <br> :memo:
Messy data with inconsistencies, mistakes and missing values <br>
:ledger: Censor dates to daily, epi- and iso-weekly, yearly and other
groupings <br>

## Installation

`{simulistabm}` is not on CRAN. Install it and `{individual}` from
GitHub:

``` r
if (!require("pak")) install.packages("pak")
pak::pak("mrc-ide/individual")
pak::pak("bquilty25/simulistabm")
```

It can be installed alongside `{simulist}`; call functions with
`simulistabm::` to be explicit about which package you are using.

## Quick start: agent-based model

``` r
library(simulistabm)
```

Define the locations, their populations and how much community mixing
stays in each location. Here three towns are placed on a line, and those
who mix outside their own town choose a destination by a gravity model
with a 10 km kernel.

``` r
towns <- data.frame(
  name = c("Ashford", "Brayton", "Cole"),
  population = c(20000, 8000, 5000),
  x = c(0, 12, 30),
  y = c(0, 0, 0)
)
geography <- create_geography(
  towns,
  mobility = mobility_gravity(kernel_km = 10),
  prob_stay = 0.85,
  seed_location = "Ashford"
)
geography
#> <simulist_geography> with 3 locations
#> Total population: 33,000
#> Mean probability of staying in a location: 0.85
#> Index case location: Ashford
```

`sim_outbreak_abm()` simulates transmission among every resident, in
households and the community, then adds delays, outcomes and names as
`sim_outbreak()` does. Hazards are per day of infectiousness. With the
default mean infectious period of about 4.9 days, these give about 1.5
community infections and 0.7 household infections per case, a
reproduction number of about 2.2. The second element of `outbreak_size`
stops the simulation early once infections exceed it.

``` r
set.seed(1)
outbreak <- sim_outbreak_abm(
  geography = geography,
  household_hazard = 0.05,
  community_hazard = 0.3,
  infectiousness = function(n) rgamma(n, shape = 0.5, rate = 0.5),
  outbreak_size = c(50, 2000)
)
#> Warning: Number of cases exceeds maximum outbreak size. 
#> Returning data early with 2114 infections.
head(outbreak$linelist)
#>    id        case_name case_type sex age location date_onset date_reporting
#> 1 137      Anna Turner confirmed   f   0  Ashford 2023-01-04     2023-01-04
#> 2  62    Trevor Garman suspected   m  22  Ashford 2023-01-09     2023-01-09
#> 3 745    Joshua Porter suspected   m  15  Ashford 2023-01-14     2023-01-14
#> 4 678 Yesenia Manchego suspected   f  74  Ashford 2023-01-15     2023-01-15
#> 5 619       Lael Swift  probable   f  22  Ashford 2023-01-18     2023-01-18
#> 6 178    David Eichner confirmed   m  88  Ashford 2023-01-18     2023-01-18
#>   date_admission   outcome date_outcome date_first_contact date_last_contact
#> 1           <NA> recovered         <NA>               <NA>              <NA>
#> 2           <NA> recovered         <NA>         2022-12-30        2023-01-04
#> 3           <NA> recovered         <NA>         2023-01-04        2023-01-14
#> 4     2023-01-19      died   2023-01-27         2023-01-07        2023-01-10
#> 5           <NA> recovered         <NA>         2023-01-12        2023-01-19
#> 6     2023-01-21      died   2023-01-30         2023-01-11        2023-01-17
#>   ct_value
#> 1     22.3
#> 2       NA
#> 3       NA
#> 4       NA
#> 5       NA
#> 6     26.2
head(outbreak$contacts)
#>             from               to age sex date_first_contact date_last_contact
#> 1 Brandon Tsosie       Ivan Moore  19   m         2023-02-26        2023-03-03
#> 2 Brandon Tsosie Bradley Petrarca  68   m         2023-03-01        2023-03-03
#> 3 Brandon Tsosie   Boyd Forrester  47   m         2023-02-28        2023-03-05
#> 4 Brandon Tsosie     Megahn Waite  55   f         2023-02-23        2023-03-05
#> 5 Brandon Tsosie  Angelito Wilson  77   m         2023-02-27        2023-03-03
#> 6 Brandon Tsosie         Al Lopez  59   m         2023-02-27        2023-03-05
#>   was_case           status from_location to_location
#> 1    FALSE   under_followup       Ashford     Ashford
#> 2    FALSE lost_to_followup       Ashford     Ashford
#> 3    FALSE   under_followup       Ashford     Ashford
#> 4    FALSE   under_followup       Ashford     Ashford
#> 5    FALSE   under_followup       Ashford     Ashford
#> 6    FALSE   under_followup       Ashford     Ashford
table(outbreak$linelist$location)
#> 
#> Ashford Brayton    Cole 
#>     799     396     119
```

Daily infections and infectious agents by location are in
`attr(outbreak, "daily")`. Without a cap the outbreak runs until it
ends, depleting each town. Interventions reduce transmission by location
from a start to an end day; here transmission is cut by 60% everywhere
from day 40:

``` r
set.seed(1)
uncontrolled <- sim_outbreak_abm(
  geography = geography,
  household_hazard = 0.05,
  community_hazard = 0.3,
  outbreak_size = c(50, 1e5)
)
set.seed(1)
controlled <- sim_outbreak_abm(
  geography = geography,
  household_hazard = 0.05,
  community_hazard = 0.3,
  interventions = data.frame(location = towns$name, start = 40, end = 200, reduction = 0.6),
  outbreak_size = c(50, 1e5)
)
c(uncontrolled = nrow(uncontrolled$linelist), controlled = nrow(controlled$linelist))
#> uncontrolled   controlled 
#>        28301         2139
```

The output works with the rest of the package, for example a real-time
snapshot that drops the last 14 days, with dates censored to epiweeks:

``` r
snapshot <- truncate_linelist(outbreak$linelist, truncation_day = 14)
head(censor_linelist(snapshot, interval = "epiweek"))
#>    id        case_name case_type sex age location date_onset date_reporting
#> 1 137      Anna Turner confirmed   f   0  Ashford   2023-W01       2023-W01
#> 2  62    Trevor Garman suspected   m  22  Ashford   2023-W02       2023-W02
#> 3 745    Joshua Porter suspected   m  15  Ashford   2023-W02       2023-W02
#> 4 678 Yesenia Manchego suspected   f  74  Ashford   2023-W03       2023-W03
#> 5 619       Lael Swift  probable   f  22  Ashford   2023-W03       2023-W03
#> 6 178    David Eichner confirmed   m  88  Ashford   2023-W03       2023-W03
#>   date_admission   outcome date_outcome date_first_contact date_last_contact
#> 1           <NA> recovered         <NA>               <NA>              <NA>
#> 2           <NA> recovered         <NA>           2022-W52          2023-W01
#> 3           <NA> recovered         <NA>           2023-W01          2023-W02
#> 4       2023-W03      died     2023-W04           2023-W01          2023-W02
#> 5           <NA> recovered         <NA>           2023-W02          2023-W03
#> 6       2023-W03      died     2023-W05           2023-W02          2023-W03
#>   ct_value
#> 1     22.3
#> 2       NA
#> 3       NA
#> 4       NA
#> 5       NA
#> 6     26.2
```

## Quick start: branching process

The rest of this section is from `{simulist}`. A line list can be
simulated by calling `sim_linelist()`. The function provides sensible
defaults to quickly generate a epidemiologically valid data set.

``` r
set.seed(1)
linelist <- sim_linelist()
head(linelist)
#>   id           case_name case_type sex age date_onset date_reporting
#> 1  1    Michael Mckenney confirmed   m  58 2023-01-01     2023-01-01
#> 2  2         Jordan Reed  probable   m  89 2023-01-01     2023-01-01
#> 3  3   Bassil el-Hossain suspected   m   3 2023-01-02     2023-01-02
#> 4  5    Joseph Dickerson suspected   m  28 2023-01-04     2023-01-04
#> 5  6 Naqiyya el-Ghaffari confirmed   f  13 2023-01-05     2023-01-05
#> 6  7     Sergio Sandoval suspected   m  84 2023-01-06     2023-01-06
#>   date_admission   outcome date_outcome date_first_contact date_last_contact
#> 1     2023-01-13      died   2023-01-21               <NA>              <NA>
#> 2           <NA> recovered         <NA>         2022-12-29        2023-01-03
#> 3     2023-01-08      died   2023-01-15         2022-12-28        2023-01-01
#> 4           <NA> recovered         <NA>         2022-12-28        2023-01-04
#> 5           <NA> recovered         <NA>         2022-12-31        2023-01-04
#> 6           <NA> recovered         <NA>         2022-12-31        2023-01-06
#>   ct_value
#> 1     24.8
#> 2       NA
#> 3       NA
#> 4       NA
#> 5     21.2
#> 6       NA
```

However, to simulate a more realistic line list using epidemiological
parameters estimated for an infectious disease outbreak we can use
previously estimated epidemiological parameters. These can be from the
`{epiparameter}` R package if available, or if these are not in the
`{epiparameter}` database yet (such as the contact distribution for
COVID-19) we can define them ourselves. Here we define a contact
distribution, period of infectiousness, onset-to-hospitalisation delay,
and onset-to-death delay.

``` r
library(epiparameter)
```

``` r
# create COVID-19 contact distribution
contact_distribution <- epiparameter::epiparameter(
  disease = "COVID-19",
  epi_name = "contact distribution",
  prob_distribution = create_prob_distribution(
    prob_distribution = "pois",
    prob_distribution_params = c(mean = 2)
  )
)
#> Citation cannot be created as author, year, journal or title is missing

# create COVID-19 infectious period
infectious_period <- epiparameter::epiparameter(
  disease = "COVID-19",
  epi_name = "infectious period",
  prob_distribution = create_prob_distribution(
    prob_distribution = "gamma",
    prob_distribution_params = c(shape = 1, scale = 1)
  )
)
#> Citation cannot be created as author, year, journal or title is missing

# create COVID-19 onset to hospital admission
onset_to_hosp <- epiparameter(
  disease = "COVID-19",
  epi_name = "onset to hospitalisation",
  prob_distribution = create_prob_distribution(
    prob_distribution = "lnorm",
    prob_distribution_params = c(meanlog = 1, sdlog = 0.5)
  )
)
#> Citation cannot be created as author, year, journal or title is missing

# get onset to death from {epiparameter} database
onset_to_death <- epiparameter::epiparameter_db(
  disease = "COVID-19",
  epi_name = "onset to death",
  single_epiparameter = TRUE
)
#> Using Linton N, Kobayashi T, Yang Y, Hayashi K, Akhmetzhanov A, Jung S, Yuan
#> B, Kinoshita R, Nishiura H (2020). "Incubation Period and Other
#> Epidemiological Characteristics of 2019 Novel Coronavirus Infections
#> with Right Truncation: A Statistical Analysis of Publicly Available
#> Case Data." _Journal of Clinical Medicine_. doi:10.3390/jcm9020538
#> <https://doi.org/10.3390/jcm9020538>.. 
#> To retrieve the citation use the 'get_citation' function
```

To simulate a line list for COVID-19 with an Poisson contact
distribution with a mean number of contacts of 2 and a probability of
infection per contact of 0.5, we use the `sim_linelist()` function. The
mean number of contacts and probability of infection determine the
outbreak reproduction number, if the resulting reproduction number is
around one it means we will likely get a reasonably sized outbreak (10 -
1,000 cases, varying due to the stochastic simulation).

***Warning***: the reproduction number of the simulation results from
the contact distribution (`contact_distribution`) and the probability of
infection (`prob_infection`); the number of infections is a binomial
sample of the number of contacts for each case with the probability of
infection (i.e. being sampled) given by `prob_infection`. If the average
number of secondary infections from each primary case is greater than 1
then this can lead to the outbreak becoming extremely large. There is
currently no depletion of susceptible individuals in the simulation
model, so the maximum outbreak size (second element of the vector
supplied to the `outbreak_size` argument) can be used to return a line
list early without producing an excessively large data set.

``` r
set.seed(1)
linelist <- sim_linelist(
  contact_distribution = contact_distribution,
  infectious_period = infectious_period,
  prob_infection = 0.5,
  onset_to_hosp = onset_to_hosp,
  onset_to_death = onset_to_death
)
head(linelist)
#>   id          case_name case_type sex age date_onset date_reporting
#> 1  1 Jennifer Pritchett confirmed   f   0 2023-01-01     2023-01-01
#> 2  2       Tyler Payson confirmed   f  28 2023-01-01     2023-01-01
#> 3  3   Faisal al-Hallal confirmed   m  77 2023-01-01     2023-01-01
#> 4  5    Collin Johnston confirmed   m  69 2023-01-01     2023-01-01
#> 5  6           Kaleb Le  probable   m  27 2023-01-01     2023-01-01
#> 6  8   Muslim el-Radwan suspected   m  60 2023-01-01     2023-01-01
#>   date_admission   outcome date_outcome date_first_contact date_last_contact
#> 1           <NA> recovered         <NA>               <NA>              <NA>
#> 2           <NA> recovered         <NA>         2022-12-30        2023-01-08
#> 3     2023-01-03      died   2023-01-10         2022-12-31        2023-01-05
#> 4           <NA> recovered         <NA>         2022-12-31        2023-01-04
#> 5           <NA> recovered         <NA>         2022-12-29        2023-01-02
#> 6     2023-01-05 recovered         <NA>         2022-12-28        2023-01-05
#>   ct_value
#> 1     27.1
#> 2     25.2
#> 3     24.8
#> 4     26.2
#> 5       NA
#> 6       NA
```

In this example, the line list is simulated using the default values
(see `?sim_linelist`). The default hospitalisation risk is assumed to be
0.2 (i.e. there is a 20% probability an infected individual becomes
hospitalised) and the start date of the outbreak is 1st January 2023. To
modify either of these, we can specify them in the function.

``` r
linelist <- sim_linelist(
  contact_distribution = contact_distribution,
  infectious_period = infectious_period,
  prob_infection = 0.5,
  onset_to_hosp = onset_to_hosp,
  onset_to_death = onset_to_death,
  hosp_risk = 0.01,
  outbreak_start_date = as.Date("2019-12-01")
)
head(linelist)
#>   id         case_name case_type sex age date_onset date_reporting
#> 1  1    Karlena Garner confirmed   f  31 2019-12-01     2019-12-01
#> 2  3 Cassandra Le Beau suspected   f  47 2019-12-01     2019-12-01
#> 3  4      Lauren Moran confirmed   f   2 2019-12-01     2019-12-01
#> 4  5       Alexis Ware confirmed   f  24 2019-12-01     2019-12-01
#> 5  6       Sophia Cruz  probable   f  56 2019-12-01     2019-12-01
#> 6  7      Amanda Quast confirmed   f   8 2019-12-01     2019-12-01
#>   date_admission   outcome date_outcome date_first_contact date_last_contact
#> 1           <NA> recovered         <NA>               <NA>              <NA>
#> 2           <NA> recovered         <NA>         2019-11-28        2019-12-03
#> 3           <NA>      died   2019-12-25         2019-11-28        2019-12-07
#> 4           <NA> recovered         <NA>         2019-11-26        2019-12-02
#> 5           <NA> recovered         <NA>         2019-11-28        2019-12-04
#> 6           <NA>      died   2019-12-16         2019-11-30        2019-12-05
#>   ct_value
#> 1     26.3
#> 2       NA
#> 3     26.2
#> 4     28.5
#> 5       NA
#> 6     28.6
```

To simulate a table of contacts of cases (i.e. to reflect a contact
tracing dataset) we can use the same parameters defined for the example
above.

``` r
contacts <- sim_contacts(
  contact_distribution = contact_distribution,
  infectious_period = infectious_period,
  prob_infection = 0.5
)
head(contacts)
#>                from                  to age sex date_first_contact
#> 1 Ivory Bing-Turner       Umar al-Edris  63   m         2022-12-30
#> 2 Ivory Bing-Turner           Scott Kim  10   m         2022-12-30
#> 3 Ivory Bing-Turner     Alyssa Harrison  24   f         2022-12-28
#> 4 Ivory Bing-Turner Nhuquynh Kodicherla  10   f         2022-12-29
#> 5 Ivory Bing-Turner      Carmen Montoya  56   f         2023-01-01
#> 6         Scott Kim          Ryan Moore   4   m         2022-12-31
#>   date_last_contact was_case         status
#> 1        2023-01-05    FALSE under_followup
#> 2        2023-01-04     TRUE           case
#> 3        2023-01-04     TRUE           case
#> 4        2023-01-02     TRUE           case
#> 5        2023-01-04     TRUE           case
#> 6        2023-01-07    FALSE under_followup
```

If both the line list and contacts table are required, they can be
jointly simulated using the `sim_outbreak()` function. This uses the
same inputs as `sim_linelist()` and `sim_contacts()` to produce a line
list and contacts table of the same outbreak (the arguments also have
the same default settings as the other functions).

``` r
outbreak <- sim_outbreak(
  contact_distribution = contact_distribution,
  infectious_period = infectious_period,
  prob_infection = 0.5,
  onset_to_hosp = onset_to_hosp,
  onset_to_death = onset_to_death
)
head(outbreak$linelist)
#>   id       case_name case_type sex age date_onset date_reporting date_admission
#> 1  1  Steffano Jones confirmed   m  87 2023-01-01     2023-01-01           <NA>
#> 2  2     Mahmoud Lee confirmed   m  17 2023-01-01     2023-01-01     2023-01-04
#> 3  4 Jaicy Collareno  probable   f  68 2023-01-01     2023-01-01           <NA>
#> 4  5  Chandler Allen confirmed   f  34 2023-01-01     2023-01-01           <NA>
#> 5  8   Destiny Crown  probable   f  47 2023-01-01     2023-01-01           <NA>
#> 6  9 Jacob Littleman  probable   m  72 2023-01-01     2023-01-01     2023-01-05
#>     outcome date_outcome date_first_contact date_last_contact ct_value
#> 1 recovered         <NA>               <NA>              <NA>     24.8
#> 2 recovered         <NA>         2022-12-29        2023-01-03     24.7
#> 3 recovered         <NA>         2022-12-29        2023-01-04       NA
#> 4 recovered         <NA>         2022-12-31        2023-01-01     25.4
#> 5 recovered         <NA>         2022-12-30        2023-01-03       NA
#> 6      died   2023-01-22         2022-12-28        2023-01-01       NA
head(outbreak$contacts)
#>              from                to age sex date_first_contact
#> 1  Steffano Jones       Mahmoud Lee  17   m         2022-12-29
#> 2  Steffano Jones   Cooper Straight  11   m         2022-12-26
#> 3     Mahmoud Lee   Jaicy Collareno  68   f         2022-12-29
#> 4 Jaicy Collareno    Chandler Allen  34   f         2022-12-31
#> 5 Jaicy Collareno Cheyenne Trujillo  58   f         2022-12-29
#> 6 Jaicy Collareno   Cameron Barajas  22   m         2022-12-29
#>   date_last_contact was_case           status
#> 1        2023-01-03     TRUE             case
#> 2        2023-01-05    FALSE   under_followup
#> 3        2023-01-04     TRUE             case
#> 4        2023-01-01     TRUE             case
#> 5        2023-01-04    FALSE lost_to_followup
#> 6        2023-01-02    FALSE   under_followup
```

## Help

To report a bug in the agent-based model or other additions, please open
an [issue](https://github.com/bquilty25/simulistabm/issues). Bugs in
functions that `{simulistabm}` shares with `{simulist}` can be reported
to
[`{simulist}`](https://github.com/epiverse-trace/simulist/issues/new/choose).

## Contribute

Contributions are welcome. Please follow the [Epiverse-TRACE
contributing
guide](https://github.com/epiverse-trace/.github/blob/main/CONTRIBUTING.md).

## Code of Conduct

Please note that `{simulistabm}`, like `{simulist}`, is released with a
[Contributor Code of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## Citing this package

`{simulistabm}` builds on `{simulist}` and `{individual}`. Please cite
`{simulist}`’s publication in the Journal of Open Source Software
(JOSS):

``` r
citation("simulistabm")
#> To cite package 'simulistabm' in publications use:
#> 
#>   Lambert J, Kucharski A, Tamayo Cuartero C (2026). "simulist: An R
#>   package to simulate disease outbreak line list and contacts data."
#>   _Journal of Open Source Software_, *11*(119), 8781.
#>   doi:10.21105/joss.08781 <https://doi.org/10.21105/joss.08781>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {simulist: An R package to simulate disease outbreak line list and contacts data},
#>     author = {Joshua W. Lambert and Adam Kucharski and Carmen {Tamayo Cuartero}},
#>     journal = {Journal of Open Source Software},
#>     year = {2026},
#>     volume = {11},
#>     number = {119},
#>     pages = {8781},
#>     doi = {10.21105/joss.08781},
#>     publisher = {The Open Journal},
#>   }
```

## Complimentary R packages

:package: :left_right_arrow: :package:
[{epiparameter}](https://epiverse-trace.github.io/epiparameter/) <br>
:package: :left_right_arrow: :package:
[{epicontacts}](https://www.repidemicsconsortium.org/epicontacts/) <br>
:package: :left_right_arrow: :package:
[{incidence2}](https://www.reconverse.org/incidence2/) <br> :package:
:left_right_arrow: :package:
[{cleanepi}](https://epiverse-trace.github.io/cleanepi/) <br>

## Related projects

This project has some overlap with other R packages. Here we list these
packages and provide a table of features and attributes that are present
for each package to help decide which package is appropriate for each
use-case.

In some cases the packages are dedicated to simulating line list and
other epidemiological data (e.g. {simulist}), in others the line list
simulation is one part of a wider R package (e.g. {EpiNow}).

- [`{LLsim}`](https://github.com/jrcpulliam/LLsim) simulates line list
  data using a stochastic SIR model with a fixed population with
  observation and reporting delays. Line list data is generated in two
  steps, 1) the SIR model simulates the outbreak (`simpleSim()`), 2) the
  outbreak data is converted into a line list (`createLineList()`).
- [`{simulacr}`](https://github.com/reconhub/simulacr) uses a branching
  process model to simulate cases and contacts for an outbreak. It
  simulates transmission of infections using other epidemiological R
  packages (`{epicontacts}` and `{distcrete}`) to parameterise and plot
  simulated data.
- [`{epidict}`](https://github.com/R4EPI/epidict) is a package that can
  be used to simulate outbreak data, including line lists, in a DHIS2
  format, and survey data that mimics the format by Kobo, using the
  function `gen_data()`. In addition, MSF outbreak data are available in
  this package as data dictionaries for Acute Jaundice Syndrome,
  Cholera, Measles and Meningitis, accessible through the function
  `msf_dict()`.
- [`{EpiNow}`](https://github.com/epiforecasts/EpiNow) - a now
  deprecated R package - includes the `simulate_cases()` and
  `generate_pseudo_linelist()` functions for generating line list data.
- [generative-nowcasting](https://github.com/adrian-lison/generative-nowcasting)
  is a set of R scripts and functions to perform epidemiological
  nowcasting. There are [functions to simulate line list
  data](https://github.com/adrian-lison/generative-nowcasting/blob/bf48e027e82ce9d42de468d0b708d010253b7475/code/utils/utils_simulate.R)
  within the repository, but the repository is not (and does not
  contain) an R package. Functions can be sourced. Cases are simulated
  with a renewal process and the simulation can incorporate
  epidemiological delays and ascertainment.

<details>

<summary>

Table of line list simulator features
</summary>

|  | {simulist} | {LLsim} | {simulacr} | {epidict} | {EpiNow} | generative-nowcasting |
|----|----|----|----|----|----|----|
| Simulates line list | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| Simulates contacts | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: |
| Parameterised with epi distributions[^1] | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: | :white_check_mark: | :white_check_mark: |
| Interoperable with {epicontacts} | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: |
| Explicit population size[^2] | :x: | :white_check_mark: | :white_check_mark: | :x: | :x: | :x: |
| R package | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: |
| Actively maintained[^3] | :white_check_mark: | :x: | :x: | :x: | :x: | :white_check_mark: |
| On CRAN | :white_check_mark: | :x: | :x: | :x: | :x: | NA |
| Unit testing[^4] | :white_check_mark: | :white_check_mark: | :x: | :white_check_mark: | :x: | NA |

</details>

If there is another package with this functionality missing from the
list that should be added, or if a package included in this list has
been updated and the table should reflect this please contribute by
making an [issue](https://github.com/epiverse-trace/simulist/issues) or
a [pull request](https://github.com/epiverse-trace/simulist/pulls).

Other R packages for simulating epidemic dynamics can be found in the
*Epidemic simulation models* section of the [Epidemiology CRAN task
view](https://CRAN.R-project.org/view=Epidemiology).

Some packages are related to {simulist} but do not simulate line list
data. These include:

- [`{outbreaks}`](https://CRAN.R-project.org/package=outbreaks) an R
  package containing a library of outbreak data sets, including line
  list data, for a variety of past and simulated outbreaks, e.g. Ebola
  and MERS.
- [`{ringbp}`](https://github.com/epiforecasts/ringbp) an R package to
  simulate cases using an individual-level transmission model with
  contact tracing.
- [`{epichains}`](https://github.com/epiverse-trace/epichains) an R
  package with functionality to simulate transmission chains using a
  branching process model.

The {outbreaks} package is useful if data from a past outbreak data or
generic line list data is required. The {ringbp} and {epichains}
packages can be used to generate case data over time which can then be
converted into a line list with some manual post-processing.

Another package for creating messy data is the
[{messy}](https://CRAN.R-project.org/package=messy) package. This can be
used, either independently or in combination with `messy_linelist()`, to
create messy line list and contacts data.

[^1]: In this context *Parameterised with epi distributions* means that
    the simulation uses epidemiological distributions (e.g. serial
    interval, infectious period) to parameterise the model and the
    parameters of these epi distributions can be modified by the user.

[^2]: *Explicit population size* refers to the simulation using a finite
    population size which is controlled by the user for the depletion of
    susceptible individuals in the model.

[^3]: We define *Actively maintained* as the repository having a commit
    to the main branch within the last 12 months.

[^4]: *Unit testing* is ticked if the package contains any form of
    testing, this can use any testing framework, for example
    [{testthat}](https://CRAN.R-project.org/package=testthat) or
    [{tinytest}](https://CRAN.R-project.org/package=tinytest).
