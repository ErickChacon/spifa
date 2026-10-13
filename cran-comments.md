## Resubmission

This is a resubmission. Per the reviewer's feedback:

* Added a `\value` tag to `dic.Rd` (the generic `dic()`, which previously
  had none).
* All vignettes now reset `options(width = ...)` to its prior value after
  changing it, via `old_options <- options(width = 100)` in the setup
  chunk and `options(old_options)` in a hidden chunk at the end.

## Test environments

* Local: Ubuntu 24.04.4 LTS, R 4.6.0 (2026-04-24), via
  `devtools::check(cran = TRUE, remote = TRUE)`
* GitHub Actions: R CMD check (Linux, macOS, Windows), via
  `.github/workflows/R-CMD-check.yaml`
* win-builder (R-devel)
* R-hub (Linux, Windows; R-devel), via `.github/workflows/rhub.yaml`

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* On win-builder, DESCRIPTION's spell-check flags "Giorgi", "Moraga", and
  "Orellana" as possibly misspelled -- these are co-authors' surnames from
  the citation, not misspellings.

## Downstream dependencies

This is a new package, so there are no downstream dependencies to check.
