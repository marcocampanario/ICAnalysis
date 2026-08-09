#' Annotate variant novelty
#'
#' Classifies split ICA allele records according to their representation in
#' population-frequency databases and dbSNP.
#'
#' A variant is classified as:
#'
#' * `"known"` when its allele frequency is greater than zero in gnomAD
#'   Genomes, gnomAD Exomes, 1000 Genomes, or ABRaOM;
#' * `"absent"` when it has no positive population frequency but has at least
#'   one dbSNP identifier;
#' * `"novel"` when it has neither a positive population frequency nor a
#'   dbSNP identifier.
#'
#' @details
#' gnomAD Genomes and gnomAD Exomes are evaluated independently. The
#' convenience column `af_gnomad` is deliberately not used for novelty
#' classification.
#'
#' Missing population frequencies, including an `af_abraom` column containing
#' only missing values do not contribute to `"known"` status. Zero frequency
#' is treated as absence from that population database.
#'
#' If `variant_novelty` already exists, it is recalculated and replaced.
#'
#' @param alleles An `S4Vectors::DataFrame` containing split ICA allele
#'   records, normally returned by [split_ica_alleles()].
#'
#' @return The input `S4Vectors::DataFrame` with an additional character
#'   column named `variant_novelty`.
#'
#' @examples
#' alleles <- S4Vectors::DataFrame(
#'     variant_id = c("chr1-100-A-C", "chr1-200-G-T", "chr1-300-C-G"),
#'     af_gnomad_genome = c(0.01, NA_real_, NA_real_),
#'     af_gnomad_exome = rep(NA_real_, 3L),
#'     af_1000g = rep(NA_real_, 3L),
#'     af_abraom = rep(NA_real_, 3L),
#'     dbsnp_id = c(NA_character_, "rs123", NA_character_)
#' )
#'
#' alleles <- annotate_variant_novelty(alleles)
#' alleles$variant_novelty
#'
#' @export
annotate_variant_novelty <- function(alleles) {
  if (!methods::is(alleles, "DataFrame")) {
    stop(
      "`alleles` must be an S4Vectors::DataFrame.",
      call. = FALSE
    )
  }

  frequency_columns <- c(
    "af_gnomad_genome",
    "af_gnomad_exome",
    "af_1000g",
    "af_abraom"
  )

  required_columns <- c(
    frequency_columns,
    "dbsnp_id"
  )

  missing_columns <- setdiff(
    required_columns,
    names(alleles)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "`alleles` is missing required column(s): ",
        paste(missing_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  numeric_frequency <- vapply(
    frequency_columns,
    function(column) {
      is.numeric(alleles[[column]])
    },
    logical(1L)
  )

  if (any(!numeric_frequency)) {
    invalid_columns <- frequency_columns[!numeric_frequency]

    stop(
      paste0(
        "Population-frequency columns must be numeric: ",
        paste(invalid_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (!is.character(alleles$dbsnp_id)) {
    stop(
      "`dbsnp_id` must be a character vector.",
      call. = FALSE
    )
  }

  has_positive_population_af <- rep(
    FALSE,
    nrow(alleles)
  )

  for (column in frequency_columns) {
    frequency <- alleles[[column]]

    has_positive_population_af <-
      has_positive_population_af |
      (
        is.finite(frequency) &
          frequency > 0
      )
  }

  dbsnp_id <- trimws(alleles$dbsnp_id)

  has_dbsnp_id <- !is.na(dbsnp_id) &
    nzchar(dbsnp_id) &
    dbsnp_id != "." &
    toupper(dbsnp_id) != "NA"

  variant_novelty <- rep(
    "novel",
    nrow(alleles)
  )

  variant_novelty[has_dbsnp_id] <- "absent"
  variant_novelty[has_positive_population_af] <- "known"

  alleles$variant_novelty <- variant_novelty
  alleles
}
