#' Validate an ICAAlleles object
#'
#' Performs shallow structural validation of an `ICAAlleles` object. It
#' verifies the required allele-level columns, their basic types, and valid
#' positive indexes and genomic positions.
#'
#' Nested annotation contents and consistency between variant identifiers,
#' coordinates, and allele indexes are deliberately not inspected.
#'
#' @param object An `ICAAlleles` object.
#'
#' @return `TRUE` when the object is valid, or a character vector describing
#'   validation problems.
#'
#' @noRd
#'
#' @keywords internal
.validate_ica_alleles <- function(object) {
  allele_data <- object@alleles

  required_types <- c(
    site_index = "integer",
    allele_index = "integer",
    chromosome = "character",
    position = "integer",
    refAllele = "character",
    altAllele = "character",
    variant_id = "character",
    af_gnomad_genome = "numeric",
    af_gnomad_exome = "numeric",
    af_gnomad = "numeric",
    af_1000g = "numeric",
    af_abraom = "numeric",
    dbsnp_id = "character"
  )

  column_occurrences <- vapply(
    names(required_types),
    function(column) {
      sum(names(allele_data) == column)
    },
    integer(1L)
  )

  invalid_occurrences <- names(
    column_occurrences[column_occurrences != 1L]
  )

  if (length(invalid_occurrences) > 0L) {
    return(
      paste0(
        "The `alleles` slot must contain exactly one column ",
        "named: ",
        paste(invalid_occurrences, collapse = ", "),
        "."
      )
    )
  }

  problems <- character()

  for (column in names(required_types)) {
    expected_type <- required_types[[column]]
    value <- allele_data[[column]]

    valid_type <- switch(
      expected_type,
      integer = is.integer(value),
      numeric = is.numeric(value),
      character = is.character(value),
      FALSE
    )

    if (!valid_type) {
      problems <- c(
        problems,
        sprintf(
          "Column `%s` must be %s.",
          column,
          expected_type
        )
      )
    }
  }

  if (is.integer(allele_data$site_index)) {
    if (anyNA(allele_data$site_index)) {
      problems <- c(
        problems,
        "Column `site_index` cannot contain missing values."
      )
    }

    if (any(allele_data$site_index < 1L, na.rm = TRUE)) {
      problems <- c(
        problems,
        "Column `site_index` must contain positive integers."
      )
    }
  }

  if (is.integer(allele_data$allele_index)) {
    if (anyNA(allele_data$allele_index)) {
      problems <- c(
        problems,
        "Column `allele_index` cannot contain missing values."
      )
    }

    if (any(allele_data$allele_index < 1L, na.rm = TRUE)) {
      problems <- c(
        problems,
        "Column `allele_index` must contain positive integers."
      )
    }
  }

  if (is.integer(allele_data$position)) {
    if (anyNA(allele_data$position)) {
      problems <- c(
        problems,
        "Column `position` cannot contain missing values."
      )
    }

    if (any(allele_data$position < 1L, na.rm = TRUE)) {
      problems <- c(
        problems,
        "Column `position` must contain positive integers."
      )
    }
  }

  if (length(problems) == 0L) {
    TRUE
  } else {
    problems
  }
}


#' ICA alternative allele records
#'
#' The `ICAAlleles` class stores allele-level records derived from an
#' `ICAPositions` object. Each row represents one alternative allele, while
#' the original ICA header and analysis metadata are retained.
#'
#' The class deliberately differs from `ICAPositions`, where each row
#' represents one unsplit VCF/ICA site.
#'
#' @name ICAAlleles-class
#'
#' @slot alleles An `S4Vectors::DataFrame` containing one row per alternative
#'   allele.
#' @slot header A list containing the original ICA header.
#' @slot metadata A list containing provenance and analysis metadata.
#'
#' @importClassesFrom S4Vectors DataFrame
#' @export
methods::setClass(
  "ICAAlleles",
  slots = c(
    alleles = "DataFrame",
    header = "list",
    metadata = "list"
  ),
  prototype = list(
    alleles = S4Vectors::DataFrame(
      site_index = integer(),
      allele_index = integer(),
      chromosome = character(),
      position = integer(),
      refAllele = character(),
      altAllele = character(),
      variant_id = character(),
      af_gnomad_genome = numeric(),
      af_gnomad_exome = numeric(),
      af_gnomad = numeric(),
      af_1000g = numeric(),
      af_abraom = numeric(),
      dbsnp_id = character()
    ),
    header = list(),
    metadata = list()
  ),
  validity = .validate_ica_alleles
)
